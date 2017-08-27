
#include "linear.hpp"
#include "random.hpp"
#include "colour.hpp"
#include "image.hpp"
#include "scene.hpp"
#include "photon.hpp"
#include "photonmapping.hpp"
//#include "directrt.hpp"
#include "ppm.hpp"

#include <thread>
#include <functional>
#include <mutex>
#include <chrono>
#include <fenv.h>

namespace RT {
  struct Options {
  //bool  indirect;
    int /*rtSamplesIndirect,
          rtSamplesDirect,
          directBranch,
          nPhotons,
          nNears,*/
          sqrtSamples,
          nPasses,
          nPhotonsPerPass,
          nThreads;
    float //nearLimit,
          initRadius,
          alpha,
          exposeK;
  };

  class SampleSet {
    int const sqrtn;
    Vector* const storage;

  public:
    SampleSet (int sqrtn) : sqrtn(sqrtn), storage(new Vector [sqrtn*sqrtn]) { }
    ~SampleSet () { delete[] storage; }

    void regen (RandBits& rng) {
      int ys[sqrtn][sqrtn],
          xs[sqrtn][sqrtn];
      for (int i = 0; i != sqrtn; i++) {
        std::iota    (ys[i], ys[i]+sqrtn, 0);
        std::shuffle (ys[i], ys[i]+sqrtn, rng);
        std::iota    (xs[i], xs[i]+sqrtn, 0);
        std::shuffle (xs[i], xs[i]+sqrtn, rng);
      }

      float const k = 1.f / float(sqrtn);
      for (int j = 0; j != sqrtn; j++) {
        for (int i = 0; i != sqrtn; i++) {
          int const dx = xs[i][j],
                    dy = ys[j][i];
          float const x = k*(float(i) + k*(float(dx) + .5f)) - .5f,
                      y = k*(float(j) + k*(float(dy) + .5f)) - .5f;
          storage[j*sqrtn + i] = {x,y,0,0};
        }
      }
    }

    Vector const* begin () const { return storage; }
    Vector const* end   () const { return storage + size(); }
    size_t size () const { return sqrtn*sqrtn; }

    Vector operator[] (int i) { return storage[i]; }
  };

/*template<typename Term>
  Colour gather
    ( Term&          term
    , Options const& opts
    , RandBits&      rng
    , Scene const&   scene
    , Ray            ray
    , Item const*
    , int            ttl
    , Colour         sourceK = {1,1,1}
    )
  {
    if (ttl == 0)
      return black;

    RayHit<Item const*> hit = scene.traceRay (ray);
    if (!hit)
      return black;
    Item const* item = hit.occ;
    Material const* mat = item->material;
    Basis const tang = Basis::fromK (hit.normal);
    Vector const incident = tang.into (ray.disp);
    Interaction const inter = mat->interact (rng, incident, sourceK);

    if (inter.type == Interaction::Type::specular) {
      Vector const bounce = tang.outOf (mat->bounceSpecular (rng, incident));
      Ray const newRay (hit.position, bounce);
      return gather
        ( term
        , opts, rng
        , scene, newRay, item
        , ttl-1
        , sourceK * inter.correctRefl
        );
    }
    else {
      return sourceK * term (hit);
    }
  }

  template<typename MkTerm>
  void eval
    ( MkTerm& mkTerm
    , Options const& opts
    , Scene const& scene
    , Camera const& camera
    , Image<float>& im
    , int nThreads
    , int sqrtNSamples
    )
  {
    float const
      xImToCam = 2.f / im.width (),
      yImToCam = 2.f / im.height (),
      aspect = (float) im.width () / im.height ();

    std::thread threads[nThreads];
    std::mutex mutex;
    int yProgress = 0;

    int const nSamples = sqrtNSamples * sqrtNSamples;

    for (int iThread = 0; iThread != nThreads; iThread++) {
      auto& thread = threads[iThread];
      thread = std::thread (
        [&, iThread] () {
          RandBits rng = threadRNG ();
          auto term = mkTerm (rng);
          SampleSet samplePosns (sqrtNSamples);

          for (;;) {
            int yPixel;
            { std::lock_guard<std::mutex> lock (mutex);
              int const progress
                = int (float(yProgress)*100.f / float(im.height()));
              fprintf (stderr, "\r%3d%%", progress);
              if (yProgress == im.height())
                return;
              yPixel = yProgress++;
            }

            float const yPixelCentre = float(yPixel) + 0.5f;

            for (int xPixel = 0; xPixel != im.width(); xPixel++) {
              float const xPixelCentre = float(xPixel) + 0.5f;

              samplePosns.regen (rng);
              Colour samples[nSamples];
              for (int i = 0; i != nSamples; i++) {
                Vector const samplePos = samplePosns[i];
                float const
                  xSample = (xPixelCentre + samplePos.x)*xImToCam - 1.f,
                  ySample = (yPixelCentre + samplePos.y)*yImToCam - 1.f;
                Ray const ray = camera.shoot (aspect * xSample, ySample);
                samples[i] = gather (term, opts, rng, scene, ray, nullptr, 1000);
              }

              Colour const
                sum = std::accumulate (samples, samples + nSamples, Colour (black)),
                mean = sum / nSamples,
                exposed = mean * opts.exposeK;
              im.at (xPixel, yPixel) = {exposed.r, exposed.g, exposed.b};
            }
          }
        }
      );
    }

    for (auto& thread : threads)
      thread.join ();

    fprintf (stderr, "\r");
  }

  void evalIndirect
    ( Options const& opts
    , Scene const& scene
    , Camera const& cam
    , PhotonMap const& map
    , Image<float>& im
    )
  {
    int const nSamples = opts.rtSamplesIndirect * opts.rtSamplesIndirect;
    fprintf (stderr, "Evaluating indirect lighting at %i spp...\n", nSamples);
    auto const mkTerm = [&] (RandBits&) {
      NearSet nears (opts.nNears, opts.nearLimit);
      return [&, ns = std::move (nears)] (auto hit) mutable {
        return indirectTermPM (ns, map, hit);
      };
    };
    eval (mkTerm, opts, scene, cam, im, opts.nThreads, opts.rtSamplesIndirect);
  }

  void evalDirect
    ( Options const& opts
    , Scene const& scene
    , Camera const& cam
    , Image<float>& im
    )
  {
    int const nSamples = opts.rtSamplesDirect * opts.rtSamplesDirect;
    fprintf (stderr, "Evaluating direct lighting at %i spp...\n", nSamples);
    auto const mkTerm = [&] (RandBits& rng) {
      return [&] (auto hit) {
        return directTermRT (opts.directBranch, scene, rng, hit);
      };
    };
    eval (mkTerm, opts, scene, cam, im, opts.nThreads, opts.rtSamplesDirect);
  }

  void render
    ( Scene const& scene
    , Camera const& cam
    , Image<uint8_t>& im
    , Options const& opts
    )
  {
    int const
    //nSamples = opts.sqrtNSamples * opts.sqrtNSamples,
      nPhotonsPerThread = opts.nPhotons/opts.nThreads;

    fprintf (stderr, "Rendering to %ix%i\n", im.width (), im.height ());

    if (opts.directBranch > 1)
      fprintf (stderr, "Direct light oversample: %i\n", opts.directBranch);

    if (opts.indirect) {
      fprintf (stderr,
        "Using %i photons, %i per estimate\n",
        opts.nPhotons, opts.nNears);
    }

    using clock = std::chrono::steady_clock;
    auto startTime = clock::now ();

    std::thread threads[opts.nThreads];

    Image<float> indirect (im.width (), im.height ());
    if (opts.indirect) {
      fprintf (stderr, "Casting photons...\n");

      std::vector<Photon> photonsRaw (opts.nPhotons);

      for (int iThread = 0; iThread != opts.nThreads; iThread++) {
        auto& thread = threads[iThread];
        thread = std::thread (
          [&, iThread] () {
            RandBits rng = threadRNG ();
            Photon* photons = photonsRaw.data () + iThread * nPhotonsPerThread;
            castPhotons (rng, scene, photons, photons + nPhotonsPerThread);
          }
        );
      }

      for (auto& thread : threads)
        thread.join ();

      fprintf (stderr, "Building photon map...\n");
      PhotonMap map (std::move (photonsRaw));

      evalIndirect (opts, scene, cam, map, indirect);
    }

    Image<float> direct (im.width (), im.height ());
    evalDirect (opts, scene, cam, direct);

    auto stopTime = clock::now ();
    auto duration = stopTime - startTime;
    auto seconds = std::chrono::duration<float> (duration);

    fprintf (stderr, "\rTook %.1f s\n", seconds.count ());

    fprintf (stderr, "Compositing...\n");
    auto& comp = indirect += direct;

    for (int y = 0; y != im.height (); y++) {
      for (int x = 0; x != im.width (); x++)
        im.at (x, y) = toSRGB (comp.at (x, y));
    }
  }*/

  void tracePath
    ( std::vector<PathNode>& nodes
    , float initR2
    , RandBits& rng
    , Scene const& scene
    , int xPixel, int yPixel
    , Ray ray, int ttl
    , Colour sourceK = {1,1,1}
    )
  {
    if (ttl == 0)
      return;

    RayHit<Item const*> hit = scene.traceRay (ray);
    if (!hit)
      return;
    Item const* item = hit.occ;
    Material const* mat = item->material;
    Basis const tang = Basis::fromK (hit.normal);
    Vector const incident = tang.into (ray.disp);
    Interaction const inter = mat->interact (rng, incident, sourceK);

    if (mat->kDiffuse != black) {
      PathNode node
        ( hit.position, hit.normal, incident
        , mat
        , xPixel, yPixel
        , sourceK
        , initR2
        );
      nodes.push_back (node);
    }

    if (inter.type == Interaction::Type::specular) {
      Vector const bounce = tang.outOf (mat->bounceSpecular (rng, incident));
      Ray const newRay (hit.position, bounce);
      tracePath
        ( nodes
        , initR2
        , rng
        , scene
        , xPixel, yPixel
        , newRay, ttl-1
        , sourceK * inter.correctRefl
        );
    }
  }

  std::vector<PathNode> tracePaths
    ( Options const& opts
    , Scene const& scene
    , Camera const& camera
    , Image<float> const& im
    )
  {
    fprintf (stderr, "Tracing paths\n");

    float const
      xImToCam = 2.f / im.width (),
      yImToCam = 2.f / im.height (),
      aspect = (float) im.width () / im.height (),
      initR2 = opts.initRadius * opts.initRadius;

    int const nSamples = opts.sqrtSamples * opts.sqrtSamples;

    RandBits rng = threadRNG ();
    SampleSet samplePosns (opts.sqrtSamples);

    std::vector<PathNode> nodes;
    nodes.reserve (nSamples * im.size ());

    for (int yPixel = 0; yPixel != im.height (); yPixel++) {
      int const progress
        = int (yPixel * 100.f / im.height ());
      fprintf (stderr, "\r%3d%%", progress);

      float const yPixelCentre = float(yPixel) + 0.5f;

      for (int xPixel = 0; xPixel != im.width(); xPixel++) {
        float const xPixelCentre = float(xPixel) + 0.5f;

        samplePosns.regen (rng);
        for (int i = 0; i != nSamples; i++) {
          Vector const samplePos = samplePosns[i];
          float const
            xSample = (xPixelCentre + samplePos.x)*xImToCam - 1.f,
            ySample = (yPixelCentre + samplePos.y)*yImToCam - 1.f;
          Ray const ray = camera.shoot (aspect * xSample, ySample);
          tracePath (nodes, initR2, rng, scene, xPixel, yPixel, ray, 1);
        }
      }
    }

    fprintf (stderr, "\r");
    return nodes;
  }

  void render
    ( Scene const& scene
    , Camera const& cam
    , Image<uint8_t>& im
    , Options const& opts
    )
  {
    Image<float> frame (im.width (), im.height ());
    auto nodes = tracePaths (opts, scene, cam, frame);

    int const
      nPhotonsPerThread = opts.nPhotonsPerPass / opts.nThreads,
      nSamples = opts.sqrtSamples * opts.sqrtSamples;

    std::thread threads[opts.nThreads];
    std::mutex mutex;

    std::vector<Photon> photonsRaw (opts.nPhotonsPerPass);

    for (int iPass = 0; iPass != opts.nPasses; iPass++) {
      fprintf (stderr, "Pass %i/%i\n", iPass+1, opts.nPasses);

      // cast
      fprintf (stderr, "Casting...");
      for (int iThread = 0; iThread != opts.nThreads; iThread++) {
        auto& thread = threads[iThread];
        thread = std::thread (
          [&, iThread] {
            RandBits rng = threadRNG ();
            Photon* photons = photonsRaw.data () + iThread * nPhotonsPerThread;
            castPhotons (rng, scene, photons, photons + nPhotonsPerThread);
          }
        );
      }

      for (auto& thread : threads)
        thread.join ();

      // build tree
      fprintf (stderr, "\rBuilding tree...");
      PhotonMap map (std::move (photonsRaw));

      // update
      fprintf (stderr, "\rUpdating estimates...");
      int const nNodes = nodes.size ();
      for (int iThread = 0; iThread != opts.nThreads; iThread++) {
        auto& thread = threads[iThread];
        thread = std::thread (
          [&, iThread] {
            int const
              begin = (nNodes *  iThread   ) / opts.nThreads,
              end   = (nNodes * (iThread+1)) / opts.nThreads;
            for (int i = begin; i != end; i++)
              nodes[i].update (map, opts.alpha);
          }
        );
      }

      for (auto& thread : threads)
        thread.join ();

      // recycle
      photonsRaw = std::move (map.array);
      fprintf (stderr, "\r");
    }

    // finalize
    for (auto const& node : nodes) {
      auto L = node.radiance ();
      //fprintf (stderr, "L: (%.3f %.3f %.3f)\n", L.r, L.g, L.b);
      frame.at (node.px, node.py) += {L.r, L.g, L.b};
    }

    // ramp
    for (int y = 0; y != im.height (); y++) {
      for (int x = 0; x != im.width (); x++)
        im.at (x, y) = toSRGB (frame.at (x, y));
    }
  }

  void writePPM (char const* path, Image<uint8_t> const& im) {
    FILE* f = fopen (path, "wb");
    fprintf (f, "P6\n%d %d\n255\n", im.width(), im.height());
    fwrite (im.raw(), sizeof(Pixel<uint8_t>), im.size(), f);
    fclose (f);
  }
}

int main () {
  using namespace RT;

  Point const
    eye    {0.5f, 0.f,-0.7f},
    centre {0.3f, 0.f, 0.0f};
  Vector const
    up     {0.f,-1.f, 0.f };
  Camera cam (eye, centre, up, 0.4f);

  Scene s;

  Lambertian const
    mWhite ({0.8f,0.8f,0.8f}),
    mRed   ({0.8f,0.0f,0.0f}),
    mGreen ({0.0f,0.8f,0.0f}),
    mBlue  ({0.0f,0.0f,0.8f});

  Plane const
    sLeftWall  ({ 1, 0, 0}, -1.5f),
    sRightWall ({-1, 0, 0}, -1.5f),
    sFloor     ({ 0,-1, 0}, -1.0f),
    sCeiling   ({ 0, 1, 0}, -1.5f),
    sRearWall  ({ 0, 0,-1}, -3.0f),
    sFrontWall ({ 0, 0, 1}, -0.5f);
  Item const
    leftWall  (s, &sLeftWall,  &mRed  ),
    rightWall (s, &sRightWall, &mGreen),
    floor     (s, &sFloor,     &mWhite),
    ceiling   (s, &sCeiling,   &mWhite),
    rearWall  (s, &sRearWall,  &mWhite),
    frontWall (s, &sFrontWall, &mWhite);

  Dielectric const mGlass ({1,1,1}, 1.458f);
  Mirror const mCu ({0.7, 0.2, 0.02});

  Sphere const sSphere ({0.3f,0.5f,2.2f}, 0.5f);
  Item   const sphere (s, &sSphere, &mWhite);

  Sphere const sCuBall ({0.65f,0.75f,1.3f}, 0.25f);
  Item   const cuBall (s, &sCuBall, &mCu);

  Sphere const sMarble ({-0.4,0.65f,1.2f}, 0.35);
  Item   const marble (s, &sMarble, &mGlass);

  Triangle const sTri ({-0.5f,-0.7f,2.9f}, {0.f,0.f,2.9f}, {0.5f,-0.7f,2.9f});
  Item     const tri (s, &sTri, &mWhite);

  Colour lamp { 150.f, 140.f, 110.f };

/*Triangle const
      uplumLeftShroudS  ({-2.0f,-0.2f,3.0f}, {-1.7f,-0.6f,3.0f}, {-2.0f,-0.6f,2.7f})
    , uplumRightShroudS ({ 2.0f,-0.2f,3.0f}, { 2.0f,-0.6f,2.7f}, { 1.7f,-0.6f,3.0f})
    , uplumLeftTopS     ({-2.0f,-0.7f,3.0f}, {-2.0f,-0.6f,2.7f}, {-1.7f,-0.6f,3.0f})
    , uplumRightTopS    ({ 2.0f,-0.7f,3.0f}, { 1.7f,-0.6f,3.0f}, { 2.0f,-0.6f,2.7f})
    ;
  Item const
    uplumLeftShroud  (s, &uplumLeftShroudS,  &mWhite),
    uplumRightShroud (s, &uplumRightShroudS, &mWhite);
  Lum const
    uplumLeftTop  (s, &uplumLeftTopS,  lamp),
    uplumRightTop (s, &uplumRightTopS, lamp);*/

  Sphere const
  /*  lum1Surf    ({ 0.7f, -1.3f, 2.3f}, 0.1f)
    , lum2Surf    ({ 0.7f, -1.3f, 1.2f}, 0.1f)
    , lum3Surf    ({-0.7f, -1.3f, 2.3f}, 0.1f)
    , lum4Surf    ({-0.7f, -1.3f, 1.2f}, 0.1f)
    , lumTriSurf  ({ 0.0f, -0.5f, 2.0f}, 0.3f)
    ,*/slBall ({ 0.1f, 0.8f, 1.5f}, 0.1f)
    ;
  Lum const
  /*  lum1    (s, &lum1Surf,    lamp)
    , lum2    (s, &lum2Surf,    lamp)
    , lum3    (s, &lum3Surf,    lamp)
    , lum4    (s, &lum4Surf,    lamp)
    , lumTri  (s, &lumTriSurf,  lamp)
    ,*/lBall (s, &slBall, lamp)
    ;

#ifndef NDEBUG
  feenableexcept (FE_DIVBYZERO | FE_INVALID);
#endif

  constexpr bool const
    wide = false;

  constexpr int const
    scale = 1,
    w0 = wide? 640 : 480,
    w = w0 * scale,
    h = 360 * scale;
  Image<uint8_t> im (w, h);

  Options opts;
  opts.sqrtSamples = 1;
  opts.initRadius = 1.0f;
  opts.alpha = 0.7f;
  opts.nPasses = 4;
  opts.nPhotonsPerPass = 100'000;
/*opts.indirect = true;
  opts.rtSamplesIndirect = 1;
  opts.rtSamplesDirect = 4;
  opts.directBranch = 1;
  opts.nPhotons = 200'000;
  opts.nNears = 50;*/
  opts.nThreads = 4;
/*opts.nearLimit = 0.01f;*/
  opts.exposeK = 0.1f;

  render (s, cam, im, opts);

  fprintf (stderr, "Done\n");

  auto outTime = std::chrono::system_clock::to_time_t (
    std::chrono::system_clock::now ());
  char path[64];
  strftime (path, sizeof (path), "out-%FT%TZ.ppm", localtime (&outTime));
  writePPM (path, im);
}

