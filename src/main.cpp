
#include "linear.hpp"
#include "random.hpp"
#include "colour.hpp"
#include "image.hpp"
#include "scene.hpp"
#include "photon.hpp"
#include "photonmapping.hpp"
//#include "directrt.hpp"
#include "ppm.hpp"
#include "monitor.hpp"

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
          directBranch,*/
          nPhotons,
        //nNears,
          sqrtSamples,
          nPasses,
        //nPhotonsPerPass,
          branch,
          nThreads;
    float //nearLimit,
          initRadius,
          alpha;
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
          storage[j*sqrtn + i] = {x,y,0};
        }
      }
    }

    Vector const* begin () const { return storage; }
    Vector const* end   () const { return storage + size(); }
    size_t size () const { return sqrtn*sqrtn; }

    Vector operator[] (int i) { return storage[i]; }
  };

  void tracePath
    ( std::vector<PathNode>& nodes
    , float initRadius
    , RandBits& rng
    , Scene const& scene
    , int xPixel, int yPixel
    , Ray ray, int ttl, int branch
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
  //Interaction const inter = mat->interact (rng, incident, sourceK);

    if (mat->kDiffuse != black) {
      PathNode node
        ( hit.position, hit.normal
        , mat
        , xPixel, yPixel
        , sourceK
        , initRadius
        );
      nodes.push_back (node);
    }

    if (mat->kSpecular != black) { //inter.type == Interaction::Type::specular) {
      Colour const newSourceK
        = (1.f/branch)
        * mat->correctRefl (Interaction::Type::specular, sourceK)
        * sourceK;

      for (int b = 0; b != branch; b++) {
        Vector const bounce = tang.outOf (mat->bounceSpecular (rng, incident));
        Ray const newRay (hit.position, bounce);
        tracePath
          ( nodes, initRadius
          , rng
          , scene
          , xPixel, yPixel
          , newRay, ttl-1, branch
          , newSourceK
          );
      }
    }
  }

  std::vector<PathNode> tracePaths
    ( Options const& opts
    , Scene const& scene
    , Camera const& camera
    , int w, int h
    )
  {
    fprintf (stderr, "Tracing paths\n");

    float const
      xImToCam = 2.f / w,
      yImToCam = 2.f / h,
      aspect = (float) w / h;

    int const nSamples = opts.sqrtSamples * opts.sqrtSamples;

    RandBits rng = threadRNG ();
    SampleSet samplePosns (opts.sqrtSamples);

    std::vector<PathNode> nodes;
    nodes.reserve (nSamples * w * h);

    for (int yPixel = 0; yPixel != h; yPixel++) {
      int const progress
        = int (yPixel * 100.f / h);
      fprintf (stderr, "\r%3d%%", progress);

      float const yPixelCentre = float(yPixel) + 0.5f;

      for (int xPixel = 0; xPixel != w; xPixel++) {
        float const xPixelCentre = float(xPixel) + 0.5f;

        samplePosns.regen (rng);
        for (int i = 0; i != nSamples; i++) {
          Vector const samplePos = samplePosns[i];
          float const
            xSample = (xPixelCentre + samplePos.x)*xImToCam - 1.f,
            ySample = (yPixelCentre + samplePos.y)*yImToCam - 1.f;
          Ray const ray = camera.shoot (aspect * xSample, ySample);
          tracePath (
            nodes, opts.initRadius,
            rng,
            scene,
            xPixel, yPixel,
            ray, opts.branch, 10
          );
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
    , Monitor& mon
    , Options const& opts
    )
  {
    using clock = std::chrono::steady_clock;
    auto startTime = clock::now ();

    auto nodes = tracePaths (opts, scene, cam, im.width (), im.height ());

    int const
      nPhotonsPerPass = opts.nPhotons / opts.nPasses,
      nPhotonsPerThread = nPhotonsPerPass / opts.nThreads,
      nSamples = opts.sqrtSamples * opts.sqrtSamples;

    std::thread threads[opts.nThreads];
    std::mutex mutex;
    int nEmitted = 0;

    std::vector<Photon> photonsRaw (nPhotonsPerPass);

    Image<uint8_t> oldIm;

    for (int iPass = 0; iPass != opts.nPasses; iPass++) {
      fprintf (stderr, "Pass %i/%i\n", iPass+1, opts.nPasses);

      // cast
      fprintf (stderr, "Casting...");
      for (int iThread = 0; iThread != opts.nThreads; iThread++) {
        auto& thread = threads[iThread];
        thread = std::thread (
          [&, iThread] {
            RandBits rng = threadRNG (iPass);
            Photon* photons = photonsRaw.data () + iThread * nPhotonsPerThread;
            int n = castPhotons (rng, scene, photons, photons + nPhotonsPerThread);
            std::lock_guard<std::mutex> lock (mutex);
            nEmitted += n;
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
      int nodeProgress = 0;

      for (int iThread = 0; iThread != opts.nThreads; iThread++) {
        auto& thread = threads[iThread];
        thread = std::thread (
          [&, iThread] {
            for (;;) {
              int begin, end;
              {
                std::lock_guard<std::mutex> lock (mutex);
                if (nodeProgress == nNodes)
                  return;
                begin = nodeProgress;
                end = nodeProgress = std::min (nodeProgress + 100, nNodes);
              }

              for (int i = begin; i != end; i++)
                nodes[i].update (map, opts.alpha);
            }
          }
        );
      }

      for (auto& thread : threads)
        thread.join ();

      // gather
      Image<float> frame (im.width (), im.height ());
      float const sampleK = 1.f / float (nSamples * nEmitted);
      for (auto const& node : nodes) {
        auto L = sampleK * node.radiance ();
        //fprintf (stderr, "L: (%.3f %.3f %.3f)\n", L.r, L.g, L.b);
        frame.at (node.px, node.py) += {L.r, L.g, L.b};
      }

      // ramp
      for (int y = 0; y != im.height (); y++) {
        for (int x = 0; x != im.width (); x++)
          im.at (x, y) = toSRGB (frame.at (x, y));
      }

      if (iPass > 0) {
        int maxDiff = 0;
        for (int y = 0; y != im.height (); y++) {
          for (int x = 0; x != im.width (); x++) {
            auto pNew = im.at (x, y);
            auto pOld = oldIm.at (x, y);
            maxDiff = std::max ({
              maxDiff,
              std::abs ((int) pNew.r - pOld.r),
              std::abs ((int) pNew.g - pOld.g),
              std::abs ((int) pNew.b - pOld.b)
            });
          }
        }

        //fprintf (stderr, "\rmax diff: %i\n", maxDiff);

        if (maxDiff == 0)
          break;
      }

      // update
      mon.update (im);
      oldIm = im;

      // recycle
      photonsRaw = std::move (map.array);
      fprintf (stderr, "\r                        \r");
    }

    auto stopTime = clock::now ();
    auto duration = stopTime - startTime;
    auto seconds = std::chrono::duration<float> (duration);

    fprintf (stderr, "\rTook %.1f s\n", seconds.count ());
    fprintf (stderr, "Emitted %i\n", nEmitted);
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
    eye    {0.4f, 0.f,-0.6f},
    centre {0.3f, 0.f, 0.0f};
  Vector const
    up     {0.f,-1.f, 0.f };
  Camera cam (eye, centre, up, 0.4f);

  Scene s;

  Lambertian const
    mWhite  ({0.8f,0.8f,0.8f}),
    mRed    ({0.8f,0.0f,0.0f}),
    mGreen  ({0.0f,0.8f,0.0f}),
    mYellow ({0.6f,0.6f,0.0f}),
    mBlue   ({0.0f,0.0f,0.8f});

  Colour lamp { 16.f, 16.f, 16.f };

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

  Dielectric const mGlass ({0.2f,0.4f,1.f}, 1.458f);
  Mirror const mCu ({0.7, 0.2, 0.02});

  Sphere const sSphere ({0.2f,0.5f,1.8f}, 0.5f);
  Item   const sphere (s, &sSphere, &mWhite);

  Sphere const sCuBall ({0.65f,0.75f,1.3f}, 0.25f);
  Item   const cuBall (s, &sCuBall, &mCu);

  Sphere const sMarble ({-0.4f,0.65f,1.2f}, 0.35f);
  Item   const marble (s, &sMarble, &mGlass);

  Sphere const sGlow ({-0.7f, 0.95f, 0.9f}, 0.05f);
  Lum    const glow (s, &sGlow, lamp);

  Triangle const sTri ({-0.5f,-0.7f,2.9f}, {0.f,0.f,2.9f}, {0.5f,-0.7f,2.9f});
  //Item     const tri (s, &sTri, &mWhite);

/*Triangle const
    uplumLeftShroudS  ({-1.5f,-0.2f,3.0f}, {-1.2f,-0.6f,3.0f}, {-1.5f,-0.6f,2.7f}),
    uplumRightShroudS ({ 1.5f,-0.2f,3.0f}, { 1.5f,-0.6f,2.7f}, { 1.2f,-0.6f,3.0f}),
    uplumLeftTopS     ({-1.5f,-0.7f,3.0f}, {-1.5f,-0.6f,2.7f}, {-1.2f,-0.6f,3.0f}),
    uplumRightTopS    ({ 1.5f,-0.7f,3.0f}, { 1.2f,-0.6f,3.0f}, { 1.5f,-0.6f,2.7f});
  Item const
    uplumLeftShroud  (s, &uplumLeftShroudS,  &mWhite),
    uplumRightShroud (s, &uplumRightShroudS, &mWhite);
  Lum const
    uplumLeftTop  (s, &uplumLeftTopS,  lamp),
    uplumRightTop (s, &uplumRightTopS, lamp);*/

/*Sphere const
    lum1Surf ({ 0.7f, -1.3f, 2.3f}, 0.1f),
    lum2Surf ({ 0.7f, -1.3f, 1.2f}, 0.1f),
    lum3Surf ({-0.7f, -1.3f, 2.3f}, 0.1f),
    lum4Surf ({-0.7f, -1.3f, 1.2f}, 0.1f);

  Lum const
    lum1 (s, &lum1Surf, lamp),
    lum2 (s, &lum2Surf, lamp),
    lum3 (s, &lum3Surf, lamp),
    lum4 (s, &lum4Surf, lamp);*/

#ifndef NDEBUG
  feenableexcept (FE_DIVBYZERO | FE_INVALID | FE_UNDERFLOW | FE_OVERFLOW);
#endif

  constexpr bool const
    wide = false;

  constexpr int const
    scale = 2,
    w0 = wide? 640 : 480,
    w = w0 * scale,
    h = 360 * scale;
  Image<uint8_t> im (w, h);

  Options opts;
  opts.sqrtSamples = 2;
  opts.branch = 80;
  opts.initRadius = 0.5f;
  opts.alpha = 0.95f;
  opts.nPhotons = 640'000;
  opts.nPasses  = 64;
/*opts.indirect = true;
  opts.rtSamplesIndirect = 1;
  opts.rtSamplesDirect = 4;
  opts.directBranch = 1;
  opts.nPhotons = 200'000;
  opts.nNears = 50;*/
  opts.nThreads = 4;
/*opts.nearLimit = 0.01f;*/

  Monitor mon (im.width (), im.height ());
  render (s, cam, im, mon, opts);

  fprintf (stderr, "Done\n");

  auto outTime = std::chrono::system_clock::to_time_t (
    std::chrono::system_clock::now ());
  char path[64];
  strftime (path, sizeof (path), "out-%FT%TZ.ppm", localtime (&outTime));
  writePPM (path, im);
}

