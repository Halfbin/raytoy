#include <cstdio>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>
#include <thread>
#include <mutex>
#include <cassert>
#include <chrono>

namespace RT {
  struct Options {
    bool indirect;
    int  sqrtNSamples,
         nPhotons,
         nNears;
    float nearLimit;
  };

  using u8 = uint8_t;

  constexpr float const infF = std::numeric_limits<float>::infinity ();

//using RandBits = std::ranlux24_base;
//using RandBits = std::knuth_b;
  using RandBits = std::shuffle_order_engine<std::minstd_rand, 256>;
  std::uniform_real_distribution<float> canon (0.f, 1.f);

  bool coin (RandBits& rng) {
    auto i = rng();
    return (i-RandBits::min()) > (RandBits::max()-RandBits::min())/2;
  }

  constexpr float const pi = M_PI;

  struct Pixel {
    uint8_t r, g, b;
    Pixel () = default;
    Pixel (Pixel const&) = default;
    constexpr Pixel (u8 r, u8 g, u8 b) : r(r), g(g), b(b) { }
    static constexpr Pixel black () { return Pixel{0,0,0}; }
  };

  struct Colour {
    float r, g, b;
    Colour () = default;
    constexpr Colour (float r, float g, float b) : r(r), g(g), b(b) { }
    Colour (Colour const&) = default;
    static constexpr Colour black () { return Colour{0,0,0}; }
    explicit operator bool () const { return *this != black(); }
    bool operator== (Colour const& other) const
      { return r==other.r && g==other.g && b==other.b; }
    bool operator!= (Colour const& other) const { return !(*this == other); }
  };

  Colour operator* (Colour a, Colour b) { return {a.r*b.r, a.g*b.g, a.b*b.b}; }
  Colour operator* (Colour c, float k) { return {c.r*k, c.g*k, c.b*k}; }
  Colour operator* (float k, Colour c) { return c*k; }
  Colour operator/ (Colour c, float k) { return {c.r/k, c.g/k, c.b/k}; }
  Colour operator+ (Colour a, Colour b) { return {a.r+b.r, a.g+b.g, a.b+b.b}; }

  Colour& operator += (Colour& a, Colour b) { return (a = a+b); }
  Colour& operator *= (Colour& a, Colour b) { return (a = a*b); }

  float max (Colour c) { return std::max (c.r, std::max (c.g, c.b)); }

  Colour clamp (Colour c) {
    return Colour
      { std::max (0.f, std::min (1.f, c.r))
      , std::max (0.f, std::min (1.f, c.g))
      , std::max (0.f, std::min (1.f, c.b))
      };
  }

  float toSRGB1 (float raw) {
    float const
      c = std::max (0.f, std::min (1.f, raw)),
      lo = 12.92f*c,
      hi = 1.055f*std::pow (c, 1.f/2.4f) - 0.055f;
    return (c <= 0.0031308)? lo : hi;
  }

  u8 toU8 (float x) {
    return u8 (255.f * x);
  }

  Pixel toSRGB (Colour c) {
    return
      { toU8 (toSRGB1 (c.r))
      , toU8 (toSRGB1 (c.g))
      , toU8 (toSRGB1 (c.b))
      };
  }

  class Image {
    const int wide, high;
    Pixel* const pixels;
  public:
    Image (int w, int h) :
      wide (w), high (h),
      pixels (new Pixel [size()])
    {
      for (size_t i = 0; i != size(); i++)
        pixels[i] = Pixel::black ();
    }
    ~Image () { delete[] pixels; }
    int width  () const { return wide; }
    int height () const { return high; }
    size_t size () const { return wide * high; }
    Pixel const* raw () const { return pixels; }
    Pixel& at (int x, int y) { return pixels [y * wide + x]; }
  };

  struct Vector {
    union {
      float components[4];
      struct { float x, y, z, w; };
    };
    constexpr Vector () : x(0), y(0), z(0), w(0) { }
    constexpr Vector (float x, float y, float z, float=0) : x(x), y(y), z(z), w(0) { }
    Vector (Vector const&) = default;
    bool operator== (Vector const& other) const
      { return x==other.x && y==other.y && z==other.z && w==other.w; }
    bool operator!= (Vector const& other) const { return !(*this == other); }
  };

  constexpr Vector const zero = Vector ();

  struct Point {
    union {
      float components[4];
      struct { float x, y, z, w; };
    };
    Point () : x(0), y(0), z(0), w(1) { }
    Point (float x, float y, float z, float=1) : x(x), y(y), z(z), w(1) { }
    Point (Point const&) = default;
  };

  Vector operator+ (Vector a, Vector b) { return { a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w }; }
  Point  operator+ (Point p,  Vector d) { return { p.x+d.x, p.y+d.y, p.z+d.z, p.w+d.w }; }
  Point  operator- (Point p,  Vector d) { return { p.x-d.x, p.y-d.y, p.z-d.z, p.w-d.w }; }
  Vector operator- (Point p,  Point  q) { return { p.x-q.x, p.y-q.y, p.z-q.z, p.w-q.w }; }
  Vector operator- (Vector a, Vector b) { return { a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w }; }
  Vector operator* (Vector a, Vector b) { return { a.x*b.x, a.y*b.y, a.z*b.z, a.w*b.w }; }

  Vector operator* (float k,  Vector v) { return { k*v.x, k*v.y, k*v.z, k*v.w }; }
  Vector operator* (Vector v, float k) { return k*v; }

  Vector operator- (Vector v) { return { -v.x, -v.y, -v.z, -v.w }; }

  float dot (Vector a, Vector b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

  float  norm2 (Vector v) { return dot(v,v); }
  float  norm  (Vector v) { return sqrt(norm2(v)); }
  Vector unit  (Vector v) { return v * (1.f / norm(v)); }

  Vector cross (Vector a, Vector b) {
    return
      { a.y*b.z - a.z*b.y
      , a.z*b.x - a.x*b.z
      , a.x*b.y - a.y*b.x
      };
  }

  float max (Vector v) { return std::max (v.x, std::max (v.y, v.z)); }

  Vector hemispherical (RandBits& rng) {
    float const
      xi1 = canon (rng),
      xi2 = canon (rng),
      cosTheta = std::sqrt(xi1),
      phi      = 2.f*pi*xi2,
      sinTheta = std::sqrt(1.f - cosTheta*cosTheta),
      x = std::cos(phi)*sinTheta,
      y = std::sin(phi)*sinTheta,
      z = cosTheta;
    return {x,y,z};
  }

  Vector spherical (RandBits& rng) {
    auto v = hemispherical (rng);
    float const sign = coin (rng)? 1.f : -1.f;
    v.z *= sign;
    return v;
  }

  struct Basis {
    Vector const i, j, k;

    Basis (Vector i, Vector j, Vector k) : i(i), j(j), k(k) { }
    Basis (Basis const&) = default;

    static Basis fromK (Vector n) {
      float const sign = std::copysignf (1.f, n.z),
                  a = -1.f / (sign + n.z),
                  b = n.x*n.y*a;
      return Basis
        ( { 1.f + sign*n.x*n.x*a, sign*b,           -sign*n.x }
        , { b,                    sign + n.y*n.y*a, -     n.y }
        , n
        );
    }

    Vector into (Vector v) const {
      return { dot (i, v), dot (j, v), dot (k, v) };
    }

    Vector outOf (Vector v) const {
      return v.x * i + v.y * j + v.z * k;
    }
  };

  struct Ray {
    Point const origin;
    Vector const disp;
    Ray () : origin{0,0,0}, disp{0,0,0} { }
    Ray (Point o, Vector d) : origin(o), disp(unit(d)) { }
    Point at (float t) const { return origin + t*disp; }
    explicit operator bool () const { return disp == Vector{0,0,0}; }
  };

  struct RayHit {
    Point  position;
    float  t;
    Vector normal;

    RayHit () : t(-1.f) { }
    RayHit (RayHit const&) = default;
    RayHit& operator= (RayHit const&) = default;
    RayHit (Point p, float t, Vector n) :
      position(p), t(t), normal(n) { }
    explicit operator bool () const { return t >= 0.f; }
  };

  class Surface {
  public:
    virtual RayHit intersect (Ray) const = 0;
  };

  class BoundedSurface : public Surface {
  public:
    virtual Ray randomEmission (RandBits&) const = 0;
    virtual Point randomPoint (RandBits&, Point towards) const = 0;
    virtual float area () const = 0;
    virtual float boundingRadius () const = 0;
    virtual Point boundingCentre () const = 0;
  };

  struct Sphere : public BoundedSurface {
    Point const centre;
    float const radius;

    Sphere (Point o, float r) : centre(o), radius(r) { }

    RayHit intersect (Ray r) const {
      Vector const s = r.origin - centre; // "separation"
      float  const ds = dot (r.disp, s), // d-dot-s
                  dd = norm2 (r.disp),  // d-dot-d
                  discrim = ds*ds - dd*(norm2(s) - radius*radius);
      if (discrim < 0.f) return { };
      float const root = std::sqrt(discrim),
                  tdd = -(ds + root);
      if (tdd < 0.0001f) return { };
      float const t = tdd / dd;
      Point const p = r.at(t);
      Vector const n = unit(p-centre);
      return RayHit (p, t, n);
    }

    Ray randomEmission (RandBits& rng) const {
      Vector const normal = spherical (rng);
      Point const origin = centre + radius * normal;
      Vector const dir = hemispherical (rng);
      Basis const tang = Basis::fromK (normal);
      return Ray (origin, tang.outOf (dir));
    }

    Point randomPoint (RandBits& rng, Point towards) const {
      Basis const facing = Basis::fromK (unit (towards - centre));
      return centre + radius * facing.outOf (hemispherical (rng));
    }

    float area () const {
      return 4.f * pi * radius * radius;
    }

    float boundingRadius () const {
      return radius;
    }

    Point boundingCentre () const {
      return centre;
    }
  };

  RayHit rayPlane (Ray r, Vector normal, Point p0) {
    float const denom = dot (r.disp, normal);
    if (std::abs(denom) < 0.000001f)
      return { };
    float const t = dot (p0 - r.origin, normal) / denom;
    Vector const n = (denom<0.f)? normal : -normal;
    return { r.at(t), t, n };
  }

  struct Plane : public Surface {
    Vector const normal;
    Point  const onPlane;
    Plane (Vector n, float d) :
      normal(unit(n)), onPlane(Point{0,0,0} + d*n) { }
    RayHit intersect (Ray r) const {
      return rayPlane (r, normal, onPlane);
    }
  };

  struct Triangle : public BoundedSurface {
    Point const p0;
    Vector const abc, def;
    Basis const basis;

    Triangle (Point p, Point a, Point b) :
      p0(p), abc(p0-a), def(p0-b),
      basis (Basis::fromK (unit (cross (abc, def))))
    {
      /*fprintf (stderr, "made tri with normal (%f %f %f)\n",
        basis.k.x, basis.k.y, basis.k.z);*/
    }

    RayHit intersect (Ray r) const {
      Vector const
        jkl = p0 - r.origin,
        defxdisp = cross (def, r.disp),
        abcxjkl  = cross (abc, jkl);
      float const
        denom = dot (abc,    defxdisp),
        beta  = dot (jkl,    defxdisp) /  denom,
        gamma = dot (r.disp, abcxjkl ) /  denom,
        t     = dot (def,    abcxjkl ) / -denom;
      if (beta <= 0.f || gamma <= 0.f || (beta+gamma) >= 1.f)
        return { };
      return { r.at(t), t, basis.k };
    }

    Ray randomEmission (RandBits& rng) const {
      Point const origin = randomPoint (rng, {0,0,0});
      Vector const dir = hemispherical (rng);
      Vector const disp = basis.outOf (dir);
      /*fprintf (stderr, "emitted photon along (%f %f %f)\n",
        disp.x, disp.y, disp.z);*/
      return Ray (origin, disp);
    }

    Point randomPoint (RandBits& rng, Point) const {
      float
        xi1 = canon (rng),
        xi2 = canon (rng);

      if (xi1 + xi2 > 1.f) {
        xi1 = (1.f - xi1);
        xi2 = (1.f - xi2);
      }

      return p0 - xi1 * abc - xi2 * def;
    }

    float area () const {
      return dot (abc, def);
    }

    float boundingRadius () const {
      norm ((abc + def) * (1.f/3.f));
    }

    Point boundingCentre () const {
      return p0 - (abc + def) * (1.f/3.f);
    }
  };

  struct Interaction {
    enum class Type { diffuse, specular, absorbed } type;
    Colour correctRefl;
  };

  class Material {
  public:
    Colour const
      kDiffuse,
      kSpecular,
      emission;

    Material (Colour kd, Colour ks, Colour e = Colour::black ()) :
      kDiffuse(kd), kSpecular(ks), emission(e) { }

    Interaction interact (RandBits& rng, Vector, Colour ci) const {
      float const
        pDiffuse  = max (kDiffuse *ci) / max (ci),
        pSpecular = max (kSpecular*ci) / max (ci),
        pAbsorb   = 1.f - pDiffuse - pSpecular,
        xi = canon (rng);

      if (xi < pDiffuse)
        return { Interaction::Type::diffuse, kDiffuse / pDiffuse };
      else if (xi - pDiffuse < pSpecular)
        return { Interaction::Type::specular, kSpecular / pSpecular };
      else
        return { Interaction::Type::absorbed, {0,0,0} };
    }

    virtual Vector bounce (RandBits& rng, Interaction::Type t, Vector incident) const {
      switch (t) {
        case Interaction::Type::diffuse: return bounceDiffuse (rng, incident);
        case Interaction::Type::specular: return bounceSpecular (rng, incident);
        case Interaction::Type::absorbed: return {0,0,0};
      }

      return {0,0,0};
    }

    virtual Vector bounceDiffuse (RandBits& rng, Vector) const {
      return hemispherical (rng);
    }

    virtual Vector bounceSpecular (RandBits&, Vector incident) const {
      return incident * Vector{1,1,-1};
    }
  };

  class Lambertian : public Material {
  public:
    Lambertian (Colour k) : Material (k, {0,0,0}) { }
  };

  float fresnel (float k, float cosine) {
    return k + (1.f - k)*pow(1 - cosine, 5.f);
  }

  class Mirror : public Material {
  public:
    Mirror (Colour k) : Material ({0,0,0}, k) { }
  /*Interaction bounce (RandBits&, Vector incident) const {
      float const
        k = 0.75,
        schlick = fresnel (k, -incident.z);
      return { incident*Vector{1,1,-1}, schlick*Colour{1,0.3,0.1} };
    }*/
  };

  class Glass : public Material {
  public:
    Glass () : Material ({0,0,0}, {1,1,1}) { }

    Vector bounceSpecular (RandBits& rng, Vector incident) const {
      float const
        c = std::abs (incident.z),
        k = 0.25f,
        pReflect = fresnel (k, c),
        xi1 = canon (rng);

      Vector const reflection = incident*Vector{1,1,-1};
      if (xi1 < pReflect)
        return reflection;

      float const
        eta = 1.7f,
        r = (incident.z>0.f)? eta : (1.f/eta),
        discrim = 1.f - r*r*(1.f - c*c);
      /*if (discrim < 0.f)
        return reflection;*/

      return { r*incident.x, r*incident.y, r*c - std::sqrt(discrim) };
    }
  };

  struct Item {
    Surface const* surface;
    Material const* material;

    Item (Surface const* surface, Material const* material) :
      surface (surface), material (material)
    { }
  };

  struct Lum : Item {
    BoundedSurface const* surface;
    Colour const power;
    Material const material;

    Lum (Lum const&) = delete;
    Lum& operator = (Lum const&) = delete;

    Lum (BoundedSurface const* surface, Colour power) :
      Item (surface, &material),
      surface (surface),
      power (power),
      material ({0,0,0}, {0,0,0}, overallRadiance ())
    { }

    Ray emit (RandBits& rng) const {
      return surface->randomEmission (rng);
    }

    Point randomPoint (RandBits& rng, Point towards) const {
      return surface->randomPoint (rng, towards);
    }

    Colour overallRadiance () const {
      return power / surface->area ();
    }

    float coverageK (Point about) const {
      float const
        r = surface->boundingRadius (),
        d = norm (about - surface->boundingCentre ()),
        cosTheta = d / std::sqrt (r*r + d*d);
      return 1.f - cosTheta;
    }
  };

  enum class Axis : char { none, X, Y, Z };

  struct RGBE {
    unsigned char r, g, b, e;

    RGBE () = default;

    explicit RGBE (Colour c) {
      int re, ge, be;
      float const
        rm = std::frexp (c.r, &re),
        gm = std::frexp (c.g, &ge),
        bm = std::frexp (c.b, &be);
      int maxExp = -128;
      if (rm != 0.f) maxExp = re;
      if (gm != 0.f) maxExp = std::max (maxExp, ge);
      if (bm != 0.f) maxExp = std::max (maxExp, be);
      if (maxExp == -128) maxExp = 0;
      r = (unsigned char) ((int) (rm * 255.f) >> (maxExp - re));
      g = (unsigned char) ((int) (gm * 255.f) >> (maxExp - ge));
      b = (unsigned char) ((int) (bm * 255.f) >> (maxExp - be));
      e = (unsigned char) (128 + maxExp);
      /*fprintf (stderr, "(%f %f %f) -> (%3i %3i %3i | %3i)\n",
        c.r, c.g, c.b, (int) r, (int) b, (int) g, (int) e);*/
    }

    explicit operator Colour () const {
      int const exp = (int) e - 128;
      Colour const c
        { std::scalbn ((int) r / 255.f, exp)
        , std::scalbn ((int) g / 255.f, exp)
        , std::scalbn ((int) b / 255.f, exp)
        };
      /*fprintf (stderr, "(%f %f %f) <- (%3i %3i %3i | %3i)\n",
        c.r, c.g, c.b, (int) r, (int) b, (int) g, (int) e);*/
      return c;
    }
  };

  struct Octo {
    signed char s, t;

    Octo () = default;

    static float signumNZ (float f) {
      return (f >= 0.f)? 1.f : -1.f;
    }

    static char toByte (float f) {
      return ((int) (f * 127.f));
    }

    explicit Octo (Vector v) {
      v = unit (v);
      using std::abs;
      float const
        k = 1.f / (abs (v.x) + abs (v.y) + abs (v.z)),
        px = v.x * k,
        py = v.y * k;
      if (v.z > 0.f) {
        s = toByte (px);
        t = toByte (py);
      }
      else {
        float const
          sf = (1.f - abs (py)) * signumNZ (px),
          tf = (1.f - abs (px)) * signumNZ (py);
        s = toByte (sf);
        t = toByte (tf);
      }
      /*fprintf (stderr, "(%f %f %f) -> (%i %i)\n",
        v.x, v.y, v.z, (int) s, (int) t);*/
    }

    explicit operator Vector () const {
      float const
        sf = ((int) s) / 127.f,
        tf = ((int) t) / 127.f;
      using std::abs;
      float const z = 1.f - abs (sf) - abs (tf);
      Vector v;
      if (z < 0.f) {
        float const
          x = (1.f - abs (tf)) * signumNZ (sf),
          y = (1.f - abs (sf)) * signumNZ (tf);
        v = unit ({ x, y, z });
      }
      else {
        v = unit ({ sf, tf, z });
      }
      assert (norm2 (v) - 1.f < 0.001f);
      /*fprintf (stderr, "(%+f %+f %+f) <- (%i %i)\n",
        v.x, v.y, v.z, (int) s, (int) t);*/
      return v;
    }
  };

  struct Photon {
    float x, y, z;
    RGBE powerPacked;
    Octo incomingPacked;
    Axis split;
    char pad[1];

    Photon () = default;
    Photon (Point p, Vector i, Colour c)
      : x (p.x), y (p.y), z (p.z)
      , powerPacked (c)
      , incomingPacked (i)
      , split (Axis::none)
      { }

    Point position () const {
      return { x, y, z };
    }

    Colour power () const {
      return static_cast<Colour> (powerPacked);
    }

    Vector incoming () const {
      return static_cast<Vector> (incomingPacked);
    }
  };

  //static_assert (sizeof (Photon) == 20, "Photon miscompiled");

  struct PhotonMap {
    static constexpr int const maxLeaf = 16;

    std::vector<Photon> array;

    PhotonMap () = default;
    PhotonMap (PhotonMap const&) = delete;
    PhotonMap& operator = (PhotonMap const&) = delete;
    PhotonMap (PhotonMap&&) = default;
    PhotonMap& operator = (PhotonMap&&) = default;
  };

  class NearSet {
  public:
    struct Element {
      float qd;
      Photon photon;
      bool operator < (Element const& other) const {
        return qd < other.qd;
      }
    };

    /*static constexpr auto const cmp
      = [] (Element const& a, Element const& b) { return a.qd < b.qd; };*/

  private:
    size_t const cap;
    float const limitQd;
    size_t n;
    Element* elems;

  public:
    NearSet (NearSet const&) = delete;

    NearSet (size_t cap, float limitQd) :
      cap (cap), limitQd (limitQd), n (0), elems (new Element[cap])
    { }

    ~NearSet () {
      delete [] elems;
    }

    bool full () const {
      return n == cap;
    }

    void clear () {
      n = 0;
    }

    float maxQd () const {
      if (!full ())
        return limitQd;
      else
        return elems[0].qd;
    }

    void insert (float qd, Photon photon) {
      if (!full ()) {
        elems[n++] = { qd, photon };
        if (full ())
          std::make_heap (elems, elems+cap);
      }
      else if (qd < maxQd ()) {
        std::pop_heap (elems, elems+cap);
        elems[cap-1] = { qd, photon };
        std::push_heap (elems, elems+cap);
      }
    }

    void recycle (Point newPos) {
      if (!full ()) {
        clear ();
        return;
      }

      float const newMaxQd = norm2 (elems[0].photon.position () - newPos);
      if (newMaxQd > maxQd ()) {
        clear ();
      }
      else {
        elems[0].qd = newMaxQd;
        for (Element* elem = elems+1; elem != end (); elem++) {
          elem->qd = norm2 (elem->photon.position () - newPos);
        }
        std::make_heap (elems, elems+cap);
      }
    }

    Element const* begin () const {
      return elems;
    }

    Element const* end () const {
      return elems+n;
    }
  };

  void findNearestPhotons
    ( Photon const* begin, Photon const* end
    , Point pos
    , NearSet& nears
    )
  {
    if (end - begin > PhotonMap::maxLeaf) {
      Photon const* node = begin + (end - begin) / 2;
      float const qd = norm2 (pos - node->position ());

      float planeDist = infF;
      switch (node->split) {
        case Axis::X: planeDist = pos.x - node->x; break;
        case Axis::Y: planeDist = pos.y - node->y; break;
        case Axis::Z: planeDist = pos.z - node->z; break;
        default: assert (false); //"tried to treat photon map leaf as node!"
      }

      float const planeQd = planeDist * planeDist;

      if (planeDist > 0.f) {
        findNearestPhotons (node+1, end, pos, nears);
        if (planeQd < nears.maxQd ())
          findNearestPhotons (begin, node, pos, nears);
      }
      else {
        findNearestPhotons (begin, node, pos, nears);
        if (planeQd < nears.maxQd ())
          findNearestPhotons (node+1, end, pos, nears);
      }

      nears.insert (qd, *node);
    }
    else {
      for (Photon const* photon = begin; photon != end; photon++) {
        float const qd = norm2 (pos - photon->position ());
        nears.insert (qd, *photon);
      }
    }
  }

  void nearestPhotons (NearSet& nears, PhotonMap const& map, Point pos) {
    nears.recycle (pos);
    auto* ptr = map.array.data ();
    findNearestPhotons (ptr, ptr + map.array.size (), pos, nears);
  }

  void buildPhotonMap (Photon* begin, Photon* end, Vector mins, Vector maxs) {
    if (end - begin <= PhotonMap::maxLeaf)
      return;

    Vector extents = maxs - mins;
    float longest = max (extents);
    Axis split;
    if      (longest == extents.x) split = Axis::X;
    else if (longest == extents.y) split = Axis::Y;
    else                           split = Axis::Z;

    switch (split) {
      case Axis::X:
        std::sort (begin, end, [] (auto const& a, auto const& b)
          { return a.x < b.x; });
        break;
      case Axis::Y:
        std::sort (begin, end, [] (auto const& a, auto const& b)
          { return a.y < b.y; });
        break;
      case Axis::Z:
        std::sort (begin, end, [] (auto const& a, auto const& b)
          { return a.z < b.z; });
        break;
      default:
        assert (false);
    }

    Photon* median = begin + (end - begin) / 2;
    median->split = split;
    Vector newMins = mins,
           newMaxs = maxs;
    switch (split) {
      case Axis::X: newMins.x = newMaxs.x = median->x; break;
      case Axis::Y: newMins.y = newMaxs.y = median->y; break;
      case Axis::Z: newMins.z = newMaxs.z = median->z; break;
      default: assert (false);
    }

    buildPhotonMap (begin, median, mins, newMaxs);
    buildPhotonMap (median+1, end, newMins, maxs);
  }

  PhotonMap mapPhotons (std::vector<Photon> in) {
    PhotonMap map { std::move (in) };

    Vector mins {infF,infF,infF}, maxs{-infF,-infF,-infF};
    for (auto& photon : map.array) {
      mins.x = std::min (mins.x, photon.x);
      mins.y = std::min (mins.y, photon.y);
      mins.z = std::min (mins.z, photon.z);
      maxs.x = std::max (maxs.x, photon.x);
      maxs.y = std::max (maxs.y, photon.y);
      maxs.z = std::max (maxs.z, photon.z);
    }

    fprintf (stderr, "mins: (%f %f %f) maxs: (%f %f %f)\n",
      mins.x, mins.y, mins.z,
      maxs.x, maxs.y, maxs.z);
    auto* p = map.array.data ();
    buildPhotonMap (p, p+map.array.size(), mins, maxs);

    return map;
  }

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

  constexpr Colour const background = Colour::black();

  Lambertian const
    whiteMaterial ({0.6f,0.6f,0.6f}),
    redMaterial   ({0.6f,0.f,0.f}),
    greenMaterial ({0.f,0.6f,0.f}),
    blueMaterial  ({0.f,0.f,0.6f});

  Plane const
    leftWallSurf  ({ 1, 0, 0}, -2.0f),
    rightWallSurf ({-1, 0, 0}, -2.0f),
    floorSurf     ({ 0,-1, 0}, -1.0f),
    ceilingSurf   ({ 0, 1, 0}, -1.5f),
    rearWallSurf  ({ 0, 0,-1}, -3.0f),
    frontWallSurf ({ 0, 0, 1},  0.0f);
  Item const
    leftWall  {&leftWallSurf,  &redMaterial  },
    rightWall {&rightWallSurf, &greenMaterial},
    floor     {&floorSurf,     &whiteMaterial},
    ceiling   {&ceilingSurf,   &whiteMaterial},
    rearWall  {&rearWallSurf,  &whiteMaterial},
    frontWall {&frontWallSurf, &whiteMaterial};

  Sphere const sph_s ({0.3f,0.5f,2.2f}, 0.5f);
  Item   const sph {&sph_s, &whiteMaterial};

  Sphere const mirrorBallSurf ({0.65f,0.75f,1.3f}, 0.25f);
  Mirror const copper ({0.7, 0.3, 0.1});
  Item   const mirrorBall {&mirrorBallSurf, &copper};

  Sphere const glassBallSurf ({-0.4f,0.6f,1.5f}, 0.4f);
  Glass  const glass;
  Item   const glassBall {&glassBallSurf, &glass};

  Triangle const triSurf ({-0.5f,-0.7f,2.9f}, {0.f,0.f,2.9f}, {0.5f,-0.7f,2.9f});
  Item     const tri {&triSurf, &whiteMaterial};

  Colour lamp { 5.0f, 4.8f, 4.6f };

  Triangle const
    uplumLeftShroudS  ({-2.0f,-0.2f,3.0f}, {-1.7f,-0.6f,3.0f}, {-2.0f,-0.6f,2.7f}),
    uplumRightShroudS ({ 2.0f,-0.2f,3.0f}, { 2.0f,-0.6f,2.7f}, { 1.7f,-0.6f,3.0f}),
    uplumLeftTopS     ({-2.0f,-0.7f,3.0f}, {-2.0f,-0.6f,2.7f}, {-1.7f,-0.6f,3.0f}),
    uplumRightTopS    ({ 2.0f,-0.7f,3.0f}, { 1.7f,-0.6f,3.0f}, { 2.0f,-0.6f,2.7f});
  Item const
    uplumLeftShroud  {&uplumLeftShroudS,  &whiteMaterial},
    uplumRightShroud {&uplumRightShroudS, &whiteMaterial};
  Lum const
    uplumLeftTop  {&uplumLeftTopS,  lamp},
    uplumRightTop {&uplumRightTopS, lamp};

  Sphere const
    lum1Surf ({ 0.7f,-1.3f,2.3f}, 0.1f),
    lum2Surf ({ 0.7f,-1.3f,1.2f}, 0.1f),
    lum3Surf ({-0.7f,-1.3f,2.3f}, 0.1f),
    lum4Surf ({-0.7f,-1.3f,1.2f}, 0.1f),
    lumTriSurf ({0.0f, -0.5f, 2.0f}, 0.3f);
  Lum const
    lum1 {&lum1Surf, lamp},
    lum2 {&lum2Surf, lamp},
    lum3 {&lum3Surf, lamp},
    lum4 {&lum4Surf, lamp},
    lumTri {&lumTriSurf, lamp};

  Item const* items[]
    { &leftWall, &rightWall
    , &floor,    &ceiling
    , &rearWall, &frontWall
    , /*&lum1,*/ &lum2, /*&lum3,*/ &lum4
    , &sph, &mirrorBall, &glassBall
    , &tri
  //, &lumTri
    , &uplumLeftShroud
    , &uplumRightShroud
    , &uplumLeftTop
    , &uplumRightTop
    };

  constexpr int const nLums = 4;
  Lum const* lums[nLums] { &lum1, &lum2, &lum3, &lum4 };

  Lum const& randomLum (RandBits& rng) {
    float const xi = canon (rng);
    return *lums[std::min (nLums-1, int(xi*nLums))];
  }

  RayHit traceRay (Ray ray, Item const*& item, Item const* source) {
    RayHit hit = { };
    item = nullptr;
    for (Item const* candidate : items) {
      auto const tryHit = candidate->surface->intersect (ray);
      bool const hitSource = tryHit.t < 0.00001f && candidate == source;
      if (!tryHit || hitSource) {
        continue;
      }
      if (!hit || tryHit.t < hit.t) {
        hit = tryHit;
        item = candidate;
      }
    }

    return hit;
  }

  Colour indirectTerm (RayHit hit, NearSet& nears, PhotonMap const& map) {
    constexpr float const
      filterK = 10.f,
      normalizeFilter = 1.f - 2.f / (3.f * filterK);

    nearestPhotons (nears, map, hit.position);
    float const radius = std::sqrt (nears.maxQd ());
    //fprintf (stderr, "radius %f\n", radius);

    Colour power = Colour::black ();
    for (auto const& elem : nears) {
      float const
        cosine = std::max (0.f, -dot (elem.photon.incoming (), hit.normal)),
        filter = 1.f - std::sqrt (elem.qd) / filterK * radius;
      power += cosine * filter * elem.photon.power ();
    }

    float const area = pi * nears.maxQd ();
    return power / (area * normalizeFilter);
  }

  Colour directTerm (RandBits& rng, Basis const& tang, RayHit hit, Item const* self) {
    Lum const& lum = randomLum (rng);
    Point const lumPoint = lum.randomPoint (rng, hit.position);
    Vector const lumDir = unit (lumPoint - hit.position);

    float const geomK = dot (hit.normal, lumDir);
    if (geomK <= 0.f)
      return Colour::black ();

    Ray shadowRay (hit.position, lumDir);
    Item const* occluder = nullptr;
    RayHit const occluderHit = traceRay (shadowRay, occluder, self);
    if (!occluderHit || occluder != &lum)
      return Colour::black ();

    float const coverageK = lum.coverageK (hit.position);
    return coverageK * geomK * lum.overallRadiance ();
  }

  Colour trace
    ( Options const& opts
    , RandBits& rng
    , PhotonMap const& map
    , NearSet& nears
    , Ray ray
    , Item const* prev
    , int ttl
    , Colour sourceK = {1,1,1}
    )
  {
    if (ttl == 0)
      return Colour::black ();

    Item const* item = nullptr;
    RayHit hit = traceRay (ray, item, prev);
    if (!hit)
      return Colour::black ();

    auto const* mat = item->material;

    Basis const tang = Basis::fromK (hit.normal);
    Vector const incident = tang.into (ray.disp);

    Interaction const inter = mat->interact (rng, incident, sourceK);
    if (inter.type == Interaction::Type::specular) {
      Vector const bounce = tang.outOf (mat->bounceSpecular (rng, incident));
      Ray const newRay (hit.position, bounce);
      return trace (
        opts, rng, map, nears, newRay, item, ttl-1, sourceK * inter.correctRefl);
    }
    else {
      Colour const
        Ll = directTerm (rng, tang, hit, item),
        Li = opts.indirect
          ? indirectTerm (hit, nears, map)
          : Colour::black (),
        Ld = mat->kDiffuse * (Ll + Li),
        Le = mat->emission,
        L  = Ld + Le;
      return sourceK * L;
    }
  }

  Photon* tracePhoton
    ( RandBits& rng
    , Ray ray
    , Item const* prev
    , Colour power
    , int ttl
    , Photon* photons, Photon* end
    , bool indirect = false
    )
  {
    if (ttl == 0 || photons == end)
      return end;

    Item const* item = nullptr;
    RayHit const hit = traceRay (ray, item, prev);
    if (!hit)
      return photons;

    auto const* mat = item->material;

    /*fprintf (stderr, "hit at (%f %f %f)\n",
      hit.position.x, hit.position.y, hit.position.z);*/

    Basis const tang = Basis::fromK (hit.normal);
    Vector const incident = tang.into (ray.disp);

    Interaction const inter = mat->interact (rng, incident, power);
    bool const store
      =   indirect
      &&  inter.type != Interaction::Type::specular
      &&  mat->kDiffuse != Colour::black ();
    if (store) {
      Photon const photon (hit.position, ray.disp, power);
      *photons++ = photon;
    }

    if (inter.type == Interaction::Type::absorbed)
      return photons;

    Vector const bounce = tang.outOf (mat->bounce (rng, inter.type, incident));
    Ray const newRay (hit.position, bounce);
    return tracePhoton (
      rng, newRay, item, power * inter.correctRefl, ttl-1, photons, end, true);
  }

  void castPhotons (RandBits& rng, Photon* photons, Photon* end) {
    float const scale = 1.f / (end - photons);
    for (Photon* ptr = photons; ptr != end;) {
      Lum const& lum = randomLum (rng);
      Ray const ray = lum.emit (rng);
      ptr = tracePhoton (rng, ray, &lum, lum.power * scale, 1000, ptr, end);
    }
  }

  RandBits threadRNG () {
    auto const id = std::this_thread::get_id ();
    std::hash<std::thread::id> threadHasher;
    auto const seed = uint_fast32_t (threadHasher (id));
    return RandBits (seed);
  }

  void render (Options const& opts, Image& im) {
    Point const
      eye    {0,0,-1.5f},
      centre {0,0, 0.5f};

    float const a = float(im.width()) / float(im.height());
    Vector const
      horizontal {a,0,0},
      vertical   {0,1,0};

    float const
      xImageToWorld = 2.f / float(im.width() ),
      yImageToWorld = 2.f / float(im.height());

    int const
      nSamples = opts.sqrtNSamples * opts.sqrtNSamples,
      nThreads = 4,
      nPhotonsPerThread = opts.nPhotons/nThreads;

    fprintf (stderr,
      "Rendering to %ix%i @ %ispp\n"
      "Using %i photons, %i per estimate\n",
      im.width (), im.height (), nSamples,
      opts.nPhotons, opts.nNears
    );

    std::mutex mutex;
    int yProgress = 0;
    std::thread threads[nThreads];

    PhotonMap map;

    if (opts.indirect) {
      fprintf (stderr, "Casting photons...\n");

      std::vector<Photon> photonsRaw (opts.nPhotons);

      for (int iThread = 0; iThread != nThreads; iThread++) {
        auto& thread = threads[iThread];
        thread = std::thread (
          [&, iThread] () {
            RandBits rng = threadRNG ();
            Photon* photons = photonsRaw.data () + iThread * nPhotonsPerThread;
            castPhotons (rng, photons, photons + nPhotonsPerThread);
          }
        );
      }

      for (auto& thread : threads)
        thread.join ();

      fprintf (stderr, "Building photon map...\n");
      map = mapPhotons (std::move (photonsRaw));
    }

    fprintf (stderr, "Gathering...\n");
    using clock = std::chrono::steady_clock;
    auto gatherStartTime = clock::now ();

    for (int iThread = 0; iThread != nThreads; iThread++) {
      auto& thread = threads[iThread];
      thread = std::thread (
        [&, iThread] () {
          RandBits rng = threadRNG ();
          SampleSet samplePosns (opts.sqrtNSamples);
          NearSet nears (opts.nNears, opts.nearLimit);

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
                  xSample = (xPixelCentre + samplePos.x)*xImageToWorld - 1.f,
                  ySample = (yPixelCentre + samplePos.y)*yImageToWorld - 1.f;
                Point const onFrame
                  = centre + xSample*horizontal + ySample*vertical;
                Ray const ray (onFrame, unit(onFrame-eye));
                samples[i] = trace (opts, rng, map, nears, ray, nullptr, 10);
              }

              Colour const
                sum = std::accumulate (
                  samples, samples + nSamples, Colour::black()),
                mean = sum / float (nSamples);
              im.at (xPixel, yPixel) = toSRGB (mean);
            }
          }
        }
      );
    }

    for (auto& thread : threads)
      thread.join ();

    auto gatherStopTime = clock::now ();
    auto gatherDuration = gatherStopTime - gatherStartTime;
    auto gatherSeconds = std::chrono::duration<float> (gatherDuration);

    fprintf (stderr, "\rTook %.1f s\n", gatherSeconds.count ());
  }

  void writePPM (char const* path, Image const& im) {
    FILE* f = fopen (path, "wb");
    fprintf (f, "P6\n%d %d\n255\n", im.width(), im.height());
    fwrite (im.raw(), sizeof(Pixel), im.size(), f);
    fclose (f);
  }
}

int main () {
  using namespace RT;

  constexpr int const
    fracHD = 3,
    w = 1440 / fracHD,
    h = 1080 / fracHD;
  Image im (w, h);

  Options opts;
  opts.indirect = true;
  opts.sqrtNSamples = 5;
  opts.nPhotons = 2'000'000;
  opts.nNears = 500;
  opts.nearLimit = 0.4f;
  render (opts, im);

  fprintf (stderr, "Done\n");

  auto outTime = std::chrono::system_clock::to_time_t (
    std::chrono::system_clock::now ());
  char path[64];
  strftime (path, sizeof (path), "out-%FT%TZ.ppm", localtime (&outTime));
  writePPM (path, im);
}

