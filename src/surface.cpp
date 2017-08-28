
#include "surface.hpp"

namespace RT {
  RayHit<Surface const*> Sphere::intersect (Ray r) const {
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
    return { p, t, n, this };
  }

  Ray Sphere::randomEmission (RandBits& rng) const {
    Vector const normal = spherical (rng);
    Point const origin = centre + radius * normal;
    Vector const dir = lambert (rng);
    Basis const tang = Basis::fromK (normal);
    return Ray (origin, tang.outOf (dir));
  }

  Point Sphere::randomPoint (RandBits& rng, Point towards) const {
    Basis const facing = Basis::fromK (unit (towards - centre));
    return centre + radius * facing.outOf (hemispherical (rng));
  }

  RayHit<Surface const*> Plane::intersect (Ray r) const {
    float const denom = dot (r.disp, normal);
    if (std::abs(denom) < 0.000001f)
      return { };
    float const t = dot (onPlane - r.origin, normal) / denom;
    Vector const n = (denom<0.f)? normal : -normal;
    return { r.at(t), t, n, this };
  }

  Triangle::Triangle (Point p, Point a, Point b)
    : p0(p), abc(p0-a), def(p0-b)
    , basis (Basis::fromK (unit (cross (abc, def))))
  {
    /*fprintf (stderr, "made tri with normal (%f %f %f)\n",
      basis.k.x, basis.k.y, basis.k.z);*/
  }

  RayHit<Surface const*> Triangle::intersect (Ray r) const {
    Vector const
      jkl = p0 - r.origin,
      defxdisp = cross (def, r.disp),
      abcxjkl  = cross (abc, jkl);
    float const
      denom = dot (abc,    defxdisp),
      beta  = dot (jkl,    defxdisp) /  denom,
      gamma = dot (r.disp, abcxjkl ) /  denom,
      t     = dot (def,    abcxjkl ) / -denom;
    if (beta < 0.f || gamma < 0.f || (beta+gamma) > 1.f)
      return { };
    Vector const n = (r.disp.z < 0.f)? basis.k : -basis.k;
    return { r.at(t), t, n, this };
  }

  Ray Triangle::randomEmission (RandBits& rng) const {
    Point const origin = randomPoint (rng, {0,0,0});
    Vector const dir = lambert (rng);
    Vector const disp = basis.outOf (dir);
    /*fprintf (stderr, "emitted photon along (%f %f %f)\n",
      disp.x, disp.y, disp.z);*/
    return Ray (origin, disp);
  }

  Point Triangle::randomPoint (RandBits& rng, Point) const {
    float
      xi1 = canon (rng),
      xi2 = canon (rng);

    if (xi1 + xi2 > 1.f) {
      xi1 = (1.f - xi1);
      xi2 = (1.f - xi2);
    }

    return p0 - xi1 * abc - xi2 * def;
  }
}

