#include "surface.hpp"

namespace RT {
  RayHit<Surface const*> Sphere::intersect (Ray r) const {
    Vector const
      s = r.origin - centre;
    float const
      ds = dot (r.disp, s),
      dd = norm2 (r.disp),
      discrim = ds*ds - dd*(norm2(s) - radius*radius);
    if (discrim < 0.f) return { };

    float const root = std::sqrt(discrim),
                s1 = -root - ds,
                s2 =  root - ds;
    float tdd = std::min (s1, s2);
    if (tdd < 0.f) tdd = std::max (s1, s2);
    if (tdd < 0.f) return { };

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
    if (std::abs (denom) < 0.000001f)
      return { };
    float const t = dot (onPlane - r.origin, normal) / denom;
    Vector const n = (denom < 0.f)? normal : -normal;
    return { r.at(t), t, n, this };
  }

  RayHit<Surface const*> Triangle::intersect (Ray r) const {
    Vector const
      dxpb   = cross (r.disp, pb),
      sep    = r.origin - p,
      sepxpa = cross (sep, pa);
    float const
      det = dot (pa, dxpb),
      rdet = 1.f / det, // == inf when |det| <= eps
      u = dot (sep,    dxpb  ) * rdet,
      v = dot (r.disp, sepxpa) * rdet,
      t = dot (pb,     sepxpa) * rdet;
    if (u < 0.f || v < 0.f || u + v > 1.f)
      return { };
    Vector const
      in = (det > 0.f)? n : -n;
    return { r.at (t), t, in, this };
  }

  Ray Triangle::randomEmission (RandBits& rng) const {
    Point const origin = randomPoint (rng, {0,0,0});
    Vector const dir = lambert (rng);
    Basis const tang = Basis::fromK (n);
    Vector const disp = tang.outOf (dir);
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
    return p + xi1 * pa + xi2 * pb;
  }
}

