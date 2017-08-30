#pragma once

#include "ray.hpp"
#include "random.hpp"
#include "basis.hpp"
#include "linear.hpp"
#include "math.hpp"

namespace RT {
  class Surface {
  public:
    virtual RayHit<Surface const*> intersect (Ray) const = 0;
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

    Sphere (Point o, float r)
      : centre(o), radius(r) { }

    RayHit<Surface const*> intersect (Ray r) const;
    Ray randomEmission (RandBits& rng) const;
    Point randomPoint (RandBits& rng, Point towards) const;

    float area () const
      { return 4.f * pi * radius * radius; }
    float boundingRadius () const
      { return radius; }
    Point boundingCentre () const
      { return centre; }
  };

  struct Plane : public Surface {
    Vector const normal;
    Point  const onPlane;

    Plane (Vector n, float d)
      : normal(unit(n))
      , onPlane(Point{0,0,0} + d*n)
      { }

    RayHit<Surface const*> intersect (Ray r) const;
  };

/*struct Paragram : public BoundedSurface {
    Point const corner;
    Vector const a, b;
  };*/

  struct Triangle : public BoundedSurface {
    Point const p;
    Vector const pa, pb, n;

    Triangle (Point p, Point a, Point b)
      : p (p), pa (a - p), pb (b - p)
      , n (unit (cross (pa, pb)))
      { }

    RayHit<Surface const*> intersect (Ray r) const;
    Ray randomEmission (RandBits& rng) const;
    Point randomPoint (RandBits& rng, Point) const;

    float area () const
      { return dot (pa, pb); }
    float boundingRadius () const
      { return norm ((pa + pb) * (1.f/3.f)); }
    Point boundingCentre () const
      { return p + (pa + pb) * (1.f/3.f); }
  };
}

