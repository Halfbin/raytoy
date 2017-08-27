#pragma once

#include "linear.hpp"
#include "math.hpp"

namespace RT {
  struct Ray {
    Point const origin;
    Vector const disp;

    Ray () : origin{0,0,0}, disp{0,0,0} { }
    Ray (Point o, Vector d)
      : origin(o), disp(unit(d)) { }

    Point at (float t) const { return origin + t * disp; }

    explicit operator bool () const
      { return disp == Vector{0,0,0}; }
  };

  template<typename Occ>
  struct RayHit {
    Point  position;
    float  t;
    Vector normal;
    Occ    occ;

    RayHit () : t(inf) { }
    RayHit (RayHit const&) = default;
    RayHit& operator= (RayHit const&) = default;
    RayHit (Point p, float t, Vector n, Occ o)
      : position(p), t(t), normal(n), occ(o) { }

    explicit operator bool () const { return std::isfinite (t); }
  };
}

