#pragma once

#include "linear.hpp"

namespace RT {
  struct Basis {
    Vector const i, j, k;

    Basis (Vector i, Vector j, Vector k)
      : i(i), j(j), k(k) { }
    Basis (Basis const&) = default;

    static Basis fromK (Vector n);

    Vector into (Vector v) const
      { return { dot (i, v), dot (j, v), dot (k, v) }; }

    Vector outOf (Vector v) const
      { return v.x * i + v.y * j + v.z * k; }
  };
}

