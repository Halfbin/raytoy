
#include "basis.hpp"

namespace RT {
  Basis Basis::fromK (Vector n) {
    float const
      s = std::copysignf (1.f, n.z),
      a = -1.f / (s + n.z),
      b = n.x * n.y * a,
      nx2a = n.x * n.x * a,
      ny2a = n.y * n.y * a;
    return Basis
      ( { 1.f + s*nx2a, s*b,      -s*n.x }
      , { b,            s + ny2a, -n.y   }
      , n
      );
  }
}

