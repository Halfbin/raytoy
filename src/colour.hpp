#pragma once

#include <algorithm>

namespace RT {
  enum Black { black };

  struct Colour {
    float r, g, b;

    Colour () = default;
    constexpr Colour (Black) : r(0), g(0), b(0) { }
    constexpr Colour (float r, float g, float b) : r(r), g(g), b(b) { }
    Colour (Colour const&) = default;

    explicit operator bool () const { return *this != black; }
    bool operator== (Colour const& other) const
      { return r==other.r && g==other.g && b==other.b; }
    bool operator!= (Colour const& other) const { return !(*this == other); }
  };

  static inline Colour operator* (Colour a, Colour b) { return {a.r*b.r, a.g*b.g, a.b*b.b}; }
  static inline Colour operator* (Colour c, float k) { return {c.r*k, c.g*k, c.b*k}; }
  static inline Colour operator* (float k, Colour c) { return c*k; }
  static inline Colour operator/ (Colour c, float k) { return {c.r/k, c.g/k, c.b/k}; }
  static inline Colour operator+ (Colour a, Colour b) { return {a.r+b.r, a.g+b.g, a.b+b.b}; }

  static inline Colour& operator += (Colour& a, Colour b) { return (a = a+b); }
  static inline Colour& operator *= (Colour& a, Colour b) { return (a = a*b); }
  static inline Colour& operator *= (Colour& c, float k) { return (c = k*c); }

  static inline float max (Colour c) { return std::max (c.r, std::max (c.g, c.b)); }

  static inline Colour clamp (Colour c) {
    return Colour
      { std::max (0.f, std::min (1.f, c.r))
      , std::max (0.f, std::min (1.f, c.g))
      , std::max (0.f, std::min (1.f, c.b))
      };
  }
}

