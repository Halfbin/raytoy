#pragma once

#include <cmath>
#include <algorithm>

namespace RT {
  struct Vector {
    union {
      float components[4];
      struct { float x, y, z, w; };
    };

    constexpr Vector ()
      : x(0), y(0), z(0), w(0) { }
    constexpr Vector (float x, float y, float z, float=0)
      : x(x), y(y), z(z), w(0) { }
    Vector (Vector const&) = default;

    bool operator== (Vector const& other) const
      { return x==other.x && y==other.y && z==other.z && w==other.w; }
    bool operator!= (Vector const& other) const
      { return !(*this == other); }
  };

  constexpr Vector const zero = Vector ();

  struct Point {
    union {
      float components[4];
      struct { float x, y, z, w; };
    };

    constexpr Point ()
      : x(0), y(0), z(0), w(1) { }
    constexpr Point (float x, float y, float z, float=1)
      : x(x), y(y), z(z), w(1) { }
    Point (Point const&) = default;
  };

  static inline Vector operator+ (Vector a, Vector b)
    { return { a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w }; }
  static inline Point  operator+ (Point p,  Vector d)
    { return { p.x+d.x, p.y+d.y, p.z+d.z, p.w+d.w }; }
  static inline Point  operator- (Point p,  Vector d)
    { return { p.x-d.x, p.y-d.y, p.z-d.z, p.w-d.w }; }
  static inline Vector operator- (Point p,  Point  q)
    { return { p.x-q.x, p.y-q.y, p.z-q.z, p.w-q.w }; }
  static inline Vector operator- (Vector a, Vector b)
    { return { a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w }; }
  static inline Vector operator* (Vector a, Vector b)
    { return { a.x*b.x, a.y*b.y, a.z*b.z, a.w*b.w }; }

  static inline Vector operator* (float k, Vector v)
    { return { k*v.x, k*v.y, k*v.z, k*v.w }; }
  static inline Vector operator* (Vector v, float k) { return k*v; }

  static inline Vector operator- (Vector v) { return { -v.x, -v.y, -v.z, -v.w }; }

  static inline float dot (Vector a, Vector b)
    { return a.x*b.x + a.y*b.y + a.z*b.z; }

  static inline float  norm2 (Vector v) { return dot (v, v); }
  static inline float  norm  (Vector v) { return std::sqrt (norm2 (v)); }
  static inline Vector unit  (Vector v) { return v * (1.f / norm (v)); }

  static inline Vector cross (Vector a, Vector b) {
    return
      { a.y*b.z - a.z*b.y
      , a.z*b.x - a.x*b.z
      , a.x*b.y - a.y*b.x
      };
  }

  static inline float max (Vector v) { return std::max (v.x, std::max (v.y, v.z)); }
}

