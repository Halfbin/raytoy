#pragma once

#include <cmath>
#include <algorithm>

namespace RT {
  struct Vector {
    union {
      float components[3];
      struct { float x, y, z; };
    };

    constexpr Vector ()
      : x(0), y(0), z(0) { }
    constexpr Vector (float x, float y, float z)
      : x(x), y(y), z(z) { }
    Vector (Vector const&) = default;

    bool operator== (Vector const& other) const
      { return x==other.x && y==other.y && z==other.z; }
    bool operator!= (Vector const& other) const
      { return !(*this == other); }
  };

  constexpr Vector const zero = Vector ();

  struct Point {
    union {
      float components[3];
      struct { float x, y, z; };
    };

    constexpr Point ()
      : x(0), y(0), z(0) { }
    constexpr Point (float x, float y, float z)
      : x(x), y(y), z(z) { }
    Point (Point const&) = default;
  };

  static inline Vector operator+ (Vector a, Vector b)
    { return { a.x+b.x, a.y+b.y, a.z+b.z }; }
  static inline Point  operator+ (Point p,  Vector d)
    { return { p.x+d.x, p.y+d.y, p.z+d.z }; }
  static inline Point  operator- (Point p,  Vector d)
    { return { p.x-d.x, p.y-d.y, p.z-d.z }; }
  static inline Vector operator- (Point p,  Point  q)
    { return { p.x-q.x, p.y-q.y, p.z-q.z }; }
  static inline Vector operator- (Vector a, Vector b)
    { return { a.x-b.x, a.y-b.y, a.z-b.z }; }
  static inline Vector operator* (Vector a, Vector b)
    { return { a.x*b.x, a.y*b.y, a.z*b.z }; }

  static inline Vector operator* (float k, Vector v)
    { return { k*v.x, k*v.y, k*v.z }; }
  static inline Vector operator* (Vector v, float k) { return k*v; }

  static inline Vector operator- (Vector v) { return { -v.x, -v.y, -v.z }; }

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

