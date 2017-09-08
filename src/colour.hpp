#pragma once

#include <algorithm>
#include <type_traits>

namespace RT {
  enum Black { black };

  template<typename X>
  struct RGB {
    X r, g, b;

    RGB () = default;
    constexpr RGB (Black) : r(0), g(0), b(0) { }
    constexpr RGB (X r, X g, X b) : r(r), g(g), b(b) { }
    RGB (RGB const&) = default;

    explicit operator bool () const { return *this != black; }
    bool operator== (RGB const& other) const
      { return r==other.r && g==other.g && b==other.b; }
    bool operator!= (RGB const& other) const { return !(*this == other); }
  };

  template<typename X>
  static inline RGB<X> operator* (RGB<X> a, RGB<X> b)
    { return {a.r*b.r, a.g*b.g, a.b*b.b}; }
  template<typename X>
  static inline RGB<X> operator* (RGB<X> c, X k)
    { return {c.r*k, c.g*k, c.b*k}; }
  template<typename X>
  static inline RGB<X> operator* (X k, RGB<X> c)
    { return c*k; }
  template<typename X>
  static inline RGB<X> operator/ (RGB<X> c, X k)
    { return {c.r/k, c.g/k, c.b/k}; }
  template<typename X>
  static inline RGB<X> operator+ (RGB<X> a, RGB<X> b)
    { return {a.r+b.r, a.g+b.g, a.b+b.b}; }

  template<typename X>
  static inline RGB<X>& operator += (RGB<X>& a, RGB<X> b)
    { return (a = a+b); }
  template<typename X>
  static inline RGB<X>& operator *= (RGB<X>& a, RGB<X> b)
    { return (a = a*b); }
  template<typename X>
  static inline RGB<X>& operator *= (RGB<X>& c, X k)
    { return (c = k*c); }

  template<typename X>
  static inline X max (RGB<X> c)
    { return std::max ({c.r, c.g, c.b}); }

  template<typename X, typename = std::enable_if_t<std::is_floating_point_v<X>>>
  static inline RGB<X> clamp (RGB<X> c) {
    return RGB<X>
      { std::max ((X) 0, std::min ((X) 1, c.r))
      , std::max ((X) 0, std::min ((X) 1, c.g))
      , std::max ((X) 0, std::min ((X) 1, c.b))
      };
  }

  using Colour = RGB<float>;
}

