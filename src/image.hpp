#pragma once

#include "colour.hpp"

#include <cstdint>
#include <cmath>

namespace RT {
  template<typename X>
  struct Pixel {
    X r, g, b;
    Pixel () = default;
    constexpr Pixel (Black) : r(0), g(0), b(0) { }
    constexpr Pixel (X r, X g, X b) : r(r), g(g), b(b) { }
    Pixel (Pixel const&) = default;
  };

  template<typename X>
  Pixel<X>& operator += (Pixel<X>& onto, Pixel<X> other) {
    onto.r += other.r;
    onto.g += other.g;
    onto.b += other.b;
    return onto;
  }

  template<typename X>
  class Image {
    int wide, high;
    std::vector<Pixel<X>> pixels;

  public:
    Image (int w, int h)
      : wide (w), high (h)
      , pixels (size (), black)
    { }
  /*Image (Image&& other)
      : wide (std::exchange (other.wide, 0))
      , high (std::exchange (other.high, 0))
      , pixels (std::move (other.pixels))
    { }
    Image& operator = (Image&& other) {
      wide = std::exchange (other.wide, 0);
      high = std::exchange (other.high, 0);
      pixels = std::move (other.pixels);
      return *this;
    }*/
    Image (Image const&) = delete;
    Image& operator = (Image const&) = delete;

    int width  () const { return wide; }
    int height () const { return high; }
    size_t size () const { return wide * high; }
    Pixel<X> const* raw () const { return pixels.data (); }
    Pixel<X>& at (int x, int y)       { return pixels [y * wide + x]; }
    Pixel<X>  at (int x, int y) const { return pixels [y * wide + x]; }
  };

  template<typename X>
  static Image<X>& operator += (Image<X>& onto, Image<X> const& other) {
    int const
      w = std::min (onto.width (), other.width ()),
      h = std::min (onto.height (), other.height ());
    for (int y = 0; y != h; y++) {
      for (int x = 0; x != w; x++)
        onto.at (x, y) += other.at (x, y);
    }
    return onto;
  }

  float toSRGB1 (float raw) {
    float const
      c = std::max (0.f, std::min (1.f, raw)),
      lo = 12.92f*c,
      hi = 1.055f*std::pow (c, 1.f/2.4f) - 0.055f;
    return (c <= 0.0031308)? lo : hi;
  }

  uint8_t toU8 (float x) {
    return uint8_t (255.f * x);
  }

  Pixel<uint8_t> toSRGB (Pixel<float> c) {
    return
      { toU8 (toSRGB1 (c.r))
      , toU8 (toSRGB1 (c.g))
      , toU8 (toSRGB1 (c.b))
      };
  }
}

