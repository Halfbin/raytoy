#pragma once

#include "colour.hpp"

#include <vector>
#include <cstdint>
#include <cmath>

namespace RT {
  template<typename Pixel>
  class Image {
    int wide, high;
    std::vector<Pixel> pixels;

  public:
    Image ()
      : wide (0), high (0)
    { }
    Image (int w, int h, Pixel init = Pixel ())
      : wide (w), high (h)
      , pixels (size (), init)
    { }
    Image (Image&& other)
      : wide (std::exchange (other.wide, 0))
      , high (std::exchange (other.high, 0))
      , pixels (std::move (other.pixels))
    { }
    Image& operator = (Image&& other) {
      wide = std::exchange (other.wide, 0);
      high = std::exchange (other.high, 0);
      pixels = std::move (other.pixels);
      return *this;
    }
    Image (Image const& other)
      : wide (other.wide), high (other.high)
      , pixels (other.pixels)
    { }
    Image& operator = (Image const& other) {
      clear ();
      wide = other.wide; high = other.high;
      pixels = other.pixels;
      return *this;
    }

    void clear () {
      wide = 0; high = 0;
      pixels.clear ();
    }

    int width  () const { return wide; }
    int height () const { return high; }
    size_t size () const { return wide * high; }

    Pixel*       raw ()       { return pixels.data (); }
    Pixel const* raw () const { return pixels.data (); }

    Pixel&       at (int x, int y)       { return pixels [y * wide + x]; }
    Pixel const& at (int x, int y) const { return pixels [y * wide + x]; }

    Pixel*       begin ()       { return raw (); }
    Pixel const* begin () const { return raw (); }
    Pixel*       end   ()       { return begin () + size (); }
    Pixel const* end   () const { return begin () + size (); }

    Pixel const* cbegin () const { return begin (); }
    Pixel const* cend   () const { return end (); }
  };

/*template<typename X>
  static Image<X>& operator += (Image<X>& onto, Image<X> const& other) {
    int const
      w = std::min (onto.width (), other.width ()),
      h = std::min (onto.height (), other.height ());
    for (int y = 0; y != h; y++) {
      for (int x = 0; x != w; x++)
        onto.at (x, y) += other.at (x, y);
    }
    return onto;
  }*/

  static inline float toSRGB1 (float raw) {
    float const
      c = std::max (0.f, std::min (1.f, raw)),
      lo = 12.92f*c,
      hi = 1.055f*std::pow (c, 1.f/2.4f) - 0.055f;
    return (c <= 0.0031308)? lo : hi;
  }

  static inline uint8_t toU8 (float x) {
    return uint8_t (255.f * x);
  }

  static inline RGB<uint8_t> toSRGB (Colour c) {
    return
      { toU8 (toSRGB1 (c.r))
      , toU8 (toSRGB1 (c.g))
      , toU8 (toSRGB1 (c.b))
      };
  }
}

