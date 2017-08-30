#pragma once

#include "linear.hpp"
#include "colour.hpp"

#include <cstdint>
#include <vector>

namespace RT {
  struct RGBE {
    using Component = uint16_t;
    Component r, g, b, e;

    RGBE () = default;
    RGBE (RGBE const&) = default;
    RGBE& operator = (RGBE const&) = default;

    explicit RGBE (Colour c);
    explicit operator Colour () const;
  };

  struct Octo {
    using Component = int16_t;
    Component s, t;

    Octo () = default;
    Octo (Octo const&) = default;
    Octo& operator = (Octo const&) = default;

    explicit Octo (Vector v);
    explicit operator Vector () const;
  };

  enum class Axis : char { none, X, Y, Z };

  struct Photon {
    Point position;
    RGBE powerPacked;
    Octo incomingPacked;
    Axis split;
    char pad[7];

    Photon () = default;
    Photon (Point p, Vector i, Colour c)
      : position (p)
      , powerPacked (c)
      , incomingPacked (i)
      , split (Axis::none)
      { }

    Colour power () const {
      return static_cast<Colour> (powerPacked);
    }

    Vector incoming () const {
      return static_cast<Vector> (incomingPacked);
    }
  };

  static_assert (sizeof (Photon) == 32, "Photon miscompiled");

  struct PhotonMap {
    std::vector<Photon> array;

    PhotonMap (std::vector<Photon> raw);
    PhotonMap (PhotonMap const&) = delete;
    PhotonMap& operator = (PhotonMap const&) = delete;
    PhotonMap (PhotonMap&&) = delete;
    PhotonMap& operator = (PhotonMap&&) = delete;

    Photon const* begin () const { return array.data (); }
    Photon const* end () const { return begin () + array.size (); }
  };
}

