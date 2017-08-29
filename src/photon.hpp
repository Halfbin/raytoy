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

  class NearSet {
  public:
    struct Element {
      float qd;
      Photon photon;
      bool operator < (Element const& other) const {
        return qd < other.qd;
      }
    };

    /*static constexpr auto const cmp
      = [] (Element const& a, Element const& b) { return a.qd < b.qd; };*/

  private:
    size_t cap;
    float limitQd;
    size_t n;
    Element* elems;
    Vector normal;

  public:
    NearSet (size_t cap, float limitQd) :
      cap (cap), limitQd (limitQd), n (0), elems (new Element[cap]),
      normal {0,0,0}
    { }

    NearSet (NearSet const&) = delete;

    NearSet (NearSet&& other)
      : cap (std::exchange (other.cap, 0))
      , limitQd (std::exchange (other.limitQd, 0))
      , n (std::exchange (other.n, 0))
      , elems (std::exchange (other.elems, nullptr))
      , normal {0,0,0}
      { }

    ~NearSet () {
      delete [] elems;
    }

    bool full () const {
      return n == cap;
    }

    void clear () {
      n = 0;
    }

    float maxQd () const {
      if (!full ())
        return limitQd;
      else
        return elems[0].qd;
    }

    void broaden (float k) {
      limitQd *= k;
    }

    void insert (float qd, Photon photon) {
      if (-dot (photon.incoming (), normal) < 0.f)
        return;
      if (!full ()) {
        elems[n++] = { qd, photon };
        if (full ())
          std::make_heap (elems, elems+cap);
      }
      else if (qd < maxQd ()) {
        std::pop_heap (elems, elems+cap);
        elems[cap-1] = { qd, photon };
        std::push_heap (elems, elems+cap);
      }
    }

    void prepare (Point, Vector newNormal) {
      normal = newNormal;
      clear ();

    /*if (!full ()) {
        clear ();
        return;
      }

      float const newMaxQd = norm2 (elems[0].photon.position () - newPos);
      if (newMaxQd > maxQd ()) {
        clear ();
      }
      else {
        elems[0].qd = newMaxQd;
        for (Element* elem = elems+1; elem != end (); elem++) {
          elem->qd = norm2 (elem->photon.position () - newPos);
        }
        std::make_heap (elems, elems+cap);
      }*/
    }

    Element const* begin () const {
      return elems;
    }

    Element const* end () const {
      return elems+n;
    }
  };

  struct PhotonMap {
  //static constexpr int const maxLeaf = 4;

    std::vector<Photon> array;

    PhotonMap (std::vector<Photon> raw);
    PhotonMap (PhotonMap const&) = delete;
    PhotonMap& operator = (PhotonMap const&) = delete;
    PhotonMap (PhotonMap&&) = delete;
    PhotonMap& operator = (PhotonMap&&) = delete;

    Photon const* begin () const { return array.data (); }
    Photon const* end () const { return begin () + array.size (); }

    void nearest (NearSet&, Point, Vector normal) const;
  };
}

