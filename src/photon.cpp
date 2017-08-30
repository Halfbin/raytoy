
#include "photon.hpp"

#include "math.hpp"

#include <cassert>
#include <limits>

namespace RT {
  constexpr int const
    rgbeScale = std::numeric_limits<RGBE::Component>::max(),
    rgbeBias = rgbeScale / 2;

  RGBE::RGBE (Colour c) {
    int re, ge, be;
    float const
      rm = std::frexp (c.r, &re),
      gm = std::frexp (c.g, &ge),
      bm = std::frexp (c.b, &be);

    int maxExp = -rgbeBias;
    if (rm != 0.f) maxExp = re;
    if (gm != 0.f) maxExp = std::max (maxExp, ge);
    if (bm != 0.f) maxExp = std::max (maxExp, be);
    if (maxExp == -rgbeBias) maxExp = 0;

    r = (RGBE::Component) ((int) (rm * (float) rgbeScale) >> (maxExp - re));
    g = (RGBE::Component) ((int) (gm * (float) rgbeScale) >> (maxExp - ge));
    b = (RGBE::Component) ((int) (bm * (float) rgbeScale) >> (maxExp - be));
    e = (RGBE::Component) (rgbeBias + maxExp);
    /*fprintf (stderr, "(%f %f %f) -> (%3i %3i %3i | %3i)\n",
      c.r, c.g, c.b, (int) r, (int) b, (int) g, (int) e);*/
  }

  RGBE::operator Colour () const {
    int const exp = (int) e - rgbeBias;
    Colour const c
      { std::scalbn ((int) r / (float) rgbeScale, exp)
      , std::scalbn ((int) g / (float) rgbeScale, exp)
      , std::scalbn ((int) b / (float) rgbeScale, exp)
      };
    /*fprintf (stderr, "(%f %f %f) <- (%3i %3i %3i | %3i)\n",
      c.r, c.g, c.b, (int) r, (int) b, (int) g, (int) e);*/
    return c;
  }

  static float signumNZ (float f) {
    return (f >= 0.f)? 1.f : -1.f;
  }

  static constexpr float const
    octoScale = (float) std::numeric_limits<Octo::Component>::max ();

  static Octo::Component toComponent (float f) {
    return (int) (f * octoScale);
  }

  Octo::Octo (Vector v) {
    using std::abs;

    v = unit (v);
    float const
      k = 1.f / (abs (v.x) + abs (v.y) + abs (v.z)),
      px = v.x * k,
      py = v.y * k;

    if (v.z > 0.f) {
      s = toComponent (px);
      t = toComponent (py);
    }
    else {
      float const
        sf = (1.f - abs (py)) * signumNZ (px),
        tf = (1.f - abs (px)) * signumNZ (py);
      s = toComponent (sf);
      t = toComponent (tf);
    }
    /*fprintf (stderr, "(%f %f %f) -> (%i %i)\n",
      v.x, v.y, v.z, (int) s, (int) t);*/
  }

  Octo::operator Vector () const {
    using std::abs;

    float const
      sf = ((int) s) / octoScale,
      tf = ((int) t) / octoScale,
      z = 1.f - abs (sf) - abs (tf);

    Vector v;
    if (z < 0.f) {
      float const
        x = (1.f - abs (tf)) * signumNZ (sf),
        y = (1.f - abs (sf)) * signumNZ (tf);
      v = unit ({ x, y, z });
    }
    else {
      v = unit ({ sf, tf, z });
    }

    assert (norm2 (v) - 1.f < 0.001f);
    /*fprintf (stderr, "(%+f %+f %+f) <- (%i %i)\n",
      v.x, v.y, v.z, (int) s, (int) t);*/
    return v;
  }

  static void buildPhotonMap
    ( Photon* begin
    , Photon* end
    , Vector mins
    , Vector maxs
    )
  {
    if (begin == end)//(end - begin <= PhotonMap::maxLeaf)
      return;

    Vector extents = maxs - mins;
    float longest = max (extents);
    Axis split;
    if      (longest == extents.x) split = Axis::X;
    else if (longest == extents.y) split = Axis::Y;
    else                           split = Axis::Z;

    switch (split) {
      case Axis::X:
        std::sort (begin, end, [] (auto const& a, auto const& b)
          { return a.position.x < b.position.x; });
        break;
      case Axis::Y:
        std::sort (begin, end, [] (auto const& a, auto const& b)
          { return a.position.y < b.position.y; });
        break;
      case Axis::Z:
        std::sort (begin, end, [] (auto const& a, auto const& b)
          { return a.position.z < b.position.z; });
        break;
      default:
        assert (false);
    }

    Photon* median = begin + (end - begin) / 2;
    median->split = split;
    Vector newMins = mins,
           newMaxs = maxs;
    switch (split) {
      case Axis::X: newMins.x = newMaxs.x = median->position.x; break;
      case Axis::Y: newMins.y = newMaxs.y = median->position.y; break;
      case Axis::Z: newMins.z = newMaxs.z = median->position.z; break;
      default: assert (false);
    }

    buildPhotonMap (begin, median, mins, newMaxs);
    buildPhotonMap (median+1, end, newMins, maxs);
  }

  PhotonMap::PhotonMap (std::vector<Photon> raw)
    : array (std::move (raw))
  {
    Vector mins {inf,inf,inf}, maxs{-inf,-inf,-inf};
    for (auto& photon : array) {
      mins.x = std::min (mins.x, photon.position.x);
      mins.y = std::min (mins.y, photon.position.y);
      mins.z = std::min (mins.z, photon.position.z);
      maxs.x = std::max (maxs.x, photon.position.x);
      maxs.y = std::max (maxs.y, photon.position.y);
      maxs.z = std::max (maxs.z, photon.position.z);
    }

  /*fprintf (stderr, "mins: (%f %f %f) maxs: (%f %f %f)\n",
      mins.x, mins.y, mins.z,
      maxs.x, maxs.y, maxs.z);*/
    auto* p = array.data ();
    buildPhotonMap (p, p + array.size(), mins, maxs);
  }
}

