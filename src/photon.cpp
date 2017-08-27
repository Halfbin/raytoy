
#include "photon.hpp"

#include "math.hpp"

#include <cassert>

namespace RT {
  RGBE::RGBE (Colour c) {
    int re, ge, be;
    float const
      rm = std::frexp (c.r, &re),
      gm = std::frexp (c.g, &ge),
      bm = std::frexp (c.b, &be);

    int maxExp = -128;
    if (rm != 0.f) maxExp = re;
    if (gm != 0.f) maxExp = std::max (maxExp, ge);
    if (bm != 0.f) maxExp = std::max (maxExp, be);
    if (maxExp == -128) maxExp = 0;

    r = (uint8_t) ((int) (rm * 255.f) >> (maxExp - re));
    g = (uint8_t) ((int) (gm * 255.f) >> (maxExp - ge));
    b = (uint8_t) ((int) (bm * 255.f) >> (maxExp - be));
    e = (uint8_t) (128 + maxExp);
    /*fprintf (stderr, "(%f %f %f) -> (%3i %3i %3i | %3i)\n",
      c.r, c.g, c.b, (int) r, (int) b, (int) g, (int) e);*/
  }

  RGBE::operator Colour () const {
    int const exp = (int) e - 128;
    Colour const c
      { std::scalbn ((int) r / 255.f, exp)
      , std::scalbn ((int) g / 255.f, exp)
      , std::scalbn ((int) b / 255.f, exp)
      };
    /*fprintf (stderr, "(%f %f %f) <- (%3i %3i %3i | %3i)\n",
      c.r, c.g, c.b, (int) r, (int) b, (int) g, (int) e);*/
    return c;
  }

  static float signumNZ (float f) {
    return (f >= 0.f)? 1.f : -1.f;
  }

  static int8_t toByte (float f) {
    return ((int) (f * 127.f));
  }

  Octo::Octo (Vector v) {
    using std::abs;

    v = unit (v);
    float const
      k = 1.f / (abs (v.x) + abs (v.y) + abs (v.z)),
      px = v.x * k,
      py = v.y * k;

    if (v.z > 0.f) {
      s = toByte (px);
      t = toByte (py);
    }
    else {
      float const
        sf = (1.f - abs (py)) * signumNZ (px),
        tf = (1.f - abs (px)) * signumNZ (py);
      s = toByte (sf);
      t = toByte (tf);
    }
    /*fprintf (stderr, "(%f %f %f) -> (%i %i)\n",
      v.x, v.y, v.z, (int) s, (int) t);*/
  }

  Octo::operator Vector () const {
    using std::abs;

    float const
      sf = ((int) s) / 127.f,
      tf = ((int) t) / 127.f,
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
          { return a.x < b.x; });
        break;
      case Axis::Y:
        std::sort (begin, end, [] (auto const& a, auto const& b)
          { return a.y < b.y; });
        break;
      case Axis::Z:
        std::sort (begin, end, [] (auto const& a, auto const& b)
          { return a.z < b.z; });
        break;
      default:
        assert (false);
    }

    Photon* median = begin + (end - begin) / 2;
    median->split = split;
    Vector newMins = mins,
           newMaxs = maxs;
    switch (split) {
      case Axis::X: newMins.x = newMaxs.x = median->x; break;
      case Axis::Y: newMins.y = newMaxs.y = median->y; break;
      case Axis::Z: newMins.z = newMaxs.z = median->z; break;
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
      mins.x = std::min (mins.x, photon.x);
      mins.y = std::min (mins.y, photon.y);
      mins.z = std::min (mins.z, photon.z);
      maxs.x = std::max (maxs.x, photon.x);
      maxs.y = std::max (maxs.y, photon.y);
      maxs.z = std::max (maxs.z, photon.z);
    }

  /*fprintf (stderr, "mins: (%f %f %f) maxs: (%f %f %f)\n",
      mins.x, mins.y, mins.z,
      maxs.x, maxs.y, maxs.z);*/
    auto* p = array.data ();
    buildPhotonMap (p, p + array.size(), mins, maxs);
  }

  static void findNearestPhotons
    ( Photon const* begin, Photon const* end
    , Point pos
    , NearSet& nears
    )
  {
    if (begin == end)
      return;

  //if (end - begin > PhotonMap::maxLeaf) {
      Photon const* node = begin + (end - begin) / 2;
      float const qd = norm2 (pos - node->position ());

      float planeDist = inf;
      switch (node->split) {
        case Axis::X: planeDist = pos.x - node->x; break;
        case Axis::Y: planeDist = pos.y - node->y; break;
        case Axis::Z: planeDist = pos.z - node->z; break;
        default: assert (false); //"tried to treat photon map leaf as node!"
      }

      float const planeQd = planeDist * planeDist;

      if (planeDist > 0.f) {
        findNearestPhotons (node+1, end, pos, nears);
        if (planeQd < nears.maxQd ())
          findNearestPhotons (begin, node, pos, nears);
      }
      else {
        findNearestPhotons (begin, node, pos, nears);
        if (planeQd < nears.maxQd ())
          findNearestPhotons (node+1, end, pos, nears);
      }

      nears.insert (qd, *node);
  /*}
    else {
      for (Photon const* photon = begin; photon != end; photon++) {
        float const qd = norm2 (pos - photon->position ());
        nears.insert (qd, *photon);
      }
    }*/
  }

  void PhotonMap::nearest (NearSet& nears, Point pos, Vector normal) const {
    nears.prepare (pos, normal);
    auto* ptr = array.data ();

    for (;;) {
      findNearestPhotons (ptr, ptr + array.size (), pos, nears);
      if (nears.full ())
        break;
      nears.broaden (1.4f);
    }
  }
}

