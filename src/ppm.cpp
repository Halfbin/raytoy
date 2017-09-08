
#include "ppm.hpp"

#include "math.hpp"
#include "basis.hpp"

#include <cassert>
#include <algorithm>
#include <cfloat>

namespace RT {
  static Colour addPhoton (PathNode const& node, Photon photon) {
    float const cosine = std::max (0.f, -dot (photon.incoming (), node.normal));
    return node.k * photon.power () * cosine * node.matDiffuse;
  }

  static float metric (PathNode const& node, Point pos, float k = 2.f) {
    Vector const sep = pos - node.position;
    float planeDist = dot (sep, node.normal);
    Point const eff = pos + (k - 1.f) * planeDist * node.normal;
    return norm2 (eff - node.position);
  }

  static void addPhotonsInRadius
    ( Gather& g, float r2
    , Photon const* begin, Photon const* end
    , PathNode const& node
    )
  {
    if (end - begin < 16) {
      for (Photon const* photon = begin; photon != end; photon++) {
        float const d2 = metric (node, photon->position);
        if (d2 >= 0.f && d2 < r2) {
          g.contrib += addPhoton (node, *photon);
          g.nPhotons++;
        }
      }
    }
    else {
      Photon const* mid = begin + (end - begin) / 2;
      float const d2 = metric (node, mid->position);

      if (d2 >= 0.f && d2 < r2) {
        g.contrib += addPhoton (node, *mid);
        g.nPhotons++;
      }

      float planeDist = inf;
      switch (mid->split) {
        case Axis::X: planeDist = node.position.x - mid->position.x; break;
        case Axis::Y: planeDist = node.position.y - mid->position.y; break;
        case Axis::Z: planeDist = node.position.z - mid->position.z; break;
        default: assert (false); //"tried to treat photon map leaf as node!"
      }
      float const planeD2 = planeDist * planeDist;

      if (planeDist > 0.f || r2 > planeD2)
        addPhotonsInRadius (g, r2, mid+1, end, node);
      if (planeDist < 0.f || r2 > planeD2)
        addPhotonsInRadius (g, r2, begin, mid, node);
    }
  }

  Gather PathNode::gather (float r2, PhotonMap const& map) const {
    Gather g (px, py);
    addPhotonsInRadius (g, r2, map.begin (), map.end (), *this);
    return g;
  }

  void Stats::add (Colour contrib, int n) {
    tau += contrib;
    dn += n;
  }

  void Stats::update (float alpha) {
    if (dn == 0)
      return;

    float const ratio
      = (nPhotons + alpha * dn)
      / (nPhotons +         dn);

    r2 *= ratio;
    if (r2 <= 1e-8f || ratio >= 1.f)
      fprintf (stderr, "ugh\n");
    tau *= ratio;
    nPhotons += alpha * dn;

    dn = 0;
  }
}

