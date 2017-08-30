
#include "ppm.hpp"

#include "math.hpp"
#include "basis.hpp"

#include <cassert>
#include <algorithm>
#include <cfloat>

namespace RT {
  static float addPhoton (PathNode& node, Photon photon) {
    float const cosine = std::max (0.f, -dot (photon.incoming (), node.normal));
    node.tau += photon.power () * cosine * node.matDiffuse;
    return 1.f;//cosine;
  }

  static float metric (PathNode const& node, Point pos, float k = 2.f) {
    Vector const sep = pos - node.position;
    float planeDist = dot (sep, node.normal);
    Point const eff = pos + (k - 1.f) * planeDist * node.normal;
    return norm2 (eff - node.position);
  }

  static float addPhotonsInRadius
    ( Photon const* begin, Photon const* end
    , PathNode& node
    )
  {
    float added = 0;

    float const r2 = node.radius * node.radius;

    if (end - begin < 16) {
      for (Photon const* photon = begin; photon != end; photon++) {
        float const d2 = metric (node, photon->position);
        if (d2 >= 0.f && d2 < r2)
          added += addPhoton (node, *photon);
      }
    }
    else {
      Photon const* mid = begin + (end - begin) / 2;
      float const d2 = metric (node, mid->position);

      float planeDist = inf;
      switch (mid->split) {
        case Axis::X: planeDist = node.position.x - mid->position.x; break;
        case Axis::Y: planeDist = node.position.y - mid->position.y; break;
        case Axis::Z: planeDist = node.position.z - mid->position.z; break;
        default: assert (false); //"tried to treat photon map leaf as node!"
      }

      if (d2 >= 0.f && d2 < r2)
        added += addPhoton (node, *mid);
      if (planeDist > 0.f || node.radius > planeDist)
        added += addPhotonsInRadius (mid+1, end, node);
      if (planeDist < 0.f || node.radius > planeDist)
        added += addPhotonsInRadius (begin, mid, node);
    }

    return added;
  }

  void PathNode::update (PhotonMap const& map, float alpha) {
    float const
      dn = addPhotonsInRadius (map.begin (), map.end (), *this);
    if (dn == 0.f)
      return;

    float const ratio
      = (nPhotons + alpha * dn)
      / (nPhotons +         dn);

    //fprintf (stderr, "r: %f, rat: %f\n", radius, ratio);
    radius *= std::sqrt (ratio);
  //radius = std::max (radius, limit);
    if (radius == 0.f || ratio >= 1.f)
      fprintf (stderr, "ugh\n");
    tau *= ratio;
    nPhotons += alpha * dn;
  }
}

