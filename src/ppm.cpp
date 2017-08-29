
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
    return cosine;
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
        float const qd = norm2 (node.position - photon->position);
        if (qd < r2)
          added += addPhoton (node, *photon);
      }
    }
    else {
      Photon const* mid = begin + (end - begin) / 2;
      float const qd = norm2 (node.position - mid->position);

      float planeDist = inf;
      switch (mid->split) {
        case Axis::X: planeDist = node.position.x - mid->position.x; break;
        case Axis::Y: planeDist = node.position.y - mid->position.y; break;
        case Axis::Z: planeDist = node.position.z - mid->position.z; break;
        default: assert (false); //"tried to treat photon map leaf as node!"
      }

      if (qd < r2)
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

