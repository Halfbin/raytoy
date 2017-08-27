
#include "ppm.hpp"

#include "math.hpp"
#include "basis.hpp"

#include <cassert>

namespace RT {
  static void addPhoton (PathNode& node, Photon photon) {
    Basis const tang = Basis::fromK (node.normal);
    node.tau += photon.power () * node.material->brdfDiffuse
      (node.incident, tang.into (photon.incoming ()));
  }

  static int addPhotonsInRadius
    ( Photon const* begin, Photon const* end
    , PathNode& node
    )
  {
    if (begin == end)
      return 0;

    Photon const* mid = begin + (end - begin) / 2;
    float const qd = norm2 (node.position - mid->position ());

    float planeDist = inf;
    switch (mid->split) {
      case Axis::X: planeDist = node.position.x - mid->x; break;
      case Axis::Y: planeDist = node.position.y - mid->y; break;
      case Axis::Z: planeDist = node.position.z - mid->z; break;
      default: assert (false); //"tried to treat photon map leaf as node!"
    }

    float const planeQd = planeDist * planeDist;

    int added = 0;
    if (node.radius2 > qd) {
      addPhoton (node, *mid);
      added++;
    }
    if (planeDist > 0.f || node.radius2 > planeQd)
      added += addPhotonsInRadius (mid+1, end, node);
    if (planeDist < 0.f || node.radius2 > planeQd)
      added += addPhotonsInRadius (begin, mid, node);
    return added;
  }

  void PathNode::update (PhotonMap const& map, float alpha) {
    int const
      dn = addPhotonsInRadius (map.begin (), map.end (), *this),
      adn = (int) (alpha * dn);

    if (dn > 0) {
      float const ratio
        = (float) (nPhotons + adn)
        / (float) (nPhotons +  dn);
      radius2 *= ratio;
      tau     *= ratio;
    }
    nPhotons += adn;
  }
}

