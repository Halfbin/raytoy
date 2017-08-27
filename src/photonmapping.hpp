
#pragma once

#include "photon.hpp"
#include "scene.hpp"

namespace RT {
  Colour indirectTermPM
    ( NearSet& nears
    , PhotonMap const& map
    , RayHit<const Item*> hit
    );

  void castPhotons
    ( RandBits& rng
    , Scene const& scene
    , Photon* photons, Photon* end
    );

  Photon* tracePhoton
    ( RandBits& rng
    , Scene const& scene
    , Ray ray
    , Item const*
    , Colour power
    , int ttl
    , Photon* photons, Photon* end
    , bool indirect = false
    );
}

