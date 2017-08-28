
#pragma once

#include "photon.hpp"
#include "scene.hpp"

namespace RT {
  Colour indirectTermPM
    ( NearSet& nears
    , PhotonMap const& map
    , RayHit<const Item*> hit
    );

  int castPhotons
    ( RandBits& rng
    , Scene const& scene
    , Photon* photons, Photon* end
    );
}

