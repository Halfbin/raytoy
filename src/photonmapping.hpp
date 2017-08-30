
#pragma once

#include "photon.hpp"
#include "scene.hpp"

namespace RT {
  int castPhotons
    ( RandBits& rng
    , Scene const& scene
    , Photon* photons, Photon* end
    );
}

