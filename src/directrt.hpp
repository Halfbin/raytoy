
#pragma once

#include "scene.hpp"
#include "random.hpp"
#include "ray.hpp"

namespace RT {
  Colour directTermRT
    ( int branch
    , Scene const& scene
    , RandBits& rng
    , RayHit<Item const*> hit
    );
}

