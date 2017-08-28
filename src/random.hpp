#pragma once

#include "linear.hpp"

#include <random>

namespace RT {
  using RandBits = std::mt19937;
  static std::uniform_real_distribution<float> canon (0.f, 1.f);

  RandBits threadRNG (int seed = 1234);

  Vector lambert       (RandBits&);
  Vector spherical     (RandBits&);
  Vector hemispherical (RandBits&);
  Vector disc          (RandBits&);
  Vector conical       (RandBits&, float halfApex);
}

