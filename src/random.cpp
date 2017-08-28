
#include "random.hpp"
#include "math.hpp"

#include <thread>

namespace RT {
  RandBits threadRNG (int extra) {
    auto const id = std::this_thread::get_id ();
    std::hash<std::thread::id> threadHasher;
    auto const seed = uint_fast32_t (threadHasher (id));
    return RandBits (seed ^ extra);
  }

  Vector lambert (RandBits& rng) {
    float const
      cosTheta = std::sqrt (canon (rng)),
      phi      = 2.f * pi * canon (rng),
      sinTheta = std::sqrt (1.f - cosTheta * cosTheta),
      x = std::cos (phi) * sinTheta,
      y = std::sin (phi) * sinTheta,
      z = cosTheta;
    return {x,y,z};
  }

  Vector spherical (RandBits& rng) {
    float const
      cosTheta = 2.f * canon (rng) - 1.f,
      phi = 2.f * pi * canon (rng),
      sinTheta = std::sqrt (1.f - cosTheta * cosTheta),
      x = std::cos (phi) * sinTheta,
      y = std::sin (phi) * sinTheta,
      z = cosTheta;
    return {x,y,z};
  }

  Vector hemispherical (RandBits& rng) {
    Vector const v = spherical (rng);
    return { v.x, v.y, std::abs (v.z) };
  }

  Vector disc (RandBits& rng) {
    float const
      rr = std::sqrt (canon (rng)),
      theta = 2.f * pi * canon (rng),
      x = rr * std::cos (theta),
      y = rr * std::sin (theta);
    return {x,y,0};
  }

  Vector conical (RandBits& rng, float halfApex) {
    float const
      lambda = 2.f * pi * canon (rng),
      theta = halfApex * std::sqrt (canon (rng)),
      sinTheta = std::sin (theta),
      x = std::cos (lambda) * sinTheta,
      y = std::sin (lambda) * sinTheta,
      z = std::cos (theta);
    return {x,y,z};
  }
}

