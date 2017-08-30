
#include "scene.hpp"

#include "math.hpp"

namespace RT {
  Lum const* Item::isLum () const
    { return nullptr; }
  Lum const* Lum::isLum () const
    { return this; }


  void Scene::add (Item const* item)
    { items.push_back (item); }
  void Scene::add (Lum const* lum)
    { lums.push_back (lum); }

  Lum const& Scene::randomLum (RandBits& rng) const {
    float const xi = canon (rng);
    int const
      n = lums.size (),
      i = std::min (n-1, int (xi*n));
    return *lums[i];
  }

  Vector Lum::randomShadow (RandBits& rng, Point from) const {
    Vector const towards = from - surface->boundingCentre ();
    float const
      r = surface->boundingRadius (),
      d = norm (towards),
      halfApex = std::asin (r / d);
    Basis const facing = Basis::fromK (-unit (towards));
    Vector v = facing.outOf (conical (rng, halfApex));
    //fprintf (stderr, "(%f %f %f)\n", v.x, v.y, v.z);
    return v;
  }

  float Lum::coverageK (Point about) const {
    float const
      r = surface->boundingRadius (),
      d = norm (about - surface->boundingCentre ()),
      cosTheta = d / std::sqrt (r*r + d*d);
    return 1.f - cosTheta;
  }

  RayHit<Item const*> Scene::traceRay (Ray ray) const {
    RayHit<Item const*> hit = { };

    for (Item const* candidate : items) {
      auto const tryHit = candidate->surface->intersect (ray);
      if (tryHit.t >= 0.001f && tryHit.t < hit.t)
        hit = { tryHit.position, tryHit.t, tryHit.normal, candidate };
    }

    return hit;
  }
}

