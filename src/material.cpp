
#include "material.hpp"

#include <cassert>
#include <algorithm>

namespace RT {
  Interaction Material::interact (RandBits& rng, Vector, Colour ci) const {
    float const
      pDiffuse  = max (kDiffuse *ci) / max (ci),
      pSpecular = max (kSpecular*ci) / max (ci),
      xi = canon (rng);

    Interaction::Type ty = Interaction::Type::absorbed;
    if (xi < pDiffuse)
      ty = Interaction::Type::diffuse;
    else if (xi - pDiffuse < pSpecular)
      ty = Interaction::Type::specular;

    return { ty, correctRefl (ty, ci) };
  }

  Colour Material::correctRefl (Interaction::Type ty, Colour ci) const {
    float const
      pDiffuse  = max (kDiffuse *ci) / max (ci),
      pSpecular = max (kSpecular*ci) / max (ci);

    switch (ty) {
      case Interaction::Type::diffuse: return kDiffuse / pDiffuse;
      case Interaction::Type::specular: return kSpecular / pSpecular;
      case Interaction::Type::absorbed:;
    }

    return {0,0,0};
  }

  Colour Material::brdfDiffuse (Vector, Vector in) const {
    return kDiffuse * std::max (0.f, -in.z);
  }

  Vector Material::bounce
    ( RandBits& rng
    , Interaction::Type t
    , Vector incident
    ) const
  {
    switch (t) {
      case Interaction::Type::diffuse: return bounceDiffuse (rng, incident);
      case Interaction::Type::specular: return bounceSpecular (rng, incident);
      case Interaction::Type::absorbed:;
    }

    return {0,0,0};
  }

  Vector Material::bounceDiffuse (RandBits& rng, Vector) const
    { return lambert (rng); }

  Vector Material::bounceSpecular (RandBits&, Vector incident) const
    { return incident * Vector{1,1,-1}; }

  float fresnel (float R0, float cosine)
    { return R0 + (1.f - R0) * std::pow (1 - cosine, 5.f); }

  Vector Dielectric::bounceSpecular (RandBits& rng, Vector incident) const {
    float const
      c = std::abs (incident.z),
      dir = (incident.z < 0.f)? 1.f : -1.f,
      sqrtR0 = dir * (1.f - eta) / (1.f + eta),
      R0 = sqrtR0 * sqrtR0,
      pReflect = fresnel (R0, c),
      xi1 = canon (rng);

    Vector const reflection = incident*Vector{1,1,-1};
    if (xi1 < pReflect)
      return reflection;

    float const
      r = (incident.z > 0.f)? eta : (1.f / eta),
      discrim = 1.f - r*r * (1.f - c*c);
    if (discrim < 0.f)
      return reflection;

    Vector const
      scaled = r * incident,
      adjust = { 0, 0, r*c - std::sqrt (discrim) },
      transmit = scaled + dir * adjust;
    assert (dot (transmit, incident) > 0.f);
    return unit (transmit);
  }
}

