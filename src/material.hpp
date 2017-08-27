#pragma once

#include "linear.hpp"
#include "colour.hpp"
#include "random.hpp"

namespace RT {
  struct Interaction {
    enum class Type { diffuse, specular, absorbed } type;
    Colour correctRefl;
  };

  struct Material {
    Colour const
      kDiffuse,
      kSpecular,
      emission;

    Material (Colour kd, Colour ks, Colour e = black)
      : kDiffuse(kd), kSpecular(ks), emission(e) { }

    Colour brdfDiffuse (Vector, Vector) const;
    Interaction interact (RandBits&, Vector, Colour ci) const;

    Vector bounce (RandBits&, Interaction::Type, Vector incident) const;

    virtual Vector bounceDiffuse (RandBits&, Vector incident) const;
    virtual Vector bounceSpecular (RandBits&, Vector incident) const;
  };

  struct Lambertian : public Material {
    Lambertian (Colour k) : Material (k, {0,0,0}) { }
  };

  struct Mirror : public Material {
    Mirror (Colour k) : Material ({0,0,0}, k) { }
  };

  struct Dielectric : public Material {
    float const eta;

    Dielectric (Colour k, float eta)
      : Material ({0,0,0}, k)
      , eta (eta)
      { }

    Vector bounceSpecular (RandBits&, Vector incident) const;
  };
}

