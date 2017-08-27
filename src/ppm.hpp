
#pragma once

#include "linear.hpp"
#include "material.hpp"
#include "colour.hpp"
#include "math.hpp"
#include "photon.hpp"

namespace RT {
  struct PathNode {
    Point position;
    Vector normal, incident;
    Material const* material;
    int px, py;
    Colour k;
    float radius2;
    int nPhotons;
    Colour tau;

    PathNode
      ( Point p, Vector n, Vector i
      , Material const* m
      , int px, int py, Colour k
      , float r2
      )
      : position (p), normal (n), incident (i)
      , material (m)
      , px (px), py (py), k (k)
      , radius2 (r2), nPhotons (0), tau (black)
      { }

    PathNode (PathNode const&) = default;
    PathNode& operator = (PathNode const&) = default;
    PathNode (PathNode&&) = default;
    PathNode& operator = (PathNode&&) = default;

    Colour radiance () const
      { return k * tau / (pi * radius2); }

    void update (PhotonMap const&, float alpha);
  };
}

