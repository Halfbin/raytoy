
#pragma once

#include "linear.hpp"
#include "material.hpp"
#include "colour.hpp"
#include "math.hpp"
#include "photon.hpp"

#include <cstdio>

namespace RT {
  struct PathNode {
    Point position;
    Vector normal;
    Colour matDiffuse;
    short px, py;
    Colour k;
    float radius;
    int nPhotons : 31;
    bool enabled : 1;
    Colour tau;

    PathNode
      ( Point p, Vector n
      , Material const* m
      , int px, int py, Colour k
      , float r
      )
      : position (p), normal (n)
      , matDiffuse (m->kDiffuse)
      , px (px), py (py), k (k)
      , radius (r), nPhotons (0), enabled (true), tau (black)
      { }

    PathNode (PathNode const&) = default;
    PathNode& operator = (PathNode const&) = default;
    PathNode (PathNode&&) = default;
    PathNode& operator = (PathNode&&) = default;

    Colour radiance () const
      { return k * tau / (pi * radius * radius); }

    void update (PhotonMap const&, float alpha);

    explicit operator bool () const { return enabled; };
  };
}

