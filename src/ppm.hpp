
#pragma once

#include "linear.hpp"
#include "material.hpp"
#include "colour.hpp"
#include "math.hpp"
#include "photon.hpp"

#include <cstdio>

namespace RT {
  struct Stats {
    float r2;
    int nPhotons, dn;
    Colour tau;

    explicit Stats (float r2)
      : r2 (r2), nPhotons (0), dn (0), tau (black)
      { }

    Colour radiance () const {
      constexpr float const twoPi2 = 2.f * pi * pi;
      return tau / (twoPi2 * r2);
    }

    void add (Colour, int dn);
    void update (float alpha);
  };

  struct Gather {
    int px, py;
    Colour contrib;
    int nPhotons;

    Gather () = default;
    Gather (int px, int py)
      : px (px), py (py)
      , contrib (black), nPhotons (0)
      { }
    Gather (Gather const&) = default;
    Gather& operator = (Gather const&) = default;
  };

  struct PathNode {
    int px, py;
    Point position;
    Vector normal;
    Colour matDiffuse;
    Colour k;

    PathNode
      ( int px, int py
      , Point p, Vector n
      , Material const* m
      , Colour k
      )
      : px (px), py (py)
      , position (p), normal (n)
      , matDiffuse (m->kDiffuse)
      , k (k)
      { }

    PathNode (PathNode const&) = default;
    PathNode& operator = (PathNode const&) = default;
    PathNode (PathNode&&) = default;
    PathNode& operator = (PathNode&&) = default;

    Gather gather (float r2, PhotonMap const&) const;
  };
}

