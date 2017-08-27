#pragma once

#include "surface.hpp"
#include "material.hpp"
#include "random.hpp"

namespace RT {
  struct Item;
  struct Lum;

  struct Camera {
    Point const eye, centre;
    Vector const vertical, horizontal;

    Camera (Point eye, Point centre, Vector up, float halfHeight)
      : eye (eye), centre (centre)
      , vertical (unit (up) * -halfHeight)
      , horizontal (unit (cross (centre - eye, up)) * halfHeight)
    { }

    Ray shoot (float nx, float ny) const {
      Point const origin = centre + horizontal * nx + vertical * ny;
      return { origin, origin - eye };
    }
  };

  struct Scene {
    std::vector<Item const*> items;
    std::vector<Lum const*>  lums;

    void add (Item const* item);
    void add (Lum const* lum);

    RayHit<const Item*> traceRay (Ray ray) const;

    Lum const& randomLum (RandBits&) const;
  };

  struct Item {
    Surface const* surface;
    Material const* material;

    Item (Scene& scene, Surface const* surface, Material const* material)
      : surface (surface), material (material)
    { scene.add (this); }

    Item (Item const&) = delete;
    Item& operator = (Item const&) = delete;
    Item (Item&&) = delete;
    Item& operator = (Item&&) = delete;

    virtual Lum const* isLum () const;
  };

  struct Lum : Item {
    BoundedSurface const* surface;
    Colour const power;
    Material const material;

    Lum (Scene& scene, BoundedSurface const* surface, Colour power)
      : Item (scene, surface, &material)
      , surface (surface)
      , power (power)
      , material ({0,0,0}, {0,0,0}, exitance ())
    { scene.add (this); }

    Lum const* isLum () const;

    Ray emit (RandBits& rng) const
      { return surface->randomEmission (rng); }

    Vector randomShadow (RandBits&, Point from) const;

    Colour exitance () const
      { return power / surface->area (); }

    float coverageK (Point about) const;
  };
}

