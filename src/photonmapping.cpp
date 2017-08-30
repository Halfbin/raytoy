
#include "photonmapping.hpp"

#include "scene.hpp"

namespace RT {
  Photon* tracePhoton
    ( RandBits& rng
    , Scene const& scene
    , Ray ray
    , Item const*
    , Colour power
    , int ttl
    , Photon* photons, Photon* end
    )
  {
    if (ttl == 0 || photons == end)
      return photons;

    RayHit<Item const*> const hit = scene.traceRay (ray);
    if (!hit)
      return photons;
    Item const* item = hit.occ;
    auto const* mat = item->material;

    /*fprintf (stderr, "hit at (%f %f %f)\n",
      hit.position.x, hit.position.y, hit.position.z);*/

    Basis const tang = Basis::fromK (hit.normal);
    Vector const incident = tang.into (ray.disp);

    Interaction const inter = mat->interact (rng, incident, power);

    bool const store
      =   inter.type != Interaction::Type::specular
      &&  mat->kDiffuse != black;
    if (store) {
      Photon const photon (hit.position, ray.disp, power);
      *photons++ = photon;
    }

    if (inter.type == Interaction::Type::absorbed)
      return photons;

    Vector const bounce = tang.outOf (mat->bounce (rng, inter.type, incident));
    Ray const newRay (hit.position, bounce);
    return tracePhoton
      ( rng
      , scene, newRay, item
      , power * inter.correctRefl
      , ttl-1
      , photons, end
      );
  }

  int castPhotons
    ( RandBits& rng
    , Scene const& scene
    , Photon* photons, Photon* end
    )
  {
    int nEmitted = 0;

    for (Photon* ptr = photons; ptr != end;) {
      Lum const& lum = scene.randomLum (rng);
      Ray const ray = lum.emit (rng);
      ptr = tracePhoton
        ( rng
        , scene, ray, &lum
        , lum.power
        , 100
        , ptr, end
        );
      nEmitted++;
    }

    return nEmitted;
  }
}

