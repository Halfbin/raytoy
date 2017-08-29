
#include "photonmapping.hpp"

#include "scene.hpp"

namespace RT {
/*Colour indirectTermPM
    ( NearSet& nears
    , PhotonMap const& map
    , RayHit<const Item*> hit
    )
  {
    map.nearest (nears, hit.position, hit.normal);
    float const radius = std::sqrt (nears.maxQd ());
    //fprintf (stderr, "radius %f\n", radius);

    Colour power = black;
    for (auto const& elem : nears) {
      float const
        cosine = std::max (0.f, -dot (elem.photon.incoming (), hit.normal)),
        filter = 1.f;//3.f * (1.f - std::sqrt (elem.qd) / radius);
      power += cosine * filter * elem.photon.power ();
    }

    float const
      angle = 2.f * pi,
      area = pi * nears.maxQd ();
    return hit.occ->material->kDiffuse * power / (angle * area);
  }*/

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

