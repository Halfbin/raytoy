
#include "directrt.hpp"

namespace RT {
  Colour directTermRT
    ( int branch
    , Scene const& scene
    , RandBits& rng
    , RayHit<Item const*> hit
    )
  {
    auto const* mat = hit.occ->material;

    Colour total = black;
    if (mat->kDiffuse != black) {
      for (int i = 0; i != branch; i++) {
        Lum const& lum = scene.randomLum (rng);
        Vector const lumDir = lum.randomShadow (rng, hit.position);

        float const geomK = dot (hit.normal, lumDir);
        if (geomK <= 0.f)
          continue;

        Ray shadowRay (hit.position, lumDir);
        auto const occHit = scene.traceRay (shadowRay);
        if (!occHit || occHit.occ != &lum)
          continue;

        float const coverageK = lum.coverageK (hit.position);
        total += coverageK * geomK * lum.exitance ();
      }
    }

    return mat->kDiffuse * (total / (float) branch)
         + mat->emission / (2.f * pi);
  }
}

