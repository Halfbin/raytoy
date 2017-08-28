
#pragma once

#include "image.hpp"

namespace RT {
  class Monitor {
    struct Impl;
    Impl* impl;

  public:
    Monitor (int w, int h);
    ~Monitor ();

    Monitor (Monitor const&) = delete;

    void update (Image<uint8_t> const&);
  };
}

