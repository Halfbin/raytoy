
#include "monitor.hpp"

#define SDL_MAIN_HANDLED
#include <SDL2/SDL.h>

#include <stdexcept>

namespace RT {
  struct Monitor::Impl {
    SDL_Window* win;
  };

  Monitor::Monitor (int w, int h)
    : impl (new Impl)
  {
    SDL_SetMainReady ();

    SDL_SetHintWithPriority (SDL_HINT_NO_SIGNAL_HANDLERS, "1", SDL_HINT_OVERRIDE);

    int fail = SDL_Init (SDL_INIT_VIDEO);
    if (fail)
      throw std::runtime_error ("sdl fail");

    impl->win = SDL_CreateWindow
      ( "render"
      , SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED
      , w, h
      , SDL_WINDOW_RESIZABLE
      );
    if (!impl->win)
      throw std::runtime_error ("SDL_CreateWindow failed");
  }

  Monitor::~Monitor () {
    SDL_Quit ();
    delete impl;
  }

  void Monitor::update (Image<uint8_t> const& image) {
    SDL_PumpEvents ();

    int w, h;
    SDL_GetWindowSize (impl->win, &w, &h);
    float
      aImage = (float) image.width () / image.height (),
      aWin   = (float) w              / h;

    //fprintf (stderr, "%ix%i\n", w, h);

    SDL_Rect dest { 0, 0, w, h };
    if (aWin < aImage) {
      dest.h = (int) (dest.h * (aWin / aImage));
      dest.y = (h - dest.h) / 2;
    }
    else {
      dest.w = (int) (dest.w * (aImage / aWin));
      dest.x = (w - dest.w) / 2;
    }

    auto src = SDL_CreateRGBSurfaceWithFormatFrom
      ( const_cast<Pixel<uint8_t>*> (image.raw ())
      , image.width (), image.height ()
      , 24
      , image.width () * 3
      , SDL_PIXELFORMAT_RGB24
      );
    if (!src)
      throw std::runtime_error ("bluh");

    auto surf = SDL_GetWindowSurface (impl->win);
    SDL_FillRect (surf, nullptr, 0);
    SDL_BlitScaled (src, nullptr, surf, &dest);
    SDL_UpdateWindowSurface (impl->win);
    SDL_FreeSurface (src);
  }
}

