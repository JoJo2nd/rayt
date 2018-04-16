#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include "SDL.h"
#include "vectormath/vectormath_aos.h"
#include "apis.h"

#define SCRN_WIDTH (960)
#define SCRN_HEIGHT (720)

int main(int argc, char** argv) {
  if (SDL_Init(SDL_INIT_VIDEO)!=0) {
    fprintf(stderr, "Error calling SDL_Init");        
  }

  SDL_Window* window = SDL_CreateWindow("Ray Tracer", 0, 0, SCRN_WIDTH, SCRN_HEIGHT, SDL_WINDOW_SHOWN);    
  SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, 0);
  SDL_Texture* main_surface = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, SCRN_WIDTH, SCRN_HEIGHT);  
  SDL_Event evt;

  for (;;) {
    while (SDL_PollEvent(&evt)) {
      if (evt.type == SDL_QUIT)
        return 0;
    }

    SDL_RenderClear(renderer);
    // void* dest_main_buffer;
    // int dest_main_bfr_pitch;
    // SDL_LockTexture(main_surface, NULL, (void**)&dest_main_buffer, &dest_main_bfr_pitch);
    // SDL_UnlockTexture(main_surface);
    
    SDL_RenderCopy(renderer, main_surface, NULL, NULL);
    SDL_RenderPresent(renderer);
  }
  return 0;
}

