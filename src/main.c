#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include "SDL.h"
#include "vectormath/vectormath_aos.h"
#include "apis.h"

#define SCRN_WIDTH (960)
#define SCRN_HEIGHT (720)

#define MAKE_U32_COLOR(r, g, b, a) (uint32_t)(((r&0xFF)  << 24) | ((g&0xFF) << 16) | ((b&0xFF) << 8) | ((a&0xFF) << 0))

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

    static int do_single_frame = 0;

    if (do_single_frame)
      continue;

    do_single_frame = 1;
    SDL_RenderClear(renderer);
    void* dest_main_buffer;
    int dest_main_bfr_pitch;
    SDL_LockTexture(main_surface, NULL, (void**)&dest_main_buffer, &dest_main_bfr_pitch);

    for (int32_t j=SCRN_HEIGHT-1; j >= 0; --j) {
      for (int32_t i=0; i < SCRN_WIDTH; ++i) {
        float r = (float)i / (float)SCRN_WIDTH;
        float g = (float)j / (float)SCRN_HEIGHT;
        float b = .2f;
        int ir = round(255.99*r);
        int ig = round(255.99*g);
        int ib = round(255.99*b);

        *(uint32_t*)(((uint8_t*)dest_main_buffer) + dest_main_bfr_pitch*j + i*4) = 
        MAKE_U32_COLOR(ir, ig, ib, 255);
      }
    }

    SDL_UnlockTexture(main_surface);
    
    SDL_RenderCopy(renderer, main_surface, NULL, NULL);
    SDL_RenderPresent(renderer);
  }
  return 0;
}

