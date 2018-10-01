#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#if USE_SDL
#include "SDL.h"
#endif
#include "apis.h"

#define SCRN_WIDTH (1280)
#define SCRN_HEIGHT (720)

#define MAKE_U32_COLOR(r, g, b, a) (uint32_t)(((r&0xFF)  << 24) | ((g&0xFF) << 16) | ((b&0xFF) << 8) | ((a&0xFF) << 0))

static struct {
  vec3_t centre;
  float radius;
} sphere = {
  .centre= { 0.f, 0.f, -1.f },
  .radius = .5f,
};


struct ray_t {
  vec3_t pt, dir;
};
typedef struct ray_t ray_t;

vec3_t ray_point_at(ray_t* ray, float t) {
  vec3_t r, vt;
  vmathV3ScalarMul(&vt, &ray->dir, t);
  vmathV3Add(&r, &ray->pt, &vt);
  return r;
}

float hit_sphere(vec3_t* centre, float radius, ray_t const* r) {
  vec3_t oc;
  vmathV3Sub(&oc, &r->pt, centre);
  float a = vmathV3Dot(&r->dir, &r->dir);
  float b = 2.f * vmathV3Dot(&oc, &r->dir);
  float c = vmathV3Dot(&oc, &oc) - radius*radius;
  float discrimiant = b*b - 4.f*a*c;
  return discrimiant;
}

vec3_t colour(const ray_t* ray) {
  static const vec3_t k1 = {1.f, 1.f, 1.f};
  static const vec3_t k2 = {.5f, .7f, 1.f};
  if (hit_sphere(&sphere.centre, sphere.radius, ray) > 0.f) {
    return vmathV3MakeFromElems_V(1.f, 0.f, 0.f);
  }
  vec3_t unit_dir, a, b, ret;
  vmathV3Normalize(&unit_dir, &ray->dir);
  float t = .5f * (unit_dir.y + 1.f);
  vmathV3ScalarMul(&a, &k1, 1.0f-t);
  vmathV3ScalarMul(&b, &k2, t);
  vmathV3Add(&ret, &a, &b);
  return ret;  
}

int main(int argc, char** argv) {
#if USE_SDL
  if (SDL_Init(SDL_INIT_VIDEO)!=0) {
    fprintf(stderr, "Error calling SDL_Init");        
  }

  SDL_Window* window = SDL_CreateWindow("Ray Tracer", 0, 0, SCRN_WIDTH, SCRN_HEIGHT, SDL_WINDOW_SHOWN);    
  SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, 0);
  SDL_Texture* main_surface = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, SCRN_WIDTH, SCRN_HEIGHT);  
  SDL_Event evt;

#endif
  for (;;) {
#if USE_SDL
    while (SDL_PollEvent(&evt)) {
      if (evt.type == SDL_QUIT)
        return 0;
    }
#endif

    static int do_single_frame = 0;

    if (do_single_frame)
      break;

    do_single_frame = 1;
#if USE_SDL
    SDL_RenderClear(renderer);
    void* dest_main_buffer;
    int dest_main_bfr_pitch;
    SDL_LockTexture(main_surface, NULL, (void**)&dest_main_buffer, &dest_main_bfr_pitch);
#else
		void* dest_main_buffer = malloc(SCRN_WIDTH*SCRN_HEIGHT*4);
		int dest_main_bfr_pitch = SCRN_WIDTH*4;
#endif

    vec3_t lower_left_corner = {-2.f, -1.f, -1.f};
    vec3_t horizontal = {4.f, 0.f, 0.f};
    vec3_t vertical = {0.f, 2.f, 0.f};
    vec3_t origin = {0.f, 0.f, 0.f};
    for (int32_t j=SCRN_HEIGHT-1; j >= 0; --j) {
      for (int32_t i=0; i < SCRN_WIDTH; ++i) {
        float u = (float)i / (float)SCRN_WIDTH;
        float v = (float)j / (float)SCRN_HEIGHT;
        vec3_t vv;
        ray_t screen_ray;
        vmathV3ScalarMul(&screen_ray.dir, &horizontal, u);
        vmathV3ScalarMul(&vv, &vertical, v);
        vmathV3Add(&screen_ray.dir, &screen_ray.dir, &vv);
        vmathV3Add(&screen_ray.dir, &screen_ray.dir, &lower_left_corner);
        screen_ray.pt = origin;
        vec3_t col = colour(&screen_ray);
        int ir = (int)round(255*col.x);
        int ig = (int)round(255*col.y);
        int ib = (int)round(255*col.z);

        *(uint32_t*)(((uint8_t*)dest_main_buffer) + dest_main_bfr_pitch*j + i*4) = MAKE_U32_COLOR(ir, ig, ib, 255);
      }
    }
#if USE_SDL
    SDL_UnlockTexture(main_surface);
    
    SDL_RenderCopy(renderer, main_surface, NULL, NULL);
    SDL_RenderPresent(renderer);
#else
		stbi_write_png("rayt.png", SCRN_WIDTH, SCRN_HEIGHT, 4, dest_main_buffer, dest_main_bfr_pitch);
		free(dest_main_buffer);
#endif
  }
  return 0;
}

