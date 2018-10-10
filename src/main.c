#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#if USE_SDL
#include "SDL.h"
#endif
#include "apis.h"

#define SCRN_HEIGHT (720)
#define SCRN_WIDTH (SCRN_HEIGHT*2)

#define MAKE_U32_COLOR(r, g, b, a) (uint32_t)(((r&0xFF)  << 0) | ((g&0xFF) << 8) | ((b&0xFF) << 16) | ((a&0xFF) << 24))

typedef struct sphere_t {
  vec3_t centre;
  float radius;
} sphere_t;

typedef struct hit_rec_t {
	float t;
	vec3_t p, normal;
} hit_rec_t;

typedef struct ray_t {
  vec3_t pt, dir;
} ray_t;

typedef struct camera_t {
  vec3_t lower_left_corner;
  vec3_t horizontal;
  vec3_t vertical;
  vec3_t origin;
} camera_t;

void cam_default(camera_t* cam) {
  vmathV3MakeFromElems(&cam->lower_left_corner, -2.f, -1.f, -1.f);
  vmathV3MakeFromElems(&cam->horizontal, 4.f, 0.f, 0.f);
  vmathV3MakeFromElems(&cam->vertical, 0.f, 2.f, 0.f);
  vmathV3MakeFromElems(&cam->origin, 0.f, 0.f, 0.f);
}

void cam_ray(ray_t* screen_ray, camera_t const* cam, float u, float v) {
  vec3_t vv;
  vmathV3ScalarMul(&screen_ray->dir, &cam->horizontal, u);
  vmathV3ScalarMul(&vv, &cam->vertical, v);
  vmathV3Add(&screen_ray->dir, &screen_ray->dir, &vv);
  vmathV3Add(&screen_ray->dir, &screen_ray->dir, &cam->lower_left_corner);
  screen_ray->pt = cam->origin;
}

enum { sphere_count = 2 };
static sphere_t sphere[sphere_count] = {
	{.centre = { 0.f, 0.f, -1.f },.radius = .5f, },
	{.centre = { 0.f, -100.5f, -1.f },.radius = 100.f, },
};

vec3_t ray_point_at(ray_t const* ray, float t) {
  vec3_t r, vt;
  vmathV3ScalarMul(&vt, &ray->dir, t);
  vmathV3Add(&r, &ray->pt, &vt);
  return r;
}

vec3_t rand_vec3() {
	vec3_t p;
	do {
		p = vmathV3Sub_V( vmathV3MakeFromElems_V(2.f*frand48(), 2.f*frand48(), 2.f*frand48()), vmathV3MakeFromScalar_V(1));
	} while (vmathV3Dot(&p, &p) >= 1.f);
	return p;
}

int hit_spheres(hit_rec_t* hit, sphere_t const* s, uint32_t s_count, float t_min, float t_max, ray_t const* r) {
	hit->t = t_max;
	for (uint32_t i = 0; i < s_count; ++i) {
		vec3_t oc;
		vmathV3Sub(&oc, &r->pt, &s[i].centre);
		float a = vmathV3Dot(&r->dir, &r->dir);
		float b = vmathV3Dot(&oc, &r->dir);
		float c = vmathV3Dot(&oc, &oc) - s[i].radius*s[i].radius;
		float discrimiant = b * b - a * c;
		if (discrimiant > 0) {
			float tmp = (-b - sqrtf(discrimiant)) / a;
			if (tmp < hit->t && tmp < t_max && tmp > t_min) {
				hit->t = tmp;
				hit->p = ray_point_at(r, tmp);
				hit->normal = vmathV3ScalarMul_V(vmathV3Sub_V(hit->p, s[i].centre), 1.f / s[i].radius);
			}
			tmp = (-b + sqrtf(discrimiant)) / a;
			if (tmp < hit->t && tmp < t_max && tmp > t_min) {
				hit->t = tmp;
				hit->p = ray_point_at(r, tmp);
				hit->normal = vmathV3ScalarMul_V(vmathV3Sub_V(hit->p, s[i].centre), 1.f / s[i].radius);
			}
		}
	}
	return hit->t < t_max;
}

vec3_t colour(ray_t const* ray) {
  static const vec3_t k1 = {1.f, 1.f, 1.f};
  static const vec3_t k2 = {.5f, .7f, 1.f};
  hit_rec_t hit;
  if (hit_spheres(&hit, sphere, sphere_count, 0.f, FLT_MAX, ray)) {
		vec3_t target = vmathV3Add_V(vmathV3Add_V(hit.p, hit.normal), rand_vec3());
		ray_t bounce = { .pt = hit.p,.dir = vmathV3Sub_V(target, hit.p) };
		return vmathV3ScalarMul_V(colour(&bounce), .5f);
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
    uint32_t samples = 100;
    camera_t cam;
    cam_default(&cam);
    for (int32_t j=0; j < SCRN_HEIGHT; ++j) {
      for (int32_t i=0; i < SCRN_WIDTH; ++i) {
        vec3_t col = {0.f, 0.f, 0.f};
        for (uint32_t s=0; s < samples; ++s) {
          ray_t screen_ray;
          float u = ((float)i + (float)drand48()) / (float)SCRN_WIDTH;
          float v = ((float)j + (float)drand48()) / (float)SCRN_HEIGHT;
          cam_ray(&screen_ray, &cam, u, v);
          col = vmathV3Add_V(col, colour(&screen_ray));

        }
        col = vmathV3ScalarDiv_V(col, (float)samples);
        int ir = (int)floorf(255*col.x);
        int ig = (int)floorf(255*col.y);
        int ib = (int)floorf(255*col.z);

				*(uint32_t*)(((uint8_t*)dest_main_buffer) + (dest_main_bfr_pitch * (SCRN_HEIGHT-(j+1))) + (i * 4)) = MAKE_U32_COLOR(ir, ig, ib, 255);
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

