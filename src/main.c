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

typedef
enum materialtype_t {
  MatLambertian,
  MatMetal,
	MatDielectric,

  MatMax,
} materialtype_t;

typedef struct matlambertian_t {
  vec3_t albedo;
} matlambertian_t;

typedef struct matmetal_t {
  vec3_t albedo;
  float roughness;
} matmetal_t;

typedef struct matdielectric_t {
	float refraction;
} matdielectric_t;

typedef struct sphere_t {
  vec3_t centre;
  float radius;
  materialtype_t mat;
  uint16_t matidx;
} sphere_t;

typedef struct hit_rec_t {
	float t;
	vec3_t p, normal;
  materialtype_t mat;
  uint16_t midx;
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

float schlick(float cosine, float ref_idx) {
  float r0 = (1.f-ref_idx)/(1.f+ref_idx);
  r0 = r0 * r0;
  return r0 + (1.f-r0)*powf(1.f-cosine, 5);
}

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

vec3_t vec_reflect(vec3_t const* v, vec3_t const* n) {
  vec3_t r = vmathV3ScalarMul_V(*n, 2.f*vmathV3Dot(v, n));
  return vmathV3Sub_V(*v, r);
}

float vec_refract(vec3_t const* v, vec3_t const* n, float idx_refraction, vec3_t* r) {
	vec3_t nv = vmathV3Normalize_V(*v);
  vec3_t nn = vmathV3Normalize_V(*n);
	float dt = vmathV3Dot(&nv, &nn);
	float discriminant = 1.f - idx_refraction*idx_refraction*(1-dt*dt);
	if (discriminant > 0.f) {
		vec3_t v1, v2;
    vmathV3ScalarMul(&v1, &nv, idx_refraction);
    float v2f = idx_refraction * dt + sqrtf(discriminant);
    vmathV3ScalarMul(&v2, &nn, v2f);
    vmathV3Sub(r, &v1, &v2);
	}
	return discriminant;
}

int lambertian_scatter(matlambertian_t const* params, ray_t const* rin, hit_rec_t const* hit, vec3_t* atten, ray_t* scatter) {
  vec3_t target = vmathV3Add_V(vmathV3Add_V(hit->p, hit->normal), rand_vec3());
  scatter->pt = hit->p;
  scatter->dir = vmathV3Sub_V(target, hit->p);
  *atten = params->albedo;
  return 1;
}

int metal_scatter(matmetal_t const* params, ray_t const* rin, hit_rec_t const* hit, vec3_t* atten, ray_t* scatter) {
  vec3_t n = vmathV3Normalize_V(rin->dir);
  vec3_t reflected = vmathV3Add_V(vec_reflect(&n, &hit->normal), vmathV3ScalarMul_V(rand_vec3(), params->roughness));
  scatter->pt = hit->p;
  scatter->dir = reflected;
  *atten = params->albedo;
  return vmathV3Dot_V(scatter->dir, hit->normal) > 0.f ? 1 : 0;
}

int dielectric_scatter_v1(matdielectric_t const* params, ray_t const* rin, hit_rec_t const* hit, vec3_t* atten, ray_t* scatter) {
  vec3_t outward_nrm;
  vec3_t reflected = vec_reflect(&rin->dir, &hit->normal);
  vec3_t refracted;
  float ni_over_ni;
  vmathV3MakeFromScalar(atten, 1.f);
  //float reflect_prob;
  //float cosine;
  if (vmathV3Dot(&rin->dir, &hit->normal) > 0) {
    vmathV3Neg(&outward_nrm, &hit->normal);
    ni_over_ni = params->refraction;
    //cosine = params->refraction * vmathV3Dot(&rin->dir, &hit->normal) / vmathV3Length(&rin->dir);
  } else {
    outward_nrm = hit->normal;
    ni_over_ni = 1.f / params->refraction;
    //cosine = -vmathV3Dot(&rin->dir, &hit->normal) / vmathV3Length(&rin->dir);
  }
  if (vec_refract(&rin->dir, &outward_nrm, ni_over_ni, &refracted) > 0) {
    scatter->pt = hit->p;
    scatter->dir = refracted;
  } else {
    scatter->pt = hit->p;
    scatter->dir = reflected;
  }

  return 1;
}

int dielectric_scatter_v2(matdielectric_t const* params, ray_t const* rin, hit_rec_t const* hit, vec3_t* atten, ray_t* scatter) {
	vec3_t outward_nrm;
  vec3_t reflected = vec_reflect(&rin->dir, &hit->normal);
  vec3_t refracted;
  float ni_over_ni;
  float reflect_prob;
  float cosine;
  vmathV3MakeFromScalar(atten, 1.f);
  if (vmathV3Dot(&rin->dir, &hit->normal) > 0) {
    vmathV3Neg(&outward_nrm, &hit->normal);
    ni_over_ni = params->refraction;
    cosine = params->refraction * vmathV3Dot(&rin->dir, &hit->normal) / vmathV3Length(&rin->dir);
  } else {
    outward_nrm = hit->normal;
    ni_over_ni = 1.f / params->refraction;
    cosine = -vmathV3Dot(&rin->dir, &hit->normal) / vmathV3Length(&rin->dir);
  }
  if (vec_refract(&rin->dir, &outward_nrm, ni_over_ni, &refracted) > 0) {
    reflect_prob = schlick(cosine, params->refraction);
  } else {
    reflect_prob = 1.f;
  }

  if (drand48() < reflect_prob) {
    scatter->pt = hit->p;
    scatter->dir = reflected;
  } else {
    scatter->pt = hit->p;
    scatter->dir = refracted;
  }
	return 1;
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
        hit->mat = s[i].mat;
        hit->midx = s[i].matidx;
			}
			tmp = (-b + sqrtf(discrimiant)) / a;
			if (tmp < hit->t && tmp < t_max && tmp > t_min) {
				hit->t = tmp;
				hit->p = ray_point_at(r, tmp);
				hit->normal = vmathV3ScalarMul_V(vmathV3Sub_V(hit->p, s[i].centre), 1.f / s[i].radius);
        hit->mat = s[i].mat;
        hit->midx = s[i].matidx;
			}
		}
	}
	return hit->t < t_max;
}

enum { 
  sphere_count = 5,
  lambertian_count = 2,
  metal_count = 3,
  dielectric_count = 1,
  call_depth_limit = 50,
};
static sphere_t sphere[sphere_count] = {
	{.centre = { 0.f, 0.f, -1.f },.radius = .5f, .mat=MatLambertian, .matidx=0},
	{.centre = { 0.f, -100.5f, -1.f },.radius = 100.f, .mat=MatLambertian,.matidx = 1},
  {.centre = { 1.f, 0.f, -1.f },.radius = .5f, .mat=MatMetal, .matidx=0},
  {.centre = { -1.f, 0.f, -1.f },.radius = .5f, .mat=MatDielectric, .matidx=0},
{ .centre = { -1.f, 0.f, -1.f },.radius = -.45f,.mat = MatDielectric,.matidx = 0 },
};
matlambertian_t mat_lamb[lambertian_count] = {
  {.albedo={.8f, .3f, .3f}},
  {.albedo={.8f, .8f, .0f}},
};
matmetal_t mat_metal[metal_count] = {
  {.albedo={.8f, .6f, .2f}, .roughness=.0f, },
  {.albedo={.8f, .6f, .2f}, .roughness=.4f, },
  {.albedo={.8f, .8f, .8f},.roughness =1.f },
};
matdielectric_t mat_dielectric[dielectric_count] = {
  {.refraction=1.6f}
};

vec3_t colour(ray_t const* ray, int call_depth) {
  static const vec3_t k1 = {1.f, 1.f, 1.f};
  static const vec3_t k2 = {.5f, .7f, 1.f};
  hit_rec_t hit;
  if (hit_spheres(&hit, sphere, sphere_count, 0.001f, FLT_MAX, ray)) {
    if (call_depth >= call_depth_limit) return vmathV3MakeFromScalar_V(0.f);
    ray_t scattered;
    vec3_t atten;
    int ret = 0;
    switch(hit.mat) {
    case MatLambertian:
      ret = lambertian_scatter(mat_lamb + hit.midx, ray, &hit, &atten, &scattered);
      break;
    case MatMetal:
      ret = metal_scatter(mat_metal + hit.midx, ray, &hit, &atten, &scattered);
      break;
    case MatDielectric:
      ret = dielectric_scatter_v2(mat_dielectric+hit.midx, ray, &hit, &atten, &scattered);
      break;
    }
		
    if (!ret) return vmathV3MakeFromScalar_V(0.f);
    return vmathV3MulPerElem_V(atten, colour(&scattered, call_depth + 1));
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
          col = vmathV3Add_V(col, colour(&screen_ray, 0));

        }
        col = vmathV3ScalarDiv_V(col, (float)samples);
				// Gamma correct the image.
				col.x = powf(col.x, 1.f / 2.2f);
				col.y = powf(col.y, 1.f / 2.2f);
				col.z = powf(col.z, 1.f / 2.2f);
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

