#pragma once

#ifdef WIN32
  #define drand48() ((double)rand() / RAND_MAX)
	#define frand48() ((float)drand48())
#endif

#include "vectormath/vectormath_aos.h"
#include "vectormath/vectormath_aos_v.h"
#include <stdint.h>
#include <stdio.h>

typedef VmathVector3 vec3_t;
typedef VmathVector4 vec4_t;
typedef VmathMatrix3 mat3_t;
typedef VmathMatrix4 mat4_t;

