
#define USE_SDL (0)
#define _CRT_SECURE_NO_WARNINGS

#include <stdint.h>

typedef uint8_t bool_t;

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#undef STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#undef STB_IMAGE_WRITE_IMPLEMENTATION
#include "main.c"