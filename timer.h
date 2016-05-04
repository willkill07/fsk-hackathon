#include <time.h>

typedef struct timespec timer;

float usTime (struct timespec start, struct timespec end) {
  return (end.tv_sec - start.tv_sec) * 1000000.0f + (end.tv_nsec - start.tv_nsec) / 1000.0f;
}

struct timespec getTime() {
	struct timespec t;
	clock_gettime (CLOCK_MONOTONIC, &t);
	return t;
}
