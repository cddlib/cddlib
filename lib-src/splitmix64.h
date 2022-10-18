/* 
 * A fast portable random number generator.
 *
 * Adapted version of splitmix64.c originally published by Sebastiano
 * Vigna to the public domain.
 * http://xoroshiro.di.unimi.it/splitmix64.c
*/

#ifndef __SPLITMIX64_H
#define __SPLITMIX64_H

#include <stdint.h>

static _Thread_local uint64_t x; /* The state can be seeded with any value. */

static void srand_splitmix64(uint64_t seed) {
	x = seed;
}

static uint64_t rand_splitmix64() {
	uint64_t z = (x += 0x9e3779b97f4a7c15);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
	return z ^ (z >> 31);
}

#endif // __SPLITMIX64_H
