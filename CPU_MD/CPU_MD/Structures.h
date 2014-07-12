/*
 * Molecular Dynamics Project.
 * Outhor: Vladimir Novikov.
 */

#ifndef PARTICLE_STRACT_H_
#define PARTICLE_STRACT_H_
#include <iostream>
using namespace std;

typedef double real;

typedef unsigned char byte;

struct cell
{
	int i;
	int j;
	int k;
};

struct real3
{
	double x;
	double y;
	double z;
};

typedef struct particleStruct
{
	byte type;
	real3 velocity;
	real3 force;
	real3 aAcc;
	real3 bAcc;
	real3 cAcc;
	real3 position;
	cell cellPos;
 	int *siNeighbors;
	int *xeNeighbors;
}particleStruct;

struct allLists
{
	particleStruct *siParticles;
	particleStruct *xeParticles;
};

#endif