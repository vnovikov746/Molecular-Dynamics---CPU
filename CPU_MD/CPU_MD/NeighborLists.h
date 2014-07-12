/*
 * Molecular Dynamics Project.
 * Outhor: Vladimir Novikov.
 */

#ifndef NEIGHBOR_LISTS_H_
#define NEOGHBOR_LISTS_H_

#include "Structures.h"
#include "Configurations.h"
#include "Constants.h"

#define H1 2971215073
#define H2 433494437
#define H3 1073807359

void divideToCells(particleStruct*, int, float);

void divideToCells(particleStruct*, particleStruct*, int, int, float);

void buildNeighbors(particleStruct*, particleStruct*, configurations*, float, float, float, float, byte, int, int);

int computeBucket(int, cell*, cell, bool, int, int, int);

#endif //NEIGHBOR_LISTS_H_