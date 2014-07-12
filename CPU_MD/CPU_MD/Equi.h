/*
 * Molecular Dynamics Project.
 * Outhor: Vladimir Novikov.
 */

#ifndef EQUI_H_
#define EQUI_H_

#include "Configurations.h"
#include "Structures.h"
#include "Constants.h"
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

const double rel_error= 1E-12;
const double Cpi      = 3.1415926535897932384626434;
const double C2sqrtPi = 1.1283791670955125738961589;
const double CsqrtPi  = 1.7724538509055160272981675;
const double Csqrt2   = 1.4142135623730950488016887;
const double C1sqrt2  = 0.7071067811865475244008444;

int writeParticlesInput(configurations*);
double erf(double);
double erfc(double);
double Phi (double);
double InvErfSmall (const double);
double lambertW2 (const double);
double InvErfcSmall (const double);
double inverf (const double);
double inverfc (const double);
void chooseVelocities(particleStruct*, int, double, configurations*, double);

#endif //EQUI_H_