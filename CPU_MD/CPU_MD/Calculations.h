/*
 *
 */

#ifndef CALCULATIONS_H_
#define CALCULATIONS_H_

#include "Structures.h"
#include "Constants.h"
#include "Configurations.h"
#include "LennardJones.h"
#include "SiPotential.h"
#include <fstream>

//---------- force between two Si particles ---------------//
double distance2(real3, real3);//r(i,j)

//-------------- calculate forces between Si --------------//
void calculateForce_Si(int, int, particleStruct*, particleStruct*, int, int, bool, bool);
//--------------------------------------------------------//

void calculateForce_Xe(int, int, particleStruct*, particleStruct*, int, int, bool);

//---------calculate total potential of Si ---------------//
double V_total_Si(int, particleStruct*);
double V_total_Si_Xe(int, particleStruct*, int, particleStruct*);
double V_total_Xe(int, particleStruct*);
void V_total(int, particleStruct*, int, particleStruct*);
//--------------------------------------------------------//

void initiateAcceleration(particleStruct*, int, double);

void predict(particleStruct*, int, float);

void correct(particleStruct*, double, int, double);

#endif	//CALCULATIONS_H_