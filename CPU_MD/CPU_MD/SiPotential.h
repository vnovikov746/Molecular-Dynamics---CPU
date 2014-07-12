
#ifndef SIPOTENTIAL_H_
#define SIPOTENTIAL_H_

#include "Structures.h"
#include "Constants.h"
#include "Configurations.h"
#include "Calculations.h"
#include "SiPotential.h"
#include <math.h>
#include <fstream>

//------------- potential between two Si ------------------// 
double f2(double);

double v2(double);
//--------------------------------------------------------//

double f2_derivative_of_rij_tag(double);

double v2_derivative_of_rix(real3, real3, double);
double v2_derivative_of_riy(real3, real3, double);
double v2_derivative_of_riz(real3, real3, double);
//--------------------------------------------------------//

//--------- force between three Si particles --------------//
double hi_derivative_of_rij_tag(double, double, double);
double hi_derivative_of_rik_tag(double, double, double);
double hj_derivative_of_rij_tag(double, double, double);
double hj_derivative_of_rik_tag(double, double, double);
double hk_derivative_of_rij_tag(double, double, double);
double hk_derivative_of_rik_tag(double, double, double);

double f3_derivative_of_rij_tag(double, double, double);
double f3_derivative_of_rik_tag(double, double, double);

double v3_derivative_of_rix(real3,real3,real3, double, double, double);
double v3_derivative_of_riy(real3,real3,real3, double, double, double);
double v3_derivative_of_riz(real3,real3,real3, double, double, double);
//--------------------------------------------------------//
//------------- potential between two Si ------------------// 
double hi(double, double, double);
double hj(double, double, double);
double hk(double, double, double);

double f3(double, double, double);

double v3(double, double, double);
//--------------------------------------------------------//

#endif //SIPOTENTIAL_H_