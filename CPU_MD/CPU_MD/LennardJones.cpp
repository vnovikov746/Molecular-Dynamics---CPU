#include "LennardJones.h"

//------------------------Lennard Jones Potential -----------------------------//
double lennardJonesForce(double dist, double sig, double eps)
{
	double sigsq = sig*sig;
	double con = 24.0*eps/sigsq;
	double dist2 = dist * dist;
	dist2 /= sigsq;
	
	double dist4 = dist2*dist2;
	double dist8 = dist4*dist4;
	double dist14 = dist2*dist4*dist8;
	double invdist8= 1.0/dist8;
	double invdist14= 1.0/dist14;
	double s = 2.0*invdist14-invdist8;
	return s * con;
}
//----------------------------------------------------------------------------//

double lennardJonesPotential(double dist, double sig, double eps)
{
	double expr = sig/dist;
	double expr2 = expr*expr;
	double expr4 = expr2*expr2;
	double expr6 = expr4*expr2;
	double expr12 = expr6*expr6;

	return 4.0*eps*(expr12-expr6);
}

