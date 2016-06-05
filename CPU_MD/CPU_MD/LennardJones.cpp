#include "LennardJones.h"

//------------------------Lennard Jones Potential -----------------------------//
double lennardJonesForce(double dist, double sig, double eps)
{
	double dist2 = dist*dist;
	double dist3 = dist2*dist;
	double dist6 = dist3*dist3;
	double dist12 = dist6*dist6;
	double dist13 = dist12*dist;
	double sig2 = sig*sig;
	double sig3 = sig2*sig;
	double sig6 = sig3*sig3;
	return ((24 * eps * (dist6 - (2 * sig6))) / dist13);
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

