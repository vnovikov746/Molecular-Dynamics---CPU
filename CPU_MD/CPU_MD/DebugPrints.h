#ifndef DEBUGPRINTS_H_
#define DEBUGPRINTS_H_

#include "Structures.h"
#include "Constants.h"
#include "Calculations.h"
#include "LennardJones.h"
#include "SiPotential.h"
#include <iostream>
#include <fstream>

#define si3resolution 600

void printSi2Potential(string);
void printSi2ForceX(string);
void printSi2ForceY(string);
void printSi2ForceZ(string);
void printSi3PotentialCos(string);
void printSi3ForceCos(string);
void printSi3PotentialDist(string);
void printSi3ForceDist(string);
void printXePotential(string);
void printXeForce(string);

////////////////////////////// call the wnted functions ///////////////////////////
void printDebug(string fileName)
{
	printSi2Potential(fileName);
	printSi2ForceX(fileName);
	printSi2ForceY(fileName);
	printSi2ForceZ(fileName);
	printSi3PotentialCos(fileName); // NOT DONE
	printSi3ForceCos(fileName); // NOT DONE
	printSi3PotentialDist(fileName); // NOT DONE
	printSi3ForceDist(fileName); // NOT DONE
	printXePotential(fileName);
	printXeForce(fileName);
}

///////////////////////////////// si2 potential ///////////////////////////////////
void printSi2Potential(string fileName)
{
	real3 i, j;
	ofstream out;
	double dist;

	i.x = i.y = i.z = 0.0;
	j.x = j.y = j.z = 0.0;

	fileName.append("\\V_SI_2");
	out.open(fileName);
	out.precision(20);

	j.x = 0.5*sigma_Si;
	for(int k = 0; k < 300; k++)
	{
		dist = distance2(i,j)/sigma_Si;

		out<<v2(dist)/epsilon_Si;

		if(k < 299)
		{
			out<<", ";
		}
		j.x += 0.01*sigma_Si;
	}
	out<<endl;
	out.close();

	fileName.append("_X");
	out.open(fileName);
	out.precision(20);
	j.x = 0.5;
	for(int k = 0; k < 300; k++)
	{
		out<<distance2(i,j);
		if(k < 299)
		{
			out<<", ";
		}
		j.x += 0.01;
	}
	out.close();
}
///////////////////////////////// si2 force on x ///////////////////////////////////
void printSi2ForceX(string fileName)
{
	real3 i, j;
	ofstream out;
	double dist;

	i.x = i.y = i.z = 0.0;
	j.x = j.y = j.z = 0.0;

	fileName.append("\\F_SI_2_X");
	out.open(fileName);
	out.precision(20);

	j.x = 0.5*sigma_Si;
	for(int k = 0; k < 300; k++)
	{
		dist = distance2(i,j);
		if(dist/sigma_Si >= a_Si)
			out<<0;
		else
			out<<v2_derivative_of_rix(i,j,dist)/epsilon_Si;
		if(k < 299)
		{
			out<<", ";
		}
		j.x += 0.01*sigma_Si;
	}
	out<<endl;
	out.close();

	fileName.append("_X");
	out.open(fileName);
	out.precision(20);
	j.x = 0.5;
	for(int k = 0; k < 300; k++)
	{
		out<<distance2(i,j);
		if(k < 299)
		{
			out<<", ";
		}
		j.x += 0.01;
	}
	out.close();
}
///////////////////////////////// si2 force on y ///////////////////////////////////
void printSi2ForceY(string fileName)
{
	real3 i, j;
	ofstream out;
	double dist;

	i.x = i.y = i.z = 0.0;
	j.x = j.y = j.z = 0.0;

	fileName.append("\\F_SI_2_Y");
	out.open(fileName);
	out.precision(20);

	j.y = 0.5*sigma_Si;
	for(int k = 0; k < 300; k++)
	{
		dist = distance2(i,j);
		if(dist/sigma_Si >= a_Si)
			out<<0;
		else
			out<<v2_derivative_of_riy(i,j,dist)/epsilon_Si;
		if(k < 299)
		{
			out<<", ";
		}
		j.y += 0.01*sigma_Si;
	}
	out<<endl;
	out.close();

	fileName.append("_X");
	out.open(fileName);
	out.precision(20);
	j.y = 0.5;
	for(int k = 0; k < 300; k++)
	{
		out<<distance2(i,j);
		if(k < 299)
		{
			out<<", ";
		}
		j.y += 0.01;
	}
	out.close();
}
///////////////////////////////// si2 force on z ///////////////////////////////////
void printSi2ForceZ(string fileName)
{
	real3 i, j;
	ofstream out;
	double dist;

	i.x = i.y = i.z = 0.0;
	j.x = j.y = j.z = 0.0;

	fileName.append("\\F_SI_2_Z");
	out.open(fileName);
	out.precision(20);

	j.z = 0.5*sigma_Si;
	for(int k = 0; k < 300; k++)
	{
		dist = distance2(i,j);
		if(dist/sigma_Si >= a_Si)
			out<<0;
		else
			out<<v2_derivative_of_riz(i,j,dist)/epsilon_Si;
		if(k < 299)
		{
			out<<", ";
		}
		j.z += 0.01*sigma_Si;
	}
	out<<endl;
	out.close();

	fileName.append("_X");
	out.open(fileName);
	out.precision(20);
	j.z = 0.5;
	for(int k = 0; k < 300; k++)
	{
		out<<distance2(i,j);
		if(k < 299)
		{
			out<<", ";
		}
		j.z += 0.01;
	}
	out.close();
}

///////////////////////////////// si3 potential Cosinus ///////////////////////////////////
void printSi3PotentialCos(string fileName)
{
	real3 i, j, k;
	ofstream out;
	double cosal, rij, rik, rjk;

	i.x = i.y = i.z = 0.0;
	j.x = j.y = j.z = 0.0;
	k.x = k.y = k.z = 0.0;

	fileName.append("\\V_SI_3_COS");
	out.open(fileName);
	out.precision(20);

	k.x = pow(2.0,(1.0/6.0))*sigma_Si;
	for(int n = -si3resolution; n < si3resolution; n++)
	{
		cosal = (1.0/3.0)+((double)n/400.0);
		j.x = -cosal*pow(2.0,(1.0/6.0))*sigma_Si;
		j.y = sqrt(1.0-(cosal*cosal))*pow(2.0,(1.0/6.0))*sigma_Si;

		rij = distance2(i,j)/sigma_Si;
		rik = distance2(i,k)/sigma_Si;
		rjk = distance2(j,k)/sigma_Si;

		out<<v3(rij,rik,rjk)/epsilon_Si;
		if(n < si3resolution-1)
		{
			out<<", ";
		}
	}
	out<<endl;
	out.close();

	fileName.append("_X");
	out.open(fileName);
	out.precision(20);
	for(int n = -si3resolution; n < si3resolution; n++)
	{
		cosal = ((1.0/3.0)+((double)n/400.0));
		out<<(cosal);
		if(n < si3resolution-1)
		{
			out<<", ";
		}
	}
	out.close();
}
///////////////////////////////// si3 force Cosinus////////////////////////////////////////////
void printSi3ForceCos(string fileName)
{
	real3 i, j, k;
	ofstream out;
	double cosal, rij, rik, rjk;

	i.x = i.y = i.z = 0.0;
	j.x = j.y = j.z = 0.0;
	k.x = k.y = k.z = 0.0;

	fileName.append("\\F_SI_3_COS");
	out.open(fileName);
	out.precision(20);

	k.x = pow(2.0,(1.0/6.0))*sigma_Si;
	for(int n = -si3resolution; n < si3resolution; n++)
	{
		cosal = (1.0/3.0)+((double)n/400.0);
		j.x = -cosal*pow(2.0,(1.0/6.0))*sigma_Si;
		j.y = sqrt(1.0-(cosal*cosal))*pow(2.0,(1.0/6.0))*sigma_Si;

		rij = distance2(i, j);
		rik = distance2(i, k);
		rjk = distance2(j, k);

		if((rik/sigma_Si < a_Si && rij/sigma_Si < a_Si) || (rjk/sigma_Si < a_Si && rij/sigma_Si < a_Si) || (rjk/sigma_Si < a_Si && rik/sigma_Si < a_Si))
		{
			out<<((v3_derivative_of_rix(i, j, k, rij, rik, rjk))/epsilon_Si) + ((v3_derivative_of_rix(i, k, j, rik, rij, rjk))/epsilon_Si);
		}
		else
		{
			out<<0;
		}

		if(n < si3resolution-1)
		{
			out<<", ";
		}
	}
	out<<endl;
	out.close();

	fileName.append("_X");
	out.open(fileName);
	out.precision(20);
	for(int n = -si3resolution; n < si3resolution; n++)
	{
		cosal = ((1.0/3.0)+((double)n/400.0));
		out<<cosal;
		if(n < si3resolution-1)
		{
			out<<", ";
		}
	}
	out.close();	
}

///////////////////////////////// si3 potential distance ///////////////////////////////////
void printSi3PotentialDist(string fileName)
{
	real3 i, j, k;
	ofstream out;
	double rij, rik, rjk;

	i.x = i.y = i.z = 0.0;
	j.x = j.y = j.z = 0.0;
	k.x = k.y = k.z = 0.0;

	fileName.append("\\V_SI_3_DIST");
	out.open(fileName);
	out.precision(20);

	k.x = pow(2.0,(1.0/6.0))*sigma_Si;
	j.y = pow(2.0,(1.0/6.0))*sigma_Si;
	j.x = -2*pow(2.0,(1.0/6.0))*sigma_Si;
	for(int n = 0; n < 1800; n++)
	{
		rij = distance2(i,j)/sigma_Si;
		rik = distance2(i,k)/sigma_Si;
		rjk = distance2(j,k)/sigma_Si;

		out<<v3(rij,rik,rjk)/epsilon_Si;
		if(n < 1799)
		{
			out<<", ";
		}
		j.x += 0.01;
	}
	out<<endl;
	out.close();

	fileName.append("_X");
	out.open(fileName);
	out.precision(20);
	j.x = -2*pow(2.0,(1.0/6.0))*sigma_Si;
	for(int n = 0; n < 1800; n++)
	{
		out<<(j.x);
		if(n < 1799)
		{
			out<<", ";
		}
		j.x += 0.01;
	}
	out.close();
}
///////////////////////////////// si3 force distance////////////////////////////////////////////
void printSi3ForceDist(string fileName)
{
	real3 i, j, k;
	ofstream out;
	double rij, rik, rjk;

	i.x = i.y = i.z = 0.0;
	j.x = j.y = j.z = 0.0;
	k.x = k.y = k.z = 0.0;

	fileName.append("\\F_SI_3_DIST");
	out.open(fileName);
	out.precision(20);

	k.x = pow(2.0,(1.0/6.0))*sigma_Si;
	j.y = pow(2.0,(1.0/6.0))*sigma_Si;
	j.x = -2*pow(2.0,(1.0/6.0))*sigma_Si;
	for(int n = 0; n < 1800; n++)
	{
		rij = distance2(i, j);
		rik = distance2(i, k);
		rjk = distance2(j, k);

		if((rik/sigma_Si < a_Si && rij/sigma_Si < a_Si) || (rjk/sigma_Si < a_Si && rij/sigma_Si < a_Si) || (rjk/sigma_Si < a_Si && rik/sigma_Si < a_Si))
		{
			out<<((((v3_derivative_of_rix(i, j, k, rij, rik, rjk)) + (v3_derivative_of_rix(i, k, j, rik, rij, rjk)))/epsilon_Si)+(((v3_derivative_of_riy(i, j, k, rij, rik, rjk)) + (v3_derivative_of_riy(i, k, j, rik, rij, rjk)))/epsilon_Si) + (((v3_derivative_of_riz(i, j, k, rij, rik, rjk)) + (v3_derivative_of_riz(i, k, j, rik, rij, rjk)))/epsilon_Si));
		}
		else
		{
			out<<0;
		}

		if(n < 1799)
		{
			out<<", ";
		}
		j.x += 0.01;
	}
	out<<endl;
	out.close();

	fileName.append("_X");
	out.open(fileName);
	out.precision(20);
	j.x = -2*pow(2.0,(1.0/6.0));
	for(int n = 0; n < 1800; n++)
	{
		out<<j.x;
		if(n < 1799)
		{
			out<<", ";
		}
		j.x += 0.01;
	}
	out.close();	
}

///////////////////////////////// xe potential ////////////////////////////////////////////
void printXePotential(string fileName)
{
	real3 i, j;
	ofstream out;
	double dist;

	i.x = i.y = i.z = 0.0;
	j.x = j.y = j.z = 0.0;

	fileName.append("\\V_XE");
	out.open(fileName);
	out.precision(20);

	j.x = 0.8*sigma_Xe_Xe;
	for(int k = 0; k < 300; k++)
	{
		dist = distance2(i,j);

		out<<lennardJonesPotential(dist,sigma_Xe_Xe,epsilon_Xe_Xe)/epsilon_Xe_Xe;

		if(k < 299)
		{
			out<<", ";
		}
		j.x += 0.01*sigma_Xe_Xe;
	}
	out<<endl;
	out.close();

	fileName.append("_X");
	out.open(fileName);
	out.precision(20);
	j.x = 0.8;
	for(int k = 0; k < 300; k++)
	{
		out<<distance2(i,j);
		if(k < 299)
		{
			out<<", ";
		}
		j.x += 0.01;
	}
	out.close();
}
///////////////////////////////// xe force //////////////////////////////////////////////////
void printXeForce(string fileName)
{
	real3 i, j;
	ofstream out;
	double dist;

	i.x = i.y = i.z = 0.0;
	j.x = j.y = j.z = 0.0;

	fileName.append("\\F_XE");
	out.open(fileName);
	out.precision(20);

	j.x = 0.8*sigma_Xe_Xe;
	for(int k = 0; k < 300; k++)
	{
		dist = distance2(i,j);

		out<<lennardJonesForce(dist,sigma_Xe_Xe,epsilon_Xe_Xe)/epsilon_Xe_Xe;

		if(k < 299)
		{
			out<<", ";
		}
		j.x += 0.01*sigma_Xe_Xe;
	}
	out<<endl;
	out.close();

	fileName.append("_X");
	out.open(fileName);
	out.precision(20);
	j.x = 0.8;
	for(int k = 0; k < 300; k++)
	{
		out<<distance2(i,j);
		if(k < 299)
		{
			out<<", ";
		}
		j.x += 0.01;
	}
	out.close();
}
#endif //DEBUGPRINTS_H_