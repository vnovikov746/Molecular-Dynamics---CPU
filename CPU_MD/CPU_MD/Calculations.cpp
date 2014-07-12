/*
 *
 */

#include "Calculations.h"

double distance2(real3 i, real3 j)
{
	return sqrt((i.x-j.x)*(i.x-j.x) + (i.y-j.y)*(i.y-j.y) + (i.z-j.z)*(i.z-j.z));
}

//-------------------------- calculate Force Si ------------------------------//
void calculateForce_Si(int MAX_SI_NEIGHBORS, int MAX_XE_NEIGHBORS, particleStruct particle, particleStruct* siParticles, particleStruct* xeParticles, int i, configurations config)
{
	real3 *siNeigborsPositions = NULL;
	real3 *xeNeigborsPositions = NULL;

	real3 iPosition = particle.position;
	int countSi = 0;
	int countXe = 0;
	double r_ij = 0.0;
	double r_ik = 0.0;
	double r_jk = 0.0;
	siParticles[i].force.x = 0.0;
	siParticles[i].force.y = 0.0;
	siParticles[i].force.z = 0.0;
	real3 jPosition;
	real3 kPosition;

	if(config.USE_NEIGHBOR_LISTS)
	{
		siNeigborsPositions = new real3[MAX_SI_NEIGHBORS];
		xeNeigborsPositions = new real3[MAX_XE_NEIGHBORS];
		for(countSi = 0; countSi < MAX_SI_NEIGHBORS && particle.siNeighbors[countSi] != -1; countSi++)
		{
			siNeigborsPositions[countSi] = siParticles[particle.siNeighbors[countSi]].position;
		}
		for(countXe = 0; countXe < MAX_XE_NEIGHBORS && particle.xeNeighbors[countXe] != -1; countXe++)
		{
			xeNeigborsPositions[countXe] = xeParticles[particle.xeNeighbors[countXe]].position;
		}
	}

	if(!config.USE_NEIGHBOR_LISTS)
	{
		countSi = config.SI_PARTICLES;
		countXe = config.XE_PARTICLES;
	}

	if(config.useLennardJonesPotentialForSi)
	{
		for(int j = 0; j < countSi; j++)
		{
			if(j != i)
			{
				if(!config.USE_NEIGHBOR_LISTS)
					jPosition = siParticles[j].position;
				else
					jPosition = siNeigborsPositions[j];

				r_ij = distance2(iPosition, jPosition);	
	//			if(r_ij/sigma_Si < 1.8)
				{
					siParticles[i].force.x += (iPosition.x-jPosition.x)*lennardJonesForce(r_ij,sigma_Si,epsilon_Si);
					siParticles[i].force.y += (iPosition.y-jPosition.y)*lennardJonesForce(r_ij,sigma_Si,epsilon_Si);
					siParticles[i].force.z += (iPosition.z-jPosition.z)*lennardJonesForce(r_ij,sigma_Si,epsilon_Si);
				}
			}
		}
	}

	else
	{	
		for(int j = 0; j < countSi; j++)
		{
			if(j != i)
			{
				if(!config.USE_NEIGHBOR_LISTS)
					jPosition = siParticles[j].position;
				else
					jPosition = siNeigborsPositions[j];
			
				r_ij = distance2(iPosition, jPosition);

				if(r_ij/sigma_Si < a_Si)
				{
					siParticles[i].force.x -= v2_derivative_of_rix(iPosition, jPosition, r_ij);
					siParticles[i].force.y -= v2_derivative_of_riy(iPosition, jPosition, r_ij);
					siParticles[i].force.z -= v2_derivative_of_riz(iPosition, jPosition, r_ij);
				}
				for(int k = 0; k < countSi; k++)
				{
					if(k != i && k != j)
					{
						if(!config.USE_NEIGHBOR_LISTS)
						{
							kPosition = siParticles[k].position;							
						}
						else
						{
							kPosition = siNeigborsPositions[k];
						}

						r_ik = distance2(iPosition, kPosition);
						r_jk = distance2(jPosition, kPosition);
						if((r_ij/sigma_Si < a_Si && r_ik/sigma_Si < a_Si) || (r_ij/sigma_Si < a_Si && r_jk/sigma_Si < a_Si) || (r_ik/sigma_Si < a_Si && r_jk/sigma_Si < a_Si))
						{
							siParticles[i].force.x -= v3_derivative_of_rix(iPosition, jPosition, kPosition, r_ij, r_ik, r_jk);
							siParticles[i].force.x -= v3_derivative_of_rix(iPosition, kPosition, jPosition, r_ik, r_ij, r_jk);
							siParticles[i].force.y -= v3_derivative_of_riy(iPosition, jPosition, kPosition, r_ij, r_ik, r_jk);
							siParticles[i].force.y -= v3_derivative_of_riy(iPosition, kPosition, jPosition, r_ik, r_ij, r_jk);
							siParticles[i].force.z -= v3_derivative_of_riz(iPosition, jPosition, kPosition, r_ij, r_ik, r_jk);
							siParticles[i].force.z -= v3_derivative_of_riz(iPosition, kPosition, jPosition, r_ik, r_ij, r_jk);
						}
					}
				}
			}
		}
	}

/*	for(int j = 0; j < countXe; j++)
	{
			if(!config.USE_NEIGHBOR_LISTS)
				jPosition = xeParticles[j].position;
			else
				jPosition = xeNeigborsPositions[j];

		jPosition.x += 0.25*config->SI_LENGTH*space_Si;
		jPosition.y += 0.25*config->SI_LENGTH*space_Si;
		jPosition.z += config->SI_HEIGHT+config->LA_SPACE;
		r_ij = distance2(iPosition, jPosition);
		if(r_ij/sigma_Si_Xe < a_Si_Xe)
		{
			siParticles[i].force.x += ((iPosition.x-jPosition.x)*(lennardJonesForce(r_ij,sigma_Si_Xe,epsilon_Si_Xe)));//-lennardJonesForce(2.5*sigma_Si,sigma_Si_Xe,epsilon_Si_Xe)));
			siParticles[i].force.y += ((iPosition.y-jPosition.y)*(lennardJonesForce(r_ij,sigma_Si_Xe,epsilon_Si_Xe)));//-lennardJonesForce(2.5*sigma_Si,sigma_Si_Xe,epsilon_Si_Xe)));
			siParticles[i].force.z += ((iPosition.z-jPosition.z)*(lennardJonesForce(r_ij,sigma_Si_Xe,epsilon_Si_Xe)));//-lennardJonesForce(2.5*sigma_Si,sigma_Si_Xe,epsilon_Si_Xe)));
		}
	}*/

	if(config.USE_NEIGHBOR_LISTS)
	{
		free(siNeigborsPositions);
		free(xeNeigborsPositions);
	}
}
//----------------------------------------------------------------------------//

void calculateForce_Xe(int MAX_SI_NEIGHBORS, int MAX_XE_NEIGHBORS, particleStruct particle, particleStruct* xeParticles, particleStruct* siParticles, int i, configurations config)
{
	real3 *siNeigborsPositions = NULL;
	real3 *xeNeigborsPositions = NULL;

	real3 iPosition = particle.position;

	int countSi = 0;
	int countXe = 0;
	double r_ij = 0.0;
	xeParticles[i].force.x = 0.0;
	xeParticles[i].force.y = 0.0;
	xeParticles[i].force.z = 0.0;
	real3 jPosition;

	if(config.USE_NEIGHBOR_LISTS)
	{
		siNeigborsPositions = new real3[MAX_SI_NEIGHBORS];
		xeNeigborsPositions = new real3[MAX_XE_NEIGHBORS];
		for(countSi = 0; countSi < MAX_SI_NEIGHBORS && particle.siNeighbors[countSi] != -1; countSi++)
		{
			siNeigborsPositions[countSi] = siParticles[particle.siNeighbors[countSi]].position;
		}
		for(countXe = 0; countXe < MAX_XE_NEIGHBORS && particle.xeNeighbors[countXe] != -1; countXe++)
		{
			xeNeigborsPositions[countXe] = xeParticles[particle.xeNeighbors[countXe]].position;
		}
	}

	if(!config.USE_NEIGHBOR_LISTS)
	{
		countSi = config.SI_PARTICLES;
		countXe = config.XE_PARTICLES;
	}

	for(int j = 0; j < countXe; j++)
	{
		if(j != i)
		{
			if(!config.USE_NEIGHBOR_LISTS)
				jPosition = xeParticles[j].position;
			else
				jPosition = xeNeigborsPositions[j];

			r_ij = distance2(iPosition, jPosition);
			if(r_ij/sigma_Xe_Xe < xe_Cluster)
			{
				xeParticles[i].force.x += (iPosition.x-jPosition.x)*lennardJonesForce(r_ij,sigma_Xe_Xe,epsilon_Xe_Xe);
				xeParticles[i].force.y += (iPosition.y-jPosition.y)*lennardJonesForce(r_ij,sigma_Xe_Xe,epsilon_Xe_Xe);
				xeParticles[i].force.z += (iPosition.z-jPosition.z)*lennardJonesForce(r_ij,sigma_Xe_Xe,epsilon_Xe_Xe);
			}
		}
	}

/*	for(int j = 0; j < countSi; j++)
	{
		if(!config.USE_NEIGHBOR_LISTS)
			jPosition = siParticles[j].position;
		else
			jPosition = siNeigborsPositions[j];

		r_ij = distance2(iPosition, jPosition);
		if(r_ij/sigma_Si_Xe < a_Si_Xe)
		{
			xeParticles[i].force.x += (iPosition.x-jPosition.x)*lennardJonesForce(r_ij,sigma_Si_Xe,epsilon_Si_Xe);
			xeParticles[i].force.y += (iPosition.y-jPosition.y)*lennardJonesForce(r_ij,sigma_Si_Xe,epsilon_Si_Xe);
			xeParticles[i].force.z += (iPosition.z-jPosition.z)*lennardJonesForce(r_ij,sigma_Si_Xe,epsilon_Si_Xe);
		}
	}*/

	if(config.USE_NEIGHBOR_LISTS)
	{
		free(siNeigborsPositions);
		free(xeNeigborsPositions);
	}
}

//----------------------calculate total potential of Si -----------------------//
double V_total_Si(int SI_PARTICLES, particleStruct* siParticles)
{
	double vTotal = 0.0;
	for(int i = 0; i < SI_PARTICLES; i++)
	{
		for(int j = i+1; j < SI_PARTICLES; j++)
		{
			vTotal += v2(distance2(siParticles[i].position, siParticles[j].position)/sigma_Si);
		}
	}
	for(int i = 0; i < SI_PARTICLES; i++)
	{
		for(int j = i+1; j < SI_PARTICLES; j++)
		{
			for(int k = j+1; k < SI_PARTICLES; k++)
			{
				vTotal += v3(distance2(siParticles[i].position,siParticles[j].position)/sigma_Si, distance2(siParticles[i].position,siParticles[k].position)/sigma_Si, distance2(siParticles[j].position,siParticles[k].position)/sigma_Si);
			}
		}
	}
	for(int i = 0; i < SI_PARTICLES; i++)
	{
//		cout<<(pow((sqrt((siParticles[i].velocity.x*siParticles[i].velocity.x)+(siParticles[i].velocity.y*siParticles[i].velocity.y)+(siParticles[i].velocity.z*siParticles[i].velocity.z))),2.0)*SiMass)/2<<endl;
//		vTotal += (pow((sqrt((siParticles[i].velocity.x*siParticles[i].velocity.x)+(siParticles[i].velocity.y*siParticles[i].velocity.y)+(siParticles[i].velocity.z*siParticles[i].velocity.z))),2.0)*SiMass)/2;
	}
	return vTotal;
}
//----------------------------------------------------------------------------//
double V_total_Si_Xe(int SI_PARTICLES, particleStruct* siParticles, int XE_PARTICLES, particleStruct* xeParticles)
{
	double vTotal = 0.0;
	for(int i = 0; i < SI_PARTICLES; i++)
	{
		for(int j = 0; j < XE_PARTICLES; j++)
		{
//			if(distance2(siParticles[i].position, xeParticles[j].position)/sigma_Si_Xe < a_Si_Xe)
				vTotal += lennardJonesPotential(distance2(siParticles[i].position, xeParticles[j].position), sigma_Si_Xe, epsilon_Si_Xe);
		}
	}
	return vTotal;
}
double V_total_Xe(int XE_PARTICLES, particleStruct* xeParticles)
{
	double vTotal = 0.0;
	for(int i = 0; i < XE_PARTICLES-1; i++)
	{
		for(int j = i+1; j < XE_PARTICLES; j++)
		{
			if(distance2(xeParticles[i].position, xeParticles[j].position)/sigma_Xe_Xe < a_Si_Xe)
				vTotal += lennardJonesPotential(distance2(xeParticles[i].position, xeParticles[j].position), sigma_Xe_Xe, epsilon_Xe_Xe);
		}
	}
	for(int i = 0; i < XE_PARTICLES; i++)
	{
		vTotal += (pow((sqrt((xeParticles[i].velocity.x*xeParticles[i].velocity.x)+(xeParticles[i].velocity.y*xeParticles[i].velocity.y)+(xeParticles[i].velocity.z*xeParticles[i].velocity.z))),2.0)*XeMass)/2;
	}
	return vTotal;
}
void V_total(int SI_PARTICLES, particleStruct* siParticles, int XE_PARTICLES, particleStruct* xeParticles)
{
	double vTotal = 0.0;
	vTotal += V_total_Si(SI_PARTICLES, siParticles);
	vTotal += V_total_Si_Xe(SI_PARTICLES, siParticles, XE_PARTICLES, xeParticles);
	vTotal += V_total_Xe(XE_PARTICLES, xeParticles);
	cout<<vTotal<<endl;
}

void initiateAcceleration(particleStruct *particles, int listSize, double mass)
{
	for(int i = 0; i < listSize; i++)
	{
		particles[i].aAcc.x = particles[i].force.x/mass;
		particles[i].aAcc.y = particles[i].force.y/mass;
		particles[i].aAcc.z = particles[i].force.z/mass;
	}
}

void predict(particleStruct *particles, int listSize, float dt)
{
	double c1 = dt;
    double c2 = c1*dt/2.0;
    double c3 = c2*dt/3.0;
    double c4 = c3*dt/4.0;

	for(int i = 0; i < listSize; i++)
	{
		particles[i].position.x += c1*particles[i].velocity.x + c2*particles[i].aAcc.x + c3*particles[i].bAcc.x + c4*particles[i].cAcc.x;
		particles[i].position.y += c1*particles[i].velocity.y + c2*particles[i].aAcc.y + c3*particles[i].bAcc.y + c4*particles[i].cAcc.y;
		particles[i].position.z += c1*particles[i].velocity.z + c2*particles[i].aAcc.z + c3*particles[i].bAcc.z + c4*particles[i].cAcc.z;
		particles[i].velocity.x += c1*particles[i].aAcc.x + c2*particles[i].bAcc.x + c3*particles[i].cAcc.x;
		particles[i].velocity.y += c1*particles[i].aAcc.y + c2*particles[i].bAcc.y + c3*particles[i].cAcc.y;
		particles[i].velocity.z += c1*particles[i].aAcc.z + c2*particles[i].bAcc.z + c3*particles[i].cAcc.z;
		particles[i].aAcc.x += c1*particles[i].bAcc.x + c2*particles[i].cAcc.x;
		particles[i].aAcc.y += c1*particles[i].bAcc.y + c2*particles[i].cAcc.y;
		particles[i].aAcc.z += c1*particles[i].bAcc.z + c2*particles[i].cAcc.z;
		particles[i].bAcc.x += c1*particles[i].cAcc.x;
		particles[i].bAcc.y += c1*particles[i].cAcc.y;
		particles[i].bAcc.z += c1*particles[i].cAcc.z;
	}
}

void correct(particleStruct *particles, double dt, int listSize, double mass)
{
	double c1 = dt ;
	double c2 = c1*dt/2.0;
	double c3 = c2*dt/3.0; 
	double c4 = c3*dt/4.0;

	double cr = GEAR1*c2;
	double cv = GEAR2*c2/c1;
	double cb = GEAR3*c2/c3;
	double cc = GEAR4*c2/c4;

	for(int i = 0; i < listSize; i++)
	{
		double axi = particles[i].force.x/mass;
      	double ayi = particles[i].force.y/mass;
      	double azi = particles[i].force.z/mass;

		double corrx = axi - particles[i].aAcc.x;
    	double corry = ayi - particles[i].aAcc.y;
    	double corrz = azi - particles[i].aAcc.z;

		particles[i].position.x += cr*corrx;
		particles[i].position.y += cr*corry;
		particles[i].position.z += cr*corrz;
		particles[i].velocity.x += cv*corrx;
		particles[i].velocity.y += cv*corry;
		particles[i].velocity.z += cv*corrz;
		particles[i].aAcc.x = axi;
		particles[i].aAcc.y = ayi;
		particles[i].aAcc.z = azi;
		particles[i].bAcc.x += cb*corrx;
		particles[i].bAcc.y += cb*corry;
		particles[i].bAcc.z += cb*corrz;
		particles[i].cAcc.x += cc*corrx;
		particles[i].cAcc.y += cc*corry;
    	particles[i].cAcc.z += cc*corrz;
	}
}