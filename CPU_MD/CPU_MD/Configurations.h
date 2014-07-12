/*
 * Molecular Dynamics Project
 * Author: Vladimir Novikov
 * This file contains the ........
 */

#ifndef CONFIGURATIONS_H_
#define CONFIGURATIONS_H_

#include <string>

struct configurations
{
	int SI_PARTICLES;
	int XE_PARTICLES;
	int MAX_SI_NEIGHBORS;
	int MAX_XE_NEIGHBORS;
	std::string* CONDITIONS_FILE;
	std::string* INPUT_FILE;
	int SI_LENGTH;
	int SI_WIDTH;
	int SI_HEIGHT;
	int XE_LENGTH;
	int XE_WIDTH;
	int XE_HEIGHT;
	double LA_SPACE;
	double TEMPERATURE;
	double TIMESTEPS;
	int STEPS;
	
	// debug parameters
	bool USE_NEIGHBOR_LISTS;
	bool DEBUG_SI;
	bool DEBUG_XE;
	bool CHOOSE_VELO_SI;
	bool CHOOSE_VELO_XE;
	bool SYS_PAUSE;
	bool ANIMATE;
	bool useLennardJonesPotentialForSi;
};

#endif	//CONFIGURATIONS_H_
