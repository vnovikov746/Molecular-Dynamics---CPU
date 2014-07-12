/*
 *
 */

#include "InitiateSystem.h"

int readConfig(string configFilename, configurations *config)
{
	string line;
	ifstream configurationFile;
	config->SI_PARTICLES = 0;
	config->XE_PARTICLES = 0;
	int pos = 0;
	int tempNumber = 0;
	configurationFile.open(configFilename);
	if(!configurationFile.good())
	{
		cout<<"Can not open configuration file"<<endl;
		return EXIT_FAILURE;
	}

	while(getline(configurationFile, line))
	{
		if(line.find("//") && !line.empty())// If comment is at the start of the line, find will return 0.
		{
			if(line.find("SI_LENGTH") == 0)
			{
				pos = line.find_first_of('=');
				config->SI_LENGTH = atoi(line.substr(pos+1, line.size()).c_str());
			}
			else if(line.find("SI_WIDTH") == 0)
			{
				pos = line.find_first_of('=');
				config->SI_WIDTH = atoi(line.substr(pos+1, line.size()).c_str());
			}
			else if(line.find("SI_HEIGHT") == 0)
			{
				pos = line.find_first_of('=');
				config->SI_HEIGHT = atoi(line.substr(pos+1, line.size()).c_str());
			}
			else if(line.find("XE_LENGTH") == 0)
			{
				pos = line.find_first_of('=');
				config->XE_LENGTH = atoi(line.substr(pos+1, line.size()).c_str());
			}
			else if(line.find("XE_WIDTH") == 0)
			{
				pos = line.find_first_of('=');
				config->XE_WIDTH = atoi(line.substr(pos+1, line.size()).c_str());
			}
			else if(line.find("XE_HEIGHT") == 0)
			{
				pos = line.find_first_of('=');
				config->XE_HEIGHT = atoi(line.substr(pos+1, line.size()).c_str());
			}
			else if(line.find("LA_SPACE") == 0)
			{
				pos = line.find_first_of('=');
				config->LA_SPACE = atof(line.substr(pos+1, line.size()).c_str());
			}
			else if(line.find("ConditionsFilename") == 0)
			{
				pos = line.find_first_of('=');
				config->CONDITIONS_FILE = new string(line.substr(pos+1, line.size()));
			}
			else if(line.find("TEMPERATURE") == 0)
			{
				pos = line.find_first_of('=');
				config->TEMPERATURE = atof(line.substr(pos+1, line.size()).c_str())*rb;
			}
			else if(line.find("USE_NEIGHBOR_LISTS") == 0)
			{
				pos = line.find_first_of('=');
				config->USE_NEIGHBOR_LISTS = atoi(line.substr(pos+1, line.size()).c_str()) > 0 ? true : false;
			}
			else if(line.find("DEBUG_SI") == 0)
			{
				pos = line.find_first_of('=');
				config->DEBUG_SI = atoi(line.substr(pos+1, line.size()).c_str()) > 0 ? true : false;
			}
			else if(line.find("DEBUG_XE") == 0)
			{
				pos = line.find_first_of('=');
				config->DEBUG_XE = atoi(line.substr(pos+1, line.size()).c_str()) > 0 ? true : false;
			}
			else if(line.find("CHOOSE_VELO_SI") == 0)
			{
				pos = line.find_first_of('=');
				config->CHOOSE_VELO_SI = atoi(line.substr(pos+1, line.size()).c_str()) > 0 ? true : false;
			}
			else if(line.find("CHOOSE_VELO_XE") == 0)
			{
				pos = line.find_first_of('=');
				config->CHOOSE_VELO_XE = atoi(line.substr(pos+1, line.size()).c_str()) > 0 ? true : false;
			}
			else if(line.find("SYS_PAUSE") == 0)
			{
				pos = line.find_first_of('=');
				config->SYS_PAUSE = atoi(line.substr(pos+1, line.size()).c_str()) > 0 ? true : false;
			}
			else if(line.find("ANIMATE") == 0)
			{
				pos = line.find_first_of('=');
				config->ANIMATE = atoi(line.substr(pos+1, line.size()).c_str()) > 0 ? true : false;
			}
			else if(line.find("useLennardJonesPotentialForSi") == 0)
			{
				pos = line.find_first_of('=');
				config->useLennardJonesPotentialForSi = atoi(line.substr(pos+1, line.size()).c_str()) > 0 ? true : false;
			}
			else if(line.find("TIMESTEPS") == 0)
			{
				pos = line.find_first_of('=');
				config->TIMESTEPS = atof(line.substr(pos+1, line.size()).c_str());
			}
			else if(line.find("STEPS") == 0)
			{
				pos = line.find_first_of('=');
				config->STEPS = atoi(line.substr(pos+1, line.size()).c_str());
			}
		}
	}

	config->INPUT_FILE = new string("particlesFile");

	configurationFile.close();

	return EXIT_SUCCESS;
}

int readInput(string inputFile, allLists *lists, configurations* config)
{
	string line;
	int posInSi = 0;
	int posInXe = 0;
	int posInLine_1 = 0;
	int posInLine_2 = 0;
	double tempNumber = 0.0;
	ifstream particlesFile;
	particlesFile.open(inputFile);
	if(!particlesFile.good())
	{
		cout<<"Can not open particles file"<<endl;
		return EXIT_FAILURE;
	}

	while(getline(particlesFile, line))
	{
		if(line.find("//") && !line.empty())// If comment is at the start of the line, find will return 0.
		{
			if(line.find("Si") == 0)
			{
				posInLine_1 = line.find_first_of('x');
				posInLine_2 = line.find_first_of('y');
				tempNumber = atof(line.substr(posInLine_1+2, posInLine_2).c_str());
				lists->siParticles[posInSi].position.x = tempNumber;
				posInLine_1 = line.find_first_of('z');
				tempNumber = atof(line.substr(posInLine_2+2, posInLine_1).c_str());
				lists->siParticles[posInSi].position.y = tempNumber;
				tempNumber = atof(line.substr(posInLine_1+2, line.size()).c_str());
				lists->siParticles[posInSi].position.z = tempNumber;

				posInSi++;
			}
			else if(line.find("Xe") == 0)
			{
				posInLine_1 = line.find_first_of('x');
				posInLine_2 = line.find_first_of('y');
				tempNumber = atof(line.substr(posInLine_1+2, posInLine_2).c_str());
				lists->xeParticles[posInXe].position.x = tempNumber;
				posInLine_1 = line.find_first_of('z');
				tempNumber = atof(line.substr(posInLine_2+2, posInLine_1).c_str());
				lists->xeParticles[posInXe].position.y = tempNumber;
				tempNumber = atof(line.substr(posInLine_1+2, line.size()).c_str());
				lists->xeParticles[posInXe].position.z = tempNumber;

				posInXe++;
			}
		}
	}

	particlesFile.close();

	return EXIT_SUCCESS;
}

int listsAlloc(configurations* config, allLists* lists)
{
	try
	{
		lists->siParticles = new particleStruct[config->SI_PARTICLES];
		lists->xeParticles = new particleStruct[config->XE_PARTICLES];
		for(int i = 0; i < config->SI_PARTICLES; i++)
		{
			lists->siParticles[i].velocity.x = 0.0;
			lists->siParticles[i].velocity.y = 0.0;
			lists->siParticles[i].velocity.z = 0.0;
			lists->siParticles[i].force.x = 0.0;
			lists->siParticles[i].force.y = 0.0;
			lists->siParticles[i].force.z = 0.0;
			lists->siParticles[i].aAcc.x = 0.0;
			lists->siParticles[i].aAcc.y = 0.0;
			lists->siParticles[i].aAcc.z = 0.0;
			lists->siParticles[i].bAcc.x = 0.0;
			lists->siParticles[i].bAcc.y = 0.0;
			lists->siParticles[i].bAcc.z = 0.0;
			lists->siParticles[i].cAcc.x = 0.0;
			lists->siParticles[i].cAcc.y = 0.0;
			lists->siParticles[i].cAcc.z = 0.0;
			lists->siParticles[i].type = Si;
			lists->siParticles[i].siNeighbors = new int[1];
			lists->siParticles[i].siNeighbors[0] = -1;
			lists->siParticles[i].xeNeighbors = new int[1];
			lists->siParticles[i].xeNeighbors[0] = -1;
		}
		for(int i = 0; i < config->XE_PARTICLES; i++)
		{
			lists->xeParticles[i].velocity.x = 0.0;
			lists->xeParticles[i].velocity.y = 0.0;
			lists->xeParticles[i].velocity.z = 0.0;
			lists->xeParticles[i].force.x = 0.0;
			lists->xeParticles[i].force.y = 0.0;
			lists->xeParticles[i].force.z = 0.0;
			lists->xeParticles[i].aAcc.x = 0.0;
			lists->xeParticles[i].aAcc.y = 0.0;
			lists->xeParticles[i].aAcc.z = 0.0;
			lists->xeParticles[i].bAcc.x = 0.0;
			lists->xeParticles[i].bAcc.y = 0.0;
			lists->xeParticles[i].bAcc.z = 0.0;
			lists->xeParticles[i].cAcc.x = 0.0;
			lists->xeParticles[i].cAcc.y = 0.0;
			lists->xeParticles[i].cAcc.z = 0.0;
			lists->xeParticles[i].type = Xe;
			lists->xeParticles[i].xeNeighbors = new int[1];
			lists->xeParticles[i].xeNeighbors[0] = -1;
			lists->xeParticles[i].siNeighbors = new int[1];
			lists->xeParticles[i].siNeighbors[0] = -1;
		}
	}
	catch(exception e)
	{
		cout<<"FAILED TO ALLOCATE MEMORY\n";
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

int listsDestroy(configurations *config, allLists *lists)
{
	try
	{
		for(int i = 0; i < config->SI_PARTICLES; i++)
		{
			free(lists->siParticles[i].siNeighbors);
			free(lists->siParticles[i].xeNeighbors);
		}
		for(int i = 0; i < config->XE_PARTICLES; i++)
		{
			free(lists->xeParticles[i].siNeighbors);
			free(lists->xeParticles[i].xeNeighbors);
		}
		free(lists->siParticles);
		free(lists->xeParticles);
	}
	catch(exception e)
	{
		cout<<"FAILED TO DEALOCATE MEMORY\n";
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
