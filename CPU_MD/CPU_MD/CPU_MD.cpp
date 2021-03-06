/*
 * Molecular Dynamics Project.
 * Outhor: Vladimir Novikov.
 */

#include "Configurations.h"
#include "InitiateSystem.h"
#include "Calculations.h"
#include "NeighborLists.h"
#include <time.h>
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#include "Equi.h"
#include <GL/glut.h> 
#include "DebugPrints.h"
#include <minmax.h>

GLfloat *xp;
GLfloat *yp;
GLfloat *zp;
configurations config;
allLists lists;
ofstream file;

int place;
int resolution = 150;
int iterationNum;

time_t startT, endT;

double prevVPsi = 0.0;
double prevVKsi = 0.0;
double prevVPxe = 0.0;
double prevVKxe = 0.0;

bool debugSi;
bool debugXe;
bool chooseVeloSi;
bool chooseVeloXe;
bool useNeighborList;
bool animation;
bool sysPause;
bool useLennardJonesPotentialForSi;

void Animate(int,char**);
void Initialize();
void Draw();
void DrawPoints();
void calculatePositions();

ofstream si_place_and_velocity; 
ofstream xe_place_and_velocity; 
ofstream si_potential; 
ofstream xe_potential; 

int main(int argc, char *argv[])
{	
	cout.precision(19);
	si_place_and_velocity.open("si_place_and_velocity");
	xe_place_and_velocity.open("xe_place_and_velocity");
	si_potential.open("si_potential"); 
	xe_potential.open("xe_potential"); 
	if(!si_place_and_velocity.good() || !xe_place_and_velocity.good() || !si_potential.good() || !xe_potential.good())
	{
		cout<<"Can not open the files for debug"<<endl;
		return EXIT_FAILURE;
	}

	si_place_and_velocity.precision(19);
	xe_place_and_velocity.precision(19);
	si_potential.precision(19);
	xe_potential.precision(19);
	iterationNum = 0;


	cout<<"read configurations\n";
	readConfig(argv[1], &config);

	///// debug outputs //////
	if(config.PRINT_GRAPHS)
	{
		printDebug(config.OUT_FOR_GRAPHS);
	}
	//////////////////////////	

	////////////////////// DEBUG CONFIGURATIONS ////////////////////
	debugSi = config.DEBUG_SI;
	debugXe = config.DEBUG_XE;
	chooseVeloSi = config.CHOOSE_VELO_SI;
	chooseVeloXe = config.CHOOSE_VELO_XE;
	useNeighborList = config.USE_NEIGHBOR_LISTS;
	animation = config.ANIMATE;
	sysPause = config.SYS_PAUSE;
	useLennardJonesPotentialForSi = config.useLennardJonesPotentialForSi;
	/////////////////////////////////////////////////////////////////

	cout<<"write input file\n";
	writeParticlesInput(&config);

	cout<<"lists allocation\n";
	if(listsAlloc(&config, &lists) == EXIT_FAILURE)
	{
		system("pause");
		exit(EXIT_FAILURE);
	}

	cout<<"read input\n";
	readInput(*(config.INPUT_FILE), &lists, &config);

	if(chooseVeloSi)
	{
		cout<<"choose Si velocities\n";
		if(config.SI_PARTICLES > 0 && chooseVeloSi)
			chooseVelocities(lists.siParticles, config.SI_PARTICLES, config.TEMPERATURE, &config, SiMass);
	}

	if(chooseVeloXe)
	{
		cout<<"choose Si velocities\n";
		if(config.XE_PARTICLES > 0 && chooseVeloXe)
			chooseVelocities(lists.xeParticles, config.XE_PARTICLES, config.TEMPERATURE, &config, XeMass);
	}
	
	if(useNeighborList)
	{
		cout<<"build neighbors\n";
		divideToCells(lists.siParticles, lists.xeParticles, config.SI_PARTICLES, config.XE_PARTICLES, si_xe_Cluster);

		buildNeighbors(lists.siParticles, lists.xeParticles, &config, si_xe_Cluster, max(config.SI_LENGTH,config.XE_LENGTH)*2, max(config.SI_WIDTH,config.XE_WIDTH)*2, (config.SI_HEIGHT+config.XE_HEIGHT+config.LA_SPACE)*2, XeSi, config.SI_PARTICLES, config.XE_PARTICLES);
		buildNeighbors(lists.xeParticles, lists.siParticles, &config, si_xe_Cluster, max(config.SI_LENGTH,config.XE_LENGTH)*2, max(config.SI_WIDTH,config.XE_WIDTH)*2, (config.SI_HEIGHT+config.XE_HEIGHT+config.LA_SPACE)*2, XeSi, config.XE_PARTICLES, config.SI_PARTICLES);

		divideToCells(lists.siParticles, config.SI_PARTICLES, si_Cluster);
		buildNeighbors(lists.siParticles, NULL, &config, si_Cluster, max(config.SI_LENGTH,config.XE_LENGTH)*2, max(config.SI_WIDTH,config.XE_WIDTH)*2, (config.XE_HEIGHT+config.SI_HEIGHT+config.LA_SPACE)*2, Si, config.SI_PARTICLES, 0);

		divideToCells(lists.xeParticles, config.XE_PARTICLES, xe_Cluster);
		buildNeighbors(lists.xeParticles, NULL, &config, xe_Cluster,  max(config.SI_LENGTH,config.XE_LENGTH)*2, max(config.SI_WIDTH,config.XE_WIDTH)*2, (config.XE_HEIGHT+config.SI_HEIGHT+config.LA_SPACE)*2, Xe, config.XE_PARTICLES, 0);
	}
	
	cout<<"calculate forces for first time\n";
	calculateForce_Si(config.MAX_SI_NEIGHBORS, config.MAX_XE_NEIGHBORS,lists.siParticles, lists.xeParticles, config.SI_PARTICLES, config.XE_PARTICLES, config.USE_NEIGHBOR_LISTS, config.useLennardJonesPotentialForSi);
	calculateForce_Xe(config.MAX_SI_NEIGHBORS, config.MAX_XE_NEIGHBORS,lists.xeParticles, lists.siParticles, config.SI_PARTICLES, config.XE_PARTICLES, config.USE_NEIGHBOR_LISTS);

	if(config.SI_PARTICLES > 0)
		initiateAcceleration(lists.siParticles, config.SI_PARTICLES, SiMass);
	if(config.XE_PARTICLES > 0)
		initiateAcceleration(lists.xeParticles, config.XE_PARTICLES, XeMass);

	if(animation)
	{
		time(&startT);
		Animate(argc, argv);
	}

	else
	{
		for(int i = 0; i < config.STEPS; i++)
		{
			cout<<"-----------------------"<<iterationNum<<" ITERATION-----------------------\n";
			calculatePositions();
		}
		si_place_and_velocity.close();
		xe_place_and_velocity.close();
		si_potential.close();
		xe_potential.close();
	}
	system("pause");
	return EXIT_SUCCESS;
}

void calculatePositions()
{
	cout<<"------------- Iteration: "<<iterationNum<<" -------------"<<endl;
	if(debugSi)
	{
//		cout<<"SI: "<<endl<<endl;
		double vP = 0.0;
		double vK = 0.0;
		double V3 = 0.0;
		real3 iPosition;
		real3 jPosition;
		real3 kPosition;

		si_place_and_velocity<<"-----------------------"<<iterationNum<<" ITERATION-----------------------\n";
		for(int i = 0; i < config.SI_PARTICLES; i++)
		{
			double potentialForOneAtom = 0;
			iPosition = lists.siParticles[i].position;
			//print place and velocity
			si_place_and_velocity<<"\n"<<i<<" particle\n"<<"place: "<<iPosition.x<<"  "<<iPosition.y<<"  "<<iPosition.z<<"\n"<<
				"velocity: "<<lists.siParticles[i].velocity.x<<"  "<<lists.siParticles[i].velocity.y<<"  "<<lists.siParticles[i].velocity.z<<"\n"<<
				"k energy: "<<(((lists.siParticles[i].velocity.x*lists.siParticles[i].velocity.x)+(lists.siParticles[i].velocity.y*lists.siParticles[i].velocity.y)+(lists.siParticles[i].velocity.z*lists.siParticles[i].velocity.z))*SiMass)/2.0<<"\n";

			for(int j = 0; j < config.SI_PARTICLES; j++)
			{
				if(i != j)
				{
					jPosition = lists.siParticles[j].position;
					double r_ij = distance2(iPosition,jPosition);
					if(useLennardJonesPotentialForSi)
					{
						potentialForOneAtom += lennardJonesPotential(r_ij,sigma_Si,epsilon_Si);
						if(i < j)
						{
							vP += lennardJonesPotential(r_ij,sigma_Si,epsilon_Si);
						}
					}
					else
					{
						potentialForOneAtom += v2(r_ij/sigma_Si);
						if(i < j)
						{
							vP += v2(r_ij/sigma_Si);
						}
						for(int k = 0; k < config.SI_PARTICLES; k++)
						{
							if(k != i && k != j)
								{
								kPosition = lists.siParticles[k].position;
								double r_ik = distance2(iPosition,kPosition);
								double r_jk = distance2(jPosition,kPosition);
								if((r_ij/sigma_Si < a_Si && r_ik/sigma_Si < a_Si) || (r_ij/sigma_Si < a_Si && r_jk/sigma_Si < a_Si) || (r_ik/sigma_Si < a_Si && r_jk/sigma_Si < a_Si))
								{
									potentialForOneAtom += v3(r_ij/sigma_Si, r_ik/sigma_Si, r_jk/sigma_Si);
									if(i < j && j < k)
									{
										vP += v3(r_ij/sigma_Si, r_ik/sigma_Si, r_jk/sigma_Si);
									}
								}
							}
						}
					}
				}
			}
			si_place_and_velocity<<"potential: "<<potentialForOneAtom<<"\n";
		}
		for(int i = 0; i < config.SI_PARTICLES; i++)
		{
			vK += (((lists.siParticles[i].velocity.x*lists.siParticles[i].velocity.x)+(lists.siParticles[i].velocity.y*lists.siParticles[i].velocity.y)+(lists.siParticles[i].velocity.z*lists.siParticles[i].velocity.z))*SiMass)/2.0;
		}

		si_potential<<"-----------------------"<<iterationNum<<" ITERATION-----------------------\n";
		si_potential<<endl<<"Ep	   	    ";
		si_potential<<"Ek		            ";
		si_potential<<"E_total"<<"\n\n";
		si_potential<<endl<<vP<<" "; // Ep
		si_potential<<vK<<" "; // Ek
		si_potential<<vP+vK<<"\n\n"; // E_total = Ep + Ek
		si_potential<<endl<<"Change in Ep    ";
		si_potential<<"Change in Ek   ";
		si_potential<<"Change in E_total"<<"\n\n";
		si_potential<<endl<<(vP-prevVPsi)<<" "; // change in Ep
		si_potential<<(vK-prevVKsi)<<" "; // change in Ek
		si_potential<<((vP-prevVPsi)+(vK-prevVKsi))<<"\n\n"; // change in E_total
		// determine which energy changed more (non if equals)
		if(abs(vK-prevVKsi) < abs(vP-prevVPsi))
			si_potential<<"vP"<<endl<<endl;
		else if(abs(vK-prevVKsi) > abs(vP-prevVPsi))
			si_potential<<"vK"<<endl<<endl;
		else
			si_potential<<"non"<<endl<<endl;
		prevVPsi = vP;
		prevVKsi = vK;
		si_place_and_velocity<<"\n\n\n\n";
	}
	/////////////////////////////////////////////////////////////////////////////////

	if(debugXe)
	{
//		cout<<"XE: "<<endl<<endl;
		double vP = 0.0;
		double vK = 0.0;
		real3 iPosition;
		real3 jPosition;

		xe_place_and_velocity<<"-----------------------"<<iterationNum<<" ITERATION-----------------------\n";
		for(int i = 0; i < config.XE_PARTICLES; i++)
		{
			double potentialForOneAtom = 0;
			iPosition = lists.xeParticles[i].position;
			//print place and velocity
			xe_place_and_velocity<<"\n"<<i<<" particle\n"<<"place: "<<iPosition.x<<"  "<<iPosition.y<<"  "<<iPosition.z<<"\n"<<
				"velocity: "<<lists.xeParticles[i].velocity.x<<"  "<<lists.xeParticles[i].velocity.y<<"  "<<lists.xeParticles[i].velocity.z<<"\n"<<
				"k energy: "<<(((lists.xeParticles[i].velocity.x*lists.xeParticles[i].velocity.x)+(lists.xeParticles[i].velocity.y*lists.xeParticles[i].velocity.y)+(lists.xeParticles[i].velocity.z*lists.xeParticles[i].velocity.z))*XeMass)/2.0<<"\n";

			for(int j = 0; j < config.XE_PARTICLES; j++)
			{
				if(i != j)
				{
					jPosition = lists.xeParticles[j].position;
					potentialForOneAtom += lennardJonesPotential(distance2(iPosition, jPosition),sigma_Xe_Xe,epsilon_Xe_Xe);
					if(i < j)
					{
						vP += lennardJonesPotential(distance2(iPosition, jPosition),sigma_Xe_Xe,epsilon_Xe_Xe);
					}
				}
			}
			xe_place_and_velocity<<"potential: "<<potentialForOneAtom<<"\n";
		}
		for(int i = 0; i < config.XE_PARTICLES; i++)
		{
			vK += (((lists.xeParticles[i].velocity.x*lists.xeParticles[i].velocity.x)+(lists.xeParticles[i].velocity.y*lists.xeParticles[i].velocity.y)+(lists.xeParticles[i].velocity.z*lists.xeParticles[i].velocity.z))*XeMass)/2.0;
		}

		vP /= 2;
		xe_potential<<"-----------------------"<<iterationNum<<" ITERATION-----------------------\n";
		xe_potential<<endl<<"Ep	   	    ";
		xe_potential<<"Ek		            ";
		xe_potential<<"E_total"<<"\n\n";
		xe_potential<<endl<<vP<<" "; // Ep
		xe_potential<<vK<<" "; // Ek
		xe_potential<<vP+vK<<"\n\n"; // E_total = Ep + Ek
		xe_potential<<endl<<(vP-prevVPxe)<<" "; // change in Ep
		xe_potential<<(vK-prevVKxe)<<" "; // change in Ek
		xe_potential<<((vP-prevVPxe)+(vK-prevVKxe))<<"\n\n"; // change in E_total
		// determine which energy changed more (non if equals)
		if(abs(vK-prevVKxe) < abs(vP-prevVPxe))
			xe_potential<<"vP"<<endl<<endl;
		else if(abs(vK-prevVKxe) > abs(vP-prevVPxe))
			xe_potential<<"vK"<<endl<<endl;
		else
			xe_potential<<"non"<<endl<<endl;
		prevVPxe = vP;
		prevVKxe = vK;
		xe_place_and_velocity<<"\n\n\n\n";
	}

	if((debugSi || debugXe) && sysPause)
	{
		system("pause");
	}

	if (config.SI_PARTICLES > 0)
	{
		predict(lists.siParticles, config.SI_PARTICLES, config.TIMESTEPS);
	}
	if (config.XE_PARTICLES > 0)
	{
		predict(lists.xeParticles, config.XE_PARTICLES, config.TIMESTEPS);
	}
	
	calculateForce_Si(config.MAX_SI_NEIGHBORS, config.MAX_XE_NEIGHBORS,lists.siParticles, lists.xeParticles, config.SI_PARTICLES, config.XE_PARTICLES, config.USE_NEIGHBOR_LISTS, config.useLennardJonesPotentialForSi);
	calculateForce_Xe(config.MAX_SI_NEIGHBORS, config.MAX_XE_NEIGHBORS,lists.xeParticles, lists.siParticles, config.SI_PARTICLES, config.XE_PARTICLES, config.USE_NEIGHBOR_LISTS);

	if (config.SI_PARTICLES > 0)
	{
		correct(lists.siParticles, config.TIMESTEPS, config.SI_PARTICLES, SiMass);
	}
	if (config.XE_PARTICLES > 0)
	{
		correct(lists.xeParticles, config.TIMESTEPS, config.XE_PARTICLES, XeMass);
	}

	iterationNum++;
}

void Animate(int argc, char *argv[])
{
	xp = new GLfloat[config.SI_PARTICLES+config.XE_PARTICLES];
	yp = new GLfloat[config.SI_PARTICLES+config.XE_PARTICLES];
	zp = new GLfloat[config.SI_PARTICLES+config.XE_PARTICLES];
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(950, 650);
	glutInitWindowPosition(50, 50);
	glutCreateWindow("Molecular Dynamics");
	Initialize();
	glutDisplayFunc(Draw);
	glutMainLoop();
	delete[] xp;
	delete[] yp;
	delete[] zp;
}

void Initialize()
{
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);
}

void Draw()
{
	glClear(GL_COLOR_BUFFER_BIT);
	DrawPoints();
	glFlush();

	calculatePositions();

	glutPostRedisplay();
}

void DrawPoints()
{
	for (place = 0; place < config.SI_PARTICLES; place++)
	{
		xp[place] = lists.siParticles[place].position.x;
		yp[place] = lists.siParticles[place].position.y;
		zp[place] = lists.siParticles[place].position.z;
	}
	for (int i = 0; i < config.XE_PARTICLES; i++, place++)
	{
		xp[place] = (lists.xeParticles[i].position.x);
		yp[place] = (lists.xeParticles[i].position.y);
		zp[place] = (lists.xeParticles[i].position.z);
	}
	for (int i = 0; i < config.SI_PARTICLES; i++)
	{
		glColor3f(1.0, 1.0, 1.0);
		glPointSize(2.5);
		glBegin(GL_POINTS);
		glVertex3f(((xp[i] / resolution) + 0.15), ((zp[i] / resolution) + 0.25), ((yp[i] / resolution)) + 0.25);
		glEnd();
	}
	for (int i = config.SI_PARTICLES; i < config.SI_PARTICLES + config.XE_PARTICLES; i++)
	{
		glColor3f(0.246, 0.554, 0.51908);
		glPointSize(2.5);
		glBegin(GL_POINTS);
		glVertex3f(((xp[i] / resolution) + 0.15), ((zp[i] / resolution) + 0.25), ((yp[i] / resolution)) + 0.25);
		glEnd();
	}
	for (int i = 0; i < config.SI_PARTICLES; i++)
	{
		glColor3f(1.0, 1.0, 1.0);
		glPointSize(2.5);
		glBegin(GL_POINTS);
		glVertex3f(((xp[i] / resolution) + 0.55), ((yp[i] / resolution) + 0.25), ((zp[i] / resolution) + 0.25));
		glEnd();
	}
	for (int i = config.SI_PARTICLES; i < config.SI_PARTICLES + config.XE_PARTICLES; i++)
	{
		glColor3f(0.246, 0.554, 0.51908);
		glPointSize(2.5);
		glBegin(GL_POINTS);
		glVertex3f(((xp[i] / resolution) + 0.55), ((yp[i] / resolution) + 0.25), ((zp[i] / resolution) + 0.25));
		glEnd();
	}
}








































/*	/////////////////////////////////////// TESTS //////////////////////////////////////////
	file.open("C:\\Users\\User\\Desktop\\aaa\\VELO_TEST");
	file.precision(20);
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		file<<lists.siParticles[i].velocity.x<<", "<<lists.siParticles[i].velocity.y<<", "<<lists.siParticles[i].velocity.z<<"\n";
	}
	file.close();
	file.open("C:\\Users\\User\\Desktop\\aaa\\NEIGHBORS_TEST");
	file.precision(20);
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		file<<"("<<lists.siParticles[i].position.x<<","<<lists.siParticles[i].position.y<<","<<lists.siParticles[i].position.z<<") - ";
		for(int j = 0; j < config.MAX_SI_NEIGHBORS && lists.siParticles[i].siNeighbors[j] != -1; j++)
		{
			if(distance2(lists.siParticles[i].position, lists.siParticles[lists.siParticles[i].siNeighbors[j]].position)/sigma_Si <= a_Si)
			{
				file<<distance2(lists.siParticles[i].position, lists.siParticles[lists.siParticles[i].siNeighbors[j]].position)/sigma_Si<<", ";
			}
		}
		file<<endl;
	}
	file.close();
	file.open("C:\\Users\\User\\Desktop\\aaa\\POTENTIAL_TEST");
	file.precision(20);
	double vTotal = 0.0;
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		for(int j = i+1; j < config.SI_PARTICLES; j++)
		{
			vTotal += v2(lists.siParticles[i].position, lists.siParticles[j].position);
		}
	}
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		for(int j = i+1; j < config.SI_PARTICLES; j++)
		{
			for(int k = j+1; k < config.SI_PARTICLES; k++)
			{
					vTotal += v3(lists.siParticles[i].position, lists.siParticles[j].position, lists.siParticles[k].position);
			}
		}
	}
	file<<vTotal<<"\n";
	file.precision(20);
	vTotal = 0;
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		vTotal += (pow((sqrt((lists.siParticles[i].velocity.x*lists.siParticles[i].velocity.x)+(lists.siParticles[i].velocity.y*lists.siParticles[i].velocity.y)+(lists.siParticles[i].velocity.z*lists.siParticles[i].velocity.z))),2.0)*SiMass)/2;
	}
	file<<vTotal<<"\n";
	file.close();
	system("pause");//////////////TESTS!!
	////////////////////////////////////// END OF TESTS ////////////////////////////////////
	*/




///////////////////////////////////////// FORCES ////////////////////////////////////////////////
/*	file.open("C:\\Users\\User\\Desktop\\aaa\\siFo");
	file.precision(20);
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		file<<lists.siParticles[i].force.x<<", ";
	}
	file<<endl;
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		file<<lists.siParticles[i].force.y<<", ";
	}
	file<<endl;
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		file<<lists.siParticles[i].force.z<<", ";
	}
	file<<endl;
	file.close();
	file.open("C:\\Users\\User\\Desktop\\aaa\\xeFo");
	file.precision(20);
	for(int i = 0; i < config.XE_PARTICLES; i++)
	{
		file<<lists.xeParticles[i].force.x<<", ";
	}
	file<<endl;
	for(int i = 0; i < config.XE_PARTICLES; i++)
	{
		file<<lists.xeParticles[i].force.y<<", ";
	}
	file<<endl;
	for(int i = 0; i < config.XE_PARTICLES; i++)
	{
		file<<lists.xeParticles[i].force.z<<", ";
	}
	file<<endl;
	file.close();
	system("pause");
	*/
//////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////// NEIGHBORS DISTANCE ////////////////////////////////////
/*	file.open("C:\\Users\\User\\Desktop\\aaa\\siNsi");
	file.precision(19);
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		for(int j = 0; j < config.MAX_SI_NEIGHBORS && lists.siParticles[i].siNeighbors[j] != -1; j++)
		{
			file<<distance2(lists.siParticles[i].position, (lists.siParticles[lists.siParticles[i].siNeighbors[j]]).position)/sigma_Si<<", ";
		}
	}
	file<<endl;
	file.close();

	file.open("C:\\Users\\User\\Desktop\\aaa\\siNxe");
	file.precision(19);
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		for(int j = 0; j < config.MAX_XE_NEIGHBORS && lists.siParticles[i].xeNeighbors[j] != -1; j++)
		{
			file<<distance2(lists.siParticles[i].position, (lists.xeParticles[lists.siParticles[i].xeNeighbors[j]]).position)/sigma_Si_Xe<<", ";
		}
	}
	file<<endl;
	file.close();

	file.open("C:\\Users\\User\\Desktop\\aaa\\xeNxe");
	file.precision(19);
	for(int i = 0; i < config.XE_PARTICLES; i++)
	{
		for(int j = 0; j < config.MAX_XE_NEIGHBORS && lists.xeParticles[i].xeNeighbors[j] != -1; j++)
		{
			file<<distance2(lists.xeParticles[i].position, (lists.xeParticles[lists.xeParticles[i].xeNeighbors[j]]).position)/sigma_Xe_Xe<<", ";
		}
	}
	file<<endl;
	file.close();

	file.open("C:\\Users\\User\\Desktop\\aaa\\xeNsi");
	file.precision(19);
	for(int i = 0; i < config.XE_PARTICLES; i++)
	{
		for(int j = 0; j < config.MAX_SI_NEIGHBORS && lists.xeParticles[i].siNeighbors[j] != -1; j++)
		{
			file<<distance2(lists.xeParticles[i].position, (lists.siParticles[lists.xeParticles[i].siNeighbors[j]]).position)/sigma_Si_Xe<<", ";
		}
	}
	file<<endl;
	file.close();
	system("pause");
	*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////// MAX NEIGHBORS ////////////////////////////////////////
/*	int count = 0;
	int maxi = 0;
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		count = 0;
		for(int j = 0; j < config.MAX_SI_NEIGHBORS && lists.siParticles[i].siNeighbors[j] != -1; j++)
		{
			count++;
		}
		if(maxi < count)
		{
			maxi = count;
		}
	}
	cout<<config.SI_PARTICLES<<endl;
	cout<<maxi<<endl;
	system("pause");
	count = 0;
	maxi = 0;
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		count = 0;
		for(int j = 0; j < config.MAX_XE_NEIGHBORS && lists.siParticles[i].xeNeighbors[j] != -1; j++)
		{
			count++;
		}
		if(maxi < count)
		{
			maxi = count;
		}
	}
	cout<<maxi<<endl;
	system("pause");
	count = 0;
	maxi = 0;
	for(int i = 0; i < config.XE_PARTICLES; i++)
	{
		count = 0;
		for(int j = 0; j < config.MAX_XE_NEIGHBORS && lists.xeParticles[i].xeNeighbors[j] != -1; j++)
		{
			count++;
		}
		if(maxi < count)
		{
			maxi = count;
		}
	}
	cout<<config.XE_PARTICLES<<endl;
	cout<<maxi<<endl;
	system("pause");
	count = 0;
	maxi = 0;
	for(int i = 0; i < config.XE_PARTICLES; i++)
	{
		count = 0;
		for(int j = 0; j < config.MAX_SI_NEIGHBORS && lists.xeParticles[i].siNeighbors[j] != -1; j++)
		{
			count++;
		}
		if(maxi < count)
		{
			maxi = count;
		}
	}
	cout<<maxi<<endl;
	system("pause");
	*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////// VELOCITIES ///////////////////////////////////////////
/*	outputFile = argv[2];
	outputFile.append("_Velocity");
	file.open(outputFile);
	double v = 0.0;
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		file<<sqrt((lists.siParticles[i].velocity.x*lists.siParticles[i].velocity.x)+(lists.siParticles[i].velocity.y*lists.siParticles[i].velocity.y)+(lists.siParticles[i].velocity.z*lists.siParticles[i].velocity.z));
		if(i < config.SI_PARTICLES-1)
		{
			file<<", ";
		}
	}
	file.close();
	system("pause");
*/

/*	outputFile = argv[2];
	outputFile.append("_Velocity_X");
	file.open(outputFile);
	double v = 0.0;
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		file<<lists.siParticles[i].velocity.x<<endl;
	}
	file.close();
	outputFile = argv[2];
	outputFile.append("_Velocity_Y");
	file.open(outputFile);
	v = 0.0;
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		file<<lists.siParticles[i].velocity.y<<endl;
	}
	file.close();
	outputFile = argv[2];
	outputFile.append("_Velocity_Z");
	file.open(outputFile);
	v = 0.0;
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		file<<lists.siParticles[i].velocity.z<<endl;
	}
	file.close();
	cout<<"ZEU"<<endl;
	system("pause");*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*

	outputFile = argv[2];
	outputFile.append("_Velocity_X");
	file.open(outputFile);
	double v = 0.0;
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		file<<lists.siParticles[i].velocity.x<<", ";
	}
	file.close();
	outputFile = argv[2];
	outputFile.append("_Velocity_Z_SI");
	file.open(outputFile);
	file.precision(19);
	v = 0.0;
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		file<<lists.siParticles[i].velocity.z<<", ";
	}
	file.close();
	outputFile = argv[2];
	outputFile.append("_Velocity_Z_XE");
	file.open(outputFile);
	file.precision(19);
	v = 0.0;
	for(int i = 0; i < config.XE_PARTICLES; i++)
	{
		file<<lists.xeParticles[i].velocity.z<<", ";
	}
	file.close();
	cout<<"ZEU"<<endl;
	system("pause");
*/

/*
	ofstream file;
	file.open("C:\\Users\\User\\Desktop\\פרוייקט\\images\\file");
	if(!file.good())
	{
		cout<<"Can not open simulations input file"<<endl;
		system("pause");
	}
	file.precision(19);
	for(int i = 0; i < config.XE_PARTICLES; i++)
	{
		file<<lists.xeParticles[i].force.x<<", ";
	}
	file<<endl;
	for(int i = 0; i < config.XE_PARTICLES; i++)
	{
		file<<lists.xeParticles[i].force.y<<", ";
	}
	file<<endl;
	for(int i = 0; i < config.XE_PARTICLES; i++)
	{
		file<<lists.xeParticles[i].force.z<<", ";
	}
	file.close();
	*/


	//////////////////////////////////// PRINT FORCES /////////////////////////////////////
	/*
		file.open("C:\\Users\\User\\Desktop\\aaa\\fSi2");
	j.x = 0.7*sigma_Si;
	for(int k = 0; k < 500; k++)
	{
		if(distance2(i,j)/sigma_Si >= a_Si)
			file<<0;
		else
			file<<v2_derivative_of_rix(i,j,distance2(i,j))/epsilon_Si;
		if(k < 499)
		{
			file<<", ";
		}
		j.x += 0.01*sigma_Si;
	}
	file<<endl;
	file.close();

	file.open("C:\\Users\\User\\Desktop\\aaa\\Si2X");
	j.x = 0.7;
	for(int k = 0; k < 500; k++)
	{
		file<<distance2(i,j);
		if(k < 499)
		{
			file<<", ";
		}
		j.x += 0.01;
	}
	file.close();
	system("pause");
*/	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*	coordinate3 m;
	coordinate3 i;
	coordinate3 j;
	m.x = 0.0;
	m.y = 0.0;
	m.z = 0.0;
	j.x = pow(2.0,(1.0/6.0))*sigma_Si;
	j.y = 0.0;
	j.z = 0.0;
	i.z = 0.0;
	outputFile = argv[2];
	outputFile.append("FX_SI_SI_SI");
	file.open(outputFile);
	double cosal;
	for(int k = -20; k < 20; k++)
	{
		cosal = (1.0/3.0)+((double)k/400.0);
		i.x = -cosal*pow(2.0,(1.0/6.0))*sigma_Si;
		i.y = sqrt(1.0-(cosal*cosal))*pow(2.0,(1.0/6.0))*sigma_Si;
//		if(distance2(i,m)/sigma_Si >= a_Si || distance2(m,j)/sigma_Si >= a_Si)
//			file<<0;
//		else
			file<<((v3_derivative_of_rix(m, j, i, distance2(m,j), distance2(m,i),distance2(i,j)))+(v3_derivative_of_rix(m, i, j, distance2(m,i), distance2(m,j),distance2(i,j))))/epsilon_Si;
		if(k < 19)
		{
			file<<", ";
		}
	}
	file<<endl;
	for(int k = -20; k < 20; k++)
	{
		cosal = ((1.0/3.0)+((double)k/400.0));
		file<<cosal;
		if(k < 19)
		{
			file<<", ";
		}
	}
	file.close();

	m.x = 0.0;
	m.y = 0.0;
	m.z = 0.0;
	j.x = pow(2.0,(1.0/6.0))*sigma_Si;
	j.y = 0.0;
	j.z = 0.0;
	i.z = 0.0;
	outputFile = argv[2];
	outputFile.append("V_SI_SI_SI");
	file.open(outputFile);
	cosal;
	for(int k = -20; k < 20; k++)
	{
		cosal = (1.0/3.0)+((double)k/400.0);
		i.x = -cosal*pow(2.0,(1.0/6.0))*sigma_Si;
		i.y = sqrt(1.0-(cosal*cosal))*pow(2.0,(1.0/6.0))*sigma_Si;
//		if(distance2(m,j)/sigma_Si >= a_Si && distance2(m,i)/sigma_Si >= a_Si)
	//		file<<0;
	//	else
		file<<v3(m,i,j)/sigma_Si;
		if(k < 19)
		{
			file<<", ";
		}
	}
	file<<endl;
	for(int k = -20; k < 20; k++)
	{
		cosal = ((1.0/3.0)+((double)k/400.0));
		file<<(cosal);
		if(k < 19)
		{
			file<<", ";
		}
	}
	file.close();*/
//////////////////////////////////////////////////////////////////////////////////////////////


/*	coordinate3 i;
	coordinate3 j;
	j.x = 0.95*sigma_Si;
	j.y = 0.0;
	j.z = 0.0;
	i.x = 0.0;
	i.y = 0.0;
	i.z = 0.0;
	outputFile = argv[2];
	outputFile.append("FX_MORSE_SI");
	file.open(outputFile);
	file.precision(20);
	for(int k = 0; k < 350; k++)
	{
		file<<morseForce(distance2(i,j))/epsilon_Si<<", ";
		j.x += 0.01;
	}
	file<<endl;
	j.x = 0.95;
	for(int k = 0; k < 350; k++)
	{
		file<<j.x<<", ";
		j.x += 0.01;
	}
	file.close();

	j.x = 0.95*sigma_Si;
	outputFile = argv[2];
	outputFile.append("V_MORSE_SI");
	file.open(outputFile);
	for(int k = 0; k < 350; k++)
	{
		file<<morsePotential(i,j)/sigma_Si<<", ";
		j.x += 0.01;
	}
	file<<endl;
	j.x = 0.95;
	for(int k = 0; k < 350; k++)
	{
		file<<j.x<<", ";
		j.x += 0.01;
	}
	file.close();
	system("pause");*/

		///////////////////////////////////// POTENTIAL AND FORCES /////////////////////////////////////////////////
/*	coordinate3 i;
	coordinate3 j;
	i.x = 0.0;
	i.y = 0.0;
	i.z = 0.0;
	j.x = 0.0;
	j.y = 0.0;
	j.z = 0.0;

	file.open("C:\\Users\\User\\Desktop\\aaa\\morseF");
	file.precision(20);
	j.x = 0.7*sigma_Si;
	for(int k = 0; k < 5000; k++)
	{
		file<<(morseForce(distance2(i,j)))/epsilon_Si<<", ";
		j.x += 0.001*sigma_Si;
	}
	file<<endl;
	j.x = 0.7;
	for(int k = 0; k < 5000; k++)
	{
		file<<j.x<<", ";
		j.x += 0.001;
	}
	file.close();
	system("pause");

	cout.precision(20);
	j.x = (pow(2.0,1.0/6.0))*sigma_Si;
	cout<<j.x<<endl;
	cout<<morsePotential(i,j)/epsilon_Si<<endl;
	j.x = 2.3516702374129673334013950923828;
	cout<<v2(i,j)/epsilon_Si<<endl;
	cout<<v2_derivative_of_rix(i,j,distance2(i,j))/epsilon_Si<<endl;
//	while(v2_derivative_of_rix(i,j,distance2(i,j))/epsilon_Si < 0.0)
//	{
//		cout<<v2_derivative_of_rix(i,j,distance2(i,j))/epsilon_Si<<endl;
//		j.x -= 0.000000000000001;
//	}
//	cout<<j.x<<endl;
//	cout<<v2_derivative_of_rix(i,j,distance2(i,j))/epsilon_Si<<endl;
	j.x = 1.1224620483093729814335330496792*sigma_Xe_Xe;
	cout<<lennardJonesPotential(distance2(i,j),sigma_Xe_Xe,epsilon_Xe_Xe)/epsilon_Xe_Xe<<endl;
	cout<<lennardJonesForce(distance2(i,j),sigma_Xe_Xe,epsilon_Xe_Xe)/epsilon_Xe_Xe<<endl;
	system("pause");

	file.open("C:\\Users\\User\\Desktop\\aaa\\vSi2");
	file.precision(20);
	j.x = 0.7*sigma_Si;
	for(int k = 0; k < 500; k++)
	{
		if(distance2(i,j)/sigma_Si >= a_Si)
			file<<0;
		else
			file<<v2(i,j)/epsilon_Si;
		if(k < 499)
		{
			file<<", ";
		}
		j.x += 0.01*sigma_Si;
	}
	file<<endl;
	file.close();*/
	/////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////// NEIGHBBORS /////////////////////////////////////////
/*	int nei = 0;
	double dist = 0.0;
	double mindist = INT_MAX;
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		for(int j = i+1; j < config.SI_PARTICLES; j++)
		{
//			dist = distance2(lists.siParticles[i].position, lists.siParticles[j].position);
			if(mindist > dist)
			{
				mindist = dist;
			}
			if(dist <= 2.4516702374129673334013950923828)
			{
//				cout<<dist/sigma_Si<<endl;
//				system("pause");
			}
		}
//		cout<<"/////////////////"<<endl;
	}
	cout<<mindist<<endl;
	system("pause");
	for(int i = 0; i < config.XE_PARTICLES; i++)
	{
		nei = 0;
		for(int j = 0; j < config.MAX_XE_NEIGHBORS && lists.xeParticles[i].xeNeighbors[j] != -1; j++)
		{
//			cout<<distance2(lists.siParticles[i].position, lists.siParticles[lists.siParticles[i].siNeighbors[j]].position)/sigma_Si<<endl;
			if(distance2(lists.xeParticles[i].position, lists.xeParticles[lists.xeParticles[i].xeNeighbors[j]].position)/sigma_Xe_Xe < a_Si_Xe)
				nei++;
		}
		cout<<nei<<endl;
	}
	system("pause");
	double four = 0;
	double  two = 0;
	double one = 0;
	double zero = 0;
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		nei = 0;
		for(int j = 0; j < min(config.MAX_SI_NEIGHBORS,config.SI_PARTICLES) && lists.siParticles[i].siNeighbors[j] != -1; j++)
		{
//			cout<<(lists.siParticles[lists.siParticles[i].siNeighbors[j]]).position.x<<" "<<(lists.siParticles[lists.siParticles[i].siNeighbors[j]]).position.y<<" "<<(lists.siParticles[lists.siParticles[i].siNeighbors[j]]).position.z<<endl;
//			cout<<distance2(lists.siParticles[i].position, lists.siParticles[lists.siParticles[i].siNeighbors[j]].position)/sigma_Si<<endl;
			if(distance2(lists.siParticles[i].position, lists.siParticles[lists.siParticles[i].siNeighbors[j]].position)/sigma_Si < a_Si)
			{
				nei++;
			}
		}
		if(nei == 4)
		{
			four++;
		}
		else if(nei == 2)
		{
			two++;
		}
		else if(nei == 1)
		{
			one++;
		}
		else if(nei == 0)
		{
			zero++;
		}
		cout<<nei<<endl;
//		system("pause");
	}
	cout<<four<<endl;
	cout<<two<<endl;
	cout<<one<<endl;
	cout<<zero<<endl;
	system("pause");

//	V_total(config.SI_PARTICLES, lists.siParticles, config.XE_PARTICLES, lists.xeParticles);
//	system("pause");*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////// VELOCITY ////////////////////////////////////////////
/*	file.open("C:\\Users\\User\\Desktop\\aaa\\_Velocity_X_50rb_6400at");
	file.precision(20);
	double v = 0.0;
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		file<<lists.siParticles[i].velocity.x<<endl;
	}
	file.close();
	file.open("C:\\Users\\User\\Desktop\\aaa\\_Velocity_Y_50rb_6400at");
	v = 0.0;
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		file<<lists.siParticles[i].velocity.y<<endl;
	}
	file.close();
	file.open("C:\\Users\\User\\Desktop\\aaa\\_Velocity_Z_50rb_6400at");
	v = 0.0;
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		file<<lists.siParticles[i].velocity.z<<endl;
	}
	file.close();
	cout<<"ZEU1"<<endl;

	temp = 0.000001*rb;
	chooseVelocities(lists.siParticles, config.SI_PARTICLES, temp, &config, SiMass);
//	chooseVelocities(lists.xeParticles, config.XE_PARTICLES, temp, &config, XeMass);
	file.open("C:\\Users\\User\\Desktop\\aaa\\_Velocity_X_0.000001rb_6400at");
	v = 0.0;
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		file<<lists.siParticles[i].velocity.x<<endl;
	}
	file.close();
	file.open("C:\\Users\\User\\Desktop\\aaa\\_Velocity_Y_0.000001rb_6400at");
	v = 0.0;
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		file<<lists.siParticles[i].velocity.y<<endl;
	}
	file.close();
	file.open("C:\\Users\\User\\Desktop\\aaa\\_Velocity_Z_0.000001rb_6400at");
	v = 0.0;
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		file<<lists.siParticles[i].velocity.z<<endl;
	}
	file.close();
	file.close();
	cout<<"ZEU2"<<endl;
	system("pause");
	*/
	//////////////////////////////////////////////////////////////////////////////////////////
/*	file.open("C:\\Users\\User\\Desktop\\aaa\\NEIGHBORS_TEST_XE");
	file.precision(20);
	for(int i = 0; i < config.XE_PARTICLES; i++)
	{
		file<<"("<<lists.xeParticles[i].position.x<<","<<lists.xeParticles[i].position.y<<","<<lists.xeParticles[i].position.z<<") - ";
		for(int j = 0; j < config.MAX_XE_NEIGHBORS && lists.xeParticles[i].xeNeighbors[j] != -1; j++)
		{
			if(distance2(lists.xeParticles[i].position, lists.xeParticles[lists.xeParticles[i].xeNeighbors[j]].position)/sigma_Xe_Xe <= a_Si_Xe)
			{
				file<<distance2(lists.xeParticles[i].position, lists.xeParticles[lists.xeParticles[i].xeNeighbors[j]].position)/sigma_Xe_Xe<<", ";
			}
		}
		file<<endl;
	}
	file.close();
	*/


/*	file.open("C:\\Users\\User\\Desktop\\aaa\\NEIGHBORS_TEST_SI");
	file.precision(20);
	int n = 0;
	for(int i = 0; i < config.SI_PARTICLES; i++)
	{
		n = 0;
		file<<"("<<lists.siParticles[i].position.x<<","<<lists.siParticles[i].position.y<<","<<lists.siParticles[i].position.z<<") - ";
		for(int j = 0; j < config.MAX_SI_NEIGHBORS && lists.siParticles[i].siNeighbors[j] != -1; j++)
		{
			if(distance2(lists.siParticles[i].position, lists.siParticles[lists.siParticles[i].siNeighbors[j]].position)/sigma_Si <= a_Si)
			{
				n++;
				file<<distance2(lists.siParticles[i].position, lists.siParticles[lists.siParticles[i].siNeighbors[j]].position)/sigma_Si<<", ";
			}
		}
		cout<<i<<": "<<n<<endl;
		file<<endl;
	}
	file.close();
	system("pause");
*/
