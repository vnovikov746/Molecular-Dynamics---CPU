/*
 * Molecular Dynamics Project.
 * Outhor: Vladimir Novikov.
 */

#include "NeighborLists.h"
#include "Configurations.h"

void divideToCells(particleStruct *particles, int numParticles, float cellSize)
{
	for(int i = 0; i < numParticles; i++)
	{
		particles[i].cellPos.i = (int)((particles[i].position.x)/cellSize);
		particles[i].cellPos.j = (int)((particles[i].position.y)/cellSize);
		particles[i].cellPos.k = (int)((particles[i].position.z)/cellSize);
	}
}

void divideToCells(particleStruct *siParticles, particleStruct *xeParticles, int siNumParticles, int xeNumParticles, float cellSize)
{
	for(int i = 0; i < siNumParticles; i++)
	{
		siParticles[i].cellPos.i = (int)((siParticles[i].position.x)/cellSize);
		siParticles[i].cellPos.j = (int)((siParticles[i].position.y)/cellSize);
		siParticles[i].cellPos.k = (int)((siParticles[i].position.z)/cellSize);
	}
	for(int i = 0; i < xeNumParticles; i++)
	{
		xeParticles[i].cellPos.i = (int)((xeParticles[i].position.x)/cellSize);
		xeParticles[i].cellPos.j = (int)((xeParticles[i].position.y)/cellSize);
		xeParticles[i].cellPos.k = (int)((xeParticles[i].position.z)/cellSize);
	}
}

/*void buildNeighbors(particleStruct *particles, configurations *config, float cluster, float length, float width, float height, byte type, int numOfParticles)
{
	int NumOfCellsX = ((int)(length/cluster))+1;
	int NumOfCellsY = ((int)(width/cluster))+1;
	int NumOfCellsZ = ((int)(height/cluster))+1;

}*/

void buildNeighbors(particleStruct *particles, particleStruct *secondParticles, configurations *config, float cluster, float length, float width, float height, byte type, int numOfFirst, int numOfSecond)
{
	int NumOfCellsX = ((int)(length/cluster))+1;
	int NumOfCellsY = ((int)(width/cluster))+1;
	int NumOfCellsZ = ((int)(height/cluster))+1;

	int numOfCells = NumOfCellsX*NumOfCellsY*NumOfCellsZ;

	//-ADD T&C-//
	int *numOfParticlesInCell = new int[numOfCells];
	cell *bucketList = new cell[numOfCells];
	//---------//

	memset(numOfParticlesInCell, 0, sizeof(int) * numOfCells);
	for(int i = 0; i < numOfCells; i++)
	{
		bucketList[i].i = INT_MIN;
		bucketList[i].j = INT_MIN;
		bucketList[i].k = INT_MIN;
	}

	int bucket;
	int maxNumInCell = 0;

	//build the hash and counter the amount of particles for each cell
	if(type == Si || type == Xe)
	{
		for(int i = 0; i < numOfFirst; i++)
		{
			//get a bucket for the cell
			bucket = computeBucket(numOfCells, bucketList, particles[i].cellPos, false, NumOfCellsX, NumOfCellsY, NumOfCellsZ);

			if(bucket >= 0)
			{
				//count particle in the cell
				numOfParticlesInCell[bucket]++;

				//get the maximum particle per cell
				if(maxNumInCell < numOfParticlesInCell[bucket])
					maxNumInCell = numOfParticlesInCell[bucket];
			}
		}
	}

	else if(type == XeSi)
	{
		for(int i = 0; i < numOfSecond; i++)
		{
			//get a bucket for the cell
			bucket = computeBucket(numOfCells, bucketList, secondParticles[i].cellPos, false, NumOfCellsX, NumOfCellsY, NumOfCellsZ);

			if(bucket >= 0)
			{
				//count particle in the cell
				numOfParticlesInCell[bucket]++;

				//get the maximum particle per cell
				if(maxNumInCell < numOfParticlesInCell[bucket])
					maxNumInCell = numOfParticlesInCell[bucket];
			}
		}
    }

	if(maxNumInCell == 0)
		maxNumInCell = 1;

	if(type == XeSi)
	{
		if(config->MAX_SI_NEIGHBORS < maxNumInCell*27)
			config->MAX_SI_NEIGHBORS = maxNumInCell*27;
		if(config->MAX_XE_NEIGHBORS < maxNumInCell*27)
			config->MAX_XE_NEIGHBORS = maxNumInCell*27;
	
		for(int i = 0; i < numOfFirst; i++)
		{
			if(particles[i].type == Si)
			{
				particles[i].xeNeighbors = new int[maxNumInCell*27];
				memset(particles[i].xeNeighbors, -1, sizeof(int) * maxNumInCell*27);
			}
			else if(particles[i].type == Xe)
			{
				particles[i].siNeighbors = new int[maxNumInCell*27];
				memset(particles[i].siNeighbors, -1, sizeof(int) * maxNumInCell*27);
			}
		}
	}
	else if(type == Si)
	{
		if(config->MAX_SI_NEIGHBORS < maxNumInCell*27)
			config->MAX_SI_NEIGHBORS = maxNumInCell*27;
		for(int i = 0; i < numOfFirst; i++)
		{
			particles[i].siNeighbors = new int[maxNumInCell*27];
			memset(particles[i].siNeighbors, -1, sizeof(int) * maxNumInCell*27);
		}
	}
	else if(type == Xe)
	{
		if(config->MAX_XE_NEIGHBORS < maxNumInCell*27)
			config->MAX_XE_NEIGHBORS = maxNumInCell*27;
		for(int i = 0; i < numOfFirst; i++)
		{
			particles[i].xeNeighbors = new int[maxNumInCell*27];
			memset(particles[i].xeNeighbors, -1, sizeof(int) * maxNumInCell*27);
		}
	}


	//-ADD T&C-//
	int *list = new int[numOfCells*maxNumInCell];
	//---------//

	memset(list, -1, sizeof(int) * numOfCells*maxNumInCell);

	//collect the particles to the buckets
	if(type  == Si || type == Xe)
	{
		for(int i = 0; i < numOfFirst; i++)
		{
			//get cell bucket
			bucket = computeBucket(numOfCells, bucketList, particles[i].cellPos, true, NumOfCellsX, NumOfCellsY, NumOfCellsZ);

			if(bucket >= 0)
			{
				//add the particle to the corresponding cell bucket
				list[bucket*maxNumInCell+numOfParticlesInCell[bucket]-1] = i;
				numOfParticlesInCell[bucket]--; // counter == how many atoms in each bucket.
			}
		}
	}

	else if(type == XeSi)
	{
		for(int i = 0; i < numOfSecond; i++)
		{
			//get cell bucket
			bucket = computeBucket(numOfCells, bucketList, secondParticles[i].cellPos, true, NumOfCellsX, NumOfCellsY, NumOfCellsZ);

			if(bucket >= 0)
			{
				//add the particle to the corresponding cell bucket
				list[bucket*maxNumInCell+numOfParticlesInCell[bucket]-1] = i;
				numOfParticlesInCell[bucket]--; // counter == how many atoms in each bucket.
			}
		}
    }
	
	int ci;
	int cj;
	int ck;
	int neighborIndex;
	cell c;
	if(type  == Si || type == Xe)
	{
		for(int i = 0; i < numOfFirst; i++)
		{
			neighborIndex = 0;
			ci = particles[i].cellPos.i;
			cj = particles[i].cellPos.j;
			ck = particles[i].cellPos.k;
			for(int ici = ci-1; ici < ci+2; ici++)
			{
				if(ici < 0 || ici > (NumOfCellsX))
					continue;
				c.i = ici;
				for(int icj = cj-1; icj < cj+2; icj++)
				{
					if(icj < 0 || icj > (NumOfCellsY))
						continue;
					c.j = icj;
					for(int ick = ck-1; ick < ck+2; ick++)
					{
						//if the cell is outside the Si area
						if(ick < 0 || ick > NumOfCellsZ)
							continue;
						c.k = ick;

						//find the bucket
						bucket = computeBucket(numOfCells, bucketList, c, true, NumOfCellsX, NumOfCellsY, NumOfCellsZ);

						//no bucket - means cell is empty
						if (bucket == INT_MIN || bucket < 0)
							continue;

						//for each particle in the cell - add it to the particle's neighbour list
						for(int k = 0; k < maxNumInCell && list[bucket*maxNumInCell+k] != -1; k++)
						{
							//if it is the same particle
							if (i == list[bucket*maxNumInCell+k])
								continue;

							//add to the list and count
							if(type == Si)
								particles[i].siNeighbors[neighborIndex] = list[bucket*maxNumInCell+k];
							else if(type == Xe)
								particles[i].xeNeighbors[neighborIndex] = list[bucket*maxNumInCell+k];
							neighborIndex++;
						}
					}
				}
			}
		}
	}

	else if(type  == XeSi)
	{
		for(int i = 0; i < numOfFirst; i++)
		{
			neighborIndex = 0;
			ci = particles[i].cellPos.i;
			cj = particles[i].cellPos.j;
			ck = particles[i].cellPos.k;
			for(int ici = ci-1; ici < ci+2; ici++)
			{
				if(ici < 0 || ici > (NumOfCellsX))
					continue;
				c.i = ici;
				for(int icj = cj-1; icj < cj+2; icj++)
				{
					if(icj < 0 || icj > (NumOfCellsY))
						continue;
					c.j = icj;
					for(int ick = ck-1; ick < ck+2; ick++)
					{
						//if the cell is outside the Si area
						if(ick < 0 || ick > NumOfCellsZ)
							continue;
						c.k = ick;

						//find the bucket
						bucket = computeBucket(numOfCells, bucketList, c, true, NumOfCellsX, NumOfCellsY, NumOfCellsZ);

						//no bucket - means cell is empty
						if (bucket == INT_MIN || bucket < 0)
							continue;

						//for each particle in the cell - add it to the particle's neighbour list
						for(int k = 0; k < maxNumInCell && list[bucket*maxNumInCell+k] != -1; k++)
						{
							//if it is the same particle
							if (i == list[bucket*maxNumInCell+k])
								continue;

							//add to the list and count
							if(particles[i].type == Si)
								particles[i].xeNeighbors[neighborIndex] = list[bucket*maxNumInCell+k];
							else if(particles[i].type == Xe)
								particles[i].siNeighbors[neighborIndex] = list[bucket*maxNumInCell+k];
							neighborIndex++;
						}
					}
				}
			}
		}
	}	

	delete[] numOfParticlesInCell;
	delete[] bucketList;
	delete[] list;
}

/*
 * Computes the bucket
 */
int computeBucket(int numOfCells, cell *bucketList, cell currentCell, bool pop, int NumOfCellsX, int NumOfCellsY, int NumOfCellsZ)
{	
	//first function
	int n = H1*currentCell.i + H2*currentCell.j + H3*currentCell.k;
	n = n%numOfCells;
	if (n<0) n+=numOfCells;

	if (pop && bucketList[n].i == INT_MIN)
		return INT_MIN;

	if (bucketList[n].i == INT_MIN || (bucketList[n].i == currentCell.i && bucketList[n].j == currentCell.j && bucketList[n].k == currentCell.k))
	{
		if (!pop)
		{
			bucketList[n].i = currentCell.i;
			bucketList[n].j = currentCell.j;
			bucketList[n].k = currentCell.k;
		}

		return n;	
	}

	//if there is a collision - try second function
	n = H1*currentCell.i ^ H2*currentCell.j ^ H3*currentCell.k;
	n = n%numOfCells;
	if (n<0) n+=numOfCells;

	if (pop && bucketList[n].i == INT_MIN)
		return INT_MIN;
	
	//if there is a collision again - stupidly jump to the next cell
	while (bucketList[n].i != INT_MIN && (bucketList[n].i != currentCell.i || bucketList[n].j != currentCell.j || bucketList[n].k != currentCell.k))
	{
		n += 1;
		n = n%numOfCells;
		if (n<0) n+=numOfCells;
	}
		
	if (!pop)
	{
		bucketList[n].i = currentCell.i;
		bucketList[n].j = currentCell.j;
		bucketList[n].k = currentCell.k;
	}
	
	return n;		
}