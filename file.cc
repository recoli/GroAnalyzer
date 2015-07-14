#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "typedef.h"

using namespace std;

// substr in c++
void substr(char *dest, const char *src, 
			unsigned int start, unsigned int cnt) 
{
   strncpy(dest, src + start, cnt + 1);
   dest[cnt] = '\0';
}


// read number of atoms and frames from gro file
void readNAtomsFrames(const char *groFileName, long int *pNAtoms, long int *pNFrames)
{
	// variables
	char line[MAX_STR_LEN];
	long int nAtoms;

	// open gro file
	FILE *groFile;
	groFile = fopen(groFileName, "r") ;
	if (NULL == groFile)
	{
		printf("Error: cannot open file traj.gro!\n");
		exit(1);
	}

	// count frames
	long int nFrames = 0;
	// first line: title
	while (fgets(line, sizeof(line), groFile) != NULL)
	{
		// second line: number of atoms
		fgets(line, sizeof(line), groFile);
		sscanf(line, "%ld", &nAtoms);

		// read atoms
		for (long int i = 0; i < nAtoms; ++ i)
		{
			fgets(line, sizeof(line), groFile);
		}

		// read box
		fgets(line, sizeof(line), groFile);

		++ nFrames;
	}

	// close gro file
	fclose(groFile);

	*pNFrames = nFrames;
	*pNAtoms  = nAtoms;
}


// read atomic information from gro file
void readGRO(FILE *groFile, long int nAtoms, Atom *atom, double *box)
{
	char line[MAX_STR_LEN], tmp[MAX_STR_LEN];
	fgets( line, sizeof( line ), groFile );
	fgets( line, sizeof( line ), groFile );
	for (long int i=0; i<nAtoms; i++) 
	{
		fgets( line, sizeof( line ), groFile );

		// residue ID
  		substr(tmp, line, 0, 5);
		sscanf(tmp, "%ld", &atom[i].resID);
		// atom ID
  		substr(tmp, line, 15, 5);
		sscanf(tmp, "%ld", &atom[i].atomID);

		// residue Name
  		substr(tmp, line, 5, 5);
		sscanf(tmp, "%s", atom[i].resName);
		// atom Name
  		substr(tmp, line, 10, 5);
		sscanf(tmp, "%s", atom[i].atomName);

		// atom coordinates
  		substr(tmp, line, 20, 120);
		sscanf(tmp, "%lf%lf%lf%lf%lf%lf",
				&atom[i].r.x, &atom[i].r.y, &atom[i].r.z,
				&atom[i].v.x, &atom[i].v.y, &atom[i].v.z);
	}
	fgets( line, sizeof( line ), groFile );
	sscanf( line, "%lf%lf%lf", &box[0], &box[1], &box[2]);
}

