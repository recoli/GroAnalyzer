void substr(char *dest, const char *src, 
			unsigned int start, unsigned int cnt);
void readNAtomsFrames(const char *groFileName, long int *pNAtoms, long int *pNFrames);
void readGRO(FILE *groFile, long int nAtoms, Atom *atom, double *box);
