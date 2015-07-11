void substr(char *dest, const char *src, 
			unsigned int start, unsigned int cnt);
void readNAtoms(FILE *groFile, long int *pNAtoms);
void readGRO(FILE *groFile, long int nAtoms, Atom *atom, double *box);
