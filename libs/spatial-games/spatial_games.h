#ifndef _SPATIAL_GAMES_H_
#define _SPATIAL_GAMES_H_

typedef struct{
	unsigned int L;
	unsigned int L2;
	double r;
	int seed;
	double tau;
	unsigned int M;
	double crit;
} Params;

typedef struct{
	unsigned int L;
	unsigned int L2;
	int *players;
	int (*neigh)[4];
	int (*links)[4];
	int *n;
} Lattice;


void getParams(Params *p, int argc, char *argv[]);

double powerlawRng(Params p);

int int_powerlawRng(Params p);

void initLattice(Lattice *net, int L);

void freeLattice(Lattice *net, int L);

#endif
