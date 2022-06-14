#ifndef _SPATIAL_GAMES_H_
#define _SPATIAL_GAMES_H_

#define k 0.1

typedef struct{
	unsigned int L;
	unsigned int L2;
	double r;
	int seed;
	double tau;
	unsigned int M;
	double crit;
	double payoffs[2][2];
} Params;

typedef struct{
	unsigned int L;
	unsigned int L2;
	int *players;
	int (*neigh)[4];
	int (*links)[4];
	int coop_count;
	int link_count;
} Lattice;

typedef struct {
	int *n;
	int **hist;
} Memory;

void getParams(Params *p, int argc, char *argv[]);

double powerlawRng(Params p);

int int_powerlawRng(Params p);

void initLattice(Lattice *net, int L);

void updateCounts(Lattice *net);

void freeLattice(Lattice *net);

void printNet(Lattice net);

void initMemory(Memory *mem, Lattice net, Params p);

void freeMemory(Memory *mem, Params p);

void printMemory(Memory mem, Params p);

void linkAvaliation(Lattice *net, Memory mem, Params p, int t, int player);

double payoffCalculation(Lattice net, Params p, int player);

void strategyUpdate(Lattice *net, double Px, double Py, int player, int neigh);

void playGame(Lattice *net, Memory *mem, Params p, int t);

void diluteLinks(Lattice *net, Memory mem, Params p, int t);

#endif
