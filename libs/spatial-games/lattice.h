#ifndef _LATTICE_H_
#define _LATTICE_H_

void openlat(int L, int opt);

/* void printlat(int *lat_pt, char *snpsht, int L, int opt, int *bl, int *lpj); */
void printlat(Lattice net, HKstuff hk, int *lpp, char *snpsht, int opt);


#endif
