#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lat2eps.h>

///////// Function headers /////////
void plotlat(int *lat_pt, char *snpsht, int L);
void openlat(int L, int opt);
void printlat(int *lat_pt, char *snpsht, int L, int opt, int *bl, int *lpj);
////////////////////////////////////

/////////// Declarations ///////////
void plotlat(int *lat_pt, char *snpsht, int L){
	int i, j;
	FILE *gp;

	gp = popen("gnuplot", "w");
	fprintf(gp, "set term epscairo enhanced color\n");
	fprintf(gp, "set output \"%s\"\n", snpsht);
	fprintf(gp, "set size square\n");
	fprintf(gp, "unset key\n");
	fprintf(gp, "unset label\n");
	fprintf(gp, "unset tics\n");
	fprintf(gp, "plot \"-\" matrix w image\n");
	// gp plots down-up
	for(i = L - 1; i >= 0; i--){
		for(j = 0; j < L; j++){
			fprintf(gp, "%d ", lat_pt[i*L + j]);
		}
		fprintf(gp, "\n");
	}
	fprintf(gp, "e\n");
	fclose(gp);
}

void openlat(int L, int opt){

	// opens lattice for snapshots
	if(opt == 1){
		lat2eps_init(L, L);
		lat2eps_set_color(0,0xFFFFFF);	// white
		/* lat2eps_set_color(1,0x000000); // black */
		/* lat2eps_set_color(1,0xff7cc5); // some pink */
		/* lat2eps_set_color(1,0xeeacf2); // perfume */
		/* lat2eps_set_color(1,0xacf2d2); // teal? */
		/* lat2eps_set_color(1,0x54e4a2); // teal?? */
		lat2eps_set_color(1, 0x5a7a7f); // slate gray
		/* lat2eps_set_color(1,0xc97424); // ochre */
		/* lat2eps_set_color(1,0x81160c); // dark burgundy */
	}

	// opens lattice for #bonds/site
	if(opt == 2){
		lat2eps_init(L, L);
		lat2eps_set_color(0,0xFFFFFF);	// white
		lat2eps_set_color(1,0xc2d6ff);	// blue
		lat2eps_set_color(2,0x8fa3ec);
		lat2eps_set_color(3,0x5c70b9);
		lat2eps_set_color(4,0x293d86);

	}

	// opens lattice for clusters
	if(opt == 3){
		lat2eps_init(L, L);
		lat2eps_set_color(0,0xFFFFFF);
		lat2eps_set_color(1,0xf4a460);
		lat2eps_set_color(2,0xf7b1b6);
	}

	// opens lattice for site state + link shade
	if(opt == 4){
		lat2eps_init(L, L);

		// site state 1
		lat2eps_set_color(0,0xbdc9cb);
		lat2eps_set_color(1,0x9cafb2);
		lat2eps_set_color(2,0x7a9498);
		lat2eps_set_color(3,0x5a7a7f);
		lat2eps_set_color(4,0x3e5558);

		// site state 0
		lat2eps_set_color(10,0xc3bacf);
		lat2eps_set_color(11,0x8875a0);
		lat2eps_set_color(12,0x604780);
		lat2eps_set_color(13,0x391a61);
		lat2eps_set_color(14,0x271243);
	}
}

void printlat(int *lat_pt, char *snpsht, int L, int opt, int *bl, int *lpj){
	// prints the lattice to eps file
	int i;
	int L2 = L*L;
	
	// lattice
	if(opt == 1){
		for(i = 0; i < L2; i++){
			lat2eps_set_site(i % L, (i / L), lat_pt[i]);
		}
		lat2eps_gen_eps(snpsht, 0, 0, L, L, 1, 3);
		lat2eps_release();
	}

	// links per site
	if(opt == 2){
		for(i = 0; i < L2; i++){
			lat2eps_set_site(i % L, (i / L), lpj[i]);
		}
		lat2eps_gen_eps(snpsht, 0, 0, L, L, 1, 3);
		lat2eps_release();
	}

	// two biggest clusters
	if(opt == 3){
		for(i = 0; i < L2; i++){
			if(lat_pt[i] == bl[0]){
				lat2eps_set_site(i % L, (i / L), 1);
			} else if(lat_pt[i] == bl[1]){
				lat2eps_set_site(i % L, (i / L), 2);
			} else{
				lat2eps_set_site(i % L, (i / L), 0);
			}
		}
		lat2eps_gen_eps(snpsht, 0, 0, L, L, 1, 3);
		lat2eps_release();
	}

	// state + links shade
	if(opt == 4){
		for(i = 0; i < L2; i++){
			if(lat_pt[i] == 1){	// Coop
				lat2eps_set_site(i % L, (i / L), lpj[i]);
			}
			else{
				lat2eps_set_site(i % L, (i / L), 10 + lpj[i]);
			}
		}
		lat2eps_gen_eps(snpsht, 0, 0, L, L, 1, 3);
		lat2eps_release();
	}

}
