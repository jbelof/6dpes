/* 

@2009, Jonathan Belof

*/


#include <6dpes.h>


void usage(char *progname) {

	fprintf(stderr, "usage: %s <PDB filename> [<b1_x> <b1_y> <b1_z> <b2_x> <b2_y> <b2_z> <b3_x> <b3_y> <b3_z>] [<xi> <yi> <zi> <ri> <xf> <yf> <zf> <rf>] [<dx> <dy> <dz> <dtheta> <dphi> <dr>] [ efield_x efield_y efield_z ]\n", progname);
	fprintf(stderr, "\t[b1,...,b3] : the cartesian basis vectors (in A) of the unit cell, each b vector is a row of the basis matrix\n");
	fprintf(stderr, "\t[xi,...,zf,rf] : the initial and final xyz c.o.m. and r coordinates for the PES generation\n");
	fprintf(stderr, "\t[dx,...,dphi,dr] : the step size for c.o.m. coordinates, r, and for the spherical polar angles (theta angle of rotation around the Y axis, phi is the angle around the Z axis)\n");
	fprintf(stderr, "\t[efield_x,...efield_z] : apply an electric field to the system, in units of e/A^2\n");
	fprintf(stderr, "\tcolumnar output is [x,y,z,theta,phi,r,E,BOND,LJ,ES,POL)]\n");
	fprintf(stderr, "MOF5 example:\n");
	fprintf(stderr, "\t$ %s MOF5+BSS.pdb 25.669 0.0 0.0 0.0 25.669 0.0 0.0 0.0 25.669 -10.0 -10.0 -10.0 0.2 10.0 10.0 10.0 2.0 0.001 0.001 0.001 0.1 0.1 0.01 0.0 0.0 0.0\n", progname);
	fprintf(stderr, "would map the PES at each point within a subcube of the unit cell spanned by [-10,-10,-10]x[10,10,10] for all angles theta=0-pi,phi=0-2*pi and bond distance 0.5 to 1.0 A in 0.1 A increments\n");
	fprintf(stderr, "for a vibrational surface only:\n");
	fprintf(stderr, "\t$ %s MOF5+BSS.pdb 25.669 0.0 0.0 0.0 25.669 0.0 0.0 0.0 25.669 -10.0 -10.0 -10.0 0.2 -10.0 -10.0 -10.0 4.0 1000.0 1000.0 1000.0 1000.0 1000.0 0.001 0.0 0.0 0.0 | awk '{print $6 \" \" $7}'\n", progname);
	fprintf(stderr, "for the vibrational surface of an H2 in a 0.5 e/A^2 field transverse to it's axis:\n");
	fprintf(stderr, "\t$ %s BSSP.pdb 1000.0 0.0 0.0 0.0 1000.0 0.0 0.0 0.0 1000.0 -10.0 -10.0 -10.0 0.2 -10.0 -10.0 -10.0 4.0 1000.0 1000.0 1000.0 1000.0 1000.0 0.001 0.0 0.0 0.5 | awk '{print $6 \" \" $7}'\n", progname);
	fprintf(stderr, "or to evaluate the vib frequency (cm^-1) by second derivative (central difference) at the minimum:\n");
	fprintf(stderr, "\t$ %s BSSP.pdb 1000.0 0.0 0.0 0.0 1000.0 0.0 0.0 0.0 1000.0 -10.0 -10.0 -10.0 0.7419 -10.0 -10.0 -10.0 0.7421 1000.0 1000.0 1000.0 1000.0 1000.0 0.0001 0.00003 0.0 0.0 | awk '{print $7}' | column | awk '{}{printf(\"%%.16f\\n\",sqrt(($3+$1-2.0*$2)/1.0e-8)*42.796);}{}'\n", progname);
	fprintf(stderr, "or to evaluate the vib frequency (cm^-1) by second derivative at global min in MOF5:\n");
	fprintf(stderr, "\t$ %s MOF5+H2.globalmin.pdb 25.669 0.0 0.0 0.0 25.669 0.0 0.0 0.0 25.669 -4.176 -4.216 8.437 0.7419 -4.176 -4.216 8.437 0.7421 1000.0 1000.0 1000.0 1000.0 1000.0 0.0001 0.0 0.0 0.0 | awk '{print $7}' | column | awk '{}{printf(\"%%.16f\\n\",sqrt(($3+$1-2.0*$2)/1.0e-8)*42.796);}{}'\n", progname);

	exit(1);

}

int main(int argc, char **argv) {

	system_t *system;
	double b1x, b1y, b1z, b2x, b2y, b2z, b3x, b3y, b3z;
	double xi, yi, zi, ri, xf, yf, zf, rf;
	double dx, dy, dz, dtheta, dphi, dr;
	double efieldx, efieldy, efieldz;

	/* check args */
	if(argc != 28) usage(argv[0]);

	if(!argv[1]) {
		fprintf(stderr, "%s", "MAIN: invalid PDB file specified");
		exit(1);
	}

	/* get the basis vectors for the unit cell */
	b1x = atof(argv[2]); b1y = atof(argv[3]); b1z = atof(argv[4]);
	b2x = atof(argv[5]); b2y = atof(argv[6]); b2z = atof(argv[7]);
	b3x = atof(argv[8]); b3y = atof(argv[9]); b3z = atof(argv[10]);

	/* read the geometry and potential parameters and setup the data structures */
	system = setup_system(argv[1], b1x, b1y, b1z, b2x, b2y, b2z, b3x, b3y, b3z);
	if(!system) {
		fprintf(stderr, "%s", "MAIN: could not initialize the data structures\n");
		exit(1);
	}

	/* get the surface boundaries and increments */
	xi = atof(argv[11]); yi = atof(argv[12]); zi = atof(argv[13]); ri = atof(argv[14]);
	xf = atof(argv[15]); yf = atof(argv[16]); zf = atof(argv[17]); rf = atof(argv[18]);
	dx = atof(argv[19]); dy = atof(argv[20]); dz = atof(argv[21]); dtheta = atof(argv[22]); dphi = atof(argv[23]); dr = atof(argv[24]);

	/* get the applied field */
	efieldx = atof(argv[25])*E2REDUCED; efieldy = atof(argv[26])*E2REDUCED; efieldz = atof(argv[27])*E2REDUCED;

	/* initialize the polarization parameters */
	system->damp_type = POLAR_DAMPING_EXPONENTIAL;
	system->field_damp = POLAR_DAMPING_FIELD_WIDTH;
	system->polar_damp = POLAR_DAMPING_EXPONENTIAL_WIDTH;
	system->polar_precision = POLAR_PRECISION;
	/* allocate the A_matrx */
	allocate_matrix(system);

	/* calculate the PES at the requested points and print the energy (in K) to stdout */
	surface_scan(system, xi, yi, zi, ri, xf, yf, zf, rf, dx, dy, dz, dtheta, dphi, dr, efieldx, efieldy, efieldz);


	/* cleanup */
	cleanup(system);
	exit(0);

}

