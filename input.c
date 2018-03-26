/* 

@2009, Jonathan Belof

*/

#include <6dpes.h>


molecule_t *read_molecules(system_t *system) {

	molecule_t *molecules, *molecule_ptr;
	atom_t *atom_ptr, *prev_atom_ptr;
	char linebuf[MAXLINE], *n;
	FILE *fp;
	char token_atom[MAXLINE], token_atomid[MAXLINE], token_atomtype[MAXLINE], token_moleculetype[MAXLINE];
	char token_frozen[MAXLINE], token_moleculeid[MAXLINE], token_x[MAXLINE], token_y[MAXLINE], token_z[MAXLINE];
	char token_mass[MAXLINE], token_charge[MAXLINE], token_alpha[MAXLINE], token_epsilon[MAXLINE], token_sigma[MAXLINE];
	int current_frozen, current_adiabatic, current_moleculeid, current_atomid;
	double current_x, current_y, current_z;
	double current_mass, current_charge, current_alpha, current_epsilon, current_sigma;
	double current_molecule_mass;
	int moveable;
	int atom_counter;

	/* allocate the start of the list */
	molecules = calloc(1, sizeof(molecule_t));
	molecule_ptr = molecules;
	molecule_ptr->id = 1;
	molecule_ptr->atoms = calloc(1, sizeof(atom_t));
	atom_ptr = molecule_ptr->atoms;
	prev_atom_ptr = atom_ptr;

	/* open the molecule input file */
	fp = fopen(system->pdb_input, "r");
	if(!fp) {
		sprintf(linebuf, "INPUT: couldn't open PDB input file %s\n", system->pdb_input);
		fprintf(stderr, "%s", linebuf);	
		return(NULL);
	}

	/* clear the linebuffer and read the tokens in */
	atom_counter = 0;
	memset(linebuf, 0, MAXLINE);
	n = fgets(linebuf, MAXLINE, fp);
	while(n) {

		/* clear the tokens */
		memset(token_atom, 0, MAXLINE);
		memset(token_atomid, 0, MAXLINE);
		memset(token_atomtype, 0, MAXLINE);
		memset(token_moleculetype, 0, MAXLINE);
		memset(token_frozen, 0, MAXLINE);
		memset(token_moleculeid, 0, MAXLINE);
		memset(token_x, 0, MAXLINE);
		memset(token_y, 0, MAXLINE);
		memset(token_z, 0, MAXLINE);
		memset(token_mass, 0, MAXLINE);
		memset(token_charge, 0, MAXLINE);
		memset(token_alpha, 0, MAXLINE);
		memset(token_epsilon, 0, MAXLINE);
		memset(token_sigma, 0, MAXLINE);

		/* parse the line */
		sscanf(linebuf, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s\n", token_atom, token_atomid, token_atomtype, token_moleculetype, token_frozen, token_moleculeid, token_x, token_y, token_z, token_mass, token_charge, token_alpha, token_epsilon, token_sigma);

		if(!strcasecmp(token_atom, "ATOM") && strcasecmp(token_moleculetype, "BOX")) {

			current_frozen = 0; current_adiabatic = 0;
			if(!strcasecmp(token_frozen, "F"))
				current_frozen = 1;
			if(!strcasecmp(token_frozen, "A"))
				current_adiabatic = 1;

			current_moleculeid = atoi(token_moleculeid);
			current_atomid = atoi(token_atomid);
			current_x = atof(token_x);
			current_y = atof(token_y);
			current_z = atof(token_z);
			current_mass = atof(token_mass);	/* mass in amu */
			current_charge = atof(token_charge);
			current_charge *= E2REDUCED;		/* convert charge into reduced units */
			current_alpha = atof(token_alpha);
			current_epsilon = atof(token_epsilon);
			current_sigma = atof(token_sigma);

			if(molecule_ptr->id != current_moleculeid) {
				molecule_ptr->next = calloc(1, sizeof(molecule_t));
				molecule_ptr = molecule_ptr->next;
				molecule_ptr->atoms = calloc(1, sizeof(atom_t));
				prev_atom_ptr->next = NULL;
				free(atom_ptr);
				atom_ptr = molecule_ptr->atoms;
			}
			strcpy(molecule_ptr->moleculetype, token_moleculetype);

			molecule_ptr->id = current_moleculeid;
			molecule_ptr->frozen = current_frozen;
			molecule_ptr->adiabatic = current_adiabatic;
			molecule_ptr->mass += current_mass;

			++atom_counter;
			atom_ptr->id = atom_counter;
			memset(atom_ptr->atomtype, 0, MAXLINE);
			strcpy(atom_ptr->atomtype, token_atomtype);
			atom_ptr->frozen = current_frozen;
			atom_ptr->adiabatic = current_adiabatic;
			atom_ptr->pos[0] = current_x;
			atom_ptr->pos[1] = current_y;
			atom_ptr->pos[2] = current_z;
			atom_ptr->mass = current_mass;
			atom_ptr->charge = current_charge;
			atom_ptr->polarizability = current_alpha;
			atom_ptr->epsilon = current_epsilon;
			atom_ptr->sigma = current_sigma;

			atom_ptr->next = calloc(1, sizeof(atom_t));
			prev_atom_ptr = atom_ptr;
			atom_ptr = atom_ptr->next;

		}

		memset(linebuf, 0, MAXLINE);
		n = fgets(linebuf, MAXLINE, fp);
	}

	/* terminate the atom list */
	prev_atom_ptr->next = NULL;
	free(atom_ptr);

	/* scan the list, make sure that there is at least one moveable molecule */
	for(molecule_ptr = molecules, moveable = 0; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		if(!molecule_ptr->frozen) ++moveable;
	}

	if(!moveable) {
		fprintf(stderr, "INPUT: no moveable molecules found, there must be at least one in your PDB file\n");
		return(NULL);
	} else
		return(molecules);

}

system_t *setup_system(char *pdb_name, double b1x, double b1y, double b1z, double b2x, double b2y, double b2z, double b3x, double b3y, double b3z) {

	system_t *system;

	/* allocate the initial system structure */
	system = calloc(1, sizeof(system_t));
	system->pbc = calloc(1, sizeof(pbc_t));
	system->pdb_input = calloc(MAXLINE, sizeof(char));
	strcpy(system->pdb_input, pdb_name);

	/* read in the input pdb and setup the data structures */
	system->molecules = read_molecules(system);
	if(!system->molecules) {
		fprintf(stderr, "INPUT: error reading in molecules from PDB\n");
		return(NULL);
	}

	/* allocate the necessary pairs */
	setup_pairs(system->molecules);

	/* setup the basis matrix of the unit cell */
	system->pbc->basis[0][0] = b1x;
	system->pbc->basis[0][1] = b1y;
	system->pbc->basis[0][2] = b1z;
	system->pbc->basis[1][0] = b2x;
	system->pbc->basis[1][1] = b2y;
	system->pbc->basis[1][2] = b2z;
	system->pbc->basis[2][0] = b3x;
	system->pbc->basis[2][1] = b3y;
	system->pbc->basis[2][2] = b3z;

	/* calculate the periodic boundary conditions */
	pbc(system->pbc);

	/* get all of the pairwise interactions, exclusions, etc. */
	pairs(system);

	/* set the ewald gaussian width appropriately */
	system->ewald_alpha = 3.5/system->pbc->cutoff;
	system->ewald_kmax = EWALD_KMAX;

	return(system);

}

#ifdef DEBUG
void test_list(molecule_t *molecules) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	pair_t *pair_ptr;

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {

		printf("moleculeid = %d\n", molecule_ptr->id);
		printf("moleculetype = %s\n", molecule_ptr->moleculetype);
		printf("molecule_frozen = %d\n", molecule_ptr->frozen);
		printf("molecule_mass = %f\n", molecule_ptr->mass);
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			printf("atomtype = %s x = %f y = %f z = %f\n", atom_ptr->atomtype, atom_ptr->pos[0], atom_ptr->pos[1], atom_ptr->pos[2]);
			printf("atom frozen = %d mass = %f, charge = %f, alpha = %f, eps = %f, sig = %f\n", atom_ptr->frozen, atom_ptr->mass, atom_ptr->charge, atom_ptr->polarizability, atom_ptr->epsilon, atom_ptr->sigma);
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next)
				if(!(pair_ptr->rd_excluded || pair_ptr->es_excluded || pair_ptr->frozen)) printf("pair = 0x%lx eps = %f sig = %f\n", pair_ptr, pair_ptr->epsilon, pair_ptr->sigma);
		}

	}

	fflush(stdout);

}

void test_molecule(molecule_t *molecule) {

	atom_t *atom_ptr;
	pair_t *pair_ptr;

	printf("moleculeid = %d\n", molecule->id);
	printf("moleculetype = %s\n", molecule->moleculetype);
	printf("molecule_frozen = %d\n", molecule->frozen);
	printf("molecule_mass = %f\n", molecule->mass);
	for(atom_ptr = molecule->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
		printf("atomtype = %s x = %f y = %f z = %f\n", atom_ptr->atomtype, atom_ptr->pos[0], atom_ptr->pos[1], atom_ptr->pos[2]);
		printf("atom frozen = %d mass = %f, charge = %f, alpha = %f, eps = %f, sig = %f\n", atom_ptr->frozen, atom_ptr->mass, atom_ptr->charge, atom_ptr->polarizability, atom_ptr->epsilon, atom_ptr->sigma);
		for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {
			printf("pair at 0x%lx\n", pair_ptr);fflush(stdout);
		}
	}

printf("...finished\n");fflush(stdout);


}


int write_molecules(system_t *system, char *filename) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	char linebuf[MAXLINE];
	FILE *fp;
	int i, j;

	fp = fopen(filename, "w");
	if(!fp) {
		sprintf(linebuf, "OUTPUT: could not write pdb to file %s\n", filename);
		fprintf(stderr, "%s", linebuf);
		return(-1);
	}

	/* write pdb */
	for(molecule_ptr = system->molecules, i = 1, j = 1; molecule_ptr; molecule_ptr = molecule_ptr->next, j++) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next, i++) {

			fprintf(fp, "ATOM  ");
			fprintf(fp, "%5d", i);		/* give each one a unique id */
			fprintf(fp, " %-4.45s", atom_ptr->atomtype);
			fprintf(fp, " %-3.3s ", molecule_ptr->moleculetype);
			if(atom_ptr->adiabatic)
				fprintf(fp, "%-1.1s", "A");
			else if(atom_ptr->frozen)
				fprintf(fp, "%-1.1s", "F");
			else
				fprintf(fp, "%-1.1s", "M");
			fprintf(fp, "%4d    ", j);		/* give each molecule a unique id */
			fprintf(fp, "%8.3f", atom_ptr->pos[0]);
			fprintf(fp, "%8.3f", atom_ptr->pos[1]);
			fprintf(fp, "%8.3f", atom_ptr->pos[2]);
			fprintf(fp, " %8.4f", atom_ptr->mass);
			fprintf(fp, " %8.4f", atom_ptr->charge/E2REDUCED);	/* convert charge back to real units */
			fprintf(fp, " %8.5f", atom_ptr->polarizability);
			fprintf(fp, " %8.5f", atom_ptr->epsilon);
			fprintf(fp, " %8.5f", atom_ptr->sigma);
			fprintf(fp, "\n");

		}
	}

	fprintf(fp, "END\n");
	fflush(fp);

	fclose(fp);
	return(0);

}


void write_states(FILE *fp, molecule_t *molecules) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	char linebuf[MAXLINE];
	int i, j;

	/* write pdb formatted states */
	for(molecule_ptr = molecules, i = 1, j = 1; molecule_ptr; molecule_ptr = molecule_ptr->next, j++) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next, i++) {

			fprintf(fp, "ATOM  ");
			fprintf(fp, "%5d", i);		/* give each one a unique id */
			fprintf(fp, " %-4.45s", atom_ptr->atomtype);
			fprintf(fp, " %-3.3s ", molecule_ptr->moleculetype);
			if(atom_ptr->adiabatic)
				fprintf(fp, "%-1.1s", "A");
			else if(atom_ptr->frozen)
				fprintf(fp, "%-1.1s", "F");
			else
				fprintf(fp, "%-1.1s", "M");
			fprintf(fp, "%4d    ", j);		/* give each molecule a unique id */
			fprintf(fp, "%8.3f", atom_ptr->pos[0]);
			fprintf(fp, "%8.3f", atom_ptr->pos[1]);
			fprintf(fp, "%8.3f", atom_ptr->pos[2]);
			fprintf(fp, " %8.4f", atom_ptr->mass);
			fprintf(fp, " %8.4f", atom_ptr->charge/E2REDUCED);	/* convert charge back to real units */
			fprintf(fp, " %8.5f", atom_ptr->polarizability);
			fprintf(fp, " %8.5f", atom_ptr->epsilon);
			fprintf(fp, " %8.5f", atom_ptr->sigma);
			fprintf(fp, "\n");

		}
	}

	fprintf(fp, "ENDMDL\n");
	fflush(fp);


}

#endif /* DEBUG */


