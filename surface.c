/* 

@2009, Jonathan Belof

*/


#include <6dpes.h>


/* generate the PES over the specified grid and output to stdout */
void surface_scan(system_t *system, double xi, double yi, double zi, double ri, double xf, double yf, double zf, double rf, double dx, double dy, double dz, double dtheta, double dphi, double dr, double efieldx, double efieldy, double efieldz) {

	int p;
	molecule_t *mol_ptr, *molecule_backup;
	atom_t *atom_ptr, *atom_backup;
	double x, y, z, theta, phi, r;
	energy_t *energy;
#ifdef DEBUG
	FILE *fp;
#endif /* DEBUG */

	/* print the comment line labeling the columns */
	//printf("# X Y Z THETA PHI R ENERGY BOND LJ ES POL\n");

#ifdef DEBUG
	fp = fopen("traj.pdb", "w");
#endif /* DEBUG */

	/* backup the molecular coordinates */
	for(mol_ptr = system->molecules; mol_ptr; mol_ptr = mol_ptr->next) {
		if(!mol_ptr->frozen) {
			molecule_backup = copy_molecule(mol_ptr);
			break;
		}
	}

	/* loop over the entire domain */
	for(x = xi; x <= xf; x += dx) {
		for(y = yi; y <= yf; y += dy) {
			for(z = zi; z <= zf; z += dz) {
				for(theta = 0.0; theta <= M_PI; theta += dtheta) {
					for(phi = 0.0; phi <= 2.0*M_PI; phi += dphi) {
						for(r = ri; r <= rf; r += dr) {

							/* restore the original state of the molecule */
							/* atomic coordinates... */
							for(atom_backup = molecule_backup->atoms, atom_ptr = mol_ptr->atoms; atom_backup; atom_backup = atom_backup->next, atom_ptr = atom_ptr->next)
								for(p = 0; p < 3; p++)
									atom_ptr->pos[p] = atom_backup->pos[p];

							/* ...and c.o.m. coordinates */
							for(p = 0; p < 3; p++) mol_ptr->com[p] = molecule_backup->com[p];

							energy = surface_point(system, x, y, z, theta, phi, r, efieldx, efieldy, efieldz);
							printf("%.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n", x, y, z, theta, phi, r, energy->total, energy->bond, energy->lj, energy->es, energy->pol);
							free(energy);

#ifdef DEBUG
							write_states(fp, system->molecules);
#endif /* DEBUG */

						} /* for r */
					} /* for phi */
				} /* for theta */
			} /* for z */
		} /* for y */
	} /* for x */

	/* free the backup structure */
	free_molecule(molecule_backup);

}


/* accepts arguments for the 5D domain point and returns the PES value there (in Kelvin) */
energy_t *surface_point(system_t *system, double x, double y, double z, double theta, double phi, double r, double efieldx, double efieldy, double efieldz) {

	molecule_t *mol_ptr, *interaction_point;
	energy_t *pes;

	/* locate the moveable interaction point */
	for(mol_ptr = system->molecules; mol_ptr; mol_ptr = mol_ptr->next)
		if(!mol_ptr->frozen) interaction_point = mol_ptr;

	/* translate the interaction site on the grid */
	translate(interaction_point, x, y, z);

	/* rotate the interaction site on the grid */
	rotate(interaction_point, theta, phi);

	/* stretch the molecule to have a bond distance of r */
	stretch(interaction_point, r);

	/* get the PES value */
	pes = energy(system, efieldx, efieldy, efieldz);

	return(pes);

}

/* translate the molecule the (x,y,z) */
void translate(molecule_t *molecule, double x, double y, double z) {

	atom_t *atom_ptr;

	/* preserve the geometry of the molecule being moved */
	for(atom_ptr = molecule->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
		atom_ptr->pos[0] += (x - molecule->com[0]);
		atom_ptr->pos[1] += (y - molecule->com[1]);
		atom_ptr->pos[2] += (z - molecule->com[2]);
	}

	/* update the c.o.m. to the new value */
	molecule->com[0] = x;
	molecule->com[1] = y;
	molecule->com[2] = z;

}

/* rotate the molecule to (theta, phi) */
/* theta is an angle about the Y axis, phi is about the Z axis */
void rotate(molecule_t *molecule, double theta, double phi) {

	atom_t *atom_ptr;
	double rotation_matrix[3][3];
	double com[3];
	int i, ii, n;
	double *new_coord_array;

	/* count the number of atoms in the molecule, and allocate new coords array */
	for(atom_ptr = molecule->atoms, n = 0; atom_ptr; atom_ptr = atom_ptr->next)
		++n;

	new_coord_array = calloc(n*3, sizeof(double));

	/* save the com coordinate */
	com[0] = molecule->com[0];
	com[1] = molecule->com[1];
	com[2] = molecule->com[2];

	/* translate the molecule to the origin */
	for(atom_ptr = molecule->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
		atom_ptr->pos[0] -= com[0];
		atom_ptr->pos[1] -= com[1];
		atom_ptr->pos[2] -= com[2];
	}

	rotation_matrix[0][0] = cos(phi)*cos(theta);
	rotation_matrix[0][1] = -sin(phi);
	rotation_matrix[0][2] = cos(phi)*(-sin(theta));
	rotation_matrix[1][0] = sin(phi)*cos(theta);
	rotation_matrix[1][1] = cos(phi);
	rotation_matrix[1][2] = sin(phi)*(-sin(theta));
	rotation_matrix[2][0] = sin(theta);
	rotation_matrix[2][1] = 0;
	rotation_matrix[2][2] = cos(theta);

	/* matrix multiply the ol' fashioned way */
	for(atom_ptr = molecule->atoms, i = 0; atom_ptr; atom_ptr = atom_ptr->next, i++) {

		ii = i*3;
		new_coord_array[ii+0] = rotation_matrix[0][0]*atom_ptr->pos[0] + rotation_matrix[0][1]*atom_ptr->pos[1] + rotation_matrix[0][2]*atom_ptr->pos[2];
		new_coord_array[ii+1] = rotation_matrix[1][0]*atom_ptr->pos[0] + rotation_matrix[1][1]*atom_ptr->pos[1] + rotation_matrix[1][2]*atom_ptr->pos[2];
		new_coord_array[ii+2] = rotation_matrix[2][0]*atom_ptr->pos[0] + rotation_matrix[2][1]*atom_ptr->pos[1] + rotation_matrix[2][2]*atom_ptr->pos[2];

	}


	/* set the new coordinates and then translate back from the origin */
	for(atom_ptr = molecule->atoms, i = 0; atom_ptr; atom_ptr = atom_ptr->next, i++) {

		ii = i*3;
		atom_ptr->pos[0] = new_coord_array[ii+0];
		atom_ptr->pos[1] = new_coord_array[ii+1];
		atom_ptr->pos[2] = new_coord_array[ii+2];

		atom_ptr->pos[0] += com[0];
		atom_ptr->pos[1] += com[1];
		atom_ptr->pos[2] += com[2];

        }

	/* free our temporary array */
	free(new_coord_array);


}


/* stretch out a molecule until it's bond distance is r */
void stretch(molecule_t *molecule, double r) {

	int p;
	double current_r, scaling_factor;
	atom_t *atom_ptr;
	double theta, phi;

	/* get our scaling factor */
	current_r = get_bond_distance(molecule);
	scaling_factor = r/current_r;

	/* scale each coordinate */
	for(atom_ptr = molecule->atoms; atom_ptr; atom_ptr = atom_ptr->next)
		for(p = 0; p < 3; p++)
			atom_ptr->pos[p] = (atom_ptr->pos[p] - molecule->com[p])*scaling_factor + molecule->com[p];

}

/* make an exact copy of src */
molecule_t *copy_molecule(molecule_t *src) {

	int i, j;
	molecule_t *dst;
	atom_t *atom_dst_ptr, *prev_atom_dst_ptr, *atom_src_ptr;
	pair_t *pair_dst_ptr, *prev_pair_dst_ptr, *pair_src_ptr;

	/* allocate the start of the new lists */
	dst = calloc(1, sizeof(molecule_t));
	/* copy molecule attributes */
	dst->id = src->id;
	strcpy(dst->moleculetype, src->moleculetype);
	dst->mass = src->mass;
	dst->frozen = src->frozen;
	dst->adiabatic = src->adiabatic;
	memcpy(dst->com, src->com, 3*sizeof(double));
	dst->next = NULL;


	/* new atoms list */
	dst->atoms = calloc(1, sizeof(atom_t));
	prev_atom_dst_ptr = dst->atoms;

	for(atom_dst_ptr = dst->atoms, atom_src_ptr = src->atoms; atom_src_ptr; atom_dst_ptr = atom_dst_ptr->next, atom_src_ptr = atom_src_ptr->next) {

		atom_dst_ptr->id = atom_src_ptr->id;
		strcpy(atom_dst_ptr->atomtype, atom_src_ptr->atomtype);
		atom_dst_ptr->frozen = atom_src_ptr->frozen;
		atom_dst_ptr->adiabatic = atom_src_ptr->adiabatic;
		atom_dst_ptr->mass = atom_src_ptr->mass;
		atom_dst_ptr->charge = atom_src_ptr->charge;
		atom_dst_ptr->polarizability = atom_src_ptr->polarizability;
		atom_dst_ptr->epsilon = atom_src_ptr->epsilon;
		atom_dst_ptr->sigma = atom_src_ptr->sigma;

		memcpy(atom_dst_ptr->pos, atom_src_ptr->pos, 3*sizeof(double));
		memcpy(atom_dst_ptr->ef_static, atom_src_ptr->ef_static, 3*sizeof(double));
		memcpy(atom_dst_ptr->ef_induced, atom_src_ptr->ef_induced, 3*sizeof(double));
		memcpy(atom_dst_ptr->mu, atom_src_ptr->mu, 3*sizeof(double));
		memcpy(atom_dst_ptr->old_mu, atom_src_ptr->old_mu, 3*sizeof(double));
		memcpy(atom_dst_ptr->new_mu, atom_src_ptr->new_mu, 3*sizeof(double));

		atom_dst_ptr->pairs = calloc(1, sizeof(pair_t));
		pair_dst_ptr = atom_dst_ptr->pairs;
		prev_pair_dst_ptr = pair_dst_ptr;
		for(pair_src_ptr = atom_src_ptr->pairs; pair_src_ptr; pair_src_ptr = pair_src_ptr->next) {

			pair_dst_ptr->frozen = pair_src_ptr->frozen;
			pair_dst_ptr->rd_excluded = pair_src_ptr->rd_excluded;
			pair_dst_ptr->es_excluded = pair_src_ptr->es_excluded;
			pair_dst_ptr->charge = pair_src_ptr->charge;
			pair_dst_ptr->epsilon = pair_src_ptr->epsilon;
			pair_dst_ptr->lrc = pair_src_ptr->lrc;
			pair_dst_ptr->sigma = pair_src_ptr->sigma;
			pair_dst_ptr->r = pair_src_ptr->r;
			pair_dst_ptr->rimg = pair_src_ptr->rimg;

			pair_dst_ptr->next = calloc(1, sizeof(pair_t));
			prev_pair_dst_ptr = pair_dst_ptr;
			pair_dst_ptr = pair_dst_ptr->next;

		}
		prev_pair_dst_ptr->next = NULL;
		free(pair_dst_ptr);
		/* handle an empty list */
		if(!atom_src_ptr->pairs) atom_dst_ptr->pairs = NULL;

		prev_atom_dst_ptr = atom_dst_ptr;
		atom_dst_ptr->next = calloc(1, sizeof(atom_t));
        }

	prev_atom_dst_ptr->next = NULL;
	free(atom_dst_ptr);


	return(dst);

}

