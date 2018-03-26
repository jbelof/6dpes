/* 

@2009, Jonathan Belof

*/

#include <6dpes.h>

/************************************************/
/* loop over all pairs and setup the following: */
/*	displacements and minimum images	*/
/*	frozen interactions			*/
/*	excluded interactions			*/
/*	neighbor charges			*/
/*	mixed LJ parameters			*/
/************************************************/

void pairs(system_t *system) {

	int i, j, n;
	int p, q;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	pair_t *pair_ptr;
	atom_t **atom_array;
	molecule_t **molecule_array;
	double img[3];
	double d[3], r, r2;
	double di[3], ri, ri2;

	/* generate an array of atom ptrs */
	for(molecule_ptr = system->molecules, n = 0, atom_array = NULL, molecule_array = NULL; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next, n++) {
			molecule_array = realloc(molecule_array, sizeof(molecule_t *)*(n + 1));
			molecule_array[n] = molecule_ptr;

			atom_array = realloc(atom_array, sizeof(atom_t *)*(n + 1));
			atom_array[n] = atom_ptr;
		}
	}


	/* loop over all atoms and pair */
	for(i = 0; i < (n - 1); i++) {
		for(j = (i + 1), pair_ptr = atom_array[i]->pairs; j < n; j++, pair_ptr = pair_ptr->next) {


			/* recalculate exclusions */
			if(molecule_array[i] == molecule_array[j]) { /* if both on same molecule, exclude all interactions */

				pair_ptr->rd_excluded = 1;
				pair_ptr->es_excluded = 1;

			} else {

				/* exclude null repulsion/dispersion interations */
				if((atom_array[i]->epsilon == 0.0) || (atom_array[i]->sigma == 0.0) || (atom_array[j]->epsilon == 0.0) || (atom_array[j]->sigma == 0.0))
					pair_ptr->rd_excluded = 1;
				else
					pair_ptr->rd_excluded = 0;

				/* exclude null electrostatic interactions */
				if((atom_array[i]->charge == 0.0) || (atom_array[j]->charge == 0.0))
					pair_ptr->es_excluded = 1;
				else
					pair_ptr->es_excluded = 0;

			}

			/* get the frozen interactions */
			pair_ptr->frozen = atom_array[i]->frozen && atom_array[j]->frozen;

			/*** minimum image START ***/
			/* get the real displacement */
			for(p = 0; p < 3; p++) {
				d[p] = atom_array[i]->pos[p] - atom_array[j]->pos[p];
				pair_ptr->d[p] = d[p];
			}

			/* matrix multiply with the inverse basis and round */
			for(p = 0; p < 3; p++) {
				for(q = 0, img[p] = 0; q < 3; q++) {
					img[p] += system->pbc->reciprocal_basis[p][q]*d[q];
				}
				img[p] = rint(img[p]);
			}

			/* matrix multiply to project back into our basis */
			for(p = 0; p < 3; p++)
				for(q = 0, di[p] = 0; q < 3; q++)
					di[p] += system->pbc->basis[p][q]*img[q];

			/* now correct the displacement */
			for(p = 0; p < 3; p++)
				di[p] = d[p] - di[p];

			/* pythagorean terms */
			for(p = 0, r2 = 0, ri2 = 0; p < 3; p++) {
				r2 += d[p]*d[p];
				ri2 += di[p]*di[p];
			}


			r = sqrt(r2);
			ri = sqrt(ri2);

			/* store the results for this pair */
			pair_ptr->r = r;
			pair_ptr->rimg = ri;
			for(p = 0; p < 3; p++)
				pair_ptr->dimg[p] = di[p];

			/*** minimum image DONE ***/

			/* get the mixed LJ parameters */
			/* Lorentz-Berthelot mixing rules */
			pair_ptr->sigma = 0.5*(atom_array[i]->sigma + atom_array[j]->sigma);
			pair_ptr->epsilon = sqrt(atom_array[i]->epsilon*atom_array[j]->epsilon);

			/* get the neighbor charge */
			pair_ptr->charge = atom_array[j]->charge;

			/* set the link */
			pair_ptr->atom = atom_array[j];
			pair_ptr->molecule = molecule_array[j];


		} /* for j */
	} /* for i */


        /* free our temporary arrays */
        free(atom_array);
	free(molecule_array);

	/* update the com of each molecule */
	update_com(system->molecules);



}

/* molecular center of mass */
void update_com(molecule_t *molecules) {

	int i;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {

		for(i = 0; i < 3; i++)
			molecule_ptr->com[i] = 0;

		for(atom_ptr = molecule_ptr->atoms, molecule_ptr->mass = 0; atom_ptr; atom_ptr = atom_ptr->next) {

			molecule_ptr->mass += atom_ptr->mass;

			for(i = 0; i < 3; i++)
				molecule_ptr->com[i] += atom_ptr->mass*atom_ptr->pos[i];
		}

		for(i = 0; i < 3; i++)
			molecule_ptr->com[i] /= molecule_ptr->mass;

	}


}


/* allocate the pair lists */
void setup_pairs(molecule_t *molecules) {

	int i, j, n;
	molecule_t *molecule_ptr, **molecule_array;
	atom_t *atom_ptr, **atom_array;
	pair_t *pair_ptr, *prev_pair_ptr;

	/* generate an array of atom ptrs */
	for(molecule_ptr = molecules, n = 0, atom_array = NULL, molecule_array = NULL; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

			atom_array = realloc(atom_array, sizeof(atom_t *)*(n + 1));
			atom_array[n] = atom_ptr;

			molecule_array = realloc(molecule_array, sizeof(molecule_t *)*(n + 1));
			molecule_array[n] = molecule_ptr;

			++n;

		}
	}

	/* setup the pairs, lower triangular */
	for(i = 0; i < (n - 1); i++) {

		atom_array[i]->pairs = calloc(1, sizeof(pair_t));
		pair_ptr = atom_array[i]->pairs;
		prev_pair_ptr = pair_ptr;

		for(j = (i + 1); j < n; j++) {

			pair_ptr->next = calloc(1, sizeof(pair_t));
			prev_pair_ptr = pair_ptr;
			pair_ptr = pair_ptr->next;

		}

		prev_pair_ptr->next = NULL;
		free(pair_ptr);

	}


}

#ifdef DEBUG
void test_pairs(molecule_t *molecules) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	pair_t *pair_ptr;

	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if(!(pair_ptr->frozen || pair_ptr->rd_excluded || pair_ptr->es_excluded)) printf("%d: charge = %f, epsilon = %f, sigma = %f, r = %f, rimg = %f\n", atom_ptr->id, pair_ptr->charge, pair_ptr->epsilon, pair_ptr->sigma, pair_ptr->r, pair_ptr->rimg);fflush(stdout);

			}
		}
	}

}
#endif

