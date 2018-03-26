/* 

@2009, Jonathan Belof

*/

#include <6dpes.h>

/* free a molecule and all of it's associated stuff */
void free_molecule(molecule_t *molecule) {

	molecule->next = NULL;

	free_pairs(molecule);
	free_atoms(molecule);
	free(molecule);

}


void free_pairs(molecule_t *molecules) {

	int i;
	pair_t **ptr_array;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	pair_t *pair_ptr;

	/* build an array of ptrs to be freed */
	for(molecule_ptr = molecules, i = 0, ptr_array = NULL; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				ptr_array = realloc(ptr_array, sizeof(pair_t *)*(i + 1));
				ptr_array[i] = pair_ptr;
				++i;

			}
		}
	}

	/* free the whole array of ptrs */
	for(--i; i >= 0; i--) free(ptr_array[i]);

	/* zero out the heads */
	for(molecule_ptr = molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
			atom_ptr->pairs = NULL;

	/* free our temporary array */
	if(ptr_array) free(ptr_array);

}

void free_atoms(molecule_t *molecules) {

	int i, n;
	atom_t **ptr_array;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;

	/* build the ptr array */
	for(molecule_ptr = molecules, i = 0, n = 0, ptr_array = NULL; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
			ptr_array = realloc(ptr_array, sizeof(atom_t *)*(i + 1));
			ptr_array[i] = atom_ptr;
			++i, ++n;
		}
	}

	/* free the whole array of ptrs */
	for(i = 0; i < n; i++) free(ptr_array[i]);

	/* free our temporary array */
	free(ptr_array);

}

void free_molecules(molecule_t *molecules) {

	int i;
	molecule_t **ptr_array = NULL;
	molecule_t *molecule_ptr;

	/* build the ptr array */
	for(molecule_ptr = molecules, i = 0; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		ptr_array = realloc(ptr_array, sizeof(molecule_t *)*(i + 1));
		ptr_array[i] = molecule_ptr;
		++i;
	}

	/* free the whole array of ptrs */
	for(--i; i >= 0; i--) free(ptr_array[i]);

	/* free our temporary array */
	free(ptr_array);

}

/* allocate the polarization matrix */
void allocate_matrix(system_t *system) {

	int i, N;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;

	/* count N */
	for(molecule_ptr = system->molecules, N = 0; molecule_ptr; molecule_ptr = molecule_ptr->next)
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
			++N;

	/* allocate */
	system->A_matrix = calloc(3*N, sizeof(double *));
	for(i = 0; i < 3*N; i++) system->A_matrix[i] = calloc(3*N, sizeof(double));

}


/* free the polarization matrix */
void free_matrix(system_t *system) {

	int i, N;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;

	/* count N */
	for(molecule_ptr = system->molecules, N = 0; molecule_ptr; molecule_ptr = molecule_ptr->next)
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
			++N;

	/* free */
	for(i = 0; i < N; i++) free(system->A_matrix[i]);
	free(system->A_matrix);

}

/* free all of our data structures */
void cleanup(system_t *system) {

	free_matrix(system);
	free_pairs(system->molecules);
	free_atoms(system->molecules);
	free_molecules(system->molecules);
	free(system->pdb_input);
	free(system);

}

