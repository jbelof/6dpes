/* 

@2009, Jonathan Belof

*/

#include <6dpes.h>


/* get the induction energy */
double polar(system_t *system, double efieldx, double efieldy, double efieldz) {

	int i, num_iterations;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr;
	double potential;


	/* get the A matrix */
	thole_amatrix(system);

	/* calculate the field vectors */
	thole_field(system, efieldx, efieldy, efieldz);

	/* solve the SCF */
	num_iterations = thole_iterative(system);

	/* calculate the polarization energy as 1/2 mu*E */
	for(molecule_ptr = system->molecules, potential = 0; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {
#ifdef POLAR_SELF
			for(i = 0; i < 3; i++) potential += atom_ptr->mu[i]*(atom_ptr->ef_static[i] + atom_ptr->ef_static_self[i]);
#else
			for(i = 0; i < 3; i++) potential += atom_ptr->mu[i]*atom_ptr->ef_static[i];
#endif /* POLAR_SELF */
		}

	}
	potential *= -0.5;

	return(potential);

}


/* calculate the dipole field tensor */
void thole_amatrix(system_t *system) {

	int i, j, ii, jj, N, p, q;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr, **atom_array;
	pair_t *pair_ptr;
	double r, damping_term1, damping_term2;
	double r3, r5;
	double v, s;

	/* generate an array of atom ptrs */
	for(molecule_ptr = system->molecules, N = 0, atom_array = NULL; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next, N++) {

			atom_array = realloc(atom_array, sizeof(atom_t *)*(N + 1));
			atom_array[N] = atom_ptr;

		}
	}

	/* zero out the matrix */
	for(i = 0; i < 3*N; i++)
		for(j = 0; j < 3*N; j++)
			system->A_matrix[i][j] = 0;

	/* set the diagonal blocks */
	for(i = 0; i < N; i++) {
		ii = i*3;

		for(p = 0; p < 3; p++) {

			if(atom_array[i]->polarizability != 0.0)
				system->A_matrix[ii+p][ii+p] = 1.0/atom_array[i]->polarizability;
			else
				system->A_matrix[ii+p][ii+p] = MAXVALUE;

		}

	}

	/* calculate each Tij tensor component for each dipole pair */
	for(i = 0; i < (N - 1); i++) {
		ii = i*3;

		for(j = (i + 1), pair_ptr = atom_array[i]->pairs; j < N; j++, pair_ptr = pair_ptr->next) {
			jj = j*3;

			/* inverse displacements */
			r3 = pow(pair_ptr->rimg, -3.0);
			r5 = pow(pair_ptr->rimg, -5.0);


			/* set the damping function */
			damping_term1 = 1.0;	/* default to no damping */
			damping_term2 = 1.0;
			if(system->damp_type == POLAR_DAMPING_LINEAR) { /* linear damping */

				s = (system->polar_damp)*pow((atom_array[i]->polarizability*atom_array[j]->polarizability), (1.0/6.0));
				v = pair_ptr->rimg/s;

				if(pair_ptr->rimg < s) {
					damping_term1 = (4.0*v*v*v - 3.0*v*v*v*v);
					damping_term2 = v*v*v*v;
				} else {
					damping_term1 = 1.0;
					damping_term2 = 1.0;
				}

			} else if(system->damp_type == POLAR_DAMPING_EXPONENTIAL) { /* exponential damping */


				damping_term1 = 1.0 - exp(-system->polar_damp*pair_ptr->rimg)*(0.5*system->polar_damp*system->polar_damp*pow(pair_ptr->rimg,2.0) + system->polar_damp*pair_ptr->rimg + 1.0);
				damping_term2 = 1.0 - exp(-system->polar_damp*pair_ptr->rimg)*(system->polar_damp*system->polar_damp*system->polar_damp*pow(pair_ptr->rimg,3.0)/6.0 + 0.5*system->polar_damp*system->polar_damp*pow(pair_ptr->rimg,2.0) + system->polar_damp*pair_ptr->rimg + 1.0);


			}


			/* build the tensor */
			for(p = 0; p < 3; p++) {
				for(q = 0; q < 3; q++) {

					system->A_matrix[ii+p][jj+q] = -3.0*pair_ptr->dimg[p]*pair_ptr->dimg[q]*damping_term2*r5;
					/* additional diagonal term */
					if(p == q)
						system->A_matrix[ii+p][jj+q] += damping_term1*r3;

				}
			}

			/* set the lower half of the tensor component */
			for(p = 0; p < 3; p++)
				for(q = 0; q < 3; q++)
					system->A_matrix[jj+p][ii+q] = system->A_matrix[ii+p][jj+q];


		} /* end j */
	} /* end i */

	free(atom_array);

}


/* calculate the field */
void thole_field(system_t *system, double efieldx, double efieldy, double efieldz) {

	molecule_t *molecule_ptr;
	atom_t *atom_ptr, **atom_array;
	pair_t *pair_ptr;
	int p;
	double r, damping;
	double efield[3];

	efield[0] = efieldx;
	efield[1] = efieldy;
	efield[2] = efieldz;

	/* clear the field */
	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {

			for(p = 0; p < 3; p++) {
				atom_ptr->ef_static[p] = 0;
				atom_ptr->ef_static_self[p] = 0;
			}

		}
	}

	/* calculate the field vectors for the system */
	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next) {


			for(pair_ptr = atom_ptr->pairs; pair_ptr; pair_ptr = pair_ptr->next) {

				if(!pair_ptr->frozen) {	/* disallow self-frozen polarization */

					r = pair_ptr->rimg;
                                        /* exponential field damping */
                                        if(system->field_damp == 0.0)
                                                damping = 1.0;
                                        else
                                                damping = 1.0 - exp(-pow(r/system->field_damp, 3.0));

#ifdef POLAR_SELF
					if(r < system->pbc->cutoff) {

						if(pair_ptr->es_excluded) {
							for(p = 0; p < 3; p++) {
								atom_ptr->ef_static_self[p] += damping*pair_ptr->charge*pair_ptr->dimg[p]/(r*r*r);
								pair_ptr->atom->ef_static_self[p] -= damping*atom_ptr->charge*pair_ptr->dimg[p]/(r*r*r);
							}
						} else {
							for(p = 0; p < 3; p++) {
								atom_ptr->ef_static[p] += damping*pair_ptr->charge*pair_ptr->dimg[p]/(r*r*r);
								pair_ptr->atom->ef_static[p] -= damping*atom_ptr->charge*pair_ptr->dimg[p]/(r*r*r);
							}
						}

					}
#else
					if(!pair_ptr->es_excluded && (r < system->pbc->cutoff)) {

						for(p = 0; p < 3; p++) {
							atom_ptr->ef_static[p] += damping*pair_ptr->charge*pair_ptr->dimg[p]/(r*r*r);
							pair_ptr->atom->ef_static[p] -= damping*atom_ptr->charge*pair_ptr->dimg[p]/(r*r*r);
						}

					}
#endif /* POLAR_SELF */

				}

			}
		}
	}

	/* apply the external field */
	for(molecule_ptr = system->molecules; molecule_ptr; molecule_ptr = molecule_ptr->next)
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next)
			for(p = 0; p < 3; p++)
				atom_ptr->ef_static[p] += efield[p];


}



/* iterative solver of the dipole field tensor */
/* returns the number of iterations required */
int thole_iterative(system_t *system) {

	int i, j, ii, jj, N, p, q;
	int iteration_counter, keep_iterating;
	double error;
	molecule_t *molecule_ptr;
	atom_t *atom_ptr, **atom_array;


	/* generate an array of atom ptrs */
	for(molecule_ptr = system->molecules, N = 0, atom_array = NULL; molecule_ptr; molecule_ptr = molecule_ptr->next) {
		for(atom_ptr = molecule_ptr->atoms; atom_ptr; atom_ptr = atom_ptr->next, N++) {

			atom_array = realloc(atom_array, sizeof(atom_t *)*(N + 1));
			atom_array[N] = atom_ptr;

		}
	}

	/* set the first guess to alpha*E */
	for(i = 0; i < N; i++)
		for(p = 0; p < 3; p++)
			atom_array[i]->mu[p] = atom_array[i]->polarizability*atom_array[i]->ef_static[p];

	/* iterative solver of the dipole field equations */
	for(iteration_counter = 0, keep_iterating = 1; keep_iterating; iteration_counter++) {

		/* divergence detection */
		/* if we fail to converge, then return dipoles as alpha*E */
		if(iteration_counter >= MAX_ITERATION_COUNT) {

			for(i = 0; i < N; i++)
				for(p = 0; p < 3; p++)
					atom_array[i]->mu[p] = atom_array[i]->polarizability*atom_array[i]->ef_static[p];

			free(atom_array);
			return(iteration_counter);

		}

		/* save the current dipole set and clear the induced field vectors */
		for(i = 0; i < N; i++) {
			for(p = 0; p < 3; p++) {

				atom_array[i]->old_mu[p] = atom_array[i]->mu[p];
				atom_array[i]->ef_induced[p] = 0;

			}
		}


		/* contract the dipoles with the field tensor */
		for(i = 0; i < N; i++) {
			ii = i*3;

			for(j = 0; j < N; j++) {
				jj = j*3;

				if(i != j) {

					for(p = 0; p < 3; p++)
						for(q = 0; q < 3; q++)
							atom_array[i]->ef_induced[p] -= system->A_matrix[ii+p][jj+q]*atom_array[j]->mu[q];

				}

			} /* end j */


			/* dipole is the sum of the static and induced parts */
			for(p = 0; p < 3; p++) {

				atom_array[i]->new_mu[p] = atom_array[i]->polarizability*(atom_array[i]->ef_static[p] + atom_array[i]->ef_static_self[p] + atom_array[i]->ef_induced[p]);

				/* Gauss-Seidel smoothing */
				atom_array[i]->mu[p] = atom_array[i]->new_mu[p];

			}


		} /* end i */


		/* immediately reiterate if any component broke tolerance, otherwise we are done */
		for(i = 0, keep_iterating = 0; i < N; i++) {
			for(p = 0; p < 3; p++) {

				error = pow((atom_array[i]->new_mu[p] - atom_array[i]->old_mu[p]), 2.0);
				if(error > pow(system->polar_precision*DEBYE2SKA, 2.0)) keep_iterating = 1;

			}
		}

		/* save the dipoles for the next pass */
		for(i = 0; i < N; i++)
			for(p = 0; p < 3; p++)
				atom_array[i]->mu[p] = atom_array[i]->new_mu[p];


	} /* end keep iterating */

	free(atom_array);


	/* return the iteration count */
	return(iteration_counter);

}


#ifdef DEBUG
void print_matrix(int N, double **matrix) {

	int i, j;

	printf("\n");
	for(i = 0; i < N; i++) {
		for(j = 0; j < N; j++) {
			printf("%.3f ", matrix[i][j]);
		}
		printf("\n");
	}
	printf("\n");

}
#endif /* DEBUG */

