
#define H			6.626068e-34		/* Planck's constant in J s */
#define HBAR			1.054571e-34		/* above divided by 2pi in J s */
#define KB			1.3806503e-23		/* Boltzmann's constant in J/K */
#define NA			6.0221415e23		/* Avogadro's number */

/* H2 specific constants */
#define H2_MASS			3.348e-27		/* mass of H2 molecule in kg */
#define H2_REDUCED_MASS		0.5*H2_MASS		/* reduced mass of H2 molecule in kg */
#define H2_De			5.074			/* dissocation energy in eV (Herzberg) */
#define H2_R0			0.74144			/* equilibrium bond distance in angstroms (Herzberg) */
#define H2_A			1.878			/* Morse alpha in 1/angstrom (Herzberg) */

#define EWALD_KMAX		7


/* conversion factors */
#define AU2ANGSTROM		0.529177249		/* convert from Bohr radii to angstroms */
#define METER2ANGSTROM		1.0e10			/* convert from meters to angstroms */
#define HARTREE2KELVIN		3.15774655e5		/* convert from Hartrees to Kelvin */
#define E2REDUCED		408.7816		/* convert from e to sqrt(K*A) */
#define ATM2REDUCED		0.0073389366		/* convert from atm to K/A^3 */
#define ATM2PASCALS		101325.0		/* convert from atm to Pascals */
#define ATM2PSI			14.6959488		/* convert from atm to psi */
#define A32CM3			1.0e-24			/* convert from A^3 to cm^3 */
#define AMU2KG			1.66053873e-27		/* convert amu's to kg */
#define DEBYE2SKA		85.10597636		/* convert from Debye to sqrt(KA)*A */
#define EV2K			1.160444e4		/* convert eV to K */


