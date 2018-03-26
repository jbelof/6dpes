
typedef struct _energy {

	double total;
	double bond;
	double lj;
	double es;
	double pol;

} energy_t;

typedef struct _pair {
	int frozen;
	int rd_excluded, es_excluded;
	double lrc;
	double charge;
	double epsilon;
	double sigma;
	double r;
	double d[3];
	double rimg;
	double dimg[3];
	struct _atom *atom;
	struct _molecule *molecule;
	struct _pair *next;
} pair_t;


typedef struct _atom {
	int id;
	char atomtype[MAXLINE];
	int frozen, adiabatic;
	double mass;
	double charge;
	double polarizability;
	double epsilon;
	double sigma;
	double pos[3];
	double ef_static[3];
	double ef_static_self[3];
	double ef_induced[3];
	double mu[3];
	double old_mu[3];
	double new_mu[3];
	pair_t *pairs;
	struct _atom *next;
} atom_t;


typedef struct _molecule {
	int id;
	char moleculetype[MAXLINE];
	double mass;
	int frozen, adiabatic;
	double com[3];
	atom_t *atoms;
	struct _molecule *next;
} molecule_t;

typedef struct _pbc {
	double basis[3][3];		/* unit cell lattice (A) */
	double reciprocal_basis[3][3];	/* reciprocal space lattice (1/A) */
	double cutoff;			/* radial cutoff (A) */
	double volume;			/* unit cell volume (A^3) */
} pbc_t;


typedef struct _system {
	double ewald_alpha;
	int ewald_kmax;
	double polar_damp, polar_precision, field_damp;
	int damp_type;
	double **A_matrix;
	char *pdb_input;
	pbc_t *pbc;
	molecule_t *molecules;
} system_t;

