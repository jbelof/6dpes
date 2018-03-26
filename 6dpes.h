
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#define MAXLINE					512
#define MAXVALUE				1.0e200

#include <math.h>

#include <structs.h>
#include <physical_constants.h>
#include <function_prototypes.h>


#define MAX_ITERATION_COUNT			32

#ifdef POLAR_NODAMP
#define POLAR_DAMPING_LINEAR			0
#define POLAR_DAMPING_EXPONENTIAL		0
#else
#define POLAR_DAMPING_LINEAR			1
#define POLAR_DAMPING_EXPONENTIAL		2
#endif /* POLAR_NODAMP */

#define POLAR_DAMPING_FIELD_WIDTH		0.0
#define POLAR_PRECISION				0.00001
/* XXX */
//#define POLAR_DAMPING_EXPONENTIAL_WIDTH		2.1304
#define POLAR_DAMPING_EXPONENTIAL_WIDTH		2.452


