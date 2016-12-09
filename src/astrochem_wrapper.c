#include <libastrochem.h>

//First, one must select a chemical network file, which is done as follows:
int verbose = 0;
char *chem_file = "osu2009.chm";
net_t network;
read_network (chem_file, &network, verbose);

//The second step is to set the physical parameters of the model. This is done by defining a structure of type phys_t, and by setting the different elements of this structure:
phys_t phys;
phys.cosmic = 1e-17;
phys.chi = 0;
phys.grain_abundance = 0;

//The third step is to allocate a work array that will contain the abundances of all species at a given time. We also need to set the initial abundances:
const char* species[]  = {"CO", "HCO(+)", "e(-)"};
const double initial_abundances[] = {1e-4, 1e-9, 1e-9};
double *abundances;
alloc_abundances (&network, &abundances);
set_initial_abundances (species, 3, initial_abundances, &network, abundances);

//The fourth step is to set the initial density, visual extinction and temperature of the gas cell:
double density = 1000;
double av = 20;
double temperature = 10;

cell_t cell;
cell.nh = &density;
cell.av = &av;
cell.tgas = &temperature;
cell.tdust = &temperature;

double abs_err = ABS_ERR_DEFAULT;
double rel_err = REL_ERR_DEFAULT;
astrochem_mem_t astrochem_mem;

//The fifth step is to initialize the solver, and to set the absolute and relative error:
solver_init (&cell, &network, &phys, abundances , density, abs_err, rel_err, &astrochem_mem );

//The sixth step is to call the solver to advance time, i.e. to compute the abundances up to a given time:
double time = 0;
time += 1e-6;
solve (&astrochem_mem, &network, abundances, time, verbose);

//This step is repeated an number of times, until the dynamical simulation finishes. Between two calls, the cell properties needs to be updated with the new values of the density, temperature, etc. that are computed in the dynamical code:

time += 1e-6;
density = 2000;
temperature = 15;
solve (&astrochem_mem, &network, abundances, time, verbose);

//The seventh and final step is to free the arrays and structures that are used by Astrochem:

solver_close (&astrochem_mem);
free_abundances (abundances);
free_network (&network);
