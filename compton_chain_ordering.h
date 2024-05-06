#ifndef compton_chain_ordering_h
#define compton_chain_ordering_h

#include "hgmt_structs.h"
#include <stdio.h>
double time_FOM(hit *hits, int *order, int num_hits);
double time_FOM_cum(hit *hits, int *order, int num_hits);
hit *initial_by_best_order(photon_path *path, double (*FOM)(hit *, int *, int));
hit *initial_by_best_time(photon_path *path);
hit *initial_by_weighted(photon_path *path, double (*FOM)(photon_path *, int *),
                         double time_weight);
#endif
