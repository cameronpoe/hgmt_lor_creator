#ifndef compton_chain_ordering_h
#define compton_chain_ordering_h

#include "hgmt_structs.h"
#include <stdio.h>
double time_FOM(photon_path path, int *order);
hit *initial_by_best_order(photon_path path, double (*FOM)(photon_path, int *));
hit *initial_by_best_time(photon_path path);
hit *initial_by_weighted(photon_path path, double (*FOM)(photon_path, int *),
                         double time_weight);
#endif
