#ifndef hgmt_event_h
#define hgmt_event_h

#include "vector_ops.h"
#include <stdio.h>

typedef unsigned int uint;

typedef struct event_ {
  vec3d location;   // 1,2,3: Position [cm]
  double dir_cos_x; // 4: Direction Cosine X
  double dir_cos_y; // 5: Direction Cosine Y
  double energy;    // 6: Energy [MeV]
  double weight;    // 7: Weight
  int particle;     // 8: Particle Type (in PDG Format)
  uint dir_cos_neg; // 9: Flag to tell if Third Direction Cosine is Negative (1
                    // means true)
  uint first_in_hist; // 10: Flag to tell if this is the First Scored Particle
                      // from this History (1 means true)
  double tof;         // 11: Time of Flight [ns]
  uint run_id;
  uint event_id;
  uint track_id;
} event;

typedef struct prim_lor_ {
  vec3d first_loc;
  vec3d second_loc;
  double first_energy;
  double second_energy;
  double first_time;
  double second_time;
} prim_lor;

typedef struct _lor {
  vec3d center;
  vec3d dir;
  double long_uncert;
  double transverse_uncert;
} lor;

#endif
