#ifndef DISTANCE_FUNCTIONS_H
#define DISTANCE_FUNCTIONS_H

double dis(CAcoord *atom, int i, int j);
double disSquared(CAcoord *atom, int i, int j);
double disAll(pdb *info, int i, int j);
double disAll_pdb_v23(pdb_v23 *info, int i, int j);
#endif

