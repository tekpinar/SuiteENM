//Purpose: nd_functions.c contains a few functions to calculate normalized distance (ND).
//         ND is used to analyze conformational transitions between
//         two protein conformations. It gives percentage of completing action per residue. 
//         ND for 'one residue' can be formulated as below:
//         ND_{i}=(d_{b,i}) / (d_{b,i}+d_{i,e})
//         where b stands for beginning structure, e stands for end structure and i stands for intermediate structure.
//         d denotes distance. ND is calculated for all residues and one can make a color plot to details of conformational 
//         transition for all amino acids. 
//         All functions related to elastic network model are in this file. 

/* Copyright (C) 2012  Mustafa Tekpinar

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; 
version 2.1 of the License.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  
02110-1301  USA

Additionally, you can obtain a copy of the licence from 
http://www.gnu.org/licenses/lgpl-2.1.html

*/
//Email: tekpinar@buffalo.edu
#ifndef  ND_FUNCTIONS_H
#define  ND_FUNCTIONS_H
void printNormalizedDistanceEachAtom(int N/*System Size*/, CAcoord *atom_beg, CAcoord *atom_cur, CAcoord *atom_end,
				     int frameNo, int printDetails, char *outFileName);
void printNormalizedDistanceEachAtom_v2(int N/*System Size*/, CAcoord *atom_beg, CAcoord *atom_cur, CAcoord *atom_end,
					int frameNo, int printDetails, char *outFileName);
void printNormalizedDistanceEachAtom_v3(int N/*System Size*/, CAcoord *atom_beg, CAcoord *atom_cur, CAcoord *atom_end,
					int frameNo, char *outFileName, int range_beg, int range_end, double coeff, int printDetails);
#endif
