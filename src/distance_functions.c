// Purpose: distance_functions.c contains source code to the object file 
//          distance_functions.o which contains functions to calculate 
//          distances between atoms.

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

#include <stdio.h>
#include <math.h>
//Personal header file 
#include <structures.h>

double dis(CAcoord *atom, int i, int j)
{
  /*Purpose: To get distance between two residues at any positions*/  
  double diff_x=(atom[i].x-atom[j].x);
  double diff_y=(atom[i].y-atom[j].y);
  double diff_z=(atom[i].z-atom[j].z);

  return (sqrt( (diff_x*diff_x) + (diff_y*diff_y) + (diff_z*diff_z) ) );
}

double disSquared(CAcoord *atom, int i, int j)
{
  /*Purpose: To get squared distance between two residues at any positions*/  
  double diff_x=(atom[i].x-atom[j].x);
  double diff_y=(atom[i].y-atom[j].y);
  double diff_z=(atom[i].z-atom[j].z);

  return ( (diff_x*diff_x) + (diff_y*diff_y) + (diff_z*diff_z)  );
}

double disAll(pdb *info, int i, int j)
{
  /*Purpose: To get distance between two atoms at any positions*/  
  double diff_x=(info[i].coordX-info[j].coordX);
  double diff_y=(info[i].coordY-info[j].coordY);
  double diff_z=(info[i].coordZ-info[j].coordZ);

  return (sqrt( (diff_x*diff_x) + (diff_y*diff_y) + (diff_z*diff_z) ) );
}

double disAll_pdb_v23(pdb_v23 *info, int i, int j)
{
  /*Purpose: To get distance between two atoms at any positions*/  
  double diff_x=(info[i].x-info[j].x);
  double diff_y=(info[i].y-info[j].y);
  double diff_z=(info[i].z-info[j].z);

  return (sqrt( (diff_x*diff_x) + (diff_y*diff_y) + (diff_z*diff_z) ) );
}
