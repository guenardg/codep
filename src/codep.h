/*************************************************************************
 
 (c) 2008-2023 Guillaume Guénard
 Université de Montréal, Montreal, Quebec, Canada
 
 **Permutation Test for Multiscale Codependence Analysis**
 
 This file is part of codep
 
 codep is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 codep is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with codep.  If not, see <https://www.gnu.org/licenses/>.
 
 C header
 
 *************************************************************************/

// Defines
#ifndef __codep_h__

#define __codep_h__

// Includes
#include<R.h>
#include<Rmath.h>

// Type declarations (empty)

// C function declarations
void mcapermute(double *phi_global0, double *tau_ind0, double *rY, int *m,
                double *rx, double *us, int *n, int *perm_global, int *perm_ind,
                int *nperm, int *ind);

#endif
