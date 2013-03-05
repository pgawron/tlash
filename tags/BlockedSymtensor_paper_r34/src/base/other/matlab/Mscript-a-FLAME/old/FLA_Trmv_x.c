/*
   libflame
   An object-based infrastructure for developing high-performance
   dense linear algebra libraries.

   Copyright (C) 2011, The University of Texas

   libflame is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as
   published by the Free Software Foundation; either version 2.1 of
   the License, or (at your option) any later version.

   libflame is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with libflame; if you did not receive a copy, see
   http://www.gnu.org/licenses/.

   For more information, please contact us at flame@cs.utexas.edu or
   send mail to:

   Field G. Van Zee and/or
   Robert A. van de Geijn
   The University of Texas at Austin
   Department of Computer Sciences
   1 University Station C0500
   Austin TX 78712
*/
/* Mex file created by interface.pl for the FLAME C function FLA_Trmvsx_external() */

#include "mex.h"
#include "FLA_Matlab2C.h"

extern double FLA_Clock();

#define NRHS 8
#define NINT 3

#define NOBJ (NRHS-NINT)


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int attr[NINT];
  FLA_Obj obj[NOBJ];
  double *dtime;

  FLA_Init();

  /* Check if the number of arguments supplied is correct */
  FLA_M2C_CheckNumArgs(NRHS, nrhs);

  /* Convert Matlab arguments into the appropriate FLAME C arguments */
  FLA_M2C_ConvertArgs(NRHS, prhs, NINT, attr, obj);

  /* If an extra argument is supplied, collect timing informaion in it. */
  if (nrhs == NRHS+1)
    dtime = FLA_M2C_ConvertDoublePtr(prhs[NRHS]);

  /* Now call the C FLAME function, timing it if the extra argument is given. */

  if (nrhs == NRHS+1)
    *dtime = FLA_Clock();

  FLA_Trmvsx_external(attr[0], attr[1], attr[2], obj[0], obj[1], obj[2], obj[3], obj[4]);

  if (nrhs == NRHS+1)
    *dtime = FLA_Clock() - *dtime;

  FLA_Finalize();

}

