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

// --- Define Fortran name-mangling macro --------------------------------------

// If the F77_FUNC name-mangling macro is undefined, then we we need to define
// it ourselves.
#ifndef F77_FUNC

  // Case 1: F77_FUNC is undefined because we're building for Windows.
  #ifdef BLIS_ENABLE_WINDOWS_BUILD

    // Check whether we need to use uppercase Fortran routine names; otherwise
    // default to lowercase.
    #ifdef BLIS_ENABLE_UPPERCASE_F77

      // Use uppercase routine names (no underscore).
      #define F77_FUNC( name_lower, name_upper ) \
              name_upper
    #else

      // Use lowercase routine names (no underscore).
      #define F77_FUNC( name_lower, name_upper ) \
              name_lower
    #endif

  // Case 2: F77_FUNC is undefined because we're in a Linux-like environment
  // that did not define it for us.
  #else

    // Check whether we need to use uppercase Fortran routine names; otherwise
    // default to lowercase.
    #ifdef BLIS_ENABLE_UPPERCASE_F77

      // Use uppercase routine names (single underscore).
      #define F77_FUNC( name_lower, name_upper ) \
              name_upper ## _
    #else

      // Use lowercase routine names (single underscore).
      #define F77_FUNC( name_lower, name_upper ) \
              name_lower ## _
    #endif

  #endif // #ifdef BLIS_ENABLE_WINDOWS_BUILD

#endif // #ifndef F77_FUNC

