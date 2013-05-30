#
#  libflame
#  An object-based infrastructure for developing high-performance
#  dense linear algebra libraries.
#
#  Copyright (C) 2011, The University of Texas
#
#  libflame is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation; either version 2.1 of
#  the License, or (at your option) any later version.
#
#  libflame is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with libflame; if you did not receive a copy, see
#  http://www.gnu.org/licenses/.
#
#  For more information, please contact us at flame@cs.utexas.edu or
#  send mail to:
#
#  Field G. Van Zee and/or
#  Robert A. van de Geijn
#  The University of Texas at Austin
#  Department of Computer Sciences
#  1 University Station C0500
#  Austin TX 78712
#

#
# config.mk
# build/config.mk.  Generated from config.mk.in by configure.
#
# Field G. Van Zee
#
# Variable definitions needed by the top-level Makefile.
# 


# Only include this block of code once
ifndef CONFIG_MK_INCLUDED
CONFIG_MK_INCLUDED := yes



#
# --- Build definitions -------------------------------------------------------
#

# This variable will identify the machine architecture (ie: x86_64, ia64, etc.)
# It will be used to construct a string that is appended to the libflame
# library name.
ARCH                           := x86_64

# The install prefix tell us where to install the libraries and header file
# directory. Notice that we support the use of DESTDIR so that advanced users
# may install to a temporary location.
INSTALL_PREFIX                 := $(DESTDIR)/Users/martin/flame

# Variables corresponding to other configure-time options.
FLA_ENABLE_VERBOSE_MAKE_OUTPUT := no
FLA_ENABLE_STATIC_BUILD        := yes
FLA_ENABLE_DYNAMIC_BUILD       := no
FLA_ENABLE_MAX_ARG_LIST_HACK   := no
FLA_ENABLE_BUILTIN_BLAS        := no
FLA_MULTITHREADING_MODEL       := none




#
# --- Utility program definitions ---------------------------------------------
#

SH         := /bin/sh
MV         := mv
CP         := cp
CAT        := cat
MKDIR      := mkdir -p
RM         := rm
RM_F       := rm -f
RM_RF      := rm -rf
SYMLINK    := ln -sf
FIND       := find
TOUCH      := touch
XARGS      := xargs
CTAGS      := ctags
ETAGS      := etags
RANLIB     := ranlib
INSTALL    := /usr/bin/install -c
MAIL       := mail
EMPTY_FILE := /dev/null




#
# --- Development tools definitions -------------------------------------------
#

# --- Determine the C compiler and related flags ---
CC           := gcc
CPPROCFLAGS  := -DHAVE_CONFIG_H -DBLIS_FROM_LIBFLAME
CMISCFLAGS   := -std=c99  
CDBGFLAGS    := -ggdb
CWARNFLAGS   := -Wall -Wno-comment
COPTFLAGS    := -O0
CVECFLAGS    := 

# Add position-independent code flag if building a dynamic library.
ifeq ($(FLA_ENABLE_DYNAMIC_BUILD),yes)
CMISCFLAGS   += -fPIC
endif

# Add explicit C99 support.
#CMISCFLAGS   += -std=c99

# Aggregate all of the flags into two groups: one for optimizable code, and
# one for code that should not be optimized.
CFLAGS       := $(CDBGFLAGS) $(COPTFLAGS) $(CVECFLAGS) $(CWARNFLAGS) $(CMISCFLAGS) $(CPPROCFLAGS)
CFLAGS_NOOPT := $(CDBGFLAGS) $(CWARNFLAGS) $(CMISCFLAGS) $(CPPROCFLAGS)

# If the user requested that he supplement configure's CFLAGS with this own,
# add them in.
# *** Notice that we ALSO modify the 'no optimization' set of flags.
#     THIS IS DANGEROUS. THE USER MUST TAKE CARE NOT TO ADD "EXTRA FLAGS"
#     THAT WOULD DISRUPT THE NUMERICAL PROPERTIES OF SLAMCH/DLAMCH. ***
ifeq (no,yes)
CFLAGS       += 
CFLAGS_NOOPT += 
endif

# If the user provided his own CFLAGS, allow them to override our own.
# *** Notice that we ALSO modify the 'no optimization' set of flags.
#     THIS IS DANGEROUS. THE USER MUST TAKE CARE NOT TO ADD "EXTRA FLAGS"
#     THAT WOULD DISRUPT THE NUMERICAL PROPERTIES OF SLAMCH/DLAMCH. ***
ifneq (,)
CFLAGS       := -ggdb -O0
CFLAGS_NOOPT := -ggdb -O0
endif


# --- Determine the C preprocessor ---
CPP          := cpp
CPPFLAGS     := -DHAVE_CONFIG_H


# --- Determine the archiver and related flags ---
AR               := ar
ARFLAGS          := cru
AR_ARG_LIST_FILE := ar_arg_list
AR_OBJ_LIST_FILE := ar_obj_list


# --- Determine the linker and related flags ---
LINKER       := $(CC)
ifeq (none, pthreads)
LDFLAGS      :=  -lpthread  @FLIBS@
else
ifeq (none, openmp)
LDFLAGS      :=    @FLIBS@
else
LDFLAGS      :=   @FLIBS@
endif
endif




#
# --- send-thanks target definitions ------------------------------------------
#

THANKS_MSG_SUBJECT := "Thanks!"
THANKS_MSG_BODY    := "Happy user."
THANKS_MSG_EMAIL   := flame_thanks@cs.utexas.edu




# end of ifndef CONFIG_MK_INCLUDED conditional block
endif
