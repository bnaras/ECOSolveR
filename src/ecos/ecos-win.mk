## Modified from original ecos. Original file saved with extension .orig
# Makefile configuration for ECOS under windows

# Whether to use Long or Int for index type
# comment it out to use ints
USE_LONG = 1

## Intel C Compiler
#CC = icc
#CFLAGS = -O3 -m64 -Wall -strict-ansi -DLDL_LONG -DDLONG
#LIBS = -lm

## GNU C Compiler
#CC = gcc

##CFLAGS += -O2 -Wall -DCTRLC=1 -Wextra -fPIC #-ansi -Werror #-ipo
ifdef USE_LONG
CFLAGS = -O2 -Wall -DCTRLC=1 -Wextra -DLDL_LONG -DDLONG -Iinclude -Iexternal/ldl/include -Iexternal/amd/include -Iexternal/SuiteSparse_config 
LDL = ldll.o
AMD = amd_l*.o amd_global.o
else
CFLAGS = -O2 -Wall -DCTRLC=1 -Wextra -Iinclude -Iexternal/ldl/include -Iexternal/amd/include -Iexternal/SuiteSparse_config 
LDL = ldl.o
AMD = amd_i*.o amd_global.o
endif

UNAME = MINGW
ISWINDOWS = 1

# we're on windows (cygwin or msys)
LDFLAGS = -lm
# shared library has extension .dll
SHAREDNAME = libecos.dll

## AR and RANLIB FOR GENERATING LIBRARIES
AR = ar
ARFLAGS = rcs
ARCHIVE = $(AR) $(ARFLAGS)
RANLIB = ranlib

## WHICH FILES TO CLEAN UP
CLEAN = *.o *.obj *.ln *.bb *.bbg *.da *.tcov *.gcov gmon.out *.bak *.d *.gcda *.gcno libecos*.a libecos*.so libecos*.dylib libecos*.dll ecos_bb_test ecostester ecostester.exe runecosexp
