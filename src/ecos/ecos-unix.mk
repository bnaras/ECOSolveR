# Makefile configuration for ECOS

# Whether to use Long or Int for index type
# comment it out to use ints
USE_LONG = 1

## Intel C Compiler
#CC = icc
#CFLAGS = -O3 -m64 -Wall -strict-ansi -DLDL_LONG -DDLONG
#LIBS = -lm

## GNU C Compiler
#CC = gcc
## Modified by bnaras (basically add prefix ECOS_ to flags)
## CFLAGS += -O2 -Wall -DCTRLC=1 -Wextra -fPIC #-ansi -Werror #-ipo
ifdef USE_LONG
## CFLAGS += -DCTRLC=1 -DLDL_LONG -DDLONG
ECOS_CFLAGS = -DCTRLC=1 -DLDL_LONG -DDLONG
LDL = ldll.o
AMD = amd_l*.o amd_global.o
else
#CFLAGS += -DCTRLC=1
ECOS_CFLAGS = -DCTRLC=1
LDL = ldl.o
AMD = amd_i*.o amd_global.o
endif

UNAME := $(shell uname)

ifeq (CYGWIN, $(findstring CYGWIN, $(UNAME)))
ISWINDOWS := 1
else ifeq (MINGW, $(findstring MINGW, $(UNAME)))
ISWINDOWS := 1
else ifeq (MSYS, $(findstring MSYS, $(UNAME)))
ISWINDOWS := 1
else
ISWINDOWS := 0
endif

ifeq ($(UNAME), Darwin)
# we're on apple, no need to link rt library
## LDFLAGS = -lm
ECOS_LDFLAGS = -lm
# shared library has extension .dylib
SHAREDNAME = libecos.dylib
else ifeq ($(ISWINDOWS), 1)
# we're on windows (cygwin or msys)
## LDFLAGS = -lm
ECOS_LDFLAGS = -lm
# shared library has extension .dll
SHAREDNAME = libecos.dll
else
# we're on a linux system, use accurate timer provided by clock_gettime()
## LDFLAGS = -lm -lrt
ECOS_LDFLAGS = -lm -lrt
# shared library has extension .so
SHAREDNAME = libecos.so
endif

## Section commented out by bnaras
## AR and RANLIB FOR GENERATING LIBRARIES
##AR = ar
##ARFLAGS = rcs
## End of Section commented out by bnaras
ARCHIVE = $(AR) $(ARFLAGS)
## Section commented out by bnaras
##RANLIB = ranlib
## End of Section commented out by bnaras

## WHICH FILES TO CLEAN UP
CLEAN = *.o *.obj *.ln *.bb *.bbg *.da *.tcov *.gcov gmon.out *.bak *.d *.gcda *.gcno libecos*.a libecos*.so libecos*.dylib libecos*.dll ecos_bb_test ecostester ecostester.exe runecosexp
