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

UNAME = MINGW
ISWINDOWS = 1

# we're on windows (cygwin or msys)
## LDFLAGS = -lm
ECOS_LDFLAGS = -lm
# shared library has extension .dll
SHAREDNAME = libecos.dll

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
