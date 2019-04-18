# -*- mode: Makefile -*-
CPLEXDIR      = .
SYSTEM        = x86-64_linux

CC  = gcc -O0

COPT  = -m64 -fno-strict-aliasing

CPLEXLIBDIR   = $(CPLEXDIR)/lib
CONCERTLIBDIR = $(CPLEXDIR)/lib

# For dynamic linking
CPLEXBINDIR   = $(CPLEXDIR)/bin
CPLEXOBJDIR   = $(CPLEXDIR)/obj
CPLEXSOLDIR   = $(CPLEXDIR)/sol
CPLEXMAKEDIR  = $(CPLEXDIR)/makedir
SRCC1 		  = $(CPLEXDIR)/src/remotec
SRCC2 		  = $(CPLEXDIR)/src/keyboard
CPLEXDATADIR  = $(CPLEXDIR)/data
CPLEXLIB      = cplex$(dynamic:yes=1290)
# run           = $(dynamic:yes=LD_LIBRARY_PATH=$(CPLEXBINDIR))
CLNDIRS   = -L$(CPLEXLIBDIR) $(dynamic:yes=-L$(CPLEXBINDIR))
CLNFLAGS  = -l$(CPLEXLIB) -lm -lpthread -ldl

CPLEXINCDIR   = $(CPLEXDIR)/include

EXINC         = $(EXDIR)/include

CFLAGS  = $(COPT)  -I$(CPLEXINCDIR)

dllpref		=	lib
dllsuff		=	.so
DLLFLAGS	=	-m64 -shared -rdynamic
RDYNAMIC	=	-rdynamic
LDLIBPATH	=	LD_LIBRARY_PATH
O		=	.o
E		=	

TRANSPORT	=	tcpip

# Include other Makefiles

include $(CPLEXMAKEDIR)/Makefile.remoteobject
include $(CPLEXMAKEDIR)/Makefile.$(TRANSPORT)
# Local Variables:
# mode: makefile
# End:
