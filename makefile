CC=gcc
OPT=-O2
GPROF=#-pg
GDB=#-g
W= -Wall -Wextra \
   -Wformat=0\
   -Warray-bounds\
   -Wempty-body\
   -Wfloat-equal\
   -Wimplicit\
   -Wmaybe-uninitialized\
   -Wmissing-braces\
   -Wparentheses\
   -Wsequence-point\
   -Wswitch\
   -Wtype-limits\
   -Wundef\
   -Wuninitialized\
   -Wunused-parameter\
   -Wunused\

STDFLAG= -std=gnu11   # -std=gnu99 can be used
VFLAGS= -DVERSION='"$(shell cat version.txt)"'
CFLAGS= -c $(STDFLAG) -MMD $(OPT) $(GPROF) $(W) $(GDB) $(VFLAGS)
OFLAGS= -lm $(GPROF)
INCL= -I$(SRCDIR)/mol -I$(SRCDIR)/math -I$(SRCDIR)/qm

OBJDIR=./build
SRCDIR=./source

molsrc=$(wildcard $(SRCDIR)/mol/*.c)
molobj=$(molsrc:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

mathsrc=$(wildcard $(SRCDIR)/math/*.c)
mathobj=$(mathsrc:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

qmsrc=$(wildcard $(SRCDIR)/qm/*.c)
qmobj=$(qmsrc:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

all : qm

qm : $(qmobj) $(molobj) $(mathobj) $(OBJDIR)/qm.o
	$(CC) $^ -o $@ $(OFLAGS)

$(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(CFLAGS) $< -o $@ $(INCL)

clean :
	rm -f $(OBJDIR)/*/*.o $(OBJDIR)/*.o  qm

cleand :
	rm -f $(OBJDIR)/*/*.d $(OBJDIR)/*.d

test :
	bash test.bash

include $(wildcard $(OBJDIR)/*.d)

