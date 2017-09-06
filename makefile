CC=gcc
OPT=-O2
GPROF=#-pg
GDB=#-g
W= -Wall -Wextra \
   -Warray-bounds\
   -Wempty-body\
   -Wfloat-equal\
   -Wimplicit\
   -Wmaybe-uninitialized\
   -Wmisleading-indentation\
   -Wmissing-braces\
   -Wparentheses\
   -Wsequence-point\
   -Wswitch\
   -Wtype-limits\
   -Wundef\
   -Wuninitialized\
   -Wunused-parameter\
   -Wunused\

CFLAGS= -c -MMD $(OPT) $(GPROF) $(W) $(GDB)
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

include $(wildcard $(OBJDIR)/*.d)

