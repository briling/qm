
CC=gcc
OPT=-O2
GPROF=#-pg
GDB=#-g
W= \
   -Warray-bounds\
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
   -Wunused\
   -Wmisleading-indentation\
   -Wempty-body\
   -Wunused-parameter\
   #-Winline\
   #-Wunsafe-loop-optimizations\
   #-W -Wall\
   -Wno-format\
   -Wconversion\
   -Wsign-compare\
   -Wjump-misses-init\
   #-Werror\

CFLAGS= -c -MMD $(OPT) $(GPROF) $(W) $(GDB)
OFLAGS= -lm $(GPROF)
INCL= -I$(SRCDIR)/mol -I$(SRCDIR)/math -I$(SRCDIR)/qm

OBJDIR=./obj
SRCDIR=./src

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
