BUILDDIR := build

SRCDIR := src

HDRS := $(wildcard $(SRCDIR)/*.h)

SRCS := $(wildcard $(SRCDIR)/*.c)

CC=mpicc
CFLAGS=-O2 -Wall -g -fopenmp

LD := $(CC)

all: bin2ppm diffbin karman

karman: $(BUILDDIR)/alloc.o $(BUILDDIR)/boundary.o $(BUILDDIR)/init.o $(BUILDDIR)/karman.o $(BUILDDIR)/simulation.o 
	$(CC) $(CFLAGS) -o $@ $^ -lm

bin2ppm: $(BUILDDIR)/bin2ppm.o $(BUILDDIR)/alloc.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

diffbin: $(SRCDIR)/diffbin.c
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(BUILDDIR)/%.o : $(SRCDIR)/%.c
	@echo compiling $<
	$(maketargetdir)
	$(CC) $(CFLAGS) -c -o $@ $<

define maketargetdir
	-@mkdir -p $(dir $@) > /dev/null 2>&1
endef

clean:
	rm -f bin2ppm karman diffbin $(OBJS)
	rm -rf $(BUILDDIR)
