BUILDDIR := build
PRODUCT := deqn

SRCDIR := src

HDRS := $(wildcard $(SRCDIR)/*.h)

CSRCS := $(wildcard $(SRCDIR)/*.C)
FSRCS += $(wildcard $(SRCDIR)/fortran/*.f90)

OBJS := $(CSRCS:$(SRCDIR)/%.C=$(BUILDDIR)/%.o)
OBJS += $(FSRCS:$(SRCDIR)/fortran/%.f90=$(BUILDDIR)/%.o)

BINARY := $(BUILDDIR)/$(PRODUCT)

# gcc flags:
CXX := g++
CXXFLAGS_DEBUG := -g -DDEBUG
CXXFLAGS_OPT := -O3
CXXFLAGS_OPENMP := -fopenmp

LD := $(CXX)

# select optimized or debug
#CXXFLAGS := $(CXXFLAGS_OPT) $(CPPFLAGS)
CXXFLAGS := $(CXXFLAGS_OPT) $(CXXFLAGS_DEBUG) $(CPPFLAGS)

# add openmp flags (comment out for serial build)
CXXFLAGS += $(CXXFLAGS_OPENMP)
LDFLAGS += $(CXXFLAGS_OPENMP)

all : $(BINARY)

$(BINARY) : $(OBJS)
	@echo linking $@
	$(maketargetdir)
	$(LD) $(LDFLAGS) -o $@ $^

$(BUILDDIR)/%.o : $(SRCDIR)/%.C
	@echo compiling $<
	$(maketargetdir)
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -c -o $@ $<

define maketargetdir
	-@mkdir -p $(dir $@) > /dev/null 2>&1
endef

clean :
	rm -f $(BINARY) $(OBJS)
	rm -rf $(BUILDDIR)
