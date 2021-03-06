ifndef ROOTSYS
$(warning *** ROOT environment not set)
endif

PROGRAMS = 
PROGRAMS_ROOT = analyze.cc

OFILES = 
OFILES_ROOT = menulib.o 

OPTIMIZE = -std=c++11 -g
INCLUDES = -I. -I./include -I$(ROOTSYS)/include
CXXFLAGS = $(OPTIMIZE) $(INCLUDES) -fPIC -Wall -W #-Werror
LINKFLAGS = $(OPTIMIZE) -fPIC

LIBS = -lm
LROOT = `root-config --glibs` -lTreePlayer
CROOT = `root-config --cflags`

%.o : %.cc
	g++ -c $(CXXFLAGS) $(CROOT) -MD $< -o $@
	#@sed -i 's,\($*\.o\)[:]*\(.*\),$@: $$\(wildcard\2\)\n\1:\2,g' $*.d

BINARIES = $(PROGRAMS:.cc=)
BINARIES_ROOT = $(PROGRAMS_ROOT:.cc=)

all: $(BINARIES) $(BINARIES_ROOT)

$(BINARIES): % : %.o $(OFILES); g++ $(LINKFLAGS) -o $@ $^ $(LIBS)
$(BINARIES_ROOT): % : %.o $(OFILES) $(OFILES_ROOT); g++ $(LINKFLAGS) -o $@ $^ $(LIBS) $(LROOT)


clean:
	rm -f $(BINARIES) core.* *.o *.d *~

-include $(OFILES:.o=.d)
-include $(PROGRAMS:.cc=.d)
