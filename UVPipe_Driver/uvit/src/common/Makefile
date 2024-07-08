CCC=g++
CC=gcc
FLAGS=-w -c -Wall
INCLUDE=-I../../include

SOURCES=uvtUtils.cpp alglibinternal.cpp  alglibmisc.cpp  ap.cpp  linalg.cpp
OBJECTS=$(SOURCES:.cpp=.o)

MODULELIBDIR=../lib
MODULELIB=$(MODULELIBDIR)/libcommon.a  
LIB=-L../../lib -lglog -lcfitsio -lpil  
RM=rm
MK=mkdir -p

all:	$(OBJECTS) $(MODULELIB)

$(MODULELIB): $(OBJECTS)
	$(MK) $(MODULELIBDIR)
	ar cr $(MODULELIB) $(OBJECTS) 

.cpp.o: 
	$(CPP) $(FLAGS) $< -o $@ $(INCLUDE) 

clean:
	$(RM) $(OBJECTS)
