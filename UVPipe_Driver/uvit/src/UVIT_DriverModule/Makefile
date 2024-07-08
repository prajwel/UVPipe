
export AS1

INCLUDE1 = -I. -I../common/ -I../../../include/  -I../uvtFullFrameAst/
CDIR=../common
SRCDIR=..


EXE=UVIT_DriverModule

SOURCES=main.cpp UVIT_DriverModule.cpp  $(CDIR)/DataInfo.cpp $(CDIR)/uvtUtils.cpp $(CDIR)/transform.cpp $(CDIR)/fft.cpp   $(CDIR)/caldb_Handler.cpp $(CDIR)/Directory.cpp $(CDIR)/spMatrix.cpp   $(CDIR)/Attitude.cpp  $(CDIR)/Database.cpp  $(SRCDIR)/uvtFullFrameAst/uvtFullFrameAst.cpp
lib_sql=-lsqlite3

include $(AS1)/uvit/makefile.common

