export AS1
INCLUDE1 = -I. -I../common/ -I../../../include/ -I $(SRCDIR)/uvtDetectStar -I $(SRCDIR)/DataIngest -I $(SRCDIR)/uvtFullFrameAst -I $(SRCDIR)/uvtComputeDrift
CDIR=../common
SRCDIR=..


EXE=uvtLevel2IM

SOURCES=main.cpp uvtLevel2IM.cpp  $(SRCDIR)/DataIngest/DataIngest.cpp $(SRCDIR)/uvtDetectStar/uvtDetectStar.cpp  $(SRCDIR)/uvtComputeDrift/uvtComputeDrift.cpp   $(CDIR)/DataInfo.cpp $(CDIR)/uvtUtils.cpp $(CDIR)/transform.cpp $(CDIR)/fft.cpp   $(CDIR)/caldb_Handler.cpp $(CDIR)/Directory.cpp $(CDIR)/spMatrix.cpp   $(CDIR)/Attitude.cpp $(CDIR)/Database.cpp $(SRCDIR)/uvtFullFrameAst/uvtFullFrameAst.cpp
lib_sql=-lsqlite3

include $(AS1)/uvit/makefile.common



