export AS1
INCLUDE1 = -I. -I../common/ -I../../../include/  -I $(SRCDIR)/DataIngest/ -I $(SRCDIR)/uvtUnitConversion -I $(SRCDIR)/uvtMaskBadPix -I $(SRCDIR)/uvtFlatFieldCorr -I $(SRCDIR)/uvtPixPadding/ -I $(SRCDIR)/uvtQEMCPCorr -I $(SRCDIR)/uvtSubDivision/ -I $(SRCDIR)/uvtAccEveryTsec -I $(SRCDIR)/uvtCosmicRayCorr/ -I $(SRCDIR)/uvtDetectStar/ -I $(SRCDIR)/uvtDetectDistCorr/ -I $(SRCDIR)/uvtOpticAssDistCorr/ -I $(SRCDIR)/uvtRefFrameCal/ -I $(SRCDIR)/uvtComputeDrift/   -I $(SRCDIR)/uvtComputeJitter/  -I $(SRCDIR)/uvtComputeThermal -I $(SRCDIR)/uvtRelAspCal
CDIR=../common
SRCDIR=..


EXE=uvtRelativeAspectIM	

SOURCES=main.cpp uvtRelativeAspectIM.cpp $(SRCDIR)/DataIngest/DataIngest.cpp $(SRCDIR)/uvtAccEveryTsec/uvtAccEveryTsec.cpp $(SRCDIR)/uvtUnitConversion/uvtUnitConversion.cpp $(SRCDIR)/uvtMaskBadPix/uvtMaskBadPix.cpp $(SRCDIR)/uvtCosmicRayCorr/uvtCosmicRayCorr.cpp $(SRCDIR)/uvtFlatFieldCorr/uvtFlatFieldCorr.cpp $(SRCDIR)/uvtPixPadding/uvtPixPadding.cpp $(SRCDIR)/uvtQEMCPCorr/uvtQEMCPCorr.cpp $(SRCDIR)/uvtSubDivision/uvtSubDivision.cpp $(SRCDIR)/uvtDetectDistCorr/uvtDetectDistCorr.cpp $(SRCDIR)/uvtOpticAssDistCorr/uvtOpticAssDistCorr.cpp $(SRCDIR)/uvtDetectStar/uvtDetectStar.cpp $(SRCDIR)/uvtRefFrameCal/uvtRefFrameCal.cpp $(SRCDIR)/uvtComputeDrift/uvtComputeDrift.cpp   $(CDIR)/DataInfo.cpp $(CDIR)/uvtUtils.cpp $(CDIR)/transform.cpp $(CDIR)/fft.cpp   $(CDIR)/caldb_Handler.cpp $(CDIR)/Directory.cpp $(CDIR)/spMatrix.cpp   $(SRCDIR)/uvtComputeJitter/uvtComputeJitter.cpp $(SRCDIR)/uvtRelAspCal/uvtRelAspCal.cpp  $(SRCDIR)/uvtComputeThermal/uvtComputeThermal.cpp
 
include $(AS1)/uvit/makefile.common

