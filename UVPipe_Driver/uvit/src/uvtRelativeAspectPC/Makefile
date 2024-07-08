export AS1

INCLUDE1 = -I. -I../common/ -I../../../include/  -I $(SRCDIR)/DataIngest/ -I $(SRCDIR)/uvtUnitConversion -I $(SRCDIR)/uvtMaskBadPix -I $(SRCDIR)/uvtFlatFieldCorr -I $(SRCDIR)/uvtPixPadding/ -I $(SRCDIR)/uvtQEMCPCorr -I $(SRCDIR)/uvtCentroidCorr/ -I $(SRCDIR)/uvtCentroidBias  -I $(SRCDIR)/uvtFrameIntegration -I $(SRCDIR)/uvtCosmicRayCorr/ -I $(SRCDIR)/uvtDetectStar/ -I $(SRCDIR)/uvtDetectDistCorr/ -I $(SRCDIR)/uvtOpticAssDistCorr/ -I $(SRCDIR)/uvtRefFrameCal/ -I $(SRCDIR)/uvtComputeDrift/ -I $(SRCDIR)/uvtComputeThermal -I $(SRCDIR)/uvtRelAspCal -I $(SRCDIR)/uvtComputeJitter/  
CDIR=../common
SRCDIR=..


EXE=uvtRelativeAspectPC	

SOURCES=main.cpp $(CDIR)/uvtUtils.cpp uvtRelativeAspectPC.cpp $(SRCDIR)/DataIngest/DataIngest.cpp $(SRCDIR)/uvtCentroidBias/uvtCentroidBias.cpp $(SRCDIR)/uvtUnitConversion/uvtUnitConversion.cpp $(SRCDIR)/uvtMaskBadPix/uvtMaskBadPix.cpp $(SRCDIR)/uvtCosmicRayCorr/uvtCosmicRayCorr.cpp $(SRCDIR)/uvtFlatFieldCorr/uvtFlatFieldCorr.cpp $(SRCDIR)/uvtPixPadding/uvtPixPadding.cpp $(SRCDIR)/uvtQEMCPCorr/uvtQEMCPCorr.cpp $(SRCDIR)/uvtCentroidCorr/uvtCentroidCorr.cpp $(SRCDIR)/uvtDetectDistCorr/uvtDetectDistCorr.cpp $(SRCDIR)/uvtOpticAssDistCorr/uvtOpticAssDistCorr.cpp $(SRCDIR)/uvtDetectStar/uvtDetectStar.cpp $(SRCDIR)/uvtRefFrameCal/uvtRefFrameCal.cpp $(SRCDIR)/uvtComputeDrift/uvtComputeDrift.cpp  $(SRCDIR)/uvtFrameIntegration/uvtFrameIntegration.cpp $(CDIR)/DataInfo.cpp  $(CDIR)/transform.cpp $(CDIR)/fft.cpp   $(CDIR)/caldb_Handler.cpp $(CDIR)/Directory.cpp  $(CDIR)/spMatrix.cpp  $(SRCDIR)/uvtComputeThermal/uvtComputeThermal.cpp  $(SRCDIR)/uvtRelAspCal/uvtRelAspCal.cpp  $(SRCDIR)/uvtComputeJitter/uvtComputeJitter.cpp 
 
include $(AS1)/uvit/makefile.common

