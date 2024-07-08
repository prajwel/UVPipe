#this make file is for the RELATIVE ASPECT SERIES FOR IM and PC mode
#!/bin/csh
SRC_PATH=src

all:
	
	make -C $(SRC_PATH)/DataIngest	
	make -C $(SRC_PATH)/uvtMaskBadPix	
	make -C $(SRC_PATH)/uvtUnitConversion	
	make -C $(SRC_PATH)/uvtFlatFieldCorr	
	make -C $(SRC_PATH)/uvtCosmicRayCorr	
	make -C $(SRC_PATH)/uvtSubDivision	
	make -C $(SRC_PATH)/uvtQEMCPCorr	
	make -C $(SRC_PATH)/uvtPixPadding
	make -C $(SRC_PATH)/uvtAccEveryTsec
	make -C $(SRC_PATH)/uvtDetectStar
	make -C $(SRC_PATH)/uvtDetectDistCorr
	make -C $(SRC_PATH)/uvtOpticAssDistCorr
	make -C $(SRC_PATH)/uvtOpticAssDistCorrL2
	make -C $(SRC_PATH)/uvtDetectDistCorrL2
	make -C $(SRC_PATH)/uvtRefFrameCal
	make -C $(SRC_PATH)/uvtComputeDrift
	make -C $(SRC_PATH)/uvtCentroidCorr
	make -C $(SRC_PATH)/uvtCentroidBias
	make -C $(SRC_PATH)/uvtFrameIntegration
	#make -C $(SRC_PATH)/uvt_IM_RA
	#make -C $(SRC_PATH)/uvt_PC_RA
	#make -C $(SRC_PATH)/uvt_IM_L2
	#make -C $(SRC_PATH)/uvt_PC_L2
	make -C $(SRC_PATH)/uvtFindWtdMean
	make -C $(SRC_PATH)/uvtRegAvg
	make -C $(SRC_PATH)/uvtShiftRot
	make -C $(SRC_PATH)/uvtComputeJitter
	make -C $(SRC_PATH)/uvtComputeThermal
	make -C $(SRC_PATH)/uvtRelAspCal
	make -C $(SRC_PATH)/uvtRelativeAspectIM
	make -C $(SRC_PATH)/uvtRelativeAspectPC
	make -C $(SRC_PATH)/uvtLevel2IM
	make -C $(SRC_PATH)/uvtLevel2PC
	make -C $(SRC_PATH)/uvtFullFrameAst
	make -C $(SRC_PATH)/ImageFrameGen
	make -C $(SRC_PATH)/UVIT_DriverModule
	make -C $(SRC_PATH)/uvtComputeDrift_manualMode
	@echo "----------------------------------------"
	@echo "Compilation of All modules for UVIT completed successfully"
	@echo "----------------------------------------"




clean:
	make -C $(SRC_PATH)/DataIngest clean
	make -C $(SRC_PATH)/uvtMaskBadPix clean;
	make -C $(SRC_PATH)/uvtUnitConversion clean
	make -C $(SRC_PATH)/uvtFlatFieldCorr clean
	make -C $(SRC_PATH)/uvtCosmicRayCorr clean
	make -C $(SRC_PATH)/uvtSubDivision clean
	make -C $(SRC_PATH)/uvtQEMCPCorr clean
	make -C $(SRC_PATH)/uvtPixPadding clean
	make -C $(SRC_PATH)/uvtAccEveryTsec clean
	make -C $(SRC_PATH)/uvtDetectStar clean
	make -C $(SRC_PATH)/uvtDetectDistCorr clean
	make -C $(SRC_PATH)/uvtOpticAssDistCorr clean
	make -C $(SRC_PATH)/uvtDetectDistCorrL2 clean
	make -C $(SRC_PATH)/uvtOpticAssDistCorrL2 clean
	make -C $(SRC_PATH)/uvtRefFrameCal clean
	make -C $(SRC_PATH)/uvtComputeDrift clean
	make -C $(SRC_PATH)/uvtCentroidCorr clean	
	make -C $(SRC_PATH)/uvtCentroidBias clean
	make -C $(SRC_PATH)/uvtFrameIntegration clean
	#make -C $(SRC_PATH)/uvt_IM_RA clean
	#make -C $(SRC_PATH)/uvt_PC_RA clean
	#make -C $(SRC_PATH)/uvt_IM_L2 clean
	#make -C $(SRC_PATH)/uvt_PC_L2 clean
	make -C $(SRC_PATH)/uvtFindWtdMean clean
	make -C $(SRC_PATH)/uvtRegAvg clean
	make -C $(SRC_PATH)/uvtShiftRot clean
	make -C $(SRC_PATH)/uvtComputeJitter clean;
	make -C $(SRC_PATH)/uvtComputeThermal clean;
	make -C $(SRC_PATH)/uvtRelAspCal clean;
	make -C $(SRC_PATH)/uvtRelativeAspectIM clean;
	make -C $(SRC_PATH)/uvtRelativeAspectPC clean;
	make -C $(SRC_PATH)/uvtLevel2IM clean
	make -C $(SRC_PATH)/uvtLevel2PC clean
	make -C $(SRC_PATH)/uvtFullFrameAst  clean
	make -C $(SRC_PATH)/ImageFrameGen clean;
	make -C $(SRC_PATH)/UVIT_DriverModule clean;
	make -C $(SRC_PATH)/uvtComputeDrift_manualMode clean;


	


distclean:
	rm  bin/*
	
		


