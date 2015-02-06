

#include "image_inpainting.h"

patchMatchParameterStruct * initialise_patch_match_parameters(
	int patchSizeX, int patchSizeY, int imgSizeX, int imgSizeY, bool verboseMode)
{
	patchMatchParameterStruct *patchMatchParams = new patchMatchParameterStruct;

	//set parameter structure
	patchMatchParams->patchSizeX = patchSizeX;
	patchMatchParams->patchSizeY = patchSizeY;
	patchMatchParams->patchSizeT = 1;
	patchMatchParams->nIters = 10; //number of propagation/random search steps in patchMatch
	patchMatchParams->w = max_int(imgSizeX,imgSizeY); //maximum search radius
	patchMatchParams->alpha = 0.5; //search radius shrinkage factor (0.5 in standard PatchMatch)
	patchMatchParams->maxShiftDistance = -1;
	patchMatchParams->partialComparison = 0;
	patchMatchParams->fullSearch = 0;
	//texture attributes
	patchMatchParams->normGradX = NULL;
	patchMatchParams->normGradY = NULL;
	patchMatchParams->verboseMode = verboseMode;
	
	return(patchMatchParams);	
}

void display_patch_match_parameters(patchMatchParameterStruct *patchMatchParams)
{

	printf("\n*************************\n");
	printf("* PATCHMATCH PARAMETERS *\n");
	printf("*************************\n");
	
	printf("Patch size X : %d\n",patchMatchParams->patchSizeX);
	printf("Patch size Y : %d\n",patchMatchParams->patchSizeY);
	printf("Number of propagation/random search iterations: %d\n",patchMatchParams->nIters);
	printf("Random search reduction factor (alpha) : %f\n",patchMatchParams->alpha);
	printf("Maximum search shift allowed (-1 for whole image) : %f\n",patchMatchParams->maxShiftDistance);
	printf("Full search (should be activated only for experimental purposes !!) : %d\n",patchMatchParams->fullSearch);
	printf("Verbose mode : %d\n",patchMatchParams->verboseMode);
}

inpaintingParameterStruct * initialise_inpainting_parameters(int nLevels, bool useFeatures,
	float residualThreshold, int maxIterations)
{
	inpaintingParameterStruct *inpaintingParams = new inpaintingParameterStruct;

	//set inpainting parameter structure
	inpaintingParams->nLevels = nLevels;
	inpaintingParams->useFeatures = useFeatures;
	inpaintingParams->residualThreshold = residualThreshold;
	inpaintingParams->maxIterations = maxIterations;
	
	return(inpaintingParams);	
}


void display_inpainting_parameters(inpaintingParameterStruct *inpaintingParams)
{
	printf("\n*************************\n");
	printf("* INPAINTING PARAMETERS *\n");
	printf("*************************\n");

	printf("Number of levels : %d\n",inpaintingParams->nLevels);
	printf("Use features : %d\n",inpaintingParams->useFeatures);
	printf("Residual threshold: %f\n",inpaintingParams->residualThreshold);
	printf("Maximum number of iterations: %d\n",inpaintingParams->maxIterations);
	
	printf("*************************\n\n");

}

void inpaint_image_wrapper(const char *fileIn,const char *fileOccIn, const char *fileOut,
			int patchSizeX, int patchSizeY, int nLevels, bool useFeatures, bool verboseMode)
{

	// *************************** //
	// ***** READ INPUTS ********* //
	// ************************** //
	//read input image
	size_t nx,ny,nc;
	size_t nOccX,nOccY,nOccC;
	
	//read input image
	printf("Reading input image\n");
	float *inputImage = read_image(fileIn,&nx,&ny,&nc);
	
	//read input occlusion
	printf("Reading input occlusion\n");
	float *inputOcc = read_image(fileOccIn,&nOccX,&nOccY,&nOccC);
	
	
	// ****************************************** //
	// **** INITIALISE PATCHMATCH PARAMETERS **** //
	// ****************************************** //
	patchMatchParameterStruct *patchMatchParams = initialise_patch_match_parameters(patchSizeX, patchSizeY, nx, ny, verboseMode);
	if (check_patch_match_parameters(patchMatchParams) == -1)
		return;
	// ****************************************** //
	// **** INITIALISE INPAINTING PARAMETERS **** //
	// ****************************************** //
	float residualThreshold = 0.1;
	int maxIterations = 10;
	inpaintingParameterStruct *inpaintingParams =
		initialise_inpainting_parameters(nLevels, useFeatures, residualThreshold, maxIterations);
	
	// ******************************** //
	// ***** CREATE IMAGE STRUCTURES*** //
	// ******************************** //
		
	nTupleVolume *imgVolIn = new nTupleVolume(nc,nx,ny,patchSizeX,patchSizeY,IMAGE_INDEXING,inputImage);
	nTupleVolume *occVolIn;
	if (nOccC == 3)		//if we need to convert the input occlusion
	{
		nTupleVolume *occVolTemp = new nTupleVolume(3,nx,ny,patchSizeX,patchSizeY,IMAGE_INDEXING,inputOcc);
		occVolIn = rgb_to_grey(occVolTemp);
		delete(occVolTemp);
	}
	else
		occVolIn = new nTupleVolume(1,nOccX,nOccY,patchSizeX,patchSizeY,IMAGE_INDEXING,inputOcc);
	
	
	occVolIn->binarise();

	// ***** CALL MAIN ROUTINE **** //
	
	nTupleVolume * imgOut = inpaint_image(imgVolIn, occVolIn, patchMatchParams, inpaintingParams);

	//write output
	//write_image(imgOut,fileOut,255);
	write_image(imgOut,fileOut);
	delete(imgOut);
}

float * inpaint_image_wrapper(float *inputImage, int nx, int ny, int nc,
	float *inputOcc, int nOccX, int nOccY, int nOccC,
	int patchSizeX, int patchSizeY, int nLevels, bool useFeatures)
{

	// ****************************************** //
	// **** INITIALISE PATCHMATCH PARAMETERS **** //
	// ****************************************** //
	patchMatchParameterStruct *patchMatchParams = initialise_patch_match_parameters(patchSizeX, patchSizeY, nx, ny);

	// ****************************************** //
	// **** INITIALISE INPAINTING PARAMETERS **** //
	// ****************************************** //
	float residualThreshold = 0.1;
	int maxIterations = 10;
	inpaintingParameterStruct *inpaintingParameters = initialise_inpainting_parameters(nLevels, useFeatures, residualThreshold, maxIterations);

	// ******************************** //
	// ***** CREATE IMAGE STRUCTURES*** //
	// ******************************** //
	
	nTupleVolume *imgVolIn = new nTupleVolume(nc,nx,ny,patchSizeX,patchSizeY,IMAGE_INDEXING,inputImage);
	nTupleVolume *occVolIn;
	if (nOccC == 3)		//if we need to convert the input occlusion
	{
		nTupleVolume *occVolTemp = new nTupleVolume(nOccC,nOccX,nOccY,patchSizeX,patchSizeY,IMAGE_INDEXING,inputOcc);
		occVolIn = rgb_to_grey(occVolTemp);
		delete(occVolTemp);
	}
	else
		occVolIn = new nTupleVolume(1,nOccX,nOccY,patchSizeX,patchSizeY,IMAGE_INDEXING,inputOcc);
	
	occVolIn->binarise();

	// ***** CALL MAIN ROUTINE **** //
	
	nTupleVolume * imgOut = inpaint_image(imgVolIn, occVolIn, patchMatchParams, inpaintingParameters);

	return(imgOut->get_data_ptr());
}

nTupleVolume * inpaint_image( nTupleVolume *imgVolIn, nTupleVolume *occVolIn,
patchMatchParameterStruct *patchMatchParams, inpaintingParameterStruct *inpaintingParams)
{
	
	// ******************************************************************** //
	// **** AUTOMATICALLY DETERMINE NUMBER OF LEVELS, IF NOT SPECIFIED **** //
	// ******************************************************************** //
	if (inpaintingParams->nLevels == -1)
	{
		inpaintingParams->nLevels =
			determine_multiscale_level_number(occVolIn,imgVolIn->patchSizeX,imgVolIn->patchSizeY);
	}
	display_inpainting_parameters(inpaintingParams);
	display_patch_match_parameters(patchMatchParams);
	
	nTupleVolume *imgVolOut;

	// ************************** //
	// **** CREATE PYRDAMIDS **** //
	// ************************** //
	nTupleVolumePyramid imgVolPyramid = create_nTupleVolume_pyramid(imgVolIn, inpaintingParams->nLevels);
	nTupleVolumePyramid occVolPyramid = create_nTupleVolume_pyramid_binary(occVolIn, inpaintingParams->nLevels);
	featurePyramid featuresVolPyramid;
	if (inpaintingParams->useFeatures == true)
	{
		double t1 = clock();
		featuresVolPyramid = create_feature_pyramid(imgVolIn, occVolIn, inpaintingParams->nLevels);
		MY_PRINTF("\n\nFeatures calculation time: %f\n",((double)(clock()-t1)) / CLOCKS_PER_SEC);
	}
	else
	{
		featuresVolPyramid.normGradX = NULL;
        featuresVolPyramid.normGradY = NULL;
        featuresVolPyramid.nLevels = -1;
	}

	//create structuring element
	nTupleVolume *structElDilate = create_structuring_element("rectangle", imgVolIn->patchSizeX, imgVolIn->patchSizeY);
	
	//show_patch_match_parameters(patchMatchParams);
	
	// ****************************************** //
	// ************* START INPAINTING *********** //
	// ****************************************** //
	
	nTupleVolume *imgVol,*normGradXvol,*normGradYvol;
	nTupleVolume *shiftVol=NULL;
	for (int level=( (inpaintingParams->nLevels)-1); level>=0; level--)
	{
		printf("Current pyramid level : %d\n",level);
		nTupleVolume *imgVolPrevious,*occVol,*occVolDilate;

		if (patchMatchParams->maxShiftDistance != -1)		
			patchMatchParams->maxShiftDistance =
			(float)( (patchMatchParams->maxShiftDistance)/( pow((float)SUBSAMPLE_FACTOR,(float)level) ));
		
		imgVol = copy_image_nTuple(imgVolPyramid[level]);
		occVol = copy_image_nTuple(occVolPyramid[level]);
		//create dilated occlusion
		occVolDilate = imdilate(occVol, structElDilate);
		
		if (featuresVolPyramid.nLevels >= 0)
		{
			normGradXvol = copy_image_nTuple((featuresVolPyramid.normGradX)[level]);
			normGradYvol = copy_image_nTuple((featuresVolPyramid.normGradY)[level]);
			//attach features to patchMatch parameters
			patchMatchParams->normGradX = normGradXvol;
			patchMatchParams->normGradY = normGradYvol;
		}
					
		//initialise solution
		if (level == ((inpaintingParams->nLevels)-1))
		{
			shiftVol = new nTupleVolume(4,imgVol->xSize,imgVol->ySize,imgVol->patchSizeX,imgVol->patchSizeY,IMAGE_INDEXING);
			shiftVol->set_all_image_values(0);
			printf("\nInitialisation started\n\n\n");
			initialise_inpainting(imgVol,occVol,featuresVolPyramid,shiftVol,patchMatchParams);
			patchMatchParams->partialComparison = 0;
			printf("\nInitialisation finished\n\n\n");
			
			if (featuresVolPyramid.nLevels >= 0)	//retrieve features from the pointers in the patchMatch parameters
			{
				normGradXvol = patchMatchParams->normGradX;
				normGradYvol = patchMatchParams->normGradY;
			}
		}
		else	//reconstruct current solution
		{
			if (featuresVolPyramid.nLevels >= 0)
			{
				reconstruct_image_and_features(imgVol, occVol,
				normGradXvol, normGradYvol,
				shiftVol, SIGMA_COLOUR);
			}
			else
			{
				reconstruct_image(imgVol,occVol,shiftVol,SIGMA_COLOUR);
				//write_shift_map(shiftVol,fileOut);
			}
		}
		
		calclulate_patch_distances(imgVol,imgVol,shiftVol,occVolDilate,patchMatchParams);
		
		//iterate ANN search and reconstruction
		int iterationNb = 0;
		imageDataType residual = FLT_MAX;
		while( (residual > (inpaintingParams->residualThreshold) ) && (iterationNb < (inpaintingParams->maxIterations) ) )
		{
			//copy current imgVol
			imgVolPrevious = copy_image_nTuple(imgVol);
			patch_match_ANN(imgVol,imgVol,shiftVol,occVolDilate,occVolDilate,patchMatchParams);
			if (featuresVolPyramid.nLevels >= 0)
			{
				reconstruct_image_and_features(imgVol, occVol,
        			normGradXvol, normGradYvol,
        			shiftVol, SIGMA_COLOUR);
			}
			else
				reconstruct_image(imgVol,occVol,shiftVol,SIGMA_COLOUR);
			residual = calculate_residual(imgVol,imgVolPrevious,occVol);
			if (patchMatchParams->verboseMode == true)
				printf("Iteration number %d, residual = %f\n",iterationNb,residual);
			iterationNb++;
		}
		//upsample shift volume, if we are not on the finest level
		if (level >0)
		{	
			nTupleVolume * shiftVolTemp = up_sample_image(shiftVol, SUBSAMPLE_FACTOR,imgVolPyramid[level-1]);
			delete(shiftVol);
			shiftVol = copy_image_nTuple(shiftVolTemp);
			shiftVol->multiply((imageDataType)SUBSAMPLE_FACTOR);
			delete(shiftVolTemp);
		}
		else
		{
			reconstruct_image(imgVol,occVol,shiftVol,SIGMA_COLOUR,3);
			imgVolOut = new nTupleVolume(imgVol);
		}
		//destroy structures
		delete(imgVol);
		delete(imgVolPrevious);
		delete(occVol);
		delete(occVolDilate);
		if (featuresVolPyramid.nLevels >= 0)
		{
			delete normGradXvol;
			delete normGradYvol;
		}
	}
	
	// ************************** //
	// **** DELETE STRUCTURES *** //
	// ************************** //
	for (int i=0; i< (inpaintingParams->nLevels); i++)
	{
		delete(imgVolPyramid[i]);
		delete(occVolPyramid[i]);
	}
	delete(imgVolPyramid);
	delete(occVolPyramid);
	delete(shiftVol);
	delete_feature_pyramid(featuresVolPyramid);
	delete(patchMatchParams);
	
	printf("Inpainting finished !\n");

	return(imgVolOut);
}


void initialise_inpainting(nTupleVolume *imgVol, nTupleVolume *occVol, featurePyramid featuresVolPyramid,
				nTupleVolume *shiftVol, patchMatchParameterStruct *patchMatchParams)
{
	int iterNb=0;
	patchMatchParams->partialComparison = 1;
	bool initialisation = true;
	nTupleVolume *occVolIter;
	occVolIter = copy_image_nTuple(occVol);
	
	seed_random_numbers((double)3);
	
	nTupleVolume *structElErode = create_structuring_element("rectangle", 3, 3);
	nTupleVolume *structElDilate = create_structuring_element("rectangle", imgVol->patchSizeX, imgVol->patchSizeY);
	
	nTupleVolume *occVolDilate = imdilate(occVol, structElDilate);
	
	//extract features images from featuresVolPyramid (coarsest level)
	nTupleVolume *normGradXvol,*normGradYvol;
	
	if (featuresVolPyramid.nLevels >= 0)
	{
		normGradXvol = copy_image_nTuple((featuresVolPyramid.normGradX)[featuresVolPyramid.nLevels-1]);
		normGradYvol = copy_image_nTuple((featuresVolPyramid.normGradY)[featuresVolPyramid.nLevels-1]);
		//attach features to patchMatch parameters
		patchMatchParams->normGradX = normGradXvol;
		patchMatchParams->normGradY = normGradYvol;
	}
	
	while ( (occVolIter->sum_nTupleVolume()) >0)
	{
		nTupleVolume *occVolErode = imerode(occVolIter, structElErode);
		nTupleVolume *occVolPatchMatch = copy_image_nTuple(occVolDilate);
		
		/***************************/
		/*******   NNSEARCH   ******/
		/***************************/		
		//indicate which pixels can be used for comparing patches
		//but we are not allowed to point to in the patchMatch (we set these pixels to 2)
		for (int x=0; x<(occVolPatchMatch->xSize); x++)
			for (int y=0; y<(occVolPatchMatch->ySize); y++)
			{
				if ( (occVolDilate->get_image_value(x,y,0) - occVolIter->get_image_value(x,y,0)) == 1)
					occVolPatchMatch->set_image_value(x,y,0,(imageDataType)2);
			}
			
		//set first guess
		nTupleVolume *firstGuess = copy_image_nTuple(shiftVol);
		//carry out patchMatch
		patch_match_ANN(imgVol,imgVol,shiftVol,occVolPatchMatch,occVolDilate,patchMatchParams,firstGuess);
		/***************************/
		/****   RECONSTRUCTION   ***/
		/***************************/
		nTupleVolume *occVolReconstruct = copy_image_nTuple(occVolIter);
		//Indicate which pixels are on the current border, and need to be inpainted
		//Also, we indicate that the pixels inside the occlusion (and not on the border) are occluded (and
		//therefore not to be used for reconstruction) but should not
		//be inpainted at the current iteration. To indicate this
		//we set them to 2
		for (int x=0; x<(occVolPatchMatch->xSize); x++)
			for (int y=0; y<(occVolPatchMatch->ySize); y++)
			{
				imageDataType valueTemp = abs(occVolIter->get_image_value(x,y,0) - occVolErode->get_image_value(x,y,0));
				occVolReconstruct->set_image_value(x,y,0,valueTemp);
				if (occVolErode->get_image_value(x,y,0) == 1)
					occVolReconstruct->set_image_value(x,y,0,(imageDataType)2);
			}

		//call reconstruction function
		if (featuresVolPyramid.nLevels >= 0)
			reconstruct_image_and_features(imgVol, occVolReconstruct,
		    	normGradXvol, normGradYvol,
		    	shiftVol, SIGMA_COLOUR, AGGREGATED_PATCHES,initialisation);
		else
		{
			reconstruct_image(imgVol,occVolReconstruct,shiftVol,SIGMA_COLOUR,AGGREGATED_PATCHES,initialisation);
		}
		
		iterNb++;
		if (patchMatchParams->verboseMode == true)
			printf("\n Initialisation iteration number : %d \n",iterNb);
		delete(occVolPatchMatch);
		delete(occVolReconstruct);
		
		//copy the information from the eroded occlusion to the current occlusion (occVolIter)
		delete(occVolIter);
		occVolIter = copy_image_nTuple(occVolErode);
		delete(occVolErode);
	}
}

