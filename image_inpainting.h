

#ifndef IMAGE_INPAINTING_H
#define IMAGE_INPAINTING_H

#include <stdlib.h>

#include "io_png.h"
#include "image_structures.h"
#include "patch_match.h"
#include "reconstruct_image.h"
#include "reconstruct_image_and_features.h"
#include "image_operations.h"
#include "morpho.h"

#ifndef SUBSAMPLE_FACTOR
#define SUBSAMPLE_FACTOR 2
#endif

    typedef struct paramInpaint
	{
		float residualThreshold;
		int maxIterations;
		int nLevels;
		bool useFeatures;
	}inpaintingParameterStruct;

patchMatchParameterStruct* initialise_patch_match_parameters(int patchSizeX, int patchSizeY, int imgSizeX, int imgSizeY, bool verboseMode=false);
inpaintingParameterStruct* initialise_inpainting_parameters(int nLevels, bool useFeatures, float residualThreshold, int maxIterations);

void display_inpainting_parameters(inpaintingParameterStruct *inpaintingParams);
void display_patch_match_parameters(patchMatchParameterStruct *patchMatchParams);

void initialise_inpainting(nTupleVolume *imgVol, nTupleVolume *occVol, featurePyramid featuresVolPyramid,
					nTupleVolume *shiftVol, patchMatchParameterStruct *patchMatchParams);

void inpaint_image_wrapper(const char *fileIn,const char *fileOccIn, const char *fileOut,
			int patchSizeX, int patchSizeY, int nLevels=-1, bool useFeatures=false, bool verboseMode=false);
void inpaint_image_wrapper(float *inputImage, int nx, int ny, int nc,
	float *inputOcc, int nOccx, int nOccy, int nOccc,
	const char *fileOut, int patchSizeX, int patchSizeY, int nLevels=-1, bool useFeatures=false, bool verboseMode=false);
						
nTupleVolume * inpaint_image( nTupleVolume *imgVolIn, nTupleVolume *occVolIn,
patchMatchParameterStruct *patchMatchParams, inpaintingParameterStruct *inpaintingParameters);


#endif
