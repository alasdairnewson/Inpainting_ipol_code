

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


parameterStruct* initialise_patch_match_parameters(int patchSizeX, int patchSizeY, int imgSizeX, int imgSizeY);
void initialise_inpainting(nTupleVolume *imgVol, nTupleVolume *occVol, featurePyramid featuresVolPyramid,
					nTupleVolume *shiftVol, parameterStruct *patchMatchParams);
void inpaint_image(const char *fileIn,const char *fileOccIn, const char *fileOut, bool useFeatures);


#endif
