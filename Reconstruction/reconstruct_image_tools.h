/*this function defines the tools necessary for the spatio-temporal patch match algorithm*/

#ifndef RECONSTRUCT_IMAGE_TOOLS
#define RECONSTRUCT_IMAGE_TOOLS
    
    #define MAXCOLOURLEVEL 255

	#include "common_reconstruct_image.h"
	#include "image_structures.h"
    
    int check_disp_field(nTupleVolume *dispField, nTupleVolume *departVolume, nTupleVolume *arrivalVolume);
    
    float get_adaptive_sigma(float *weights, int weightsLength, float sigmaPercentile);
    
    int estimate_best_colour(nTupleVolume *imgVol, float *weights, int weightsLength,
                                float *colours, int i, int j, int k);
    
#endif
