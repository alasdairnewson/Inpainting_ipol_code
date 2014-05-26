/*this is the include function for the colour estimation*/

#ifndef IMAGE_AND_FEATURE_RECONSTRUCTION_H
#define IMAGE_AND_FEATURE_RECONSTRUCTION_H


    #include "image_structures.h"
   	#include "common_reconstruct_image.h"
    #include "reconstruct_image_tools.h"

    int check_is_occluded( nTupleVolume *imgVolOcc, int x, int y, int t);

    int check_disp_field(nTupleVolume *dispField, nTupleVolume *departVolume, nTupleVolume *arrivalVolume, nTupleVolume *occVol);
	
    void reconstruct_image_and_features(nTupleVolume* imgVol, nTupleVolume* occVol,
        nTupleVolume *normGradXvol, nTupleVolume *normGradYvol,
        nTupleVolume* dispField, float sigmaColour, int reconstructionType=0, bool initialisation=false);

#endif
