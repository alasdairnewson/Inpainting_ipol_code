/*this is the include function for the colour estimation*/

#ifndef RECONSTRUCT_IMAGE_H
#define RECONSTRUCT_IMAGE_H


	#include "common_reconstruct_image.h"
    #include "image_structures.h"
    #include "reconstruct_image_tools.h"

    int check_is_occluded( nTupleVolume *imgVolOcc, int x, int y, int t);

    int check_disp_field(nTupleVolume *dispField, nTupleVolume *departVolume, nTupleVolume *arrivalVolume, nTupleVolume *occVol);
	
    void reconstruct_image(nTupleVolume* imgVol, nTupleVolume* imgVolModified, nTupleVolume* occVol,
            nTupleVolume* dispField, float sigmaColour, int reconstructionType=0, bool initialisation=false);

#endif
