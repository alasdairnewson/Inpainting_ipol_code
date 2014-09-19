
//this function declares the patch match measure with which we compare patches


#ifndef PATCH_MATCH_MEASURE_H
#define PATCH_MATCH_MEASURE_H

    #include "common_patch_match.h"

    float ssd_patch_measure(nTupleVolume *imgVolA, nTupleVolume *imgVolB, nTupleVolume *dispField,
    nTupleVolume *occVol, int xA, int yA, int tA,int xB, int yB, int tB, float minVal, const patchMatchParameterStruct *params);

#endif
