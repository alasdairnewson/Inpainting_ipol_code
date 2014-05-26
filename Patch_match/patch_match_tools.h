//this function defines the tools necessary for the spatio-temporal patch match algorithm

#ifndef PATCH_MATCH_TOOLS
#define PATCH_MATCH_TOOLS

	#include "common_patch_match.h"
    #include "patch_match_measure.h"
    
    bool check_max_shift_distance(int xShift, int yShift, int tShift, const parameterStruct *params);

    int check_is_occluded( nTupleVolume *imgVolOcc, int x, int y, int t);
    
    void calclulate_patch_distances(nTupleVolume *departVolume,nTupleVolume *arrivalVolume,nTupleVolume *dispField, nTupleVolume *occVol,
		const parameterStruct *params);

    bool check_already_used_patch( nTupleVolume *dispField, int x, int y, int t, int dispX, int dispY, int dispT);

    int check_in_boundaries( nTupleVolume *imgVol, int x, int y, int t, const parameterStruct *params);

    int check_in_inner_boundaries( nTupleVolume *imgVol, int x, int y, int t, const parameterStruct *params);
    
    //full search
    void patch_match_full_search(nTupleVolume *dispField, nTupleVolume *imgVolA,nTupleVolume *imgVolB,
            nTupleVolume *occVol,nTupleVolume *modVol, const parameterStruct *params);
    
    //shift volume initialisation
    void initialise_displacement_field(nTupleVolume *dispField, nTupleVolume *departVolume, nTupleVolume *arrivalVolume, nTupleVolume *firstGuessVolume, nTupleVolume *occVol, const parameterStruct *params);
	
    //Random search
	int patch_match_random_search(nTupleVolume *dispField, nTupleVolume *imgVolA, nTupleVolume *imgVolB,
            nTupleVolume *occVol,  nTupleVolume *modVol, const parameterStruct *params);
    
    //propagation functions
    int patch_match_propagation(nTupleVolume *dispField, nTupleVolume *departVolume, nTupleVolume *arrivalVolume,
            nTupleVolume *occVol,  nTupleVolume *modVol,
		const parameterStruct *params, int iterationNb);

    int patch_match_long_propagation(nTupleVolume *dispField, nTupleVolume *imgVolA, nTupleVolume *imgVolB,
            nTupleVolume *occVol,  nTupleVolume *modVol, const parameterStruct *params);

    float calclulate_patch_error(nTupleVolume *departVolume,nTupleVolume *arrivalVolume,nTupleVolume *dispField, nTupleVolume *occVol,
		int xA, int yA, int tA, float minError, const parameterStruct *params);

	float get_min_correct_error(nTupleVolume *dispField,nTupleVolume *departVol,nTupleVolume *arrivalVol, nTupleVolume *occVol,
							int x, int y, int t, int beforeOrAfter, int *correctInd, float *minVector, float minError,
                            const parameterStruct *params);

	float ssd_minimum_value(nTupleVolume *imgVolA, nTupleVolume *imgVolB, nTupleVolume *occVol, int xA, int yA, int tA,
						int xB, int yB, int tB, float minVal, const parameterStruct *params);

    //utility functions
	int check_disp_field(nTupleVolume *dispField, nTupleVolume *departVolume, nTupleVolume *arrivalVolume,
            nTupleVolume *occVol, const parameterStruct *params);
#endif
