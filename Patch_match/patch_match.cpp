//this function calculates the approximate nearest neighbour field for patch match
//in a volume, with multivalued pixels

#include "patch_match.h"

//this function calculates a nearest neighbour field, from imgVolA to imgVolB
void patch_match_ANN(nTupleVolume *imgVolA, nTupleVolume *imgVolB, 
        nTupleVolume *dispField, nTupleVolume *imgVolOcc,nTupleVolume *imgVolMod,
        const patchMatchParameterStruct *params,nTupleVolume *firstGuessVol)
{
	//decalarations
    long startTime;//clock_t startTime;
    long propagationTime,randomSearchTime;

	//check certain parameters
	if((imgVolA->nTupleSize) != (imgVolB->nTupleSize) )
	{
		MY_PRINTF("Error in patch_match_ANN, the size of the vector associated to each pixel is different for the two image volumes.");
		return;
	}
	if( (imgVolA->patchSizeX != (imgVolB->patchSizeX)) || (imgVolA->patchSizeY != (imgVolB->patchSizeY)) ||
		(imgVolA->patchSizeT != (imgVolB->patchSizeT))  )	//check that the patch sizes are equal
	{
		MY_PRINTF("Error in patch_match_ANN, the size of the patches are not equal in the two image volumes.");
		return;
	}
	if ( ( imgVolA->patchSizeX > imgVolA->xSize) || ( imgVolA->patchSizeY > imgVolA->ySize) || ( imgVolA->patchSizeT > imgVolA->tSize) ||
		( imgVolA->patchSizeX > imgVolB->xSize) || ( imgVolA->patchSizeY > imgVolB->ySize) || ( imgVolA->patchSizeT > imgVolB->tSize)
		)	//check that the patch size is less or equal to each dimension in the images
	{
		MY_PRINTF("Error in patch_match_ANN, the patch size is to large for one or more of the dimensions of the image volumes.");
		return;
	}
    
    //#pragma omp parallel
    //mexPrintf("Hello\n");
    
    //cuda_full_search(dispField,imgVolA,imgVolB,imgVolOcc,imgVolMod);
    //mexPrintf("dispField[1] : %f\n",dispField->values[1]);
    //return(dispField);
    propagationTime = 0.0;
    randomSearchTime = 0.0;
	long startTimeTotalPatchMatch = getMilliSecs();
    if (params->fullSearch == 1)
    {
        //memset(dispField->values,0,(size_t)(dispField->xSize)*(dispField->ySize)*(dispField->tSize)*(dispField->nTupleSize)*sizeof(float));
        if (firstGuessVol!= NULL)
    	{
        	initialise_displacement_field(dispField, imgVolA,imgVolB, firstGuessVol, imgVolOcc,params);
    	}
        //memset(dispField->values,0,(size_t)(dispField->xSize)*(dispField->ySize)*(dispField->tSize)*(dispField->nTupleSize)*sizeof(float));
        startTime = clock();
        patch_match_full_search(dispField, imgVolA,imgVolB, imgVolOcc,imgVolMod,params);
        propagationTime = propagationTime + double(clock() - startTime);
    }
    else    //normal patchMatch
    {
    	if (firstGuessVol!= NULL)
    	{
    		if ( (params->verboseMode) == true)
		    	MY_PRINTF("Initialisation\n");
		    long startTimeInitialisation = getMilliSecs();
		    initialise_displacement_field(dispField, imgVolA, imgVolB, firstGuessVol, imgVolOcc,params);
		    if ( (params->verboseMode) == true)
		    	MY_PRINTF("Initialisation time in s: %f\n",fabs(startTimeInitialisation-getMilliSecs())/1000);
	    }
        //show_nTuple_volume(dispField);
        if (check_disp_field(dispField, imgVolA, imgVolB,imgVolOcc,params) == -1)
            return;
        for (int i=0; i<(params->nIters); i++)
        {
        	patch_match_one_iteration_patch_level(dispField, imgVolA, imgVolB,
        	imgVolOcc, imgVolMod, params, i);
        }
    }
    if ( (params->verboseMode) == true)
    {
		MY_PRINTF("Total PatchMatch execution time in s: %f\n",fabs(startTimeTotalPatchMatch-getMilliSecs())/1000);
	}

}
