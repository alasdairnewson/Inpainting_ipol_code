/*this function defines the functions which are tools used for the
colour estimation*/


#include "reconstruct_image_tools.h"

int check_disp_field(nTupleVolume *dispField, nTupleVolume *departVolume, nTupleVolume *arrivalVolume)
{
	int dispValX,dispValY,dispValT,hPatchSizeX,hPatchSizeY,hPatchSizeT;
	int xB,yB,tB;
	int i,j,k,returnVal;

	hPatchSizeX = (int)floor((float)((departVolume->patchSizeX)/2));	/*half the patch size*/
	hPatchSizeY = (int)floor((float)((departVolume->patchSizeY)/2));	/*half the patch size*/
	hPatchSizeT = (int)floor((float)((departVolume->patchSizeT)/2));	/*half the patch size*/

	returnVal = 0;
	for (k=hPatchSizeT; k< ((dispField->tSize) -hPatchSizeT); k++)
		for (j=hPatchSizeY; j< ((dispField->ySize) -hPatchSizeY); j++)
			for (i=hPatchSizeX; i< ((dispField->xSize) -hPatchSizeX); i++)
			{
				dispValX = (int)dispField->get_value(i,j,k,0);
				dispValY = (int)dispField->get_value(i,j,k,1);
				dispValT = (int)dispField->get_value(i,j,k,2);

				/*if ( (fabs(dispValX) > w) || (fabs(dispValY) > w) || (fabs(dispValT) > w))
				{
					MY_PRINTF("Error, the displacement is greater than the minimum value w : %d.\n",w);
					MY_PRINTF(" dispValX : %d\n dispValY : %d\n dispValT : %d\n",dispValX,dispValY,dispValT);
					returnVal= -1;
				}*/

				xB = dispValX + i;
				yB = dispValY + j;
				tB = dispValT + k;

				if ( (xB <hPatchSizeX) || (yB <hPatchSizeY) || (tB <hPatchSizeT) || 
					(xB >= (arrivalVolume->xSize - hPatchSizeX)) || (yB >= (arrivalVolume->ySize - hPatchSizeY)) || (tB >= (arrivalVolume->tSize - hPatchSizeT)))
				{
					MY_PRINTF("Error, the displacement is incorrect.\n");
					MY_PRINTF("xA : %d\n yA : %d\n tA : %d\n",i,j,k);
					MY_PRINTF(" dispValX : %d\n dispValY : %d\n dispValT : %d\n",dispValX,dispValY,dispValT);
					MY_PRINTF(" xB : %d\n yB : %d\n tB : %d\n",xB,yB,tB);
					returnVal= -1;
				}
				/*else if (check_is_occluded(occVol,xB,yB,tB) == 1)
				{
					MY_PRINTF("Error, the displacement leads to an occluded pixel.\n");
					MY_PRINTF(" xB : %d\n yB : %d\n tB : %d\n",xB,yB,tB);
					returnVal= -1;
				}*/
			}
	return(returnVal);
}

/*this function gets the nth percentile of the current weights*/
float get_adaptive_sigma(float *weights, int weightsLength, float sigmaPercentile)
{
    int i, weightsInd, percentileInd;
    float *weightsTemp, adaptiveSigmaOut;
    float percentile = (float)(sigmaPercentile)/((float)100);
    
    weightsTemp = (float*)malloc((size_t)weightsLength*sizeof(float));
    weightsInd = 0;
    for (i=0; i<weightsLength;i++)
    {
        if (weights[i] != -1)   /*we want to use this patch */
        {
            weightsTemp[weightsInd] = weights[i];
            weightsInd = weightsInd+1;
        }
    }
    weightsInd = weightsInd-1;
    std::sort(weightsTemp,weightsTemp+(weightsInd));
    
    percentileInd = (int)floor((float)percentile*weightsInd);
    
    adaptiveSigmaOut = sqrt(weightsTemp[percentileInd]);
    free(weightsTemp);
    return(adaptiveSigmaOut);
}

/*this function retieves the highest mode in the colour space of the
 different colours available for reconstructing a pixel*/
int estimate_best_colour(nTupleVolume *imgVol, float *weights, int weightsLength,
                                float *colours, int i, int j, int k)
{
    int ii;
    int minWeightInd;
    float minWeight;
    
    
    minWeight = INT_MAX;
    minWeightInd = 0;
    for (ii=0; ii<weightsLength; ii++)
    {
        if (weights[ii] != -1)
        {
			if (weights[ii] < minWeight)
            {
                minWeight = weights[ii];
                minWeightInd = ii;
            }
        }
    }

    imgVol->set_value(i,j,k,0,(imageDataType)(colours[minWeightInd]));
    imgVol->set_value(i,j,k,1,(imageDataType)(colours[minWeightInd+weightsLength]));
    imgVol->set_value(i,j,k,2,(imageDataType)(colours[minWeightInd+2*weightsLength]));
    
    return(1);
}
