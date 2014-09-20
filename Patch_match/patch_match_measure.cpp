
//this function defines the patch match measure with which we compare patches

#include "patch_match_measure.h"

#define LAMDA_MEASURE 0

float ssd_patch_measure(nTupleVolume *imgVolA, nTupleVolume *imgVolB, nTupleVolume *occVol, int xA, int yA, int tA,
int xB, int yB, int tB, float minVal, const patchMatchParameterStruct *params)
{
	//declarations
	int i,j,k, p,xAtemp, yAtemp, tAtemp, xBtemp, yBtemp, tBtemp;
	int xMinA,yMinA,tMinA;
	int xMinB,yMinB,tMinB;
    int sumOcc, occA;
	imageDataType tempVal;
    float beta = 50.0;

	imageDataType ssd = 0;

	if ( ((imgVolA->patchSizeX) != (imgVolB->patchSizeX)) || ((imgVolA->patchSizeY) != (imgVolB->patchSizeY)) ||
		((imgVolA->patchSizeT) != (imgVolB->patchSizeT)) )
	{
		MY_PRINTF("Error in ssd_minimum_value, the patch sizes are not equal.\n");
		return -1;
	}

    xMinA = xA-imgVolA->hPatchSizeX;
    yMinA = yA-imgVolA->hPatchSizeY;
    tMinA = tA-imgVolA->hPatchSizeT;
    xMinB = xB-imgVolB->hPatchSizeX;
    yMinB = yB-imgVolB->hPatchSizeY;
    tMinB = tB-imgVolB->hPatchSizeT;
    
    sumOcc = 0;
    
    if (params->partialComparison)
    {
        for (k=0; k<imgVolA->patchSizeT; k++)
            for (j=0; j<imgVolA->patchSizeY; j++)
                for (i=0; i<imgVolA->patchSizeX; i++)
                {
                    xAtemp = xMinA + i;
                    yAtemp = yMinA + j;
                    tAtemp = tMinA + k;
                    xBtemp = xMinB + i;
                    yBtemp = yMinB + j;
                    tBtemp = tMinB + k;
                    //do not compare the edges in any case
                    if ((!check_in_boundaries( imgVolA, xAtemp, yAtemp, tAtemp)))
                        continue;   //if we are not in the boundaries, do not compare

                    sumOcc = sumOcc + (int)( (!(occVol->get_value(xAtemp,yAtemp,tAtemp,0))) == 1);

                }
    }
    else    //calculate patch size
    {
        int patchSizeXtemp = min_int(xA + imgVolA->hPatchSizeX,imgVolA->xSize-1) - max_int(xA - imgVolA->hPatchSizeX,0) + 1;
        int patchSizeYtemp = min_int(yA + imgVolA->hPatchSizeY,imgVolA->ySize-1) - max_int(yA - imgVolA->hPatchSizeY,0) + 1;
        int patchSizeTtemp = min_int(tA + imgVolA->hPatchSizeT,imgVolA->tSize-1) - max_int(tA - imgVolA->hPatchSizeT,0) + 1;
        
        sumOcc = patchSizeXtemp * patchSizeYtemp * patchSizeTtemp;
    }
	sumOcc = max_int(sumOcc,1);
    
	for (k=0; k<imgVolA->patchSizeT; k++)
		for (j=0; j<imgVolA->patchSizeY; j++)
			for (i=0; i<imgVolA->patchSizeX; i++)
			{
				xAtemp = xMinA + i;
				yAtemp = yMinA + j;
				tAtemp = tMinA + k;
                
				xBtemp = xMinB + i;
				yBtemp = yMinB + j;
				tBtemp = tMinB + k;
                

                //do not compare if we are not in the boundaries
                if ((!check_in_boundaries( imgVolA, xAtemp, yAtemp, tAtemp)))
                    occA = 1;
                else
                    occA = 0;
                if (occA == 1)
                    continue;   //we do not wish to compare this pixel
                /*if we want partial patch comparison*/
                if (params->partialComparison && occVol->xSize >0)
                    occA = (int)(*(occVol->get_value_ptr(xAtemp, yAtemp, tAtemp,0)) == 1);
                if (occA == 1)
                    continue;   //we do not wish to compare this pixel
                
                /* similarity */
				for (p=0; p<imgVolA->nTupleSize; p++)
				{
					imageDataType imgVolAval = imgVolA->get_image_value(xAtemp, yAtemp,p);
					imageDataType imgVolBval = imgVolB->get_image_value(xBtemp, yBtemp,p);
					tempVal = (imgVolAval) - (imgVolBval);
					ssd = ssd + (imageDataType)(((tempVal)*(tempVal))/sumOcc);
                    //ssd = ssd + (abs(tempFloat))/sumOcc;
				}
                
                if( params->normGradX != NULL)
                {            
                    imageDataType normGradXtemp = (params->normGradX)->get_value(xAtemp,yAtemp,tAtemp,0) 
                    							- (params->normGradX)->get_value(xBtemp,yBtemp,tBtemp,0);
                                        
                    imageDataType normGradYtemp = (params->normGradY)->get_value(xAtemp,yAtemp,tAtemp,0) 
                    							- (params->normGradY)->get_value(xBtemp,yBtemp,tBtemp,0);

                    ssd = ssd + beta*normGradXtemp*normGradXtemp/sumOcc;
                    ssd = ssd + beta*normGradYtemp*normGradYtemp/sumOcc;
                }
                
		if ((minVal != -1) && (ssd > minVal))
                {
			return(-1);
                }
			}
    
	return(ssd);
}
