//this function defines the functions which are tools used for the
//spatio-temporal patch_match algorithm


#include "patch_match_tools.h"

//check if the displacement values have already been used
bool check_already_used_patch( nTupleVolume *dispField, int x, int y, int t, int dispX, int dispY, int dispT)
{
    
    if ( (((int)dispField->get_value(x,y,t,0)) == dispX) && 
            (((int)dispField->get_value(x,y,t,1)) == dispY) &&
            (((int)dispField->get_value(x,y,t,2)) == dispT)
            )
        return 1;
    else
        return 0;
    
}

//see if the maximum shift distance is respected
bool check_max_shift_distance(int xShift, int yShift, int tShift, const patchMatchParameterStruct *params)
{
	int distance;
	
	distance = (int)floor( (float)sqrt( (float)xShift*xShift + yShift*yShift + tShift*tShift) );
	
	if (params->maxShiftDistance != -1)
		return(distance <= params->maxShiftDistance);
	else
		return(true);
}


//check if the pixel is occluded
int check_is_occluded( nTupleVolume *imgVolOcc, int x, int y, int t)
{
	if (imgVolOcc->xSize == 0)
		return 0;
	if ( (imgVolOcc->get_value(x,y,t,0)) > 0)
		return 1;
	else
		return 0;
}

void calclulate_patch_distances(nTupleVolume *departVolume,nTupleVolume *arrivalVolume,nTupleVolume *dispField, nTupleVolume *occVol,
		const patchMatchParameterStruct *params)
{
	for (int k=0; k< (dispField->tSize); k++)
		for (int j=0; j< (dispField->ySize); j++)
			for (int i=0; i< (dispField->xSize); i++)
			{
				if (check_in_inner_boundaries(departVolume, i, j, k,params) == 1)
				{
					int xShift,yShift,tShift;
					xShift = dispField->get_value(i,j,k,0);
					yShift = dispField->get_value(i,j,k,1);
					tShift = dispField->get_value(i,j,k,2);
					float ssdTemp = ssd_patch_measure(departVolume, arrivalVolume,dispField,occVol, i, j, k, i+xShift, j+yShift, k+tShift, -1, params);
					dispField->set_value(i,j,k,3,ssdTemp);
				}
				else
					dispField->set_value(i,j,k,3,FLT_MAX);
			
			}
			

}

float calclulate_patch_error(nTupleVolume *departVolume,nTupleVolume *arrivalVolume,nTupleVolume *dispField, nTupleVolume *occVol,
		int xA, int yA, int tA, float minError, const patchMatchParameterStruct *params)
{
	int xB, yB, tB;
	float errorOut;

	xB = (xA) + (int)dispField->get_value(xA,yA,tA,0);
	yB = (yA) + (int)dispField->get_value(xA,yA,tA,1);
	tB = (tA) + (int)dispField->get_value(xA,yA,tA,2);

	errorOut = ssd_patch_measure(departVolume, arrivalVolume,dispField,occVol, xA, yA, tA, xB, yB, tB, minError, params);
	return(errorOut);
}

int check_disp_field(nTupleVolume *dispField, nTupleVolume *departVolume, nTupleVolume *arrivalVolume, nTupleVolume *occVol, const patchMatchParameterStruct *params)
{
	int dispValX,dispValY,dispValT,hPatchSizeX,hPatchSizeY,hPatchSizeT;
	int xB,yB,tB;
	int i,j,k,returnVal;

	hPatchSizeX = (int)floor((float)((departVolume->patchSizeX)/2));	//half the patch size
	hPatchSizeY = (int)floor((float)((departVolume->patchSizeY)/2));	//half the patch size
	hPatchSizeT = (int)floor((float)((departVolume->patchSizeT)/2));	//half the patch size

	returnVal = 0;
	for (k=hPatchSizeT; k< ((dispField->tSize) -hPatchSizeT); k++)
		for (j=hPatchSizeY; j< ((dispField->ySize) -hPatchSizeY); j++)
			for (i=hPatchSizeX; i< ((dispField->xSize) -hPatchSizeX); i++)
			{
				dispValX = (int)dispField->get_value(i,j,k,0);
				dispValY = (int)dispField->get_value(i,j,k,1);
				dispValT = (int)dispField->get_value(i,j,k,2);


				xB = dispValX + i;
				yB = dispValY + j;
				tB = dispValT + k;

				if ( check_in_inner_boundaries(arrivalVolume, xB, yB, tB,params) == 0)
				{
					MY_PRINTF("Error, the displacement is incorrect.\n");
					MY_PRINTF("xA : %d\n yA : %d\n tA : %d\n",i,j,k);
					MY_PRINTF(" dispValX : %d\n dispValY : %d\n dispValT : %d\n",dispValX,dispValY,dispValT);
					MY_PRINTF(" xB : %d\n yB : %d\n tB : %d\n",xB,yB,tB);
					returnVal= -1;
				}
				else if (check_is_occluded(occVol,xB,yB,tB) == 1)
				{
					MY_PRINTF("Error, the displacement leads to an occluded pixel.\n");
					MY_PRINTF(" xB : %d\n yB : %d\n tB : %d\n",xB,yB,tB);
					returnVal= -1;
				}
			}
	return(returnVal);

}

void patch_match_full_search(nTupleVolume *dispField, nTupleVolume *imgVolA,nTupleVolume *imgVolB,
        nTupleVolume *occVol, nTupleVolume *modVol, const patchMatchParameterStruct *params)
{
    float minSSD,ssdTemp;
    int hPatchSizeX, hPatchSizeY, hPatchSizeT;
    int bestX, bestY, bestT;
    
    hPatchSizeX = (int)floor((float)((dispField->patchSizeX)/2));	//half the patch size
	hPatchSizeY = (int)floor((float)((dispField->patchSizeY)/2));	//half the patch size
	hPatchSizeT = (int)floor((float)((dispField->patchSizeT)/2));	//half the patch size
            
    //#pragma omp parallel for shared(dispField, occVol, imgVolA, imgVolB) private(kk,jj,ii,minSSD,ssdTemp,bestX,bestY,bestT,i,j,k)
	for (int k=0; k< ((dispField->tSize)); k++)
		for (int j=0; j< ((dispField->ySize)); j++)
        {
			for (int i=0; i< ((dispField->xSize)); i++)
            {
                minSSD = FLT_MAX;
                bestX = INT_MAX;
                bestY = INT_MAX;
                bestT = INT_MAX;
                if (modVol->xSize >0)
                    if (modVol->get_value(i,j,k,0) == 0)   //if we don't want to modify this match
                       continue;
                //search for the best match
                for (int kk=hPatchSizeT; kk< ((dispField->tSize) -hPatchSizeT); kk++)
                    for (int jj=hPatchSizeY; jj< ((dispField->ySize) -hPatchSizeY); jj++)
                        for (int ii=hPatchSizeX; ii< ((dispField->xSize) -hPatchSizeX); ii++)
                        {
                            if (occVol->xSize >0)
                                if (check_is_occluded(occVol,ii,jj,kk))   //if this pixel is occluded, continue
                                    continue;
                                    
                            if (check_max_shift_distance(ii-i,jj-j,kk-k,params) == false)
                            	continue;

                            ssdTemp = ssd_patch_measure(imgVolA, imgVolB, dispField,occVol, i, j, k,ii, jj, kk, minSSD,params);
                            if ( (ssdTemp != -1) && (ssdTemp <= minSSD))   //we have a new best match
                            {
                                minSSD = ssdTemp;
                                bestX = ii - i;
                                bestY = jj - j;
                                bestT = kk - k;
                            }
                        }
                if (bestX==FLT_MAX || bestY==FLT_MAX || bestT==FLT_MAX)
                {
                	MY_PRINTF("Here patch_match_full_search. Error, a correct nearest neighbour was not found for the patch centred at :\nx:%d\ny:%d\nt:%d\n\n",i,j,k);
                }
                dispField->set_value(i,j,k,0,(imageDataType)bestX);
                dispField->set_value(i,j,k,1,(imageDataType)bestY);
                dispField->set_value(i,j,k,2,(imageDataType)bestT);
                dispField->set_value(i,j,k,3,minSSD);
            }
        }

}

void initialise_displacement_field(nTupleVolume *dispField, nTupleVolume *departVolume, 
            nTupleVolume *arrivalVolume, nTupleVolume *firstGuessVolume, nTupleVolume *occVol, const patchMatchParameterStruct *params)
{
	//declarations
	int xDisp, yDisp, tDisp;
	int xMin,xMax,yMin,yMax,tMin,tMax;
	int xFirst,yFirst,tFirst;
	int isNotOcc;
    float ssdTemp;

	int hPatchSizeX = (int)floor(((float)arrivalVolume->patchSizeX)/2);
	int hPatchSizeY = (int)floor(((float)arrivalVolume->patchSizeY)/2);
	int hPatchSizeT = (int)floor(((float)arrivalVolume->patchSizeT)/2);

	int hPatchSizeCeilX = (int)ceil(((float)arrivalVolume->patchSizeX)/2);
	int hPatchSizeCeilY = (int)ceil(((float)arrivalVolume->patchSizeY)/2);
	int hPatchSizeCeilT = (int)ceil(((float)arrivalVolume->patchSizeT)/2);

	for (int i=0; i< (dispField->xSize); i++)
		for (int j=0; j< (dispField->ySize); j++)
			for (int k=0; k< (dispField->tSize); k++)
			{
				isNotOcc = 0;
                //if there is a valid first guess
                while(isNotOcc == 0)
                {
                    //if there is a first guess, and it is in the inner boundaries, and respects the minimum shift distance
                    if ( (firstGuessVolume->xSize >0) && (check_in_inner_boundaries(arrivalVolume,i+(int)firstGuessVolume->get_value(i,j,k,0),j+
                        (int)firstGuessVolume->get_value(i,j,k,1),k+(int)firstGuessVolume->get_value(i,j,k,2),params )) &&
                        (check_max_shift_distance((int)firstGuessVolume->get_value(i,j,k,0),
                        (int)firstGuessVolume->get_value(i,j,k,1),(int)firstGuessVolume->get_value(i,j,k,2),params ))
                        )
                    {
                        //if it is not occluded, we take the initial first guess and continue
                        if (!check_is_occluded(occVol,i+(int)firstGuessVolume->get_value(i,j,k,0),j+(int)firstGuessVolume->get_value(i,j,k,1),
                            k+(int)firstGuessVolume->get_value(i,j,k,2)) )
                        {
                            xDisp = (int)firstGuessVolume->get_value(i,j,k,0);
                            yDisp = (int)firstGuessVolume->get_value(i,j,k,1);
                            tDisp = (int)firstGuessVolume->get_value(i,j,k,2);
                            isNotOcc = 1;
                            continue;
                        }
                        else    //otherwise, we set up the calculation of a random initial starting point, centred on the initial guess
                        {
                            xFirst = i+(int)firstGuessVolume->get_value(i,j,k,0);
                            yFirst = j+(int)firstGuessVolume->get_value(i,j,k,1);
                            tFirst = k+(int)firstGuessVolume->get_value(i,j,k,2);
                            xMin = max_int(xFirst-params->w,hPatchSizeX);
                            xMax = min_int(xFirst+params->w,arrivalVolume->xSize - hPatchSizeX -1);
                            yMin = max_int(yFirst-params->w,hPatchSizeY);
                            yMax = min_int(yFirst+params->w,arrivalVolume->ySize - hPatchSizeY -1);
                            tMin = max_int(tFirst-params->w,hPatchSizeT);
                            tMax = min_int(tFirst+params->w,arrivalVolume->tSize - hPatchSizeT -1);
                            dispField->set_value(i,j,k,0,-1);
                            dispField->set_value(i,j,k,1,-1);
                            dispField->set_value(i,j,k,2,-1);
                            dispField->set_value(i,j,k,3,-1);
                            
                        }
                    }
                    else    //by default, set the displacement to float_max
                    {
                        dispField->set_value(i,j,k,0,(imageDataType)FLT_MAX);
                        dispField->set_value(i,j,k,1,(imageDataType)FLT_MAX);
                        dispField->set_value(i,j,k,2,(imageDataType)FLT_MAX);
                    }   
                    if (arrivalVolume->xSize == arrivalVolume->patchSizeX)	//special case where the patch size is the size of the dimension
                    {
                        xDisp = 0;
                    }
                    else{
                        if ( (dispField->get_value(i,j,k,0) == FLT_MAX) || (firstGuessVolume->xSize == 0))  //default behaviour
                        {
                            xDisp = ((rand()%( (arrivalVolume->xSize) -2*hPatchSizeCeilX-1)) + hPatchSizeX)-i;
                        }
                        else    //based on an initial guess
                        {
                            xDisp = (int)((int)round_float(rand_float_range((float)(xMin),(float)(xMax))) - i);
                        }
                    }
                    if (arrivalVolume->ySize == arrivalVolume->patchSizeY)	//special case where the patch size is the size of the dimension
                    {
                        yDisp = 0;
                    }
                    else{
                        if ( (dispField->get_value(i,j,k,1) == FLT_MAX) || (firstGuessVolume->xSize == 0))  //default behaviour
                        {
                            yDisp = ((rand()%( (arrivalVolume->ySize) -2*hPatchSizeCeilY-1)) + hPatchSizeY)-j;
                        }
                        else    //based on an initial guess
                        {
                            yDisp = (int)(round_float(rand_float_range((float)(yMin),(float)(yMax))) - j);
                        }
                    }
                    if (arrivalVolume->tSize == arrivalVolume->patchSizeT)	//special case where the patch size is the size of the dimension
                    {
                        tDisp = 0;
                    }
                    else{
                        if ( (dispField->get_value(i,j,k,2) == FLT_MAX) || (firstGuessVolume->xSize == 0))  //default behaviour
                        {
                            tDisp = ((rand()%( (arrivalVolume->tSize) -2*hPatchSizeCeilT-1)) + hPatchSizeT)-k;
                        }
                        else    //based on an initial guess
                        {
                            tDisp = (int)(round_float(rand_float_range((float)(tMin),(float)(tMax))) - k);
                        }
                    }

                    isNotOcc = (!(check_is_occluded(occVol,xDisp+i,yDisp+j,tDisp+k))
                             &&(check_in_inner_boundaries(arrivalVolume,xDisp+i,yDisp+j,tDisp+k,params))
                             &&(check_max_shift_distance(xDisp,yDisp,tDisp+k,params))
                             );
                }
                //if everything is all right, set the displacements
                dispField->set_value(i,j,k,0, (imageDataType)(xDisp));
                dispField->set_value(i,j,k,1, (imageDataType)(yDisp));
                dispField->set_value(i,j,k,2, (imageDataType)(tDisp));

                if (check_in_inner_boundaries(departVolume,i,j,k,params))
                {
                    ssdTemp = ssd_patch_measure(departVolume, arrivalVolume, dispField,occVol, i, j, k, i+xDisp, j+yDisp, k+tDisp, -1,params);
                    if(ssdTemp ==-1)
                        ssdTemp = FLT_MAX;
                }
                else
                    ssdTemp = FLT_MAX;
                dispField->set_value(i,j,k,3,(imageDataType)ssdTemp); //set the ssd error
            }
}

int patch_match_random_search(nTupleVolume *dispField, nTupleVolume *imgVolA, nTupleVolume *imgVolB,
        nTupleVolume *occVol, nTupleVolume *modVol, const patchMatchParameterStruct *params)
{
	//create random number seed
	int xRand,yRand,tRand;
	int randMinX,randMaxX,randMinY,randMaxY,randMinT,randMaxT;
	int hPatchSizeX,hPatchSizeY,hPatchSizeT, zMax;
	int xTemp,yTemp,tTemp,wTemp,wMax;
    int *xDisp, *yDisp, *tDisp;
	float ssdTemp;
    int nbModified = 0;

	hPatchSizeX = (int)floor((float)((dispField->patchSizeX)/2));	//half the patch size
	hPatchSizeY = (int)floor((float)((dispField->patchSizeY)/2));	//half the patch size
	hPatchSizeT = (int)floor((float)((dispField->patchSizeT)/2));	//half the patch size
    
    xDisp = new int[1];
    yDisp = new int[1];
    tDisp = new int[1];

	//calculate the maximum z (patch search index)
    wMax = min_int(params->w, max_int(max_int(imgVolB->xSize,imgVolB->ySize),imgVolB->tSize));
	zMax = (int)ceil((float) (- (log((float)(wMax)))/(log((float)(params->alpha)))) );
    
    nTupleVolume *wValues = new nTupleVolume(1,zMax,1,1,imgVolA->indexing);
    //store the values of the maximum search parameters
    for (int z=0; z<zMax; z++)
    {
        wValues->set_value(z,0,0,0,
                (imageDataType)round_float((params->w)*((float)pow((float)params->alpha,z)))
                );
    }

	for (int k=0; k< ((dispField->tSize) ); k++)
	{
		for (int j=0; j< ((dispField->ySize) ); j++)
			for (int i=0; i< ((dispField->xSize) ); i++)
			{
                if (modVol->xSize >0)
                    if (modVol->get_value(i,j,k,0) == 0)   //if we don't want to modify this match
                        continue;
				ssdTemp = dispField->get_value(i,j,k,3); //get the saved ssd value
                
                for (int z=0; z<zMax; z++)	//test for different search indices
                {
                    xTemp = i+(int)dispField->get_value(i,j,k,0);	//get the arrival position of the current offset
                    yTemp = j+(int)dispField->get_value(i,j,k,1);	//get the arrival position of the current offset
                    tTemp = k+(int)dispField->get_value(i,j,k,2);	//get the arrival position of the current offset

                    wTemp = wValues->get_value(z,0,0,0);
                    // X values
                    randMinX = max_int(xTemp - wTemp,hPatchSizeX);
                    randMaxX = min_int(xTemp + wTemp,imgVolB->xSize - hPatchSizeX - 1);
                    // Y values
                    randMinY = max_int(yTemp - wTemp,hPatchSizeY);
                    randMaxY = min_int(yTemp + wTemp,imgVolB->ySize - hPatchSizeY - 1);
                    // T values
                    randMinT = max_int(tTemp - wTemp,hPatchSizeT);
                    randMaxT = min_int(tTemp + wTemp,imgVolB->tSize - hPatchSizeT - 1);

                    //new positions in the image imgB
                    xRand = rand_int_range(randMinX, randMaxX);	//random values between xMin and xMax, clamped to the sizes of the image B
                    yRand = rand_int_range(randMinY, randMaxY);	//random values between yMin and yMax, clamped to the sizes of the image B
                    tRand = rand_int_range(randMinT, randMaxT);	//random values between tMin and tMax, clamped to the sizes of the image B

                    if (check_is_occluded(occVol,xRand,yRand,tRand))
                        continue;	//the new position is occluded
                    if (check_in_inner_boundaries(imgVolB,xRand,yRand,tRand,params) == 0)
                        continue;	//the new position is not in the inner boundaries
                    if (check_max_shift_distance( (xRand-i),(yRand-j),(tRand-k),params) == false)
                    	continue;	//the new position is too far away

                    ssdTemp =  ssd_patch_measure(imgVolA, imgVolB, dispField,occVol, i, j, k, xRand, yRand, tRand, ssdTemp,params);

                    if (ssdTemp != -1)	//we have a better match
                    {
                        dispField->set_value(i,j,k,0, (imageDataType)(xRand-i));
                        dispField->set_value(i,j,k,1, (imageDataType)(yRand-j));
                        dispField->set_value(i,j,k,2, (imageDataType)(tRand-k));
                        dispField->set_value(i,j,k,3, (imageDataType)(ssdTemp));

                        nbModified = nbModified+1;
                    }
                    else
                        ssdTemp = dispField->get_value(i,j,k,3); //set the saved ssd value bakc to its proper (not -1) value
                }
			}
	}

    delete xDisp;
    delete yDisp;
    delete tDisp;
    delete wValues;
    return(nbModified);
}

//one iteration of the propagation of the patch match algorithm
int patch_match_propagation(nTupleVolume *dispField, nTupleVolume *departVolume, nTupleVolume *arrivalVolume, nTupleVolume *occVol,  
        nTupleVolume *modVol, const patchMatchParameterStruct *params, int iterationNb)
{
	//declarations
	int *correctInd;
    int nbModified;
	float currentError, *minVector;
	
	correctInd = (int*)malloc((size_t)sizeof(int));

	minVector = (float*)malloc((size_t)3*sizeof(float));
	
    
	//iterate over all displacements (do not go to edge, where the patches may not be defined)
    nbModified = 0;
	if (iterationNb&1)	//if we are on an odd iteration
	{
        for (int k=((dispField->tSize)-1); k>= 0; k--)
			for (int j=((dispField->ySize) -1); j>= 0; j--)
                for (int i=((dispField->xSize) -1); i>= 0; i--)
				{
                    if (modVol->xSize >0)
                        if (modVol->get_value(i,j,k,0) == 0)   //if we don't want to modify this match
                            continue;
					
                    //calculate the error of the current displacement
					currentError = dispField->get_value(i,j,k,3);
                    
					get_min_correct_error(dispField,departVolume,arrivalVolume,occVol,
							i, j, k, iterationNb&1, correctInd,minVector,currentError,params);
					//if the best displacement is the current one. Note : we have taken into account the case
					//where none of the diplacements around the current pixel are valid
					if (*correctInd == -1)	//if the best displacement is the current one
					{
						dispField->set_value(i,j,k,3,currentError);
					}
					else	//we copy the displacement from another better one
					{
						if ((*correctInd) == 0){
							copy_pixel_values_nTuple_volume(dispField,dispField, min_int(i+1,((int)dispField->xSize)-1), j, k, i, j, k);
                            nbModified++;
                        }

						else if((*correctInd) == 1){
							copy_pixel_values_nTuple_volume(dispField,dispField, i, min_int(j+1,((int)dispField->ySize)-1), k, i, j, k);
                            nbModified++;
                        }
						else if( (*correctInd) == 2){
							copy_pixel_values_nTuple_volume(dispField,dispField, i, j, min_int(k+1,((int)dispField->tSize)-1), i, j, k);
                            nbModified++;
                        }
                        else
                            MY_PRINTF("Error, correct ind not chosen\n.");
						//now calculate the error of the patch matching
						currentError = calclulate_patch_error(departVolume,arrivalVolume,dispField,occVol,i,j,k, -1,params);
						dispField->set_value(i,j,k,3,currentError);
					}
				}
	}
	else 	//if we are on an even iteration
	{
		for (int k=0; k< ((dispField->tSize) ); k++)
			for (int j=0; j< ((dispField->ySize) ); j++)
				for (int i=0; i< ((dispField->xSize) ); i++)
				{
                    if (modVol->xSize >0)
                        if (modVol->get_value(i,j,k,0) == 0)   //if we don't want to modify this match
                            continue;
                    
					//calculate the error of the current displacement
					currentError = dispField->get_value(i,j,k,3);//calclulate_patch_error(departVolume,arrivalVolume,dispField,occVol,i,j,k, -1,params);
					//get the upper, left and before patch distances
					get_min_correct_error(dispField,departVolume,arrivalVolume,occVol,
							i, j, k, iterationNb&1, correctInd, minVector,currentError,params);

					//if the best displacement is the current one. Note : we have taken into account the case
					//where none of the diplacements around the current pixel are valid
					if (*correctInd == -1)
					{
						dispField->set_value(i,j,k,3,currentError);
                    }
					else	//we copy the displacement from another better one
					{
						if ( (*correctInd) == 0){
							copy_pixel_values_nTuple_volume(dispField,dispField, max_int(i-1,0), j, k, i, j, k);
                            nbModified++;
                        }

						else if( (*correctInd) == 1){
							copy_pixel_values_nTuple_volume(dispField,dispField, i, max_int(j-1,0), k, i, j, k);
                            nbModified++;
                        }

						else if( (*correctInd) == 2){
							copy_pixel_values_nTuple_volume(dispField,dispField, i, j, max_int(k-1,0), i, j, k);
                            nbModified++;
                        }
                        else
                            MY_PRINTF("Error, correct ind not chosen\n.");
						//now calculate the error of the patch matching
						currentError = calclulate_patch_error(departVolume,arrivalVolume,dispField,occVol,i,j,k, -1,params);
						dispField->set_value(i,j,k,3,currentError);
					}
				}
	}
	free(correctInd);
	free(minVector);
    
    return(nbModified);
}

/******************************************/
/******************************************/
/******   PATCH LEVEL INTERLEAVING   ******/
/******************************************/
/******************************************/

int patch_match_one_iteration_patch_level(nTupleVolume *dispField, nTupleVolume *departVolume, nTupleVolume *arrivalVolume,
        nTupleVolume *occVol, nTupleVolume *modVol, const patchMatchParameterStruct *params, int iterationNb)
{
	int wMax, zMax;
	//calculate the maximum z (patch search index)
    wMax = min_int(params->w, max_int(max_int(arrivalVolume->xSize,arrivalVolume->ySize),arrivalVolume->tSize));
	zMax = (int)ceil((float) (- (log((float)(wMax)))/(log((float)(params->alpha)))) );
    
    nTupleVolume *wValues = new nTupleVolume(1,zMax,1,1,departVolume->indexing);
    //store the values of the maximum search parameters
    for (int z=0; z<zMax; z++)
    {
        wValues->set_value(z,0,0,0,
                (imageDataType)round_float((params->w)*((float)pow((float)params->alpha,z)))
                );
    }

	for (int k=0; k< ((dispField->tSize) ); k++)
		for (int j=0; j< ((dispField->ySize) ); j++)
			for (int i=0; i< ((dispField->xSize) ); i++)
			{
				//propagation
				patch_match_propagation_patch_level(dispField, departVolume, arrivalVolume, occVol,  
        		modVol, params, iterationNb, i, j, k);
				
				//random search
				patch_match_random_search_patch_level(dispField, departVolume, arrivalVolume,
        		occVol, modVol, params, i, j, k, wValues);
			}
	delete wValues;
}

int patch_match_random_search_patch_level(nTupleVolume *dispField, nTupleVolume *imgVolA, nTupleVolume *imgVolB,
        nTupleVolume *occVol, nTupleVolume *modVol, const patchMatchParameterStruct *params, int i, int j, int k,
        nTupleVolume *wValues)
{
	//create random number seed
	int xRand,yRand,tRand;
	int randMinX,randMaxX,randMinY,randMaxY,randMinT,randMaxT;
	int hPatchSizeX,hPatchSizeY,hPatchSizeT;
	int xTemp,yTemp,tTemp,wTemp;
    int *xDisp, *yDisp, *tDisp;
	float ssdTemp;
    int nbModified = 0;

	hPatchSizeX = (int)floor((float)((dispField->patchSizeX)/2));	//half the patch size
	hPatchSizeY = (int)floor((float)((dispField->patchSizeY)/2));	//half the patch size
	hPatchSizeT = (int)floor((float)((dispField->patchSizeT)/2));	//half the patch size
    
    xDisp = new int[1];
    yDisp = new int[1];
    tDisp = new int[1];

	if (modVol->xSize >0)
		if (modVol->get_value(i,j,k,0) == 0)   //if we don't want to modify this match
			return(0);
	ssdTemp = dispField->get_value(i,j,k,3); //get the saved ssd value
	
	for (int z=0; z<(wValues->xSize); z++)	//test for different search indices
	{
		xTemp = i+(int)dispField->get_value(i,j,k,0);	//get the arrival position of the current offset
		yTemp = j+(int)dispField->get_value(i,j,k,1);	//get the arrival position of the current offset
		tTemp = k+(int)dispField->get_value(i,j,k,2);	//get the arrival position of the current offset

		wTemp = wValues->get_value(z,0,0,0);
		// X values
		randMinX = max_int(xTemp - wTemp,hPatchSizeX);
		randMaxX = min_int(xTemp + wTemp,imgVolB->xSize - hPatchSizeX - 1);
		// Y values
		randMinY = max_int(yTemp - wTemp,hPatchSizeY);
		randMaxY = min_int(yTemp + wTemp,imgVolB->ySize - hPatchSizeY - 1);
		// T values
		randMinT = max_int(tTemp - wTemp,hPatchSizeT);
		randMaxT = min_int(tTemp + wTemp,imgVolB->tSize - hPatchSizeT - 1);

		//new positions in the image imgB
		xRand = rand_int_range(randMinX, randMaxX);	//random values between xMin and xMax, clamped to the sizes of the image B
		yRand = rand_int_range(randMinY, randMaxY);	//random values between yMin and yMax, clamped to the sizes of the image B
		tRand = rand_int_range(randMinT, randMaxT);	//random values between tMin and tMax, clamped to the sizes of the image B

		if (check_is_occluded(occVol,xRand,yRand,tRand))
			continue;	//the new position is occluded
		if (check_in_inner_boundaries(imgVolB,xRand,yRand,tRand,params) == 0)
			continue;	//the new position is not in the inner boundaries
		if (check_max_shift_distance( (xRand-i),(yRand-j),(tRand-k),params) == false)
			continue;	//the new position is too far away

		ssdTemp =  ssd_patch_measure(imgVolA, imgVolB, dispField,occVol, i, j, k, xRand, yRand, tRand, ssdTemp,params);

		if (ssdTemp != -1)	//we have a better match
		{
			dispField->set_value(i,j,k,0, (imageDataType)(xRand-i));
			dispField->set_value(i,j,k,1, (imageDataType)(yRand-j));
			dispField->set_value(i,j,k,2, (imageDataType)(tRand-k));
			dispField->set_value(i,j,k,3, (imageDataType)(ssdTemp));

			nbModified = nbModified+1;
		}
		else
			ssdTemp = dispField->get_value(i,j,k,3); //set the saved ssd value bakc to its proper (not -1) value
	}
	
    delete xDisp;
    delete yDisp;
    delete tDisp;
    return(nbModified);
}

//one iteration of the propagation of the patch match algorithm, for a SINGLE patch
int patch_match_propagation_patch_level(nTupleVolume *dispField, nTupleVolume *departVolume, nTupleVolume *arrivalVolume, nTupleVolume *occVol,  
        nTupleVolume *modVol, const patchMatchParameterStruct *params, int iterationNb, int i, int j, int k)
{
	//declarations
	int *correctInd;
    int nbModified=0;
	float currentError, *minVector;
	
	correctInd = new int;

	minVector = (float*)malloc((size_t)3*sizeof(float));

	//calculate the error of the current displacement
	currentError = dispField->get_value(i,j,k,3);
                    
	get_min_correct_error(dispField,departVolume,arrivalVolume,occVol,
	i, j, k, iterationNb&1, correctInd,minVector,currentError,params);
	
	//if the best displacement is the current one. Note : we have taken into account the case
	//where none of the diplacements around the current pixel are valid
	if (*correctInd == -1)	//if the best displacement is the current one
	{
		dispField->set_value(i,j,k,3,currentError);
		return(0);
	}
	if (iterationNb&1)	//if we are on an odd iteration
	{
		if ((*correctInd) == 0){
			copy_pixel_values_nTuple_volume(dispField,dispField, min_int(i+1,((int)dispField->xSize)-1), j, k, i, j, k);
			nbModified++;
		}

		else if((*correctInd) == 1){
			copy_pixel_values_nTuple_volume(dispField,dispField, i, min_int(j+1,((int)dispField->ySize)-1), k, i, j, k);
			nbModified++;
		}
		else if( (*correctInd) == 2){
			copy_pixel_values_nTuple_volume(dispField,dispField, i, j, min_int(k+1,((int)dispField->tSize)-1), i, j, k);
			nbModified++;
		}
		else
			MY_PRINTF("Error, correct ind not chosen\n.");
		//now calculate the error of the patch matching
		currentError = calclulate_patch_error(departVolume,arrivalVolume,dispField,occVol,i,j,k, -1,params);
		dispField->set_value(i,j,k,3,currentError);
	}
	else		//even iteration
	{
		if ( (*correctInd) == 0){
			copy_pixel_values_nTuple_volume(dispField,dispField, max_int(i-1,0), j, k, i, j, k);
			nbModified++;
		}

		else if( (*correctInd) == 1){
			copy_pixel_values_nTuple_volume(dispField,dispField, i, max_int(j-1,0), k, i, j, k);
			nbModified++;
		}

		else if( (*correctInd) == 2){
			copy_pixel_values_nTuple_volume(dispField,dispField, i, j, max_int(k-1,0), i, j, k);
			nbModified++;
		}
		else
			MY_PRINTF("Error, correct ind not chosen\n.");
		//now calculate the error of the patch matching
		currentError = calclulate_patch_error(departVolume,arrivalVolume,dispField,occVol,i,j,k, -1,params);
		dispField->set_value(i,j,k,3,currentError);
	}
	
	delete correctInd;
	delete minVector;
	
	return(nbModified);
}


//this function returns the minimum error of the patch differences around the pixel at (i,j,k)
//and returns the index of the best position in correctInd :
// -1 : none are correct
// 0 : left/right
// 1 : upper/lower
// 2 : before/after
float get_min_correct_error(nTupleVolume *dispField,nTupleVolume *departVol,nTupleVolume *arrivalVol, nTupleVolume *occVol,
							int x, int y, int t, int beforeOrAfter, int *correctInd, float *minVector, float minError,
                            const patchMatchParameterStruct *params)
{
	float minVal;
	int i;
    int dispX, dispY, dispT;

	minVal = minError;

	*correctInd = -1;	//initialise the correctInd vector to -1
    for (i=0;i<NDIMS;i++)
		minVector[i] = -1;

	if (beforeOrAfter == 0)	//we are looking left, upper, before : we are on an even iteration
	{
        dispX = (int)dispField->get_value((int)max_int(x-1,0),y,t,0); dispY = (int)dispField->get_value((int)max_int(x-1,0),y,t,1); dispT = (int)dispField->get_value((int)max_int(x-1,0),y,t,2); 
        if ( check_in_inner_boundaries(arrivalVol,x+dispX,y+dispY,t+dispT,params) && (!check_is_occluded(occVol, x+dispX, y+dispY, t+dispT)) &&
                (!check_already_used_patch( dispField, x, y, t, dispX, dispY, dispT)))
            minVector[0] = ssd_patch_measure(departVol, arrivalVol, dispField,occVol,x,y,t,x+dispX,y+dispY,t+dispT,minError,params);
        
        dispX = (int)dispField->get_value(x,(int)max_int(y-1,0),t,0); dispY = (int)dispField->get_value(x,(int)max_int(y-1,0),t,1); dispT = (int)dispField->get_value(x,(int)max_int(y-1,0),t,2); 
        if ( check_in_inner_boundaries(arrivalVol,x+dispX,y+dispY,t+dispT,params) && (!check_is_occluded(occVol, x+dispX, y+dispY, t+dispT)) && 
                (!check_already_used_patch( dispField, x, y, t, dispX, dispY, dispT)))
            minVector[1] = ssd_patch_measure(departVol, arrivalVol, dispField,occVol,x,y,t,x+dispX,y+dispY,t+dispT,minError,params);
        
        dispX = (int)dispField->get_value(x,y,(int)max_int(t-1,0),0); dispY = (int)dispField->get_value(x,y,(int)max_int(t-1,0),1); dispT = (int)dispField->get_value(x,y,(int)max_int(t-1,0),2); 
        if ( check_in_inner_boundaries(arrivalVol,x+dispX,y+dispY,t+dispT,params) && (!check_is_occluded(occVol, x+dispX, y+dispY, t+dispT)) &&
                (!check_already_used_patch( dispField, x, y, t, dispX, dispY, dispT)))
            minVector[2] = ssd_patch_measure(departVol, arrivalVol, dispField,occVol,x,y,t,x+dispX,y+dispY,t+dispT,minError,params);//
        

		for (i=0;i<NDIMS;i++)
		{
            if (minVector[i] == -1)
                continue;
            if ( minVector[i] < minVal)
			{
				switch (i){
					case 0 :
							minVal = minVector[i];
							*correctInd = 0;
						break;
					case 1 :
							minVal = minVector[i];
							*correctInd = 1;
						break;
					case 2 :
							minVal = minVector[i];
							*correctInd = 2;
						break;
					default :
						MY_PRINTF("Error in get_min_correct_error. The index i : %d is above 3.\n",i);
				}
			}

		}


	}
	else	//we are looking right, lower, after
	{
        
        dispX = (int)dispField->get_value((int)min_int(x+1,(departVol->xSize)-1),y,t,0); dispY = (int)dispField->get_value((int)min_int(x+1,(departVol->xSize)-1),y,t,1); dispT = (int)dispField->get_value((int)min_int(x+1,(departVol->xSize)-1),y,t,2); 
        if ( check_in_inner_boundaries(arrivalVol,x+dispX,y+dispY,t+dispT,params) && (!check_is_occluded(occVol, x+dispX, y+dispY, t+dispT)) && 
                (!check_already_used_patch( dispField, x, y, t, dispX, dispY, dispT)))
            minVector[0] = ssd_patch_measure(departVol, arrivalVol, dispField,occVol,x,y,t,x+dispX,y+dispY,t+dispT,minError,params);
        
        dispX = (int)dispField->get_value(x,(int)min_int(y+1,(departVol->ySize)-1),t,0); dispY = (int)dispField->get_value(x,(int)min_int(y+1,(departVol->ySize)-1),t,1); dispT = (int)dispField->get_value(x,(int)min_int(y+1,(departVol->ySize)-1),t,2); 
        if ( check_in_inner_boundaries(arrivalVol,x+dispX,y+dispY,t+dispT,params) && (!check_is_occluded(occVol, x+dispX, y+dispY, t+dispT)) && 
                (!check_already_used_patch( dispField, x, y, t, dispX, dispY, dispT)))
            minVector[1] = ssd_patch_measure(departVol, arrivalVol, dispField,occVol,x,y,t,x+dispX,y+dispY,t+dispT,minError,params);
        
        dispX = (int)dispField->get_value(x,y,(int)min_int(t+1,(departVol->tSize)-1),0); dispY = (int)dispField->get_value(x,y,(int)min_int(t+1,(departVol->tSize)-1),1); dispT = (int)dispField->get_value(x,y,(int)min_int(t+1,(departVol->tSize)-1),2); 
        if ( check_in_inner_boundaries(arrivalVol,x+dispX,y+dispY,t+dispT,params) && (!check_is_occluded(occVol, x+dispX, y+dispY, t+dispT)) &&
                (!check_already_used_patch( dispField, x, y, t, dispX, dispY, dispT)))
            minVector[2] = ssd_patch_measure(departVol, arrivalVol,dispField,occVol,x,y,t,x+dispX,y+dispY,t+dispT,minError,params);

		for (i=0;i<NDIMS;i++)
		{
            if (minVector[i] == -1)
                continue;
			if ( minVector[i] < minVal)
			{
				switch (i){
					case 0 :
							minVal = minVector[i];
							*correctInd = 0;
						break;
					case 1 :
							minVal = minVector[i];
							*correctInd = 1;
						break;
					case 2 :
							minVal = minVector[i];
							*correctInd = 2;
						break;
					default :
						MY_PRINTF("Error in get_min_correct_error. The index i : %d is above 3.\n",i);
				}
			}

		}
	}

	if ( (*correctInd) == -1)	//if none of the displacements are valid 
	{
		minVal = -1;
	}
	return(minVal);
}
