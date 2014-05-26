/*this function estimates the colours of pixels, given a shift map
showing th*/

#include "reconstruct_image.h"

/*this function calculates a nearest neighbour field, from imgVolA to imgVolB*/
void reconstruct_image(nTupleVolume* imgVol, nTupleVolume* imgVolModified, nTupleVolume* occVol,
        nTupleVolume* dispField, float sigmaColour, int reconstructionType, bool initialisation)
{
	int useAllPatches;
	if (initialisation==true)
		useAllPatches = 0;
	else
		useAllPatches = 1;
	
	/*decalarations*/
    int iMin,iMax,jMin,jMax,kMin,kMax;
    int weightInd;
    int xDisp, yDisp, tDisp,xDispShift,yDispShift,tDispShift;
    int hPatchSizeX,hPatchSizeY,hPatchSizeT;
    int nbNeighbours;
    int correctInfo;
    float alpha, adaptiveSigma;
    float *weights,sumWeights, avgColourR, avgColourG, avgColourB, *colours;

    hPatchSizeX = imgVol->hPatchSizeX;
    hPatchSizeY = imgVol->hPatchSizeY;
    hPatchSizeT = imgVol->hPatchSizeT;
    
    /*allocate the (maximum) memory for the weights*/
    nbNeighbours = (imgVol->patchSizeX)*(imgVol->patchSizeY)*(imgVol->patchSizeT);
    weights = (float*)malloc((size_t)(nbNeighbours*sizeof(float)));
    colours = (float*)malloc((size_t)(NCHANNELS*nbNeighbours*sizeof(float)));
    
	/*check certain parameters*/
	if( (imgVol->patchSizeX != (imgVol->patchSizeX)) || (imgVol->patchSizeY != (imgVol->patchSizeY)) ||
		(imgVol->patchSizeT != (imgVol->patchSizeT))  )	/*check that the patch sizes are equal*/
	{
		MY_PRINTF("Error in estimate_colour, the size of the patches are not equal in the two image volumes.");
		return;
	}
	if ( ( imgVol->patchSizeX > imgVol->xSize) || ( imgVol->patchSizeY > imgVol->ySize) || ( imgVol->patchSizeT > imgVol->tSize) )	/*check that the patch size is less or equal to each dimension in the images*/
	{
		MY_PRINTF("Error in estimate_colour, the patch size is to large for one or more of the dimensions of the image volume.");
		return;
	}

    for (int k=0; k<(occVol->tSize); k++)
        for (int j=0; j<(occVol->ySize); j++)
            for (int i=0; i<(occVol->xSize); i++)
            {    
                if ( ((occVol->get_value(i,j,k,0)) == 0) || ((occVol->get_value(i,j,k,0) == 2) )  )
                    continue;
                else    /*an occluded pixel (therefore to be modified)*/
                {
                    if (reconstructionType == 1 )
                    {
                        xDisp = i + (int)dispField->get_value(i,j,k,0);
                        yDisp = j + (int)dispField->get_value(i,j,k,1);
                        tDisp = k + (int)dispField->get_value(i,j,k,2);

                        ////if pure replacing of pixels
                        copy_pixel_values_nTuple_volume(imgVol, imgVolModified,xDisp, yDisp, tDisp, i, j, k);
                        ///set_value_nTuple_volume(occVol,i,j,k,2,0);
                        ///set_value_nTuple_volume(imgVol,i,j,k,0,0);
                        continue;
                    }
                     
                    
                    for (int ii=0;ii<(imgVol->patchSizeX)*(imgVol->patchSizeY)*(imgVol->patchSizeT); ii++)
					{
						weights[ii] = (float)-1;
						colours[ii] = (float)-1;
						colours[ii + nbNeighbours] = (float)-1;
						colours[ii + 2*nbNeighbours] = (float)-1;
					}

                    sumWeights = 0.0;
                    alpha = FLT_MAX;
                    correctInfo = 0;
                    avgColourR = 0.0;
                    avgColourG = 0.0;
                    avgColourB = 0.0;
                    
                    iMin = max_int(i - hPatchSizeX,0);
                    iMax = min_int(i + hPatchSizeX,(imgVol->xSize)-1 );
                    jMin = max_int(j - hPatchSizeY,0);
                    jMax = min_int(j + hPatchSizeY,(imgVol->ySize)-1 );
                    kMin = max_int(k - hPatchSizeT,0);
                    kMax = min_int(k + hPatchSizeT,(imgVol->tSize)-1 );
                    
                    /*
                    MY_PRINTF("iMin : %d, iMax : %d\n",iMin,iMax);
                    MY_PRINTF("jMin : %d, jMax : %d\n",jMin,jMax);
                    MY_PRINTF("kMin : %d, kMax : %d\n",kMin,kMax);*/
                    /*first calculate the weights*/
                    for (int kk=kMin; kk<=kMax;kk++)
                        for (int jj=jMin; jj<=jMax;jj++)
                            for (int ii=iMin; ii<=iMax;ii++)
                            {
                                /*get ssd similarity*/
                                xDisp = ii + (int)dispField->get_value(ii,jj,kk,0);
                                yDisp = jj + (int)dispField->get_value(ii,jj,kk,1);
                                tDisp = kk + (int)dispField->get_value(ii,jj,kk,2);
                                /*(spatio-temporally) shifted values of the covering patches*/
                                xDispShift = xDisp - (ii-i);
                                yDispShift = yDisp - (jj-j);
                                tDispShift = tDisp - (kk-k);
                        
                                 if (useAllPatches == 1)
                                 {
                                     
                                    alpha = (float)min_float(dispField->get_value(ii,jj,kk,3),alpha); 
                                    weightInd = (int)((kk-kMin)*(imgVol->patchSizeX)*(imgVol->patchSizeY) + (jj-jMin)*(imgVol->patchSizeX) + ii-iMin);
                                    weights[weightInd] = dispField->get_value(ii,jj,kk,3);
                                    
                                    
                                    colours[weightInd] = (float)(imgVol->get_value(xDispShift,yDispShift,tDispShift,0));
                                    colours[weightInd + nbNeighbours] = (float)(imgVol->get_value(xDispShift,yDispShift,tDispShift,1));
                                    colours[weightInd + 2*nbNeighbours] = (float)(imgVol->get_value(xDispShift,yDispShift,tDispShift,2));
                                    correctInfo = 1;
                                 }
                                 else   /*only use some of the patches*/
                                 {
                                     if (((occVol->get_value(ii,jj,kk,0)) == 0) || (occVol->get_value(ii,jj,kk,0) ==-1))
                                     {
                                        alpha = (float)min_float(dispField->get_value(ii,jj,kk,3),alpha); 
                                        weightInd = (int)((kk-kMin)*(imgVol->patchSizeX)*(imgVol->patchSizeY) + (jj-jMin)*(imgVol->patchSizeX) + ii-iMin);
                                        weights[weightInd] = dispField->get_value(ii,jj,kk,3);
                                        
                                        colours[weightInd] = (float)(imgVol->get_value(xDispShift,yDispShift,tDispShift,0));
                                        colours[weightInd + nbNeighbours] = (float)(imgVol->get_value(xDispShift,yDispShift,tDispShift,1));
                                        colours[weightInd + 2*nbNeighbours] = (float)(imgVol->get_value(xDispShift,yDispShift,tDispShift,2));
                                        correctInfo = 1;
                                     }
                                     else
                                     {
                                        weightInd = (int)((kk-kMin)*(imgVol->patchSizeX)*(imgVol->patchSizeY) + (jj-jMin)*(imgVol->patchSizeX) + ii-iMin);
                                        weights[weightInd] = -1;
                                         
                                        colours[weightInd] = -1;
                                        colours[weightInd + nbNeighbours] = -1;
                                        colours[weightInd + 2*nbNeighbours] = -1;
                                        
                                        continue;
                                     }
                                 }
                            }
                    
                    alpha = max_float(alpha,1);
                    if (correctInfo == 0)
                        continue;

                    if (reconstructionType == 3)
                    {
                        estimate_best_colour(imgVol,imgVolModified, weights, nbNeighbours, colours, i, j, k);
                        continue;
                    }
                    //get the 75th percentile of the distances for setting the adaptive sigma
                    adaptiveSigma = get_adaptive_sigma(weights,(imgVol->patchSizeX)*(imgVol->patchSizeY)*(imgVol->patchSizeT),sigmaColour);
					adaptiveSigma = max_float(adaptiveSigma,(float)0.1);
                    
                    /* ///MY_PRINTF("alpha : %f\n",alpha);
                    //adjust the weights : note, the indices which are outside the image boundaries
                    //will have no influence on the final weights (they are initialised to 0)  */
                    for (int kk=kMin; kk<=kMax;kk++)
                        for (int jj=jMin; jj<=jMax;jj++)
                            for (int ii=iMin; ii<=iMax;ii++)
                            {
                                if (useAllPatches)
                                {
                                    /*weights = exp( -weights/(2*sigma*alpha))*/
                                    weightInd = (int)((kk-kMin)*(imgVol->patchSizeX)*(imgVol->patchSizeY) + (jj-jMin)*(imgVol->patchSizeX) + ii-iMin);
                                    weights[weightInd] = (float)(exp( - ((weights[weightInd])/(2*adaptiveSigma*adaptiveSigma)) ));/*exp( - ((weights[ii])/(2*sigmaColour*sigmaColour*alpha)) );*/
                                    //
                                    sumWeights = (float)(sumWeights+weights[weightInd]);
                                }
                                else   /*only use some of the patches*/
                                {
                                     if (((occVol->get_value(ii,jj,kk,0)) == 0) || (occVol->get_value(ii,jj,kk,0) ==-1))
                                     {
                                        /*weights = exp( -weights/(2*sigma*alpha))*/
                                        weightInd = (int)((kk-kMin)*(imgVol->patchSizeX)*(imgVol->patchSizeY) + (jj-jMin)*(imgVol->patchSizeX) + ii-iMin);
                                        weights[weightInd] = (float)(exp( - ((weights[weightInd])/(2*adaptiveSigma*adaptiveSigma)) ));/*exp( - ((weights[ii])/(2*sigmaColour*sigmaColour*alpha)) );*/
                                        //
                                        sumWeights = (float)(sumWeights+weights[weightInd]);
                                     }
                                     else
                                         continue;
                                }
                            }

                    /*now calculate the pixel value(s)*/
                    for (int kk=kMin; kk<=kMax;kk++)
                        for (int jj=jMin; jj<=jMax;jj++)
                            for (int ii=iMin; ii<=iMax;ii++)
                            {
                                if (useAllPatches)
                                {
                                    weightInd = (int)((kk-kMin)*(imgVol->patchSizeX)*(imgVol->patchSizeY) + (jj-jMin)*(imgVol->patchSizeX) + ii-iMin);
                                    /*get ssd similarity*/
                                    xDisp = ii + (int)dispField->get_value(ii,jj,kk,0);
                                    yDisp = jj + (int)dispField->get_value(ii,jj,kk,1);
                                    tDisp = kk + (int)dispField->get_value(ii,jj,kk,2);
                                    /*(spatio-temporally) shifted values of the covering patches*/
                                    xDispShift = xDisp - (ii-i);
                                    yDispShift = yDisp - (jj-j);
                                    tDispShift = tDisp - (kk-k);
                                    avgColourR = avgColourR + (float)(weights[weightInd])*(imgVol->get_value(xDispShift,yDispShift,tDispShift,0));
                                    avgColourG = avgColourG + (float)(weights[weightInd])*(imgVol->get_value(xDispShift,yDispShift,tDispShift,1));
                                    avgColourB = avgColourB + (float)(weights[weightInd])*(imgVol->get_value(xDispShift,yDispShift,tDispShift,2));
                                }
                                else
                                {
                                     if (((occVol->get_value(ii,jj,kk,0)) == 0) || (occVol->get_value(ii,jj,kk,0) ==-1))
                                     {
                                        weightInd = (int)((kk-kMin)*(imgVol->patchSizeX)*(imgVol->patchSizeY) + (jj-jMin)*(imgVol->patchSizeX) + ii-iMin);
                                        /*get ssd similarity*/
                                        xDisp = ii + (int)dispField->get_value(ii,jj,kk,0);
                                        yDisp = jj + (int)dispField->get_value(ii,jj,kk,1);
                                        tDisp = kk + (int)dispField->get_value(ii,jj,kk,2);
                                        /*(spatio-temporally) shifted values of the covering patches*/
                                        xDispShift = xDisp - (ii-i);
                                        yDispShift = yDisp - (jj-j);
                                        tDispShift = tDisp - (kk-k);
                                        avgColourR = avgColourR + (float)(weights[weightInd])*(imgVol->get_value(xDispShift,yDispShift,tDispShift,0));
                                        avgColourG = avgColourG + (float)(weights[weightInd])*(imgVol->get_value(xDispShift,yDispShift,tDispShift,1));
                                        avgColourB = avgColourB + (float)(weights[weightInd])*(imgVol->get_value(xDispShift,yDispShift,tDispShift,2));
                                     }
                                     else
                                         continue;
                                }
                            }
                         /*MY_PRINTF("SumWeights : %f\n",sumWeights);*/
                    imgVolModified->set_value(i,j,k,0,(imageDataType)(avgColourR/(sumWeights)));
                    imgVolModified->set_value(i,j,k,1,(imageDataType)(avgColourG/(sumWeights)));
                    imgVolModified->set_value(i,j,k,2,(imageDataType)(avgColourB/(sumWeights)));
                    /*set_value_nTuple_volume(occVol,i,j,k,0,0);*/
                }
            }

        free(weights);
        free(colours);
        return;
}
