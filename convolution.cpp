

#include "convolution.h"


nTupleVolume * create_convolution_kernel(const char * kernelType, int xSizeKernel, int ySizeKernel, int tSizeKernel, float stdDev)
{
	nTupleVolume *convKernel = new nTupleVolume(1,xSizeKernel,ySizeKernel,tSizeKernel,0,0,0,IMAGE_INDEXING);

	
	if (strcmp(kernelType,"gaussian") == 0)	//gaussian filter
	{
		for (int x=(int)-(floor(convKernel->xSize)/2); x<=(int)(floor(convKernel->xSize)/2); x++)
			for (int y=(int)-(floor(convKernel->ySize)/2); y<=(int)(floor(convKernel->ySize)/2); y++)
				for (int t=(int)-(floor(convKernel->tSize)/2); t<=(int)(floor(convKernel->tSize)/2); t++)
				{
					float distance = (x)*(x) + (y)*(y) + (t)*(t);
					convKernel->set_value(x+(int)(floor(convKernel->xSize)/2),
					y+(int)(floor(convKernel->ySize)/2),t+(int)(floor(convKernel->tSize)/2),0, exp( -(distance)/(2*stdDev*stdDev) ) );
				}
	}
	else if(strcmp(kernelType,"average") == 0)
	{
		//set convolution kernel
		for (int x=0; x<(convKernel->xSize); x++)
			for (int y=0; y<(convKernel->ySize); y++)
				for (int t=0; t<(convKernel->tSize); t++)
					convKernel->set_value(x,y,t,0,(imageDataType)1 );
	}
	else
	{
		MY_PRINTF("Here create_convolution_kernel. The convolution kernel type is unrecognised.\n");
	}
	
	return(convKernel);
}

nTupleVolume * normalised_convolution_masked(nTupleVolume *imgVol, nTupleVolume *convKernel, nTupleVolume *occlusionMask)
{
	nTupleVolume * imgVolConvolved = new nTupleVolume(imgVol->nTupleSize,imgVol->xSize, imgVol->ySize,imgVol->tSize,
	imgVol->patchSizeX, imgVol->patchSizeY, imgVol->patchSizeT,imgVol->indexing);

	for (int p=0; p<imgVol->nTupleSize; p++)
		for (int x=0; x<imgVol->xSize; x++)
			for (int y=0; y<imgVol->ySize; y++)
				for (int t=0; t<imgVol->tSize; t++)
				{
					imageDataType sumConv = 0;
					imageDataType convTemp = 0;
					for (int xKernel=0; xKernel<(convKernel->xSize); xKernel++)
						for (int yKernel=0; yKernel<(convKernel->ySize); yKernel++)
							for (int tKernel=0; tKernel<(convKernel->tSize); tKernel++)
							{
								int xShift = x - (int)floor((convKernel->xSize)/2) + xKernel;
								int yShift = y - (int)floor((convKernel->ySize)/2) + yKernel;
								int tShift = t - (int)floor((convKernel->tSize)/2) + tKernel;
						
								if ( xShift>=0 && xShift<(imgVol->xSize) && yShift>=0 && yShift<(imgVol->ySize) 
								&& tShift>=0 && tShift<(imgVol->tSize))
								{
									if ( occlusionMask != NULL)
									{
										if (occlusionMask->get_value(xShift,yShift,tShift,0) ==0)
										{
											convTemp = convTemp + ( imgVol->get_value(xShift,yShift,tShift,p) )*
											(convKernel->get_value(xKernel,yKernel,tKernel,0) );
											sumConv = sumConv + convKernel->get_value(xKernel,yKernel,tKernel,0);
										}
									}
									else
									{
										convTemp = convTemp + ( imgVol->get_value(xShift,yShift,tShift,p) )*
											(convKernel->get_value(xKernel,yKernel,tKernel,0) );
										sumConv = sumConv + convKernel->get_value(xKernel,yKernel,tKernel,0);
									}
								}
							}
					if (sumConv == 0)	//if no unmasked pixels are available here
						imgVolConvolved->set_value(x,y,t,p,(imageDataType)0);
					else
						imgVolConvolved->set_value(x,y,t,p,(imageDataType)convTemp/sumConv);
				}

	return imgVolConvolved;
}

nTupleVolume * normalised_convolution_masked_separable(nTupleVolume *imgVol, nTupleVolume *convKernelX, nTupleVolume *convKernelY, nTupleVolume *occlusionMask)
{
	nTupleVolume * imgVolConvolved = new nTupleVolume(imgVol->nTupleSize,imgVol->xSize, imgVol->ySize,imgVol->tSize,
	imgVol->patchSizeX, imgVol->patchSizeY, imgVol->patchSizeT,imgVol->indexing);

	for (int p=0; p<imgVol->nTupleSize; p++)
		for (int x=0; x<imgVol->xSize; x++)
			for (int y=0; y<imgVol->ySize; y++)
				for (int t=0; t<imgVol->tSize; t++)
				{
					imageDataType sumConv = 0;
					imageDataType convTemp = 0;
					for (int xKernel=0; xKernel<(convKernelX->xSize); xKernel++)
						for (int yKernel=0; yKernel<(convKernelX->ySize); yKernel++)
							for (int tKernel=0; tKernel<(convKernelX->tSize); tKernel++)
							{
								int xShift = x - (int)floor((convKernelX->xSize)/2) + xKernel;
								int yShift = y - (int)floor((convKernelX->ySize)/2) + yKernel;
								int tShift = t - (int)floor((convKernelX->tSize)/2) + tKernel;
						
								if ( xShift>=0 && xShift<(imgVol->xSize) && yShift>=0 && yShift<(imgVol->ySize) 
								&& tShift>=0 && tShift<(imgVol->tSize))
								{
									if ( occlusionMask != NULL)
									{
										if (occlusionMask->get_value(xShift,yShift,tShift,0) ==0)
										{
											convTemp = convTemp + ( imgVol->get_value(xShift,yShift,tShift,p) )*
											(convKernelX->get_value(xKernel,yKernel,tKernel,0) );
											sumConv = sumConv + convKernelX->get_value(xKernel,yKernel,tKernel,0);
										}
									}
									else
									{
										convTemp = convTemp + ( imgVol->get_value(xShift,yShift,tShift,p) )*
											(convKernelX->get_value(xKernel,yKernel,tKernel,0) );
										sumConv = sumConv + convKernelX->get_value(xKernel,yKernel,tKernel,0);
									}
								}
							}
					if (sumConv == 0)	//if no unmasked pixels are available here
						imgVolConvolved->set_value(x,y,t,p,(imageDataType)0);
					else
						imgVolConvolved->set_value(x,y,t,p,(imageDataType)convTemp/sumConv);
				}

	//now convolve with the other filter
	for (int p=0; p<imgVolConvolved->nTupleSize; p++)
		for (int x=0; x<imgVolConvolved->xSize; x++)
			for (int y=0; y<imgVolConvolved->ySize; y++)
				for (int t=0; t<imgVolConvolved->tSize; t++)
				{
					imageDataType sumConv = 0;
					imageDataType convTemp = 0;
					for (int xKernel=0; xKernel<(convKernelY->xSize); xKernel++)
						for (int yKernel=0; yKernel<(convKernelY->ySize); yKernel++)
							for (int tKernel=0; tKernel<(convKernelY->tSize); tKernel++)
							{
								int xShift = x - (int)floor((convKernelY->xSize)/2) + xKernel;
								int yShift = y - (int)floor((convKernelY->ySize)/2) + yKernel;
								int tShift = t - (int)floor((convKernelY->tSize)/2) + tKernel;
						
								if ( xShift>=0 && xShift<(imgVol->xSize) && yShift>=0 && yShift<(imgVol->ySize) 
								&& tShift>=0 && tShift<(imgVol->tSize))
								{
									if ( occlusionMask != NULL)
									{
										if (occlusionMask->get_value(xShift,yShift,tShift,0) ==0)
										{
											convTemp = convTemp + ( imgVolConvolved->get_value(xShift,yShift,tShift,p) )*
											(convKernelY->get_value(xKernel,yKernel,tKernel,0) );
											sumConv = sumConv + convKernelY->get_value(xKernel,yKernel,tKernel,0);
										}
									}
									else
									{
										convTemp = convTemp + ( imgVolConvolved->get_value(xShift,yShift,tShift,p) )*
											(convKernelY->get_value(xKernel,yKernel,tKernel,0) );
										sumConv = sumConv + convKernelY->get_value(xKernel,yKernel,tKernel,0);
									}
								}
							}
					if (sumConv == 0)	//if no unmasked pixels are available here
						imgVolConvolved->set_value(x,y,t,p,(imageDataType)0);
					else
						imgVolConvolved->set_value(x,y,t,p,(imageDataType)convTemp/sumConv);
				}
	return imgVolConvolved;
}
