

//definitions for the morphological operations

#include "morpho.h"

nTupleVolume* create_structuring_element(const char * structType, int xSize, int ySize)
{
	
	if (strcmp(structType,"rectangle") == 0)
	{
		nTupleVolume *structElOut =
		new nTupleVolume(2,xSize,ySize,0,0,IMAGE_INDEXING);
	
		int xMin = -floor((xSize)/2);
		int yMin = -floor((ySize)/2);
		
		for (int x=0; x<xSize; x++)
			for (int y=0; y<ySize; y++)
			{
				int xTemp = xMin+x;
				int yTemp = yMin+y;
				structElOut->set_image_value(x,y,0,(imageDataType)xTemp);
				structElOut->set_image_value(x,y,1,(imageDataType)yTemp);
			}
		
		return(structElOut);
	}
	else
	{
		MY_PRINTF("Error in create_structuring_element. The structuring element type is not recognised.\n");
		return NULL;
	}

}

nTupleVolume* imerode(nTupleVolume *imgVol, nTupleVolume *structEl)
{
	//create output (eroded) image
	nTupleVolume *imgVolEroded =
	new nTupleVolume(1,imgVol->xSize,imgVol->ySize,imgVol->patchSizeX,imgVol->patchSizeY,imgVol->indexing);
	
	//erosion
	for (int x=0;x<imgVol->xSize;x++)
		for (int y=0;y<imgVol->ySize;y++)
		{
			imageDataType newValue = imgVol->get_image_value(x,y,0);
			for (int xStructEl=0;xStructEl<structEl->xSize;xStructEl++)
			{
				for (int yStructEl=0;yStructEl<structEl->ySize;yStructEl++)
				{
					int xShift = x+structEl->get_image_value(xStructEl,yStructEl,0);
					int yShift = y+structEl->get_image_value(xStructEl,yStructEl,1);
					if (xShift>=0 && xShift<(imgVol->xSize) && yShift>=0 && yShift<(imgVol->ySize))
					{
						imageDataType currValue = imgVol->get_image_value(xShift,yShift,0);
				  		if (currValue<newValue)
				  			newValue=currValue;	  
					}
				}
			}
			//set eroded value
			imgVolEroded->set_image_value(x,y,0,newValue);
		}
	return(imgVolEroded);
}

nTupleVolume* imdilate(nTupleVolume *imgVol, nTupleVolume *structEl)
{

	//create output (eroded) image
	nTupleVolume *imgVolDilated =
	new nTupleVolume(1,imgVol->xSize,imgVol->ySize,imgVol->patchSizeX,imgVol->patchSizeY,imgVol->indexing);
	
	//dilation
	for (int x=0;x<imgVol->xSize;x++)
		for (int y=0;y<imgVol->ySize;y++)
		{
			imageDataType newValue = imgVol->get_image_value(x,y,0);
			for (int xStructEl=0;xStructEl<structEl->xSize;xStructEl++)
			{
				for (int yStructEl=0;yStructEl<structEl->ySize;yStructEl++)
				{
					int xShift = x+structEl->get_image_value(xStructEl,yStructEl,0);
					int yShift = y+structEl->get_image_value(xStructEl,yStructEl,1);
					if (xShift>=0 && xShift<(imgVol->xSize) && yShift>=0 && yShift<(imgVol->ySize))
					{
						imageDataType currValue = imgVol->get_image_value(xShift,yShift,0);
				  		if (currValue>newValue)
				  			newValue=currValue;	  
					}
				}
			}
			//set eroded value
			imgVolDilated->set_image_value(x,y,0,newValue);
		}
	return(imgVolDilated);

}

