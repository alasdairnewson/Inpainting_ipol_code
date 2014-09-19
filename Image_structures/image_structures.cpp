//this functions holds the definitions for the image_structures used in the spatio-temporal patchMatch

#include "image_structures.h"

float min_float(float a, float b)
{
    if (a<b)
        return a;
    else
        return b;
}

float max_float(float a, float b)
{
    if (a>b)
        return a;
    else
        return b;
}

int min_int(int a, int b)
{
    if (a<b)
        return a;
    else
        return b;
}

int max_int(int a, int b)
{
    if (a>b)
        return a;
    else
        return b;
}

float rand_float_range(float a, float b)
{
	if (a == b)
		return a;
	else
		return ((b-a)*((float)rand()/RAND_MAX))+a;
}

int rand_int_range(int a, int b)
{
	if (a == b)
		return a;
	else
		return ( rand()%(b-a+1) + a);
}

float round_float(float a)
{
	float aFloor,aCeil;
	aFloor = (float)floor((float)a);
	aCeil = (float)ceil((float)a);
	if (a<0)	//less than 0
	{
		if ( (a-aCeil) < -0.5)
			return(aFloor);	//round up
		else
			return(aCeil);
	}
	else	//greater or equal to 0
	{
		if ( (a-aFloor) < 0.5)	//round down
			return(aFloor);
		else
			return(aCeil);
	}
}

long getMilliSecs()
{ 
    timeval t;    
    gettimeofday(&t, NULL);
    return t.tv_sec*1000 + t.tv_usec/1000;
}

char* int_to_string(int value)
{
	std::ostringstream os ;
	os << value ;
	std::string stringTemp = os.str();
	char * stringOut = new char[stringTemp.length() + 1];
	std::strcpy(stringOut,stringTemp.c_str());
	
	stringOut[stringTemp.length()] = '\0';
	return (stringOut);
}

int sub_to_ind(nTupleVolume* imgVol, int x, int y, int t)
{
	if (imgVol->indexing == ROW_FIRST)
    	return( t*(imgVol->xSize)*(imgVol->ySize) + y*(imgVol->xSize) + x);
	else if(imgVol->indexing == COLUMN_FIRST)
    	return( t*(imgVol->xSize)*(imgVol->ySize) + x*(imgVol->ySize) + y);
	else
	{
		MY_PRINTF("Here sub_to_ind (image_structures.cpp). Error: unknown image indexing type\n.");
		return(-1);
	}
}

void ind_to_sub(nTupleVolume* imgVol, int linearIndex, int *x, int *y, int *t)
{
	if (imgVol->indexing == ROW_FIRST)
	{
		*t = (int)floor((float) linearIndex/((imgVol->xSize)*(imgVol->ySize)));
		
		*y = (int) floor((float) (linearIndex - (*t)*(imgVol->xSize)*(imgVol->ySize)) /
		        ( imgVol->xSize ) );
		*x = (int) (linearIndex - (*t)*(imgVol->xSize)*(imgVol->ySize) - (*y)*(imgVol->xSize));
	}
	else if (imgVol->indexing == COLUMN_FIRST)
	{
		*t = (int)floor((float) linearIndex/((imgVol->xSize)*(imgVol->ySize)));
		
		*x = (int) floor((float) (linearIndex - (*t)*(imgVol->xSize)*(imgVol->ySize)) /
		        ( imgVol->ySize ) );
		*y = (int) (linearIndex - (*t)*(imgVol->xSize)*(imgVol->ySize) - (*x)*(imgVol->ySize));
	}
	else
	{
		MY_PRINTF("Here ind_to_sub (image_structures.cpp). Error: unknown image indexing type\n.");
	}
}

void patch_index_to_sub(nTupleVolume *imgVol, int patchIndex, int *colourInd,int *xInd, int *yInd, int *tInd)
{
    ASSERT(patchIndex>=0 && (patchIndex<((imgVol->patchSizeT)*(imgVol->patchSizeX)*(imgVol->patchSizeY)*(imgVol->nTupleSize)) ));
        *tInd = (int)floor((float) patchIndex/((imgVol->patchSizeX)*(imgVol->patchSizeY)*(imgVol->nTupleSize)));
	    
		*xInd = (int) floor((float) (patchIndex - (*tInd)*(imgVol->patchSizeX)*(imgVol->patchSizeY)*(imgVol->nTupleSize)) /
				(( imgVol->patchSizeY )*(imgVol->nTupleSize)) );
		*yInd = (int) floor((float)(patchIndex - (*tInd)*(imgVol->patchSizeX)*(imgVol->patchSizeY)*(imgVol->nTupleSize)
			- (*xInd)*(imgVol->patchSizeY)*(imgVol->nTupleSize))/
                (imgVol->nTupleSize));
        *colourInd = (int) (patchIndex - (*tInd)*(imgVol->patchSizeX)*(imgVol->patchSizeY)*(imgVol->nTupleSize) 
                            - (*xInd)*(imgVol->patchSizeY)*(imgVol->nTupleSize) 
							- (*yInd)*(imgVol->nTupleSize) );
}

void show_patch_match_parameters(patchMatchParameterStruct *patchMatchParams)
{
	MY_PRINTF("Patch size x : %d\n", patchMatchParams->patchSizeX);
	MY_PRINTF("Patch size y : %d\n", patchMatchParams->patchSizeY);
	MY_PRINTF("Patch size t : %d\n", patchMatchParams->patchSizeT);
	MY_PRINTF("nIters : %d\n", patchMatchParams->nIters);
	MY_PRINTF("w : %d\n", patchMatchParams->w);
	MY_PRINTF("alpha : %f\n", patchMatchParams->alpha);
	MY_PRINTF("partialComparison : %d\n", patchMatchParams->partialComparison);
    MY_PRINTF("fullSearch : %d\n", patchMatchParams->fullSearch);
    
    MY_PRINTF("\n");
}

int check_in_boundaries( nTupleVolume *imgVol, int x, int y, int t)
{
    if ( (x>=0) && (y>=0) && (t>=0) &&
    (x < ( (imgVol->xSize) )) && (y < ( (imgVol->ySize) )) && (t < ( (imgVol->tSize) )))
    {
            return 1;
    }
    else
        return 0;
}

//check if the pixel is in the inner boundary : that is, none of the pixels in the patch centred on this pixel are outside of the image boundary 
int check_in_inner_boundaries( nTupleVolume *imgVol, int x, int y, int t, const patchMatchParameterStruct *params)
{

	int hPatchSizeX = imgVol->hPatchSizeX;
	int hPatchSizeY = imgVol->hPatchSizeY;
	int hPatchSizeT = imgVol->hPatchSizeT;
    
    if ( (x>=hPatchSizeX) && (y>=hPatchSizeY) && (t>=hPatchSizeT) &&
	(x < ( (imgVol->xSize)-hPatchSizeX )) && (y < ( (imgVol->ySize)-hPatchSizeY )) && (t < ( (imgVol->tSize)-hPatchSizeT )))
    {
        return 1;
    }
    else
        return 0;
}

nTupleVolume::nTupleVolume()    //create empty volume
{
    xSize = 0;
    ySize = 0;
    tSize = 0;
	nDims = 4;
    patchSizeX = 0;
	patchSizeY = 0;
	patchSizeT = 0;
	hPatchSizeX = 0;
	hPatchSizeY = 0;
	hPatchSizeT = 0;
    nT = 0;
    nY = 0;
    nX = 0;
    nC = 0;
	nElsTotal = 0;
	values = NULL;
    indexing = -1;
    destroyValues = 0;
}

nTupleVolume::nTupleVolume(nTupleVolume *imgVolIn)
{
    //copy information
	nTupleSize = imgVolIn->nTupleSize;

	xSize = imgVolIn->xSize;
	ySize = imgVolIn->ySize;
	tSize = imgVolIn->tSize;
	patchSizeX = imgVolIn->patchSizeX;
	patchSizeY = imgVolIn->patchSizeY;
	patchSizeT = imgVolIn->patchSizeT;
	hPatchSizeX = imgVolIn->hPatchSizeX;
	hPatchSizeY = imgVolIn->hPatchSizeY;
	hPatchSizeT = imgVolIn->hPatchSizeT;

	indexing = imgVolIn->indexing;
	nT = imgVolIn->nT;
	nY = imgVolIn->nY;
	nX = imgVolIn->nX;
	nC = imgVolIn->nC;

	nElsTotal = imgVolIn->nElsTotal;
	//copy the image info
	values = (imageDataType*)malloc( (size_t) nElsTotal*nTupleSize*sizeof(imageDataType));
    memcpy(values,imgVolIn->get_value_ptr(0, 0, 0, 0),nElsTotal*nTupleSize*sizeof(imageDataType));

    destroyValues = 1;
}

nTupleVolume::nTupleVolume(int nTupleSizeIn, int xSizeIn, int ySizeIn, int tSizeIn, int indexingIn)
{
    //declarations
	int i;
	nTupleSize = nTupleSizeIn;

	xSize = xSizeIn;
	ySize = ySizeIn;
	tSize = tSizeIn;
	nDims = 4;
	patchSizeX = 0;
	patchSizeY = 0;
	patchSizeT = 0;
	hPatchSizeX = 0;
	hPatchSizeY = 0;
	hPatchSizeT = 0;

    if (indexingIn == 0)    //row first
    {
	nC = (xSize)*(ySize)*(tSize);
        nT = (xSize)*(ySize);
        nY = (xSize);
        nX = 1;
    }
    else if (indexingIn == 1)   //column first
    {
        nC = (ySize)*(xSize)*(tSize);
        nT = (ySize)*(xSize);
        nX = (ySize);
        nY = 1;
    }
    else
    {
        MY_PRINTF("Unknown indexing : %d\n", indexingIn);
    }

	nElsTotal = (xSize)*(ySize)*(tSize);
	//get dimensions and total number of elements
	values = (imageDataType*)malloc( (size_t) nElsTotal*nTupleSize*sizeof(imageDataType));

	for (i=0; i<(nElsTotal*nTupleSize); i++)
		values[i] = 0;
    
    indexing = indexingIn;  //row first
    destroyValues = 1;
}

//special case for images
nTupleVolume::nTupleVolume(int nTupleSizeIn, int xSizeIn, int ySizeIn, int indexingIn)
{
    //declarations
	int i;
	nTupleSize = nTupleSizeIn;

	xSize = xSizeIn;
	ySize = ySizeIn;
	tSize = 1;
	nDims = 4;
	patchSizeX = 0;
	patchSizeY = 0;
	patchSizeT = 0;
	hPatchSizeX = 0;
	hPatchSizeY = 0;
	hPatchSizeT = 0;

    if (indexingIn == 0)    //row first
    {
		nC = (xSize)*(ySize)*(tSize);
        nT = (xSize)*(ySize);
        nY = (xSize);
        nX = 1;
    }
    else if (indexingIn == 1)   //column first
    {
        nC = (ySize)*(xSize)*(tSize);
        nT = (ySize)*(xSize);
        nX = (ySize);
        nY = 1;
    }
    else
    {
        MY_PRINTF("Unknown indexing : %d\n", indexingIn);
    }

	nElsTotal = (xSize)*(ySize)*(tSize);
	//get dimensions and total number of elements
	values = (imageDataType*)malloc( (size_t) nElsTotal*nTupleSize*sizeof(imageDataType));

	for (i=0; i<(nElsTotal*nTupleSize); i++)
		values[i] = 0;
    
    indexing = indexingIn;  //row first
    destroyValues = 1;
}

nTupleVolume::nTupleVolume(int nTupleSizeIn, int xSizeIn, int ySizeIn, int tSizeIn, 
            int patchSizeXIn, int patchSizeYIn, int patchSizeTIn, int indexingIn)
{
	//declarations
	int i;
	nTupleSize = nTupleSizeIn;

	xSize = xSizeIn;
	ySize = ySizeIn;
	tSize = tSizeIn;
	nDims = 4;
	patchSizeX = patchSizeXIn;
	patchSizeY = patchSizeYIn;
	patchSizeT = patchSizeTIn;
	hPatchSizeX = (int)floor((float)patchSizeX/2);
	hPatchSizeY = (int)floor((float)patchSizeY/2);
	hPatchSizeT = (int)floor((float)patchSizeT/2);


    if (indexingIn == 0)    //row first
    {
	nC = (xSize)*(ySize)*(tSize);
        nT = (xSize)*(ySize);
        nY = (xSize);
        nX = 1;
    }
    else if (indexingIn == 1)   //column first
    {
        nC = (ySize)*(xSize)*(tSize);
        nT = (ySize)*(xSize);
        nX = (ySize);
        nY = 1;
    }
    else
    {
        MY_PRINTF("Unknown indexing : %d\n", indexingIn);
    }

	nElsTotal = (xSize)*(ySize)*(tSize);
	//get dimensions and total number of elements
	values = (imageDataType*)malloc( (size_t) nElsTotal*nTupleSize*sizeof(imageDataType));

	for (i=0; i<(nElsTotal*nTupleSize); i++)
		values[i] = 0;
    
    indexing = indexingIn;  //row first
    destroyValues = 1;
}

//create image volume with an already existing array for the values
nTupleVolume::nTupleVolume(int nTupleSizeIn, int xSizeIn, int ySizeIn, int tSizeIn, 
            int patchSizeXIn, int patchSizeYIn, int patchSizeTIn, int indexingIn, imageDataType* valuesIn)
{
	//declarations
	nTupleSize = nTupleSizeIn;

	xSize = xSizeIn;
	ySize = ySizeIn;
	tSize = tSizeIn;
	nDims = 4;
	patchSizeX = patchSizeXIn;
	patchSizeY = patchSizeYIn;
	patchSizeT = patchSizeTIn;
	hPatchSizeX = (int)floor((float)patchSizeX/2);
	hPatchSizeY = (int)floor((float)patchSizeY/2);
	hPatchSizeT = (int)floor((float)patchSizeT/2);


    if (indexingIn == 0)    //row first
    {
	nC = (xSize)*(ySize)*(tSize);
        nT = (xSize)*(ySize);
        nY = (xSize);
        nX = 1;
    }
    else if (indexingIn == 1)   //column first
    {
        nC = (ySize)*(xSize)*(tSize);
        nT = (ySize)*(xSize);
        nX = (ySize);
        nY = 1;
    }
    else
    {
        MY_PRINTF("Unknown indexing : %d\n", indexingIn);
    }

	nElsTotal = (xSize)*(ySize)*(tSize);
	//get dimensions and total number of elements
	values = valuesIn;
    
    indexing = indexingIn;
    destroyValues = 0;
}

//create image volume whose temporal dimension is equal to 1 (an image)
nTupleVolume::nTupleVolume(int nTupleSizeIn, int xSizeIn, int ySizeIn, 
            int patchSizeXIn, int patchSizeYIn, int indexingIn)
{
	//declarations
	nTupleSize = nTupleSizeIn;

	xSize = xSizeIn;
	ySize = ySizeIn;
	tSize = 1;
	nDims = 4;
	patchSizeX = patchSizeXIn;
	patchSizeY = patchSizeYIn;
	patchSizeT = 1;
	hPatchSizeX = (int)floor((float)patchSizeX/2);
	hPatchSizeY = (int)floor((float)patchSizeY/2);
	hPatchSizeT = (int)0;


    if (indexingIn == 0)    //row first
    {
		nC = (xSize)*(ySize)*(tSize);
        nT = (xSize)*(ySize);
        nY = (xSize);
        nX = 1;
    }
    else if (indexingIn == 1)   //column first
    {
        nC = (ySize)*(xSize)*(tSize);
        nT = (ySize)*(xSize);
        nX = (ySize);
        nY = 1;
    }
    else
    {
        MY_PRINTF("Unknown indexing : %d\n", indexingIn);
    }

	nElsTotal = (xSize)*(ySize)*(tSize);
	//get dimensions and total number of elements
	values = (imageDataType*)malloc( (size_t) nElsTotal*nTupleSize*sizeof(imageDataType));
    
    indexing = indexingIn;
    destroyValues = 1;
}


//create image volume with an already existing array for the values, and whose temporal dimension is equal to 1 (an image)
nTupleVolume::nTupleVolume(int nTupleSizeIn, int xSizeIn, int ySizeIn, 
            int patchSizeXIn, int patchSizeYIn, int indexingIn, imageDataType* valuesIn)
{
	//declarations
	nTupleSize = nTupleSizeIn;

	xSize = xSizeIn;
	ySize = ySizeIn;
	tSize = 1;
	nDims = 4;
	patchSizeX = patchSizeXIn;
	patchSizeY = patchSizeYIn;
	patchSizeT = 1;
	hPatchSizeX = (int)floor((float)patchSizeX/2);
	hPatchSizeY = (int)floor((float)patchSizeY/2);
	hPatchSizeT = (int)0;


    if (indexingIn == 0)    //row first
    {
	nC = (xSize)*(ySize)*(tSize);
        nT = (xSize)*(ySize);
        nY = (xSize);
        nX = 1;
    }
    else if (indexingIn == 1)   //column first
    {
        nC = (ySize)*(xSize)*(tSize);
        nT = (ySize)*(xSize);
        nX = (ySize);
        nY = 1;
    }
    else
    {
        MY_PRINTF("Unknown indexing : %d\n", indexingIn);
    }

	nElsTotal = (xSize)*(ySize)*(tSize);
	//get dimensions and total number of elements
	values = valuesIn;
    
    indexing = indexingIn;
    destroyValues = 0;
}

nTupleVolume::~nTupleVolume()
{
	if((xSize)> 0)
    {
        if (destroyValues  == 1)
            free(values);
    }
}

imageDataType nTupleVolume::get_value(int x, int y, int t, int z)
{
	//check parameters
	if( (x<0) || (y<0) || (t<0) || (z<0) || (x>=xSize) || (y>=ySize) || (t>=tSize) || (z>=nTupleSize))
	{
		MY_PRINTF("Error, in get_value_nTuple_volume. At least one of the indices is incorrect.\n");
		MY_PRINTF(" x = %d \n y = %d\n t = %d\n z = %d\n",x,y,t,z);
		return(-1);
	}
    return( values[ (t*(nT)) + (y*(nY)) + (x*(nX)) + z*(nC)] );
}

//special get function for an image
imageDataType nTupleVolume::get_image_value(int x, int y, int z)
{
	//check parameters
	if( (x<0) || (y<0) || (z<0) || (x>=xSize) || (y>=ySize) || (z>=nTupleSize))
	{
		MY_PRINTF("Error, in get_image_value. At least one of the indices is incorrect.\n");
		MY_PRINTF(" x = %d \n y = %d\n z = %d\n",x,y,z);
		return(-1);
	}
    return( values[ (y*(nY)) + (x*(nX)) + z*(nC)] );
}

imageDataType* nTupleVolume::get_value_ptr(int x, int y, int t, int z)
{
	//check parameters
	if( (x<0) || (y<0) || (t<0) || (x>=xSize) || (y>=ySize) || (t>=tSize))
	{
		MY_PRINTF("Error, in get_value_nTuple_volume. At least one of the indices is incorrect.\n");
		MY_PRINTF(" x = %d \n y = %d\n t = %d\n",x,y,t);
		return(NULL);
	}
	return( (values) + (  (t*(nT)) + (y*(nY)) + (x*(nX)) + z*(nC) ) );
}

imageDataType* nTupleVolume::get_data_ptr()
{
	return(values);
}

void nTupleVolume::set_value(int x, int y, int t, int z, imageDataType value)
{
	//check parameters
	if( (x<0) || (y<0) || (t<0) || (z<0) || (x>=xSize) || (y>=ySize) || (t>=tSize) || (z>=nTupleSize))
	{
		MY_PRINTF("Error, at least one of the indices is incorrect.\n");
		MY_PRINTF("x = %d \n y = %d\n t = %d\n z = %d\n",x,y,t,z);
	}
	values[ (t*(nT)) + (y*(nY)) + (x*(nX)) + z*(nC)] = value;
}

//special set function for an image
void nTupleVolume::set_image_value(int x, int y, int z, imageDataType value)
{
	//check parameters
	if( (x<0) || (y<0) || (z<0) || (x>=xSize) || (y>=ySize) || (z>=nTupleSize))
	{
		MY_PRINTF("Error, at least one of the indices is incorrect.\n");
		MY_PRINTF("x = %d \n y = %d\n z = %d\n",x,y,z);
	}
	values[ (y*(nY)) + (x*(nX)) + z*(nC)] = value;
}

void nTupleVolume::set_all_image_values(imageDataType value)
{
	for (int x=0; x<xSize; x++)
		for(int y=0; y<ySize; y++)
			for(int t=0; t<tSize; t++)
				for (int p=0; p<nTupleSize; p++)
					set_value(x,y,t,p,value);
}

void nTupleVolume::add(imageDataType addScalar)
{
	for (int x=0; x<xSize; x++)
		for(int y=0; y<ySize; y++)
			for(int t=0; t<tSize; t++)
				for (int p=0; p<nTupleSize; p++)
					set_value(x,y,t,p,(imageDataType) (addScalar+get_value(x,y,t,p)));
}

imageDataType nTupleVolume::sum_nTupleVolume()
{
	int sumVolume = 0;
	for (int x=0; x<(xSize); x++)
		for(int y=0; y<(ySize); y++)
			for(int t=0; t<tSize; t++)
				for (int p=0; p<nTupleSize; p++)//nTupleSize
					sumVolume = sumVolume + (int)round_float(get_value(x,y,t,p));
	return((imageDataType)sumVolume);
}

void nTupleVolume::multiply(imageDataType multiplyFactor)
{
	for (int x=0; x<xSize; x++)
		for(int y=0; y<ySize; y++)
			for(int t=0; t<tSize; t++)
				for (int p=0; p<nTupleSize; p++)
					set_value(x,y,t,p,(imageDataType) multiplyFactor*get_value(x,y,t,p));
}

imageDataType nTupleVolume::max_value()
{
	imageDataType maxVal = 0;
	for (int x=0; x<xSize; x++)
		for(int y=0; y<ySize; y++)
			for(int t=0; t<tSize; t++)
				for (int p=0; p<nTupleSize; p++)
				{
					imageDataType currVal = get_value(x,y,t,p);
					if (currVal>maxVal)
						maxVal = currVal;
				}
	return(maxVal);
}

imageDataType nTupleVolume::min_value()
{
	imageDataType minVal = 100000;
	for (int x=0; x<xSize; x++)
		for(int y=0; y<ySize; y++)
			for(int t=0; t<tSize; t++)
				for (int p=0; p<nTupleSize; p++)
				{
					imageDataType currVal = get_value(x,y,t,p);
					if (currVal < minVal)
						minVal = currVal;
				}
	return(minVal);
}

void nTupleVolume::absolute_value()
{
	for (int x=0; x<xSize; x++)
		for(int y=0; y<ySize; y++)
			for(int t=0; t<tSize; t++)
				for (int p=0; p<nTupleSize; p++)
				{
					imageDataType currVal = get_value(x,y,t,p);
					set_value(x,y,t,p,(imageDataType)fabs(currVal));
				}
}

imageDataType nTupleVolume::mean_value()
{
	imageDataType meanVal = 0;
	int numberOfElements = 0;
	for (int x=0; x<xSize; x++)
		for(int y=0; y<ySize; y++)
			for(int t=0; t<tSize; t++)
				for (int p=0; p<nTupleSize; p++)
				{
					meanVal = meanVal + get_value(x,y,t,p);
					numberOfElements++;
				}
	return(meanVal/numberOfElements);
}

void nTupleVolume::binarise()
{
	for (int x=0; x<xSize; x++)
		for(int y=0; y<ySize; y++)
			for(int t=0; t<tSize; t++)
				for (int p=0; p<nTupleSize; p++)
				{
					imageDataType currValue = get_value(x,y,t,p);
					if (currValue>0)	//set this value to 1
					{
						set_value(x,y,t,p,(imageDataType)1);
					}
					else
					{
						set_value(x,y,t,p,(imageDataType)0);
					}
				}
}


void clamp_coordinates(nTupleVolume* imgVolA, int *x, int *y, int *t)
{
    //x coordinates
    *x = max_int(min_int( *x, (imgVolA->xSize)-1),0);
    *y = max_int(min_int( *y, (imgVolA->ySize)-1),0);
    *t = max_int(min_int( *t, (imgVolA->tSize)-1),0);
}

void copy_pixel_values_nTuple_volume(nTupleVolume *imgVolA, nTupleVolume *imgVolB, int xA, int yA, int tA, int xB, int yB, int tB)
{
	if(imgVolA->nTupleSize != imgVolB->nTupleSize)
	{
		MY_PRINTF("Here copy_pixel_values_nTuple_volume. Error, the image volumes do not have the same number of channels.\n");
	}
	
	for (int p=0; p<imgVolA->nTupleSize; p++)
	{
		imageDataType pixelValueTemp = imgVolA->get_value(xA,yA,tA,p);
		imgVolB->set_value(xB,yB,tB,p,pixelValueTemp);
	}

}


void copy_pixel_values_nTuple_image(nTupleVolume *imgVolA, nTupleVolume *imgVolB, int x1, int y1, int x2, int y2)
{
	if(imgVolA->nTupleSize != imgVolB->nTupleSize)
	{
		MY_PRINTF("Here copy_pixel_values_nTuple_image. Error, the image volumes do not have the same number of channels.\n");
	}

	for (int p=0; p<imgVolA->nTupleSize; p++)
	{
		imageDataType pixelValueTemp = imgVolA->get_image_value(x1,y1,p);
		imgVolB->set_image_value(x2,y2,p,pixelValueTemp);
	}

}


nTupleVolume* copy_image_nTuple(nTupleVolume *imgVol)
{
	nTupleVolume *imgVolOut = new nTupleVolume(imgVol->nTupleSize, imgVol->xSize, imgVol->ySize, imgVol->patchSizeX, imgVol->patchSizeY, imgVol->indexing);
	
	for (int x=0; x< (int)imgVolOut->xSize; x++)
		for (int y=0; y< (int)imgVolOut->ySize; y++)
			for (int p=0; p< (int)imgVolOut->nTupleSize; p++)
			{
				imageDataType valTemp = imgVol->get_image_value(x,y,p);
				imgVolOut->set_image_value(x,y,p,valTemp);
			}
	return(imgVolOut);
}

imageDataType calculate_residual(nTupleVolume *imgVol, nTupleVolume *imgVolPrevious, nTupleVolume *occVol)
{
	imageDataType residual = 0.0;
	int sumOcc = 0;
	
	for (int x=0; x<(int)imgVol->xSize; x++)
		for (int y=0; y<(int)imgVol->ySize; y++)
			for (int t=0; t<(int)imgVol->tSize; t++)
				for (int p=0; p<(int)imgVol->nTupleSize; p++)
					if (occVol->get_value(x,y,t,0) > 0)
					{
						residual = residual + (imageDataType)fabs( (float)(imgVol->get_value(x,y,t,p)) - (float)(imgVolPrevious->get_value(x,y,t,p))); 
						sumOcc++;
					}
	return( residual/( (imageDataType)sumOcc));
}





