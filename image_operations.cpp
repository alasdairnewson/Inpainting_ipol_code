

#include "image_operations.h"

void seed_random_numbers( double inputSeed)
{
 	srand ( (unsigned int)inputSeed );   //create random number seed
}

float * read_image(const char *fileIn, size_t *nx, size_t *ny, size_t *nc)
{

	float * pixel_stream = NULL;
	pixel_stream = read_png_f32(fileIn,nx,ny,nc);

	if (pixel_stream == NULL)
	{
		printf("Unable to get the image\n");
		return(NULL);
	}
	
	printf("xSize : %d\n",(int)*nx);
	printf("ySize : %d\n",(int)*ny);
	printf("nChannels : %d\n",(int)*nc);
	
	return(pixel_stream);
}

nTupleVolume *make_colour_wheel()
{

	int RY = 15;
	int YG = 6;
	int GC = 4;
	int CB = 11;
	int BM = 13;
	int MR = 6;
	
	int nCols = RY + YG + GC + CB + BM + MR;
	
	nTupleVolume * colourWheel = new nTupleVolume(3,nCols,1,IMAGE_INDEXING);
	colourWheel->set_all_image_values(0);
	
	int col = 0;
	//RY
	for (int i=0;i<RY; i++)
	{
		colourWheel->set_image_value(i,0,0,(imageDataType)255);
		colourWheel->set_image_value(i,0,1,(imageDataType)floor(255*((float)i) / ( (float)RY)) );
	}
	col = col+RY;

	//YG
	for (int i=col;i<(col+YG); i++)
	{
		colourWheel->set_image_value(i,0,0,(imageDataType)(255 - (float)floor( 255 * ( ((float)i-col) / ((float)YG) ) ) ) );
		colourWheel->set_image_value(i,0,1,(imageDataType)255);
	}
	col = col+YG;

	//GC
	for (int i=col;i<(col+GC); i++)
	{
		colourWheel->set_image_value(i,0,1,(imageDataType)255);
		colourWheel->set_image_value(i,0,2,(imageDataType) floor( (float) 255 * ( (float) (i-col) )/  ((float)GC)) );
	}
	col = col+GC;

	//CB
	for (int i=col;i<(col+CB); i++)
	{
		colourWheel->set_image_value(i,0,1,(imageDataType) 255 - floor(255*( (float) (i-col) / ( (float) CB ) ) ) );
		colourWheel->set_image_value(i,0,2,(imageDataType)255);
	}
	col = col+CB;

	//BM
	for (int i=col;i<(col+BM); i++)
	{
		colourWheel->set_image_value(i,0, 2, (imageDataType) 255);
		colourWheel->set_image_value(i,0, 0, (imageDataType) floor(255* (float) ( (float)(i-col))/( (float)BM)) );
	}
	col = col+BM;

	//MR
	for (int i=col;i<(col+MR); i++)
	{
		colourWheel->set_image_value(i,0, 2, (imageDataType) 255 - floor( (float) 255* ( (float) (i-col) )/ ((float) MR)) );
		colourWheel->set_image_value(i,0, 0, (imageDataType)255);
	}
	
	return(colourWheel);

}

void write_image(nTupleVolume *imgVol, const char *fileName, imageDataType normalisationScalar)
{
	int readWriteSuccess;

		std::ostringstream os;
		os << fileName << ".png";
		const char* stringOut = (os.str()).c_str();
		
		//write image
		if (normalisationScalar>0)	//for better visulisation of results
		{
			//set minimum value to 0
			imageDataType minValue = imgVol->min_value();
			imgVol->add(-minValue);
			//normalise values
			imageDataType maxValue = imgVol->max_value();
			imgVol->multiply((imageDataType)1/maxValue);
			//set maximum value to normalisationScalar
			imgVol->multiply((imageDataType)normalisationScalar);
		}
		
	  	readWriteSuccess = write_png_f32((char*)stringOut, imgVol->get_data_ptr(),
		              imgVol->xSize, imgVol->ySize, imgVol->nTupleSize);
		if (readWriteSuccess == -1)
		{
			printf("Unable to write the image\n");
			return;
		}
}

void write_shift_map(nTupleVolume *shiftVol, const char *fileName)
{
	int readWriteSuccess;
	
	imageDataType maxNorm;
	maxNorm = sqrt( (shiftVol->xSize)*(shiftVol->xSize) + (shiftVol->ySize)*(shiftVol->ySize));
	
	nTupleVolume *shiftVolColour;
	shiftVolColour = new nTupleVolume(3, shiftVol->xSize, shiftVol->ySize, IMAGE_INDEXING);
	
	//create colour wheel
	nTupleVolume *colourWheel = make_colour_wheel();
	
	//convert shift volume to colour
	for (int x=0; x<(shiftVol->xSize); x++)
		for (int y=0; y<(shiftVol->ySize); y++)
		{
			imageDataType u,v;
			u = shiftVol->get_image_value(x,y,0);
			v = shiftVol->get_image_value(x,y,1);
			
			//calculate norm and normalise this with respect to the maximum
			//distance in the image (the diagonal)
			float vectorNorm = sqrt(u*u + v*v);
			vectorNorm = vectorNorm/maxNorm;
			
			//calculate angle
			float a = atan2(-v, -u)/(M_PI);
			
			float fk = ((a+1) /2) * ((float) ((colourWheel->xSize)-1) ) + 1;  // -1~1 maped to 1~ncols
			int k0 = (int)floor(fk);
			int k1 = k0+1;
			if (k1>= (colourWheel->xSize) )
				k1 = (colourWheel->xSize)-1;
			
			float f = fk -k0;
			
			for (int i=0; i< (colourWheel->nTupleSize); i++)
			{
				float col0 = (float) ((float)colourWheel->get_image_value(k0,0,i))/((float)255);
				float col1 = (float) ((float)colourWheel->get_image_value(k1,0,i))/((float)255);
				float col = (1-f)*col0 + f * col1;
				
				col = 1-vectorNorm * (1 - col);
				
				shiftVolColour->set_image_value(x,y,i,(imageDataType)col);
			}
		}
	
	//set minimum value to 0
	imageDataType minValue = shiftVolColour->min_value();
	shiftVolColour->add(-minValue);
	//normalise values
	imageDataType maxValue = shiftVolColour->max_value();
	shiftVolColour->multiply((imageDataType)1/maxValue);
	//set maximum value to normalisationScalar
	shiftVolColour->multiply((imageDataType)255);
	
	std::ostringstream os;
	os << fileName << "_shift_map.png";
	const char* stringOut = (os.str()).c_str();
	
  	readWriteSuccess = write_png_f32((char*)stringOut, shiftVolColour->get_data_ptr(),
	              shiftVolColour->xSize, shiftVolColour->ySize, shiftVolColour->nTupleSize);
	if (readWriteSuccess == -1)
	{
		printf("Unable to write the image\n");
		delete(colourWheel);
		return;
	}
	delete(colourWheel);

}

void write_image_pyramid(nTupleVolumePyramid imgVolPyramid, int nLevels, char *fileName, imageDataType normalisationScalar)
{
	int readWriteSuccess;
	for (int i=0; i<nLevels; i++)
	{
		nTupleVolume * currImgVol = imgVolPyramid[i];
		
		std::ostringstream os;
		os << fileName << "_level_" << int_to_string(i) << ".png";
		const char* stringOut = (os.str()).c_str();
		
		//write image
		if (normalisationScalar>0)	//for better visulisation of results
		{
			//set minimum value to 0
			imageDataType minValue = currImgVol->min_value();
			currImgVol->add(-minValue);
			//normalise values
			imageDataType maxValue = currImgVol->max_value();
			currImgVol->multiply((imageDataType)1/maxValue);
			//set maximum value to normalisationScalar
			currImgVol->multiply((imageDataType)normalisationScalar);
		}
		
	  	readWriteSuccess = write_png_f32((char*)stringOut, currImgVol->get_data_ptr(),
		              currImgVol->xSize, currImgVol->ySize, currImgVol->nTupleSize);
		if (readWriteSuccess == -1)
		{
			printf("Unable to write the image\n");
			return;
		}
	}
}

//image subsampling : the image is subsampled by 1/subSampleFactor
// that is to say : the size of imgVolOut is (1/subSampleFactor)*size(imgVol)
//Careful : this function does not deal with aliasing problems
nTupleVolume * sub_sample_image(nTupleVolume *imgVol, float subSampleFactor)
{
	int xSizeOut,ySizeOut;
	
	xSizeOut = (int) ceil( (imgVol->xSize)/subSampleFactor);
	ySizeOut = (int) ceil( (imgVol->ySize)/subSampleFactor);
	
	nTupleVolume * imgVolOut = new nTupleVolume(imgVol->nTupleSize,xSizeOut,ySizeOut,imgVol->patchSizeX,imgVol->patchSizeY,imgVol->indexing);
	
	for (int x=0; x< (imgVolOut->xSize); x++)
		for (int y=0; y<(imgVolOut->ySize); y++)
			for (int p=0; p<imgVolOut->nTupleSize; p++)
				imgVolOut->set_image_value((int)x,(int)y,p,imgVol->get_image_value(subSampleFactor*x,subSampleFactor*y,p));
				
	return(imgVolOut);
}

//image upsampling : the image is upsampled by upSampleFactor
//that is to say : the size of imgVolOut is upSampleFactor*size(imgVol)
//the image is upsamled by nearest neighbour interpolation
nTupleVolume * up_sample_image(nTupleVolume *imgVol, float upSampleFactor, nTupleVolume *imgVolFine)
{
	int xSizeOut,ySizeOut;
	
	if (imgVolFine == NULL)
	{
		xSizeOut = (int) ceil( (imgVol->xSize)*upSampleFactor);
		ySizeOut = (int) ceil( (imgVol->ySize)*upSampleFactor);
	}
	else
	{
		xSizeOut = (int) (imgVolFine->xSize);
		ySizeOut = (int) (imgVolFine->ySize);
	}
	
	nTupleVolume * imgVolOut = new nTupleVolume(imgVol->nTupleSize,xSizeOut,ySizeOut,imgVol->patchSizeX,imgVol->patchSizeY,imgVol->indexing);
	
	for (int x=0; x< (imgVolOut->xSize); x++)
		for (int y=0; y<(imgVolOut->ySize); y++)
			for (int p=0; p<imgVolOut->nTupleSize; p++)
			{
				imgVolOut->set_image_value((int)x,(int)y,p,
				imgVol->get_image_value((int)(floor((1/upSampleFactor)*x)), (int)(floor((1/upSampleFactor)*y)),p));
			}
				
	return(imgVolOut);
}

nTupleVolume * rgb_to_grey(nTupleVolume * imgVolIn)
{
	nTupleVolume *imgGreyOut = new nTupleVolume(1,imgVolIn->xSize,imgVolIn->ySize,imgVolIn->patchSizeX,imgVolIn->patchSizeY,imgVolIn->indexing);
	
	for (int x=0; x<(imgVolIn->xSize); x++)
		for (int y=0; y<(imgVolIn->ySize); y++)
		{
			float greyValue = 0.2126*imgVolIn->get_image_value(x,y,0) + 0.7152*imgVolIn->get_image_value(x,y,1) + 0.0722*imgVolIn->get_image_value(x,y,2);
			imgGreyOut->set_image_value(x,y,0,(imageDataType)greyValue);
		}
	
	return(imgGreyOut);
}

nTupleVolume * image_gradient_x(nTupleVolume * imgVolIn)
{
	int destroyGreyImg = 0;
	nTupleVolume *imgVolGrey;
	nTupleVolume *gradX = new nTupleVolume(1,imgVolIn->xSize,imgVolIn->ySize,imgVolIn->patchSizeX,imgVolIn->patchSizeY,imgVolIn->indexing);
	
	//if we need to convert the image to greyscale
	if (imgVolIn->nTupleSize ==3)
	{
		imgVolGrey = rgb_to_grey(imgVolIn);
		destroyGreyImg = 1;
	}
	else
		imgVolGrey = imgVolIn;


	for (int x=0; x<(imgVolIn->xSize); x++)
		for (int y=0; y<(imgVolIn->ySize); y++)
		{
			int xMin = max_int(x-1,0);
			int xMax = min_int(x+1,(imgVolIn->xSize)-1);
			
			//gradient calculation
			imageDataType gradTemp = ( imgVolGrey->get_image_value(xMax,y,0) - imgVolGrey->get_image_value(xMin,y,0) )
			/((imageDataType)xMax-xMin);
			gradX->set_image_value(x,y,0,gradTemp);
		}
	
	if (destroyGreyImg>0)
		delete imgVolGrey;
	return gradX;

}

nTupleVolume * image_gradient_y(nTupleVolume * imgVolIn)
{
	int destroyGreyImg = 0;
	nTupleVolume *imgVolGrey;
	nTupleVolume *gradY = new nTupleVolume(1,imgVolIn->xSize,imgVolIn->ySize,imgVolIn->patchSizeX,imgVolIn->patchSizeY,imgVolIn->indexing);
	
	//if we need to convert the image to greyscale
	if (imgVolIn->nTupleSize ==3)
	{
		imgVolGrey = rgb_to_grey(imgVolIn);
		destroyGreyImg = 1;
	}
	else
		imgVolGrey = imgVolIn;


	for (int x=0; x<(imgVolIn->xSize); x++)
		for (int y=0; y<(imgVolIn->ySize); y++)
		{
			int yMin = max_int(y-1,0);
			int yMax = min_int(y+1,(imgVolIn->ySize)-1);
			
			//gradient calculation
			imageDataType gradTemp = ( imgVolGrey->get_image_value(x,yMax,0) - imgVolGrey->get_image_value(x,yMin,0) )
			/((imageDataType)yMax-yMin);
			gradY->set_image_value(x,y,0,gradTemp);
		}
	
	if (destroyGreyImg>0)
		delete imgVolGrey;
	return gradY;
}

nTupleVolumePyramid create_nTupleVolume_pyramid_binary(nTupleVolume * imgVol, int nLevels)
{
	nTupleVolumePyramid pyramidOut = (nTupleVolume**)malloc( (size_t)nLevels*sizeof(nTupleVolume*));
	
	float sigmaFilter = 1.5;
	int convSize = 3;
	nTupleVolume *convKernel = create_convolution_kernel("gaussian", convSize, convSize, 1, sigmaFilter);
	
	
	for (int i=0; i<nLevels; i++)
	{
		if(i==0)	//first level
		{
			nTupleVolume* imgVolTemp = copy_image_nTuple(imgVol);
			imgVolTemp->binarise();
			pyramidOut[i] = imgVolTemp;
		}
		else
		{
			nTupleVolume * imgVolTemp = normalised_convolution_masked(pyramidOut[i-1],convKernel);
			pyramidOut[i] = sub_sample_image(imgVolTemp, 2);
			(pyramidOut[i])->binarise();
			delete(imgVolTemp);
		}
	}
	
	//delete subsampling convolution filter
	delete(convKernel);
	
	return(pyramidOut);

}

nTupleVolumePyramid create_nTupleVolume_pyramid(nTupleVolume * imgVol, int nLevels)
{
	nTupleVolumePyramid pyramidOut = (nTupleVolume**)malloc( (size_t)nLevels*sizeof(nTupleVolume*));
	
	float sigmaFilter = 1.5;
	int convSize = 3;
	nTupleVolume *convKernel = create_convolution_kernel("gaussian", convSize, convSize, 1, sigmaFilter);
	
	
	for (int i=0; i<nLevels; i++)
	{
		if(i==0)	//first level
		{
			nTupleVolume* imgVolTemp = copy_image_nTuple(imgVol);
			pyramidOut[i] = imgVolTemp;
		}
		else
		{
			nTupleVolume * imgVolTemp = normalised_convolution_masked(pyramidOut[i-1],convKernel);
			pyramidOut[i] = sub_sample_image(imgVolTemp, 2);
			delete(imgVolTemp);
		}
	}
	
	//delete subsampling convolution filter
	delete(convKernel);
	
	return(pyramidOut);

}

featurePyramid create_feature_pyramid(nTupleVolume * imgVol, nTupleVolume * occVol, int nLevels)
{
	featurePyramid featurePyramidOut;
	nTupleVolume* imgVolTemp;
	//create the smoothing filter
	nTupleVolume *convKernelX = create_convolution_kernel("average", (int)pow(2,(float)(nLevels)), 1);
	nTupleVolume *convKernelY = create_convolution_kernel("average", 1, (int)pow(2,(float)(nLevels)));
	//nTupleVolume *convKernel = create_convolution_kernel("average", (int)pow(2,(float)(nLevels)), (int)pow(2,(float)(nLevels)));
	
	//calculate the image gradient
	nTupleVolume * imgVolGradX = image_gradient_x(imgVol);
	nTupleVolume * imgVolGradY = image_gradient_y(imgVol);
	
	//take absolute value of gradients and then average these values
	imgVolGradX->absolute_value();
	imgVolGradY->absolute_value();
	nTupleVolume * imgVolGradXavg = normalised_convolution_masked_separable(imgVolGradX, convKernelX,convKernelY, occVol);//normalised_convolution_masked(imgVolGradX, convKernel, occVol);//
	nTupleVolume * imgVolGradYavg = normalised_convolution_masked_separable(imgVolGradY, convKernelX,convKernelY, occVol);//normalised_convolution_masked(imgVolGradY, convKernel, occVol);//
	
	write_image(imgVolGradXavg,"American_norm_grad_x",255);
	write_image(imgVolGradYavg,"American_norm_grad_y",255);
	
	nTupleVolumePyramid normGradXPyramid = (nTupleVolume**)malloc( (size_t)nLevels*sizeof(nTupleVolume*));
	nTupleVolumePyramid normGradYPyramid = (nTupleVolume**)malloc( (size_t)nLevels*sizeof(nTupleVolume*));
	
	//subsample images
	for (int i=0; i<nLevels; i++)
	{
		if(i==0)	//first level
		{
			//x gradient
			imgVolTemp = copy_image_nTuple(imgVolGradXavg);
			normGradXPyramid[i] = imgVolTemp;
			//y gradient
			imgVolTemp = copy_image_nTuple(imgVolGradYavg);
			normGradYPyramid[i] = imgVolTemp;
		}
		else
		{
			//x gradient
			imgVolTemp = normGradXPyramid[i-1];
			normGradXPyramid[i] = sub_sample_image(imgVolTemp, 2);
			//x gradient
			imgVolTemp = normGradYPyramid[i-1];
			normGradYPyramid[i] = sub_sample_image(imgVolTemp, 2);
		}
	}
	
	featurePyramidOut.normGradX = normGradXPyramid;
	featurePyramidOut.normGradY = normGradYPyramid;
	featurePyramidOut.nLevels = nLevels;
	
	delete imgVolGradX;
	delete imgVolGradY;
	delete imgVolGradXavg;
	delete imgVolGradYavg;
	delete convKernelX;
	delete convKernelY;
	//delete convKernel;
	
	return(featurePyramidOut);
}


void delete_feature_pyramid(featurePyramid featurePyramidIn)
{
	 nTupleVolumePyramid normGradXpyramid = featurePyramidIn.normGradX;
	 nTupleVolumePyramid normGradYpyramid = featurePyramidIn.normGradY;
	 
	 for (int i=0; i<(featurePyramidIn.nLevels); i++)
	 {
	 	delete(normGradXpyramid[i]);
	 	delete(normGradYpyramid[i]);
	 }
	delete(normGradXpyramid);
	delete(normGradYpyramid);
}

int determine_multiscale_level_number(nTupleVolume *occVolIn, int patchSizeX, int patchSizeY)
{

	int nLevels;
	int maxOccDistance=0;
	int maxPatchSize = (int) max_int(patchSizeX,patchSizeY);
	nTupleVolume *structElErode = create_structuring_element("rectangle", 3, 3);
	
	nTupleVolume *occVol = copy_image_nTuple(occVolIn);
	
	while(occVol->sum_nTupleVolume() >0)
	{
		nTupleVolume *occVolTemp = imerode(occVol,structElErode);
		delete(occVol);
		occVol = copy_image_nTuple(occVolTemp);
		delete(occVolTemp);
		maxOccDistance++;
	}
	
	maxOccDistance = 2*maxOccDistance;
	
	nLevels = (int) floor( (float)
				(
				log( ((float) maxOccDistance) /((float)maxPatchSize) )
					) /
					( (float) log(2) )
					);
	
	delete(occVol);
	return(nLevels);

}



