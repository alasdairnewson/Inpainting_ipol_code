

//definitions for the convolution operations

#ifndef CONVOLUTIONS_H
#define CONVOLUTIONS_H

#include "image_structures.h"

nTupleVolume * create_convolution_kernel(const char * kernelType, int xSizeKernel, int ySizeKernel, int tSizeKernel=1, float stdDev = 1.5);
nTupleVolume * normalised_convolution_masked(nTupleVolume *imgVol, nTupleVolume *convKernel, nTupleVolume *occlusionMask=NULL);
nTupleVolume * normalised_convolution_masked_separable(nTupleVolume *imgVol, nTupleVolume *convKernelX, nTupleVolume *convKernelY, nTupleVolume *occlusionMask=NULL);

#endif
