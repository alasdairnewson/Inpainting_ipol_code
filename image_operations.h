

//definitions for certain basic image operations

#ifndef IMAGE_OPERATIONS_H
#define IMAGE_OPERATIONS_H

#include <stdlib.h>

#include "image_structures.h"
#include "io_png.h"
#include "convolution.h"
#include "morpho.h"

//random seed
void seed_random_numbers( double inputSeed=0.0);

//visualisation tool
nTupleVolume *make_colour_wheel();

//reading and writing functions
float * read_image(const char *fileIn, size_t *nx, size_t *ny, size_t *nc);
void write_image(nTupleVolume *imgVol, const char *fileName, imageDataType normalisationScalar=0);
void write_image_pyramid(nTupleVolumePyramid imgVolPyramid, int nLevels, const char *fileName, imageDataType normalisationScalar=0);
void write_shift_map(nTupleVolume *shiftVol, const char *fileName);

nTupleVolume * sub_sample_image(nTupleVolume *imgVol, float subSampleFactor);
nTupleVolume * up_sample_image(nTupleVolume *imgVol, float upSampleFactor, nTupleVolume *imgVolFine=NULL);

nTupleVolume * rgb_to_grey(nTupleVolume * imgVolIn);

nTupleVolume * image_gradient_x(nTupleVolume * imgVolIn);
nTupleVolume * image_gradient_y(nTupleVolume * imgVolIn);

nTupleVolumePyramid create_nTupleVolume_pyramid_binary(nTupleVolume * imgVol, int nLevels);
nTupleVolumePyramid create_nTupleVolume_pyramid(nTupleVolume * imgVol, int nLevels);

featurePyramid create_feature_pyramid(nTupleVolume * imgVol, nTupleVolume * occVol, int nLevels);
void delete_feature_pyramid(featurePyramid featurePyramidIn);

int determine_multiscale_level_number(nTupleVolume *occVolIn, int patchSizeX, int patchSizeY);

#endif
