

//decalarations for morphological operations

#ifndef MORPHO_H
#define MORPHO_H

#include "image_structures.h"

nTupleVolume* create_structuring_element(const char * structType, int xSize, int ySize);

nTupleVolume* imerode(nTupleVolume* imgVol, nTupleVolume* structEl);

nTupleVolume* imdilate(nTupleVolume* imgVol, nTupleVolume* structEl);

#endif
