
#makefile for the image inpainting routine

CC=g++
PATCHMATCHOBJS = patch_match.o patch_match_tools.o patch_match_measure.o
RECONSTRUCTIONOBJS = reconstruct_image_and_features.o reconstruct_image.o reconstruct_image_tools.o
OBJS = inpaint_image_main.o image_inpainting.o $(RECONSTRUCTIONOBJS) $(PATCHMATCHOBJS) image_operations.o convolution.o morpho.o image_structures.o io_png.o
PATCHMATCHSOURCE = Patch_match/patch_match.cpp Patch_match/patch_match_tools.cpp Patch_match/patch_match_measure.cpp
PATCHMATCHHEADERS = Patch_match/patch_match.h Patch_match/patch_match_tools.h Patch_match/patch_match_measure.h Patch_match/common_patch_match.h
RECONSTRUCTIONSOURCE = Reconstruction/reconstruct_image_and_features.cpp Reconstruction/reconstruct_image.cpp Reconstruction/reconstruct_image_tools.cpp
RECONSTRUCTIONHEADERS = Reconstruction/reconstruct_image_and_features.h Reconstruction/reconstruct_image.h Reconstruction/reconstruct_image_tools.h Reconstruction/common_reconstruct_image.h
DEBUG = -g #-pg
OPENMPFLAG = -fopenmp
CXXFLAGS= -Wall -Wextra -O2 -stdlib=libstdc++ # $(DEBUG) 
LIB=-L/usr/informix/lib/c++
INC= -IReconstruction/ -IPatch_match/ -IImage_structures/

inpaint_image: $(OBJS)
	$(CC) $(CXXFLAGS) $(OBJS) -o inpaint_image -lpng
inpaint_image_main.o: inpaint_image_main.cpp image_inpainting.o $(PATCHMATCHHEADERS) morpho.h image_operations.h convolution.h Image_structures/image_structures.h io_png.h
	$(CC) $(CXXFLAGS) $(INC) -c inpaint_image_main.cpp
image_inpainting.o: image_inpainting.cpp image_inpainting.h $(PATCHMATCHHEADERS) morpho.h image_operations.h convolution.h Image_structures/image_structures.h io_png.h
	$(CC) $(CXXFLAGS) $(INC) -c image_inpainting.cpp
reconstruct_image_and_features.o: $(RECONSTRUCTIONSOURCE) $(RECONSTRUCTIONHEADERS) Image_structures/image_structures.h
	$(CC) $(CXXFLAGS) $(INC) -c $(RECONSTRUCTIONSOURCE)
patch_match.o: $(PATCHMATCHSOURCE) $(PATCHMATCHHEADERS) Image_structures/image_structures.h
	$(CC) $(CXXFLAGS) $(INC) -c $(PATCHMATCHSOURCE)
image_operations.o: image_operations.cpp image_operations.h convolution.h morpho.h Image_structures/image_structures.h
	$(CC) $(CXXFLAGS) $(INC) -c image_operations.cpp
convolution.o: convolution.cpp convolution.h Image_structures/image_structures.h
	$(CC) $(CXXFLAGS) $(INC) -c convolution.cpp
io_png.o: io_png.c io_png.h
	$(CC) $(CXXFLAGS) -c io_png.c
morpho.o: morpho.cpp morpho.h Image_structures/image_structures.h
	$(CC) $(CXXFLAGS) $(INC) -c morpho.cpp
image_structures.o: Image_structures/image_structures.cpp Image_structures/image_structures.h
	$(CC) $(CXXFLAGS) $(INC) -c Image_structures/image_structures.cpp
all: inpaint_image	

clean: inpaint_image_main.cpp
	rm -f *.o
