Author : Alasdair Newson
Date : 14 April 2014
Subject : Image inpainting code guide

Alasdair Newson, <alasdairnewson@gmail.com>, Duke University, USA
Andres Almansa, <andres.almansa@telecom-paristech.fr>, France
Yann Gousseau, <yann.gousseau@telecom-paristech.fr>, France
Patrick Perez, <patrick.perez@technicolor.com>, France

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This file explains how to use the image inpainting code in this directory.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% INSTALLATION %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

To install the software, you need to have a c++ compiler.
The code has been tested with the gcc compiler.

Extract the source code to a directory and cd to this directory
in a shell. Now simply execute the following command :

$ make

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% USING THE CODE %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

If this is done successfully, you can execute the image inpainting
code. To do this, type:

./inpaint_image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   PARAMETERS   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Adaptive weights parameters:
    -patchSizeX : patch column size (default : 7)
    -patchSizeY : patch row size (default : 7)
    -nLevels : number of pyramid levels (by default -1, which means that it
    is determined automatically by the algorithm)
    -useFeatures : use texture features, 0 = false, 1 = true (default, 1)
    -v : verbose, 0 = false, 1 = true (default, 0)

For the moment, the inpupt file name are hard coded into the program.
You can change this by renaming the "fileIn" and "fileOcc" strings in the
"inpaint_image_main.cpp" file.

The main body of the inpainting code may be found in "image_inpainting.cpp".

The output file will be written as "input_name_output.png".

Have fun !
