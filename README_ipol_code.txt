    Author : Alasdair Newson
Date : 14 April 2014
Subject : Image inpainting code guide

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

make

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% USING THE CODE %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

If this is done successfully, you can execute the image inpainting
code. To do this, type:

./inpaint_image

For the moment, the inpupt file name are hard coded into the program.
You can change this by renaming the "fileIn" and "fileOcc" strings in the
"inpaint_image_main.cpp" file.

The main body of the inpainting code may be found in "image_inpainting.cpp".

The output file will be written as "input_name_output.png".

Have fun !
