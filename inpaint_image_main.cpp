
//main function for the image inpainting routine


#include <cstring>
#include <cstdio>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "image_inpainting.h"

int main(int argc, char* argv[])
{

	time_t startTime,stopTime;

	const char *fileIn = (argc >= 2) ? argv[1] : "American.png";
	const char *fileInOcc = (argc >= 3) ? argv[2] : "American_occlusion.png";
	const char *fileOut = (argc >= 4) ? argv[3] : "American_inpainted.png";
	const char * patchSizeX = (argc >= 5) ? argv[4] : "7";
	const char * patchSizeY = (argc >= 6) ? argv[5] : "7";
	const char * nLevels = (argc >= 7) ? argv[6] : "-1";
	const char * useFeatures = (argc >= 8) ? argv[7] : "1";
	const char * verboseMode = (argc >= 9) ? argv[8] : "0";
	
	time(&startTime);//startTime = clock();
	
	inpaint_image_wrapper(fileIn,fileInOcc,fileOut,
		atoi(patchSizeX), atoi(patchSizeY), atoi(nLevels), (bool)atoi(useFeatures), (bool)atoi(verboseMode));
	
	time(&stopTime);
	printf("\n\nTotal execution time: %f\n",fabs(difftime(startTime,stopTime)));
	//MY_PRINTF("\n\nTotal execution time: %f\n",((double)(clock()-startTime)) / CLOCKS_PER_SEC);

	
	return(0);

}
