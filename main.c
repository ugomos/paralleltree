/**                       main.c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Code of the algorithm described in the paper "A Hybrid Shared-Memory Parallel Max-Tree Algorithm for Extreme Dynamic-Range Images" by Ugo Moschini, Arnold Meijster and Michael H.F. Wilkinson, University of Groningen, The Netherlands

November 2015, Ugo Moschini (u.moschini@rug.nl)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Compilation of the sources
	- The code can be compiled and built using the script ./compile.sh and ./build.sh.
	- Required libraries: lpthread, lcfitsio, lfreeimage
	* Compilation: gcc -O2 -Wall -std=gnu99 -pedantic -fexpensive-optimizations -funroll-loops -c main.c quanttree.c handleimages.c quantizedimage.c refinetree.c 	                 radixsort.c filter.c -lpthread -lcfitsio -lfreeimage
	* Building: g++ -o main -L$HOME/lib -I$HOME/lib  main.o quanttree.o handleimages.o quantizedimage.o refinetree.o radixsort.o filter.o -lpthread -lcfitsio -lfreeimage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Running the code:

Usage: ./main <nthreads> <input image> <lambda> <bits per pixel> <output image> <is3D>.
	<nthreads> : number of threads of the parallel program.
	<input image> : the image whose tree is built.
	<lamba> : the algorithm performs area filtering: all the connected components with area smaller than lambda will be deleted in the output image.
        <bits per pixel> : number of bits per pixel of the data type carried by the input image (8,16,32,...).
        <output image> : the filtered output image (components smaller than lamba are removed).
        <is3D> : '0' means the input image is a 2D image; '1' means the input image is a 3D volume.

The present code was tested with:  - 2D .tif images (8 - 16 bits integers per pixel)
	                    	       - 2D .fits images (32 bits floating point)
                   		           - 3D .fits volumes (32 bits floating point)
We included 4 test images:
1. img_smallcube_float.fits: 32-bit floating point 3D. Contains a section of the Westerbork radio astronomy cube, courtesy of P. Serra.
2. img_galaxy.fit: 32-bit floating point 2D. 
3. img_galaxy.tif: 16-bit integers 2D.
4. img_mountain.tif: 8-bit integer 2D.

To test the four images, the code must be run with the following parameters (e.g., on 2 threads, lambda=0), respectively:
(1) ./main 2 img_smallcube_float.fits 0 32 out.fits 1
(2) ./main 2 img_galaxy.fit 0 32 out.fits 0
(3) ./main 2 img_galaxy.tif 0 16 out.tif 0
(4) ./main 2 img_mountain.tif 0 8 out.tif 0

Note: the file 'common.h' contains two important lines, initialized like:
#define USEFLOATPOINT 1  // '1' input image has floating point data; '0' input image has integer data.
typedef float greyval_t; // data type carried by the pixel (unsigned short, int, float, double, ...)

- Such inizialitation is valid for image (1) and (2).
It means that the input image is expected to be a floating point image. 

- For image (3) and (4) the lines would become:
#define USEFLOATPOINT 0 // '1' input image has floating point data; '0' input image has integer data.
typedef unsigned short greyval_t; // data type carried by the pixel (unsigned short, int, float, double, ...)

 
 * */

#include "common.h"
#include "quanttree.h"
#include "handleimages.h"
#include "radixsort.h"
#include "quantizedimage.h"
#include "refinetree.h"
#include "filter.h"

char *getExtension(const char *filename)
{
    char *dot = strrchr(filename, '.'); // A pointer to the last occurrence of character in str.
    return (char *)dot + 1;
}

void InitBuckets()
{
    NUMBUCKETS = 65536;
}

void printTimings()
{
    musec = (float)(end_sort - start_sort)/((float)tickspersec);
    printf("Sorting: %f s.\n", musec);

    musec = (float)(end_quimg - start_quimg)/((float)tickspersec);
    printf("Create Quantized Image: %f s.\n", musec);

    musec = (float)(end_qutree - start_qutree)/((float)tickspersec);
    printf("Quantized Tree: %f s.\n", musec);

    musec = (float)(end_ref - start_ref)/((float)tickspersec);
    printf("Refinement Tree: %f s.\n", musec);

    musec = (float)(end_filt - start_filt)/((float)tickspersec);
    printf("Filtering: %f s.\n", musec);

    musec = (float)(end_filt - start_sort)/((float)tickspersec);
    printf("Wall-Clock time: %f s.\n", musec);
}


int main (int argc, char *argv[])
{
    char *imgfname, *outfname;
    is3D = false;
    bitsPerPixel = 8;
    int r = 0;
    tickspersec = sysconf(_SC_CLK_TCK);

    if (argc<7)
    {
        printf("Usage: %s <nthreads> <input image> <lambda> <bits per pixel> <output image> <is3D>.\n", argv[0]);
        exit(0);
    }

    nthreads = MIN(atoi(argv[1]), MAXTHREADS); // number of threads used for the parallel sorting algorithm and for the parallel building of the quantized max tree
    inputQTZLEVELS = nthreads;
    imgfname = argv[2];
    lambda = atof(argv[3]);
    bitsPerPixel = atoi(argv[4]);
    outfname = argv[5];
    is3D = (atoi(argv[6]));
    width = height = depth = size = size2D = 0;

    char *exto = getExtension(imgfname);
    if( (strcmp(exto, "tif") == 0))
        if(!ReadTIFF(imgfname)) exit(0);

    if(strcmp(exto, "fits") == 0 || strcmp(exto, "fit") == 0)
    {
        if(is3D)
        {
            if(!ReadFITS3D(imgfname))
            {
                printf("fits image failed to read\n");
                return(-1);
            }
        }
        else
            ReadFITS(imgfname);
    }

    size = width*height*depth;
    if(is3D)
        size2D = width*height;
    else
        size2D = width;


    //printf("Command: %s\n", argv[0]);
    printf("FreeImage version: %s\n", FreeImage_GetVersion());
    printf("Filtering image '%s' using attribute area with threshold lambda=%f.\n", imgfname, lambda);
    printf("Image %s: Width=%ld Height=%ld Depth=%ld Size=%ld Size2D=%ld. ", imgfname, width, height, depth, size, size2D);
    printf("number of threads of the parallel algorithm: %d.\n", nthreads);
    printf("On this machine: Size of unsigned short=%lu. Size of int=%lu. Size of long=%ld. Size of float=%ld. Size of double=%ld. (bits per pixel)\n", sizeof(unsigned short)*8, sizeof(int)*8, sizeof(long)*8, sizeof(float)*8, sizeof(double)*8);
    //hmin = FindMin(gval); hmax = FindMax(gval);
    //printf("Min=%d. Max=%d.\n", hmin, hmax);

    /**************************************************************/
    printf("/*** Sort the pixels ***/\n");
    InitBuckets();
    start_sort = times(&tstruct);
    // run radix sort
    RunRadixSortUnsigned((int) ceil((double)bitsPerPixel/(16)));
    end_sort = times(&tstruct);

    /******************************************************************/
    printf("/*** Calculate the quantized image ***/\n");
    start_quimg = times(&tstruct);
    CalculateQuantizedImage();
    end_quimg = times(&tstruct);
    musec = (float)(end_quimg - start_quimg)/((float)tickspersec);

    /******************************************************************/
    printf("/*** Build the max tree of the quantized image. (threads %d)***/\n", nthreads);
    start_qutree = times(&tstruct);
    BuildMaxTreeOfQuantizedImage();
    end_qutree = times(&tstruct);
    printf("Pilot max-tree built.\n");

    /******************************************************************/
    nthreadsRef = numQTZLEVELS;
    printf("/*** Refinement phase (threads %d) ***/\n", nthreadsRef);
    outRef =  malloc(size*sizeof(greyval_t));
    RefineTreeBerger(nthreadsRef);
    free(gval_qu);
    printf("Refined max-tree built.\n");

    /******************************************************************/
    printf("Init filtering\n");
    start_filt = times(&tstruct);
    Filter();
    end_filt = times(&tstruct);
    printf("End filtering\n");
    /******************************************************************/
    printTimings();

    /**** Write output image   ****/
    exto = getExtension(outfname);
    if((strcmp(exto, "tif") == 0))
    {
        r = WriteTIFF(outfname, outRef);
        if (!r)  printf("Image written to '%s'\n", outfname);
    }

    if(strcmp(exto, "fit") == 0 || strcmp(exto, "fits") == 0  )
    {
        if(is3D)
        {
            WriteFITS3D(outfname, imgfname, outRef);
            printf("Image written to '%s'\n", outfname);
        }
        else //2D
        {
            WriteFITS(outfname, imgfname, outRef);
            printf("Image written to '%s'\n", outfname);
        }
    }

    printf("\n");
    free(gval);
    free(SORTED);
    free(zpar);
    free(node_qu);
    free(node_ref);
    free(outRef);

    FreeImage_DeInitialise();
    return(0);
}
