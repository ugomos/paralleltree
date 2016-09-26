%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Code of the algorithm described in the paper "A Hybrid Shared-Memory Parallel Max-Tree Algorithm for Extreme Dynamic-Range Images" by Ugo Moschini, Arnold Meijster and Michael H.F. Wilkinson, University of Groningen, The Netherlands

November 2015, Ugo Moschini (u.moschini@rug.nl , ugomoschini@gmail.com)

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

The present code was tested with:      - 2D .tif images (8 - 16 bits integers per pixel)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Example of output: 
./main 2 img_smallcube_float.fits 0 32 outcube.fits 1

Statistics of 60 x 60 x 70 image. Number of axis = 3. Bits per pixel=-32.
FreeImage version: 3.15.4
Filtering image 'img_smallcube_float.fits' using attribute area with threshold lambda=0.000000.
Image img_smallcube_float.fits: Width=60 Height=60 Depth=70 Size=252000 Size2D=3600. number of threads of the parallel algorithm: 2.
On this machine: Size of unsigned short=16. Size of int=32. Size of long=64. Size of float=32. Size of double=64. (bits per pixel)
/*** Sort the pixels ***/
Radix Sort (steps=2)
/*** Calculate the quantized image ***/
/*** Build the max tree of the quantized image. (threads 2)***/
Pilot max-tree built.
/*** Refinement phase (threads 2) ***/
Refined max-tree built.
Init filtering
0) lwb=0; upb=126000.
1) lwb=126000; upb=252000.
End filtering
Sorting: 0.020000 s.
Create Quantized Image: 0.000000 s.
Quantized Tree: 0.020000 s.
Refinement Tree: 0.060000 s.
Filtering: 0.000000 s.
Wall-Clock time: 0.100000 s.
Image written to 'outcube.fits'
