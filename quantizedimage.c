/** quantizedimage.c
 * The code computes the quantized image, quantizing the input image into a number of intensities
 * equal to the number of threads. 
 * 
 */
 
#include "common.h"
#include "quantizedimage.h"
#include "quanttree.h"

ThreadQntImgData *MakeThreadQntImgData(int numthreads)
{
    ThreadQntImgData *data = malloc(numthreads *sizeof(ThreadQntImgData));
    int i;

    for (i=0; i<numthreads; i++)
    {
        data[i].self=i;
    }
    return(data);
}

void FreeThreadQntImgData(ThreadQntImgData *data, int numthreads)
{
    free(data);
}

void RunThreadsWritingQuantizedImage(ThreadQntImgData *thdata, int nthreads)
{
    int thread;
    for (thread=0; thread<nthreads; ++thread)
    {
        pthread_create(threadID+thread, NULL, wqi, (void *) (thdata + thread));
    }
    for (thread=0; thread<nthreads; ++thread)
    {
        pthread_join(threadID[thread], NULL);
    }
}

/* CreateQuantizedImage */
int CreateQuantizedImage()
{
    pixel_t numRemPixels = size;
    pixel_t expectedNumPixelPerQtzLev = numRemPixels/(inputQTZLEVELS);

    int currentQtzLev = 0;
    pixel_t px=0, count = 0;
    greyval_t prevBoundGval = 0;
    int idxThread = 0;

    while(px<size)
    {
        pxStartPosition[idxThread]	= px;

        // iterate until the optimal number of pixels for the current thread
        while(px<size && (count < expectedNumPixelPerQtzLev || currentQtzLev == inputQTZLEVELS-1))
        {
            gval_qu[SORTED[px]] = currentQtzLev;
            count++;
            prevBoundGval = gval[SORTED[px]];
            px++;
        }

        // iterate until the next intensity change to map different intensities on different threads
        while(px<size && prevBoundGval == gval[SORTED[px]])
        {
            gval_qu[SORTED[px]] = currentQtzLev;
            count++;
            px++;
        }

        // update thread info
        pxEndPosition[idxThread] = px-1;

        // work out the optimal number of pixels for the next thread
        numRemPixels = numRemPixels - count;
        expectedNumPixelPerQtzLev = numRemPixels/(inputQTZLEVELS-currentQtzLev);

        currentQtzLev++;
        idxThread++;
        count=0;
    }

    return currentQtzLev; // return the number of quantization levels used
}

int WorkOutBoundaries()
{
    pixel_t numRemPixels = size;
    pixel_t expectedNumPixelPerQtzLev = numRemPixels/(inputQTZLEVELS);
    //printf("Optimal number of pixels per thread = %ld\n", expectedNumPixelPerQtzLev);

    int currentQtzLev = 0;
    pixel_t px=0;
    greyval_t prevBoundGval = 0;
    int idxThread = 0;

    while(px<size)
    {
        pxStartPosition[idxThread]	= px;

        px += (expectedNumPixelPerQtzLev - 1); // point to the last gvalue of that thread.

        if(px<size)
        {
            prevBoundGval = gval[SORTED[px]];
            px++;

            while( (px<size) && ( (prevBoundGval == gval[SORTED[px]]) || (currentQtzLev == numQTZLEVELS-1) ) )
            {
                px++;
            }

            pxEndPosition[idxThread] = px - 1;
        }
        else
            pxEndPosition[idxThread] = size-1;

        currentQtzLev++;
        idxThread++;
    }

    return currentQtzLev; // return the number of quantization levels used
}

void *wqi(void *arg)
{
    ThreadQntImgData *threfdata = (ThreadQntImgData *) arg;
    int self = threfdata->self;

    pixel_t i;
    for(i=pxStartPosition[self]; i <= pxEndPosition[self]; i++)
    {
        gval_qu[SORTED[i]] = self;
    }
    return NULL;
}

void CreateQuantizedImageInParallel(int numQtzLevsUsed)
{
    ThreadQntImgData *thqntimgdata = MakeThreadQntImgData(numQtzLevsUsed);
    RunThreadsWritingQuantizedImage(thqntimgdata, numQtzLevsUsed);
    FreeThreadQntImgData(thqntimgdata, numQtzLevsUsed);
}

void CalculateQuantizedImage()
{
    if(is3D && nthreads > depth)
    {
        printf("Exit. Number of threads must be equal to or larger than the third dimension of the cube.\n");
        exit(0);
    }

    /* create the quantized image */
    gval_qu = malloc(size*sizeof(int));
    if (gval_qu==NULL)
    {
        fprintf(stderr, "out of memory! \n");
        free(gval_qu);
        exit(0);
    }

    pxStartPosition = calloc(inputQTZLEVELS, sizeof(pixel_t));
    pxEndPosition = calloc(inputQTZLEVELS, sizeof(pixel_t));

    // First Way: in parallel
    numQTZLEVELS = WorkOutBoundaries();

    if(numQTZLEVELS != inputQTZLEVELS)
    {
        printf("Number of computed quantized levels %d differs from the input parameter %d. Trying an other method (sacrifice load balance).\n", numQTZLEVELS, inputQTZLEVELS);
        // Second Way: use the old sequential (worst load balance)
        numQTZLEVELS = CreateQuantizedImage();
        if(numQTZLEVELS != inputQTZLEVELS)
        {
            printf("ERROR: input quantized levels = %d, but actually used %d ! ", inputQTZLEVELS, numQTZLEVELS);
            printf("Exiting: please reduce the number of threads!\n");
            //printf("Sacrifice load balance and maintain the input number of quantization levels using e.g. \"CreateQuantizedImage()\".\n");
            exit(0);
        }
        else {printf("Number of threads now equal to %d.\n", numQTZLEVELS);}
        return;
    }
    // if everything goes good:
    CreateQuantizedImageInParallel(numQTZLEVELS);
    return;
}
