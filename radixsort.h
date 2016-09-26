#ifndef _RADIXSORT_H
#define _RADIXSORT_H

typedef struct
{
    int self;
    pixel_t *hist;
} ThreadSortingDataRS;

typedef struct
{
    int self;
} ThreadFlipData;


ThreadSortingDataRS *thsortingdataRS;
pixel_t *histogram;
int numStepsRS;
pixel_t *HISTOGRAMRS[2];

void RunRadixSortUnsigned();
void RunRSCountingSort(ThreadSortingDataRS *allThdata, int nthreads);
void *csortRS(void *arg);
void LocalHistRS(int self, pixel_t *SORTED, pixel_t *hist, int step);
ThreadSortingDataRS *MakeThreadRadixSortData(int numthreads);
void FreeThreadSortingDataRS(ThreadSortingDataRS *data, int numthreads);
void CreateSortedArrayRS(int self, pixel_t *sortedNew, pixel_t *sortedOld, pixel_t *hist, int step);

#endif
