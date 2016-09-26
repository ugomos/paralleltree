#ifndef _QUANTIZEDIMAGE_H
#define _QUANTIZEDIMAGE_H

typedef struct
{
    int self;
} ThreadQntImgData;

int CreateQuantizedImage();
void CalculateQuantizedImage();
int WorkOutBoundaries();

void CreateQuantizedImageInParallel(int numQtzLevsUsed);
ThreadQntImgData *MakeThreadQntImgData(int numthreads);
void FreeThreadQntImgData(ThreadQntImgData *data, int numthreads);
void RunThreadsWritingQuantizedImage(ThreadQntImgData *thdata, int nthreads);
void *wqi(void *arg);

#endif
