#ifndef _REFINETREE_H
#define _REFINETREE_H

#include "common.h"

typedef struct
{
    int self;
} ThreadRefData;

void RefineTreeBerger(int nthreads);
void TreeAndCollectionsCreate(pixel_t size);
void RunRefinementThreads(ThreadRefData *thdata, int nthreads);

// Manage parallel
ThreadRefData *threfdata;
ThreadRefData *MakeThreadRefData(int numthreads);
void FreeThreadRefData(ThreadRefData *data, int numthreads);
void *rnc(void *arg);

//int GetNeighborsBerger(pixel_t p, pixel_t *neighbors);
int GetNeighborsBerger(pixel_t p, pixel_t *neighbors, pixel_t lwb, pixel_t upb);
void RefTreeAreaFilterBerger(double lambda, greyval_t *out);

pixel_t FINDROOT(pixel_t p);
#endif
