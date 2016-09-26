/* It builds a max tree of the quantized image + support functions for
 * threads: Barrier, semaphores, ...*/

#ifndef _MAXTREE_H
#define _MAXTREE_H

#define QueueFirst(hq,h)  (hq[h].Pixels[hq[h].Head++])
#define QueueAdd(hq,h,p)  hq[h].Pixels[hq[h].Tail++] = p
#define QueueNotEmpty(hq,h)  (hq[h].Head != hq[h].Tail)

typedef struct
{
    int self;
} ThreadAllocQntTree;

typedef struct Queue
{
    pixel_t *Pixels;
    pixel_t Head;
    pixel_t Tail; /* First free place in queue, empty if Head=Tail */
} Queue;

typedef struct
{
    int self;
    Queue *thisqueue;
} ThreadData;

ThreadData *thdata;

void *SafeMalloc(int n);
void *SafeCalloc(int nmemb, pixel_t size);
void Psa(int p);
void Vsa(int p);
void Barrier(int self, int nthreads);
pixel_t levroot(pixel_t x, int *gvalues);
pixel_t Par(pixel_t x, int *gvalues);
void levrootfix(pixel_t lwb, pixel_t upb, int *gvalues);

Queue *QueueCreate(pixel_t imgsize);
void FillPixelsInQueue(Queue *hq, pixel_t *numpixelsperlevel);
void QueueDelete(Queue *hq);

int GetNeighbors(pixel_t p, pixel_t x, pixel_t y, pixel_t z, pixel_t *neighbors, pixel_t lwb, pixel_t upb);
int LocalTreeFlood(int self,  Queue *set, pixel_t *lero, int lev, long *thisarea, MaxNode *nodes, int *gvalues);
void Connect(pixel_t x, pixel_t y, MaxNode *node, int *gvalues);
void Fuse(int self, int i, MaxNode *node, int *gvalues);
ThreadData *MakeThreadData(int numthreads);
void FreeThreadData(ThreadData *data, int numthreads);
void BuildQuantizedTree(ThreadData *thdata, int nthreads);
void *ccaf(void *arg);

int BuildMaxTreeOfQuantizedImage();

#endif
