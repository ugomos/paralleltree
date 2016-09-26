/***         quanttree.c
 * The code implements the building of the pilot max-tree.
 * It is based on the algorithm presented in:
 * "Concurrent computation of attribute filters on shared memory parallel machines",
 * by Wilkinson, M. H. F., Gao, H., Hesselink, W. H., Jonker, J-E. & Meijster, A.
 * Oct-2008 In : Ieee transactions on pattern analysis and machine intelligence. 30, 10, p. 1800-1813 14 p.
 * 
 */ 

#include "common.h"
#include "quanttree.h"
#include "handleimages.h"

pthread_mutex_t samut[MAXTHREADS];
pthread_cond_t  sacv[MAXTHREADS];
int             saval[MAXTHREADS];
pthread_mutex_t barriermutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t barriercv = PTHREAD_COND_INITIALIZER;

/************************** safe malloc and calloc ************************************/

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

void *SafeMalloc(int n)
{
    void *ptr;
    pthread_mutex_lock(&mutex);
    ptr = malloc(n);
    if (ptr==NULL)
    {
        fprintf (stderr, "Error: out of memory. "
                 "Could not allocate %d bytes\n", n);
    }
    pthread_mutex_unlock(&mutex);
    return ptr;
}

void *SafeCalloc(int nmemb, pixel_t size)
{
    void *ptr;
    pthread_mutex_lock(&mutex);
    ptr = calloc(nmemb, size);
    if (ptr==NULL)
    {
        fprintf (stderr, "Error: out of memory. Could not allocate %ld bytes\n", nmemb*size);
    }
    pthread_mutex_unlock(&mutex);
    return ptr;
}

/************************** semaphore ************************************/

void Psa(int p)
{
    pthread_mutex_lock(&samut[p]);
    while (saval[p] <= 0)
        pthread_cond_wait(&sacv[p], &samut[p]);
    saval[p]--;
    pthread_mutex_unlock(&samut[p]);
}

void Vsa(int p)
{
    pthread_mutex_lock(&samut[p]);
    saval[p]++;
    pthread_mutex_unlock(&samut[p]);
    pthread_cond_broadcast(&sacv[p]);
}

/************************** barrier ************************************/
int barcnt = 0;

void Barrier(int self, int nthreads)
{
    pthread_mutex_lock(&barriermutex);
    barcnt++;
    if (barcnt == nthreads)
    {
        barcnt = 0;  /* for reuse of routine */
        pthread_cond_broadcast (&barriercv);
    }
    else
    {
        pthread_cond_wait (&barriercv, &barriermutex);
    }
    pthread_mutex_unlock(&barriermutex);
}

/************************** level root ************************************/

pixel_t levroot(pixel_t x, int *gvalues)
{
    pixel_t r=x, y ;
    greyval_t gv=gvalues[x];
    if (r==bottom) return bottom;
    while ((node_qu[r].parent!=bottom) && (gv==gvalues[node_qu[r].parent]))
        r = node_qu[r].parent;

    while (x!=r)
    {
        y=node_qu[x].parent;
        node_qu[x].parent=r;
        x=y;
    }
    return r;
}

pixel_t Par(pixel_t x, int *gvalues)
{
    return (pixel_t) levroot(node_qu[x].parent, gvalues);
}

void levrootfix(pixel_t lwb, pixel_t upb, int *gvalues)
{
    pixel_t z, u, x;

    for (x=lwb; x<upb; x++)
    {
        u = levroot(x, gvalues);
        if (x!=u) node_qu[x].parent=u;
        else
        {
            z = Par(x, gvalues);
            node_qu[x].parent=z;
        }
    }
}

/************************** queue ************************************/
Queue *QueueCreate(pixel_t imgsize)
{
    Queue *hq;

    hq = calloc(numQTZLEVELS, sizeof(Queue));
    if (hq==NULL)
    {
        fprintf (stderr, "Out of memory!");
        return(NULL);
    }
    hq->Pixels = calloc(imgsize, sizeof(pixel_t));
    if (hq->Pixels==NULL)
    {
        free(hq);
        fprintf (stderr, "Out of memory!");
        return(NULL);
    }
    return(hq);
} /* QueueCreate */

void FillPixelsInQueue(Queue *hq, pixel_t *numpixelsperlevel)
{
    int i;
    hq->Head = hq->Tail = 0;

    for (i=1; i<numQTZLEVELS; i++)
    {
        hq[i].Pixels = hq[i-1].Pixels + numpixelsperlevel[i-1];
        hq[i].Head = hq[i].Tail = 0;
    }
} /* FillPixelsInQueue */

void QueueDelete(Queue *hq)
{
    free(hq->Pixels);
    free(hq);
} /* QueueDelete */

/****************** Construction of local tree using algorithm flood****************************/
int GetNeighbors(pixel_t p, pixel_t x, pixel_t y, pixel_t z,
                 pixel_t *neighbors, pixel_t lwb, pixel_t upb)
{
    int n=0;
    x = p % width;
    if (x<(width-1))       neighbors[n++] = p+1;
    if (y>0)               neighbors[n++] = p-width; // never exec in 2D
    if (x>0)               neighbors[n++] = p-1;
    if (y<(height-1))      neighbors[n++] = p+width; // never exec in 2D
    if (depth>1)
    {
        if (z>lwb)          neighbors[n++] = p-size2D;
        if (z<(upb-1))           neighbors[n++] = p+size2D;
    }
    return(n);
}

/* parameter MaxNode *nodes must be the par array of the quantized image: 'node_qu';
 * parameter greyval_t *gvalues must be the quantized image: 'gval_qu' */
int LocalTreeFlood(int self,  Queue *set, pixel_t *lero, int lev, long *thisarea, MaxNode *nodes, int *gvalues)
{
    pixel_t neighbors[CONNECTIVITY];

    long area = *thisarea;
    long childarea;

    pixel_t p, q, x, y, z, lwb, upb;
    int fq;
    int numneighbors, i, m;

    lwb = LWB(self, nthreads);
    upb = UPB(self, nthreads);

    while(QueueNotEmpty(set, lev))
    {
        area++;
        p = QueueFirst(set, lev);

        x = p % width;
        y = (p % size2D) / width;
        z = p / size2D;
        numneighbors = GetNeighbors(p, x, y, z, neighbors, lwb/size2D, upb/size2D);
        for (i=0; i<numneighbors; i++)
        {
            q = neighbors[i];
            if (!reached_qu[q])
            {
                fq = gvalues[q];
                reached_qu[q] = true;

                if (lero[fq]==bottom)
                {
                    lero[fq] = q;
                }
                else
                {
                    if( ( (gval[q] < gval[lero[fq]]) || ( (gval[q] == gval[lero[fq]]) && (q < lero[fq]) )) )
                    {
                        nodes[lero[fq]].parent = q; //new lero
                        lero[fq] = q; //new lero
                    }
                    nodes[q].parent = lero[fq];
                }

                QueueAdd(set, fq, q);

                if (fq > lev)
                {
                    childarea = 0;
                    do
                    {
                        fq = LocalTreeFlood(self, set, lero, fq, &childarea, nodes, gvalues);
                        if(fq>=numQTZLEVELS)
                        {
                            return(fq);
                        }
                    }
                    while (fq!=lev);
                    area += childarea;
                }
            }
        }
    }
    m = lev-1;
    while ((m>=0) && (lero[m]==bottom))
    {
        m--;
    }

    if (m>=0) nodes[lero[lev]].parent = lero[m];
    nodes[lero[lev]].Area = area;
    lero[lev] = bottom;
    *thisarea = area;
    return(m);
} /* MaxTreeFlood */

/************************ Merge trees *************************************/

void Connect(pixel_t x, pixel_t y, MaxNode *node, int *gvalues)
{
    long area = 0, area1 = 0;

    pixel_t h, z;

    x = levroot(x, gvalues);
    y = levroot(y, gvalues);
    //if (gvalues[x] < gvalues[y]) // changed to keep the levelroot of each component corresponding to the pixel with the minimum intensity value in the original image.
    if ( gval[x] < gval[y] || ((gval[x] == gval[y]) && (x < y)) )
    {
        h=x;
        x=y;
        y=h;
    }
    while ((x!=y) && (y!=bottom))
    {
        z = Par(x, gvalues);
        //if ((z!=bottom) && (gvalues[z]>=gvalues[y]))
        //if ( (z!=bottom) && ((gval[z]>=gval[y])))
        if ( (z!=bottom) && ( (gval[z]>gval[y]) || ((gval[z]==gval[y]) && (z>=y)) ) ) //z>=y
        {
            node[x].Area += area;
            x = z;
        }
        else
        {
            area1 = node[x].Area + area;
            area = node[x].Area;
            node[x].Area = area1;

            node[x].parent = y;
            x = y;
            y = z;
        }
    }
    if (y==bottom)
    {
        while(x!=bottom)
        {
            node[x].Area += area;
            x = Par(x, gvalues);
        }
    }
}


void Fuse(int self, int i, MaxNode *node, int *gvalues) /* fuse regions [LWB(self), UPB(self+i-1)) and  [LWB(self+i), UPB(self+2i-1)) vertically */
{
    pixel_t p, q, x, y;

    greyval_t prevmin, curmin;

    /* get the horizontal boundary */

    p  = LWB(self+i, nthreads);
    q = p - size2D;

    x = p % width;
    y = (p /width) % height;

    /*  printf("Region %d merger with %d: (%d,%d,%d)\n",self, self+i, x,y,z);*/

    for ( y = 0 ; y < height ; y++ )
    {
        Connect(p, q, node, gvalues);
        //prevmin = MIN(gvalues[p], gvalues[q]);
        prevmin = MIN(gval[p], gval[q]);

        p++;
        q++;
        for ( x = 1 ; x < width ; x++, p++, q++ )
        {
            //curmin = MIN(gvalues[p], gvalues[q]); // changed from GVALUES to GVAL to keep the levelroot of each component corresponding to the pixel with the minimum intensity value in the original image.
            curmin = MIN(gval[p], gval[q]);
            if (curmin > prevmin)
            {
                Connect(p, q, node, gvalues);
            }
            prevmin = curmin;
        }
    }
}

ThreadData *MakeThreadData(int numthreads)
{
    ThreadData *data = malloc(numthreads *sizeof(ThreadData));
    int i;

    for (i=0; i<numthreads; i++)
    {
        data[i].self=i;
        data[i].thisqueue= QueueCreate( UPB(i, numthreads)-LWB(i, numthreads));
    }
    return(data);
}


void FreeThreadData(ThreadData *data, int numthreads)
{
    int i;
    for (i=0; i<numthreads; i++)
        QueueDelete(data[i].thisqueue);
    free(data);
}


void *ccaf(void *arg)
{
    ThreadData *thdata = (ThreadData *) arg;
    int self = thdata->self, q, i;
    pixel_t x;
    long area=0;
    pixel_t numpixelsperlevel[numQTZLEVELS], lero[numQTZLEVELS], xm;
    Queue *set = thdata->thisqueue;

    for (i=0; i<numQTZLEVELS; i++)
    {
        numpixelsperlevel[i] = 0;
        lero[i] = bottom;
    }

    xm = LWB(self, nthreads);

    for (x=xm; x<UPB(self, nthreads); x++)
    {
        numpixelsperlevel[gval_qu[x]]++;
        node_qu[x].parent = bottom;
        //if (gval_qu[xm]>gval_qu[x]) xm = x;
        if (gval[xm]>gval[x]) xm = x;
    }
    FillPixelsInQueue(set, numpixelsperlevel);

    QueueAdd(set, gval_qu[xm], xm);
    reached_qu[xm] = true;

    // the level root for gval_qu[] for this thread is equals to
    // the minimum original value for the quantized component.
    lero[gval_qu[xm]] = xm;

    LocalTreeFlood(self, set, lero, gval_qu[xm], &area, node_qu, gval_qu);

    i = 1;
    q = self;
    while ((self+i<nthreads) && (q%2 == 0))
    {
        Psa(self+i);  /* wait to glue with righthand neighbor */
        Fuse(self, i, node_qu, gval_qu);
        i = 2*i;
        q = q/2;
    }
    if (self != 0)
    {
        Vsa(self);  /* signal lefthand neighbor */
    }
    return NULL;
}

void BuildQuantizedTree(ThreadData *thdata, int nthreads)
{
    int thread;
    for (thread=0; thread<nthreads; ++thread)
    {
        pthread_create(threadID+thread, NULL, ccaf, (void *) (thdata + thread));
    }
    for (thread=0; thread<nthreads; ++thread)
    {
        pthread_join(threadID[thread], NULL);
    }
}

int BuildMaxTreeOfQuantizedImage()
{
    // Build the max tree of the quantized image
    node_qu = calloc((size_t)size, sizeof(MaxNode));
    if (node_qu==NULL)
    {
        fprintf(stderr, "out of memory! \n");
        free(gval);
        free(gval_qu);
        return(-1);
    }

    reached_qu = calloc((size_t)size, sizeof(bool));
    if (reached_qu==NULL)
    {
        fprintf(stderr, "out of memory!\n");
        free(node_qu);
        free(gval_qu);
        return(-1);
    }

    thdata = MakeThreadData(nthreads);
    BuildQuantizedTree(thdata,nthreads);

    FreeThreadData(thdata, nthreads);
    free(reached_qu);
    return 0;
}
