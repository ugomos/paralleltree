/** filter.c
 * The code implements connected area opening,
 * filtering the final refined hybrid tree by removing the connected components (node)
 * with are smaller than lambda.
 * 
 * 
 * 
 * 
 */

#include "common.h"
#include "stdlib.h"
#include "filter.h"

void RefTreeAreaFilterBerger(double lambda, greyval_t *out)
{
    pixel_t v, u, w, parent, lwb, upb;
    greyval_t val;

    lwb = 0;
    upb = size;
    for (v=lwb; v<upb; v++)
    {
        if (node_ref[v].isFiltered == false)
        {
            w = v;
            parent = node_ref[w].parent;

            while ((parent != w) && (node_ref[w].isFiltered == false) &&
                    //((gval[w] == gval[parent]) || (node_ref[w].Area <= lambda)))
                    ((gval[w] == gval[parent]) || (node_ref[w].Area <= lambda) ))
            {
                w = parent;
                parent = node_ref[w].parent;
            }
            if (node_ref[w].isFiltered == true)
            {
                val = node_ref[w].filter;
            }
            else if (node_ref[w].Area > lambda)
                //else if (node_ref[w].attributes->area > lambda)
            {
                val = gval[w];
            }
            else
            {
                val = 0; // criterion cannot be satisfied
            }

            u = v;
            while (u!=w)
            {
                if ((lwb<=u) && (u<upb))
                {
                    node_ref[u].filter = val;
                    node_ref[u].isFiltered = true;
                }
                u = node_ref[u].parent;
            }
            if ((lwb<=w) && (w<upb))
            {
                node_ref[w].filter = val;
                node_ref[w].isFiltered = true;
            }
        }

        out[v] = node_ref[v].filter;
    }
}


void *runfilt(void *arg)
{
    ThreadFiltData *thfiltdata = (ThreadFiltData *) arg;
    int self = thfiltdata->self;
    greyval_t *out = outRef;

    pixel_t v, u, w, parent, lwb, upb;
    greyval_t val;
    //lwb = 0;
    //upb = size;
    lwb = LWB(self, nthreads);
    upb = UPB(self, nthreads);

    printf("%d) lwb=%ld; upb=%ld.\n", self, lwb, upb);

    for (v=lwb; v<upb; v++)
    {
        if (node_ref[v].isFiltered == false)
        {
            w = v;
            parent = node_ref[w].parent;

            while ((parent != w) && (node_ref[w].isFiltered == false) &&
                    //((gval[w] == gval[parent]) || (node_ref[w].Area <= lambda)))
                    ((gval[w] == gval[parent]) || (node_ref[w].Area <= lambda) ))
            {
                w = parent;
                parent = node_ref[w].parent;
            }
            if (node_ref[w].isFiltered == true)
            {
                val = node_ref[w].filter;
            }
            //else if (node_ref[w].Area > lambda)
            else if (node_ref[w].Area > lambda)
            {
                val = gval[w];
            }
            else
            {
                val = 0; // criterion cannot be satisfied
            }

            u = v;
            while (u!=w)
            {
                if ((lwb<=u) && (u<upb))
                {
                    node_ref[u].filter = val;
                    node_ref[u].isFiltered = true;
                }
                u = node_ref[u].parent;
            }
            if ((lwb<=w) && (w<upb))
            {
                node_ref[w].filter = val;
                node_ref[w].isFiltered = true;
            }
        }

        out[v] = node_ref[v].filter;
    }
    return NULL;
}

void FreeThreadFiltData(ThreadFiltData *data, int numthreads)
{
    free(data);
}

ThreadFiltData *MakeThreadFiltData(int numthreads)
{
    ThreadFiltData *data = malloc(numthreads *sizeof(ThreadFiltData));
    int i;

    for (i=0; i<numthreads; i++)
    {
        data[i].self=i;
    }
    return(data);
}

void RunFilter(ThreadFiltData *thdata, int nthreads)
{
    int thread;
    for (thread=0; thread<nthreads; ++thread)
    {
        pthread_create(threadID+thread, NULL, runfilt, (void *) (thdata + thread));
    }
    for (thread=0; thread<nthreads; ++thread)
    {
        pthread_join(threadID[thread], NULL);
    }
}

void ParallelFilter(int nthreads)
{
    ThreadFiltData *thfiltdata = MakeThreadFiltData(nthreadsRef);
    RunFilter(thfiltdata, nthreadsRef);
    FreeThreadFiltData(thfiltdata, nthreadsRef);
}

void Filter()
{
    //RefTreeAreaFilterBerger(lambda, outRef);
    ParallelFilter(nthreadsRef);
}
