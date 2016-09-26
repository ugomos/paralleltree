/**        refinetree.c
 *  The code implements the refinement phase of the algorithm. The effect is 
 *  that the array node_ref[] contains the refined "hybrid" max-tree of the original image.
 * 
 *  It is based on the Berger et al. algorithm ("Effective Component Tree Computation with Application to
		  Pattern Recognition in Astronomical Imaging", Christophe Berger and Thierry Geraud and Roland
		  Levillain and Nicolas Widynski and Anthony Baillard and Emmanuel Bertin
 * 
 */

#include "common.h"
#include "stdlib.h"
#include "refinetree.h"
#include "quanttree.h"
#include "handleimages.h"

void TreeAndCollectionsCreate(pixel_t size)
{
    zpar = malloc(size*sizeof(pixel_t));

    pixel_t i;
    node_ref = calloc(size, sizeof(Node));

    if (node_ref)
    {
        for (i=0; i<size; i++)
        {
            node_ref[i].Area = 1;
            node_ref[i].parent = -1;
            zpar[i] = -1;
        }
    }
    else
    {
        printf("ERROR in TreeCollectionsCreate.\n");
        exit(0);
    }
    return;
}

ThreadRefData *MakeThreadRefData(int numthreads)
{
    ThreadRefData *data = malloc(numthreads *sizeof(ThreadRefData));
    int i;

    for (i=0; i<numthreads; i++)
    {
        data[i].self=i;
        //printf("Thread %d from pxStartPos=%ld to pxEndPos=%ld \n", i, pxStartPosition[i], pxEndPosition[i]);
    }
    return(data);
}


void FreeThreadRefData(ThreadRefData *data, int numthreads)
{
    free(data);
}

void RunRefinementThreads(ThreadRefData *thdata, int nthreadsRef)
{
    int thread;
    for (thread=0; thread<nthreadsRef; ++thread)
    {
        pthread_create(threadID+thread, NULL, rnc, (void *) (thdata + thread));
    }
    for (thread=0; thread<nthreadsRef; ++thread)
    {
        pthread_join(threadID[thread], NULL);
    }
}

// it works for both 3D and 2D if width = size2D
int GetNeighborsBerger(pixel_t p, pixel_t *neighbors, pixel_t lwb, pixel_t upb)
{
    pixel_t x, y, z;
    int n=0;
    x = p % width;
    y = (p % size2D) / width;
    z = p / size2D;

    if (x<(width-1))       neighbors[n++] = p+1;
    if (y>0)               neighbors[n++] = p-width; // never exec in 2D
    if (x>0)               neighbors[n++] = p-1;
    if (y<(height-1))      neighbors[n++] = p+width; // never exec in 2D
    if (depth>1)
    {
        if (z>0)            neighbors[n++] = p-size2D;
        if (z<(depth-1))    neighbors[n++] = p+size2D;
    }
    return(n);
}

pixel_t DescendRoots(pixel_t q, int myLev)
{
    pixel_t curr = q;
    //while((node_qu[curr].parent != bottom) && (gval_qu[ node_qu[curr].parent ] > myLev))
    while(gval_qu[ node_qu[curr].parent ] > myLev)
    {
        curr = node_qu[curr].parent;
    }
    return curr;
}

// path compression without recursion
/*pixel_t FINDROOT(pixel_t p)
{
	pixel_t first = p;
	while(p != zpar[p])
	{
		p = zpar[p];
	}

	pixel_t nextEl;
	pixel_t next = first;
	while(next != zpar[p])
	{
		nextEl = zpar[next];
		zpar[next] = zpar[p];
		next = nextEl;
	}

	return zpar[p];
}*/


pixel_t FINDROOT(pixel_t p)
{
    if( (zpar[p] != p) )
        zpar[p] = FINDROOT(zpar[p]);

    return zpar[p];
}

void *rnc(void *arg)
{
    ThreadRefData *threfdata = (ThreadRefData *) arg;
    int self = threfdata->self;

    // analyse the pixels in decreasing order of intensity
    int numneighbors;
    pixel_t neighbors[CONNECTIVITY];
    pixel_t i, j, p, q, r, ancest_quFound, tobezipped;
    pixel_t lwb, upb;

    lwb = pxStartPosition[self];
    upb = pxEndPosition[self];
    //printf("Thread %d: lwb=%ld, upb=%ld\n", self, lwb, upb);

    // go through the sorted pixels of your partition from the highest to the lowest intensity
    for(i = upb; i>=lwb; i--) // i>=lwb is correct, otherwise a node could have a null parent if it is reachable only through the pixel of position lwb
    {
        p = SORTED[i];
        zpar[p] = p;

        numneighbors = GetNeighborsBerger(p, neighbors, lwb, upb);

        for (j=0; j<numneighbors; j++)
        {
            q = neighbors[j];

            if(gval_qu[q] > self)
            {
                // descend level roots in node_qu[] and stop when a level root pointing to the quantization level of my thread is found
                ancest_quFound = DescendRoots(q, self);

                if (node_ref[ancest_quFound].parent == bottom)
                {
                    node_ref[ancest_quFound].parent = p;
                    node_ref[p].Area += node_qu[ancest_quFound].Area;
                }
                else //perform zipping phase
                {
                    tobezipped = FINDROOT(node_ref[ancest_quFound].parent);

                    if(tobezipped != p)
                    {
                        node_ref[tobezipped].parent = p;
                        node_ref[p].Area += node_ref[tobezipped].Area;
                        zpar[tobezipped] = p;
                    }
                }
            }
            else if(gval_qu[q] == self)
            {
                // this is not needed if you check when filtering that the root node as bottom parent pointer.
                // if( ((node_qu[p].parent == bottom) )) // || (gval_qu[p] == gval_qu[node_qu[p].parent]) ))
                if( ((node_qu[p].parent == bottom) /*))*/ || (gval_qu[p] == gval_qu[node_qu[p].parent]) ))
                {
                    node_ref[p].parent = p;
                }

                if(zpar[q] != -1)
                    //if( (gval[q] > gval[p]) || ( (gval[q] == gval[p]) && (q > p)) )
                {
                    r = FINDROOT(q);
                    if(r != p)
                    {
                        node_ref[r].parent = p;
                        zpar[r] = p;

                        node_ref[p].Area += node_ref[r].Area;
                    }
                }
            }
        }
    }

    //Barrier(self, nthreadsRef);
    return NULL;
}


void RefineTreeBerger(int nthreads)
{
    start_ref = times(&tstruct);
    threfdata = MakeThreadRefData(nthreads);
    TreeAndCollectionsCreate(size);
    RunRefinementThreads(threfdata, nthreads);
    end_ref = times(&tstruct);
}
