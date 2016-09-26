#ifndef _FILTER_H
#define _FILTER_H

#include "common.h"

typedef struct
{
    int self;
} ThreadFiltData;

void RefTreeAreaFilterBerger(double lambda, greyval_t *out);
void Filter();

#endif
