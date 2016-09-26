#!/bin/bash
#./compile.sh
g++ -o main -L$HOME/lib -I$HOME/lib  main.o quanttree.o handleimages.o quantizedimage.o refinetree.o radixsort.o filter.o -lpthread -lcfitsio -lfreeimage
rm *.o
