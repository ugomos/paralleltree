#include "common.h"
#include "handleimages.h"

/****** Image read/write functions ******************************/
/** Generic image loader
@param lpszPathName Pointer to the full file name
@param flag Optional load flag constant
@return Returns the loaded dib if successful, returns NULL otherwise
*/
FIBITMAP* GenericLoader(const char* lpszPathName, int flag)
{
    FREE_IMAGE_FORMAT fif = FIF_UNKNOWN;
//  check the file signature and deduce its format
// (the second argument is currently not used by FreeImage)
    fif = FreeImage_GetFileType(lpszPathName, 0);
    if(fif == FIF_UNKNOWN)
    {
        // no signature ?
        // try to guess the file format from the file extension
        fif = FreeImage_GetFIFFromFilename(lpszPathName);
    }
    // check that the plugin has reading capabilities ...
    if((fif != FIF_UNKNOWN) && FreeImage_FIFSupportsReading(fif))
    {
        // ok, let's load the file
        FIBITMAP *dib = FreeImage_Load(fif, lpszPathName, flag);
        // unless a bad file format, we are done !
        return dib;
    }
    return NULL;
}


void ReadFITS(char *filename)
{
    fitsfile *fptr;   /* FITS file pointer, defined in fitsio.h */
    int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
    int bitpix, naxis;
    long naxes[2] = {1,1}, fpixel[2] = {1,1};
    greyval_t *CURRENT;

    if (!fits_open_file(&fptr, filename, READONLY, &status))
    {
        if (!fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status) )
        {
            if (naxis > 2 || naxis == 0)
                printf("Error: only 1D or 2D images are supported\n");
            else
            {
                width = naxes[0];
                depth = naxes[1];
                height = 1;
                size = naxes[0]*naxes[1];

                // change to greyval_t format
                gval = (greyval_t *) malloc(naxes[0]*naxes[1] * sizeof(greyval_t));
                CURRENT = gval;

                if (CURRENT == NULL)
                {
                    printf("Memory allocation error\n");
                    return;
                }

                /* TUSHORT here is correct since greyval_t has type 'unsigned short', otherwise it must be changed */
                fits_read_pix(fptr, TUSHORT, fpixel, naxes[0]*naxes[1], NULL, CURRENT, NULL, &status);
                //fits_read_pix(fptr, Tlong, fpixel, naxes[0]*naxes[1], NULL, CURRENT, NULL, &status);

                /*int ii; for (ii = 0; ii < naxes[0]*naxes[1]; ii++)
                   printf("%d ", CURRENT[ii]);
                printf("\n");
                */
            }
        }
        fits_close_file(fptr, &status);
    }
    printf("\nHere. Statistics of %ld x %ld  image. Bits per pixel=%d.\n", naxes[0], naxes[1], bitpix);

    /**
    bitpix variable can assume:
    BYTE_IMG   =   8   ( 8-bit byte pixels, 0 - 255)
    SHORT_IMG  =  16   (16 bit integer pixels)
    LONG_IMG   =  32   (32-bit integer pixels)
    FLOAT_IMG  = -32   (32-bit floating point pixels)
    DOUBLE_IMG = -64   (64-bit floating point pixels)
    **/
    return;
}

greyval_t FindMax(greyval_t *gval)
{
    greyval_t max = gval[0];
    pixel_t i;
    for (i = 1; i<size; i++)
    {
        max = (gval[i]>max) ? gval[i] : max;
    }
    return max;
}

greyval_t FindMin(greyval_t *gvalues)
{
    greyval_t min = gvalues[0];
    pixel_t i;

    for (i = 1; i<size; i++)
    {
        min = (gval[i]<min) ? gval[i] : min;
    }
    return min;
}

short ReadTIFF(char *fnm)
{
    FIBITMAP *dib = GenericLoader(fnm,0);
    unsigned long  bitsperpixel;
    //greyval_t *im;
    //unsigned int x,y,i;
    pixel_t x,y,i;
    if (dib == NULL) return 0;

    bitsperpixel =  FreeImage_GetBPP(dib);
    depth = FreeImage_GetHeight(dib),
    width = FreeImage_GetWidth(dib);
    size = width*depth;
    //printf("size = %ld\n", size);
    height = 1;
    gval = malloc(size*sizeof(greyval_t));
    if (gval==NULL)
    {
        fprintf (stderr, "Out of memory!");
        return(0);
    }
    switch(bitsperpixel)
    {
    case 8:
        i=0;
        for(y = 0; y <depth; y++)
        {
            BYTE *bits = (BYTE *)FreeImage_GetScanLine(dib, depth - y -1);
            //printf("y=%d\n",y);
            for(x = 0; x < width; x++,i++)
            {
                gval[i] = bits[x];
            }
        }

        FreeImage_Unload(dib);
        return 1;
    case 16:
        i=0;
        for(y = 0; y < depth; y++)
        {
            unsigned short *bits = (unsigned short *)FreeImage_GetScanLine(dib,depth - y -1);
            for(x = 0; x < width; x++,i++)
            {
                gval[i] = bits[x];
            }
        }
        FreeImage_Unload(dib);
        return 1;
    default :
        FreeImage_Unload(dib);
        fprintf(stderr, "unsupported format\n");
        exit(-1);
        return 0;
    }
}

int WriteFITS(char *filename, char *inputimagefilename, greyval_t *out)
{
    fitsfile *outfptr, *infptr;  /* FITS file pointers */
    int status = 0;  /* CFITSIO status value MUST be initialized to zero! */

    char str1[100];
    strcpy(str1, "!"); // '!' symbol makes the output image, if existing, to be overwritten

    // create the new empty output file
    if (!fits_create_file(&outfptr, strcat(str1, filename), &status) )
    {
        // copy all the header keywords from original image to new output file
        fits_open_file(&infptr, inputimagefilename, READONLY, &status); // open input images
        fits_copy_header(infptr, outfptr, &status);
        fits_close_file(infptr, &status);

        long fpixel[2];
        fpixel[0] = fpixel[1] = 1; // start to copy from the first element in the array gval (row=1, col 1)
        fits_write_pix(outfptr, TUSHORT, fpixel, width*depth, out, &status); // write new values to output image
        //fits_write_pix(outfptr, TUSHORT, fpixel, 2048*1489, gval, &status); // write new values to output image

        printf("File should be written.\n");
        /**
        bitpix variable can assume:
        BYTE_IMG   =   8   ( 8-bit byte pixels, 0 - 255)
        SHORT_IMG  =  16   (16 bit integer pixels)
        LONG_IMG   =  32   (32-bit integer pixels)
        FLOAT_IMG  = -32   (32-bit floating point pixels)
        DOUBLE_IMG = -64   (64-bit floating point pixels)
        **/
    }
    fits_close_file(outfptr, &status);
    return 0;
} /* WriteFITS */


/* WriteTIFF_2 */
void WriteTIFF_2(char *fname, greyval_t *img, pixel_t width, pixel_t height, int bitspp)
{
    FIBITMAP *outmap;
    pixel_t i, y, x;
    //FREE_IMAGE_FORMAT fif = FreeImage_GetFIFFromFilename(fname);

    if (bitspp == 8)
    {
        ubyte *imagebuf;
        RGBQUAD *pal;
        outmap = FreeImage_AllocateT(FIT_BITMAP,width,depth,bitspp,0xFF,0xFF,0xFF);
        pal = FreeImage_GetPalette(outmap);
        for (i = 0; i < 256; i++)
        {
            pal[i].rgbRed = i;
            pal[i].rgbGreen = i;
            pal[i].rgbBlue = i;
        }
        i = 0;
        //for (y=0; y<depth; y++)
        for (y=depth-1; y>=0; y--)
        {
            imagebuf = FreeImage_GetScanLine(outmap,y);
            for (x=0; x<width; x++,i++)
            {
                //imagebuf[x]=out[i];
                imagebuf[x]=img[i];
            }
        }
    }
    else
    {
        unsigned short *imagebuf;
        outmap = FreeImage_AllocateT(FIT_UINT16,width,depth,16,0xFFFF,0xFFFF,0xFFFF);
        i = 0;
        //for (y=0; y<depth; y++)
        for (y=depth-1; y>=0; y--)
        {
            imagebuf = (unsigned short *)FreeImage_GetScanLine(outmap,y);
            for (x=0; x<width; x++,i++)
                imagebuf[x]=img[i];

        }
    }
    //FreeImage_Save(fif,outmap,fname,0);
    FreeImage_Save(FIF_TIFF,outmap,fname, TIFF_NONE);
    FreeImage_Unload(outmap);
}


int WriteTIFF(char *fname, greyval_t *img)
{
    if (FindMax(img)>255)
        WriteTIFF_2(fname, img, width, depth, 16);
    else WriteTIFF_2(fname, img, width, depth, 8);

    return 0;
}


int ReadFITS3D(char *filename)
{
    fitsfile *fptr;   /* FITS file pointer, defined in fitsio.h */
    int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
    int bitpix, naxis;
    long naxes[3] = {1,1,1}, fpixel[3] = {1,1,1};
    greyval_t *CURRENT;

    if (!fits_open_file(&fptr, filename, READONLY, &status))
    {
        if (!fits_get_img_param(fptr, 3, &bitpix, &naxis, naxes, &status) )
        {
            if (naxis > 4 || naxis == 0)
            {
                printf("Error: only 1D - 2D - 3D images are supported\n");
                return 0;
            }
            else
            {
                width = naxes[0];
                height = naxes[1];
                depth = naxes[2];
                size2D = naxes[0]*naxes[1];
                size = naxes[0]*naxes[1]*naxes[2];

                // change to greyval_t format
                gval = (greyval_t *) malloc(naxes[0]*naxes[1]*naxes[2] * sizeof(greyval_t));
                CURRENT = gval;

                if (CURRENT == NULL)
                {
                    printf("Memory allocation error.\n");
                    return 0;
                }

                fits_read_pix(fptr, TFLOAT, fpixel, naxes[0]*naxes[1]*naxes[2], NULL, CURRENT, NULL, &status);

                /*long ii;
                for (ii = 0; ii < size; ii++)
                   printf("%ld) %f \n", ii, gval[ii]);
                 */

            }
        }
        fits_close_file(fptr, &status);
    }
    printf("\nStatistics of %ld x %ld x %ld image. Number of axis = %d. Bits per pixel=%d.\n", naxes[0], naxes[1], naxes[2], naxis, bitpix);

    /**
    bitpix variable can assume:
    BYTE_IMG   =   8   ( 8-bit byte pixels, 0 - 255)
    SHORT_IMG  =  16   (16 bit integer pixels)
    LONG_IMG   =  32   (32-bit integer pixels)
    FLOAT_IMG  = -32   (32-bit floating point pixels)
    DOUBLE_IMG = -64   (64-bit floating point pixels)
    **/
    return 1;
}

// it reads the datacube provided by Kapteyn
int ReadFITS3D_CUBE(char *filename)
{
    fitsfile *fptr;   /* FITS file pointer, defined in fitsio.h */
    int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
    int bitpix, naxis;
    long naxes[4] = {1,1,1}, fpixel[4] = {1,1,1,1};
    greyval_t *CURRENT;

    if (!fits_open_file(&fptr, filename, READONLY, &status))
    {
        if (!fits_get_img_param(fptr, 4, &bitpix, &naxis, naxes, &status) )
        {
            if (naxis > 4 || naxis == 0)
            {
                printf("Error: only 1D - 2D - 3D images are supported\n");
                return 0;
            }
            else
            {
                printf("\nParameters: %ld x %ld x %ld x %ld image. Number of axis = %d. Bits per pixel=%d.\n", naxes[0], naxes[1], naxes[2], naxes[3], naxis, bitpix);
                width = naxes[0];
                height = naxes[1];
                depth = naxes[3];
                //depth = naxes[2];
                size2D = naxes[0]*naxes[1];
                size = naxes[0]*naxes[1]*naxes[3];
                //size = naxes[0]*naxes[1]*naxes[2];

                // change to greyval_t format
                gval = (greyval_t *) malloc(naxes[0]*naxes[1]*naxes[3] * sizeof(greyval_t));
                //gval = (greyval_t *) malloc(naxes[0]*naxes[1]*naxes[2] * sizeof(greyval_t));
                CURRENT = gval;

                if (CURRENT == NULL)
                {
                    printf("Memory allocation error.\n");
                    return 0;
                }

                fits_read_pix(fptr, TFLOAT, fpixel, naxes[0]*naxes[1]*naxes[3], NULL, CURRENT, NULL, &status);
                //fits_read_pix(fptr, TFLOAT, fpixel, naxes[0]*naxes[1]*naxes[2], NULL, CURRENT, NULL, &status);

                /*long ii;
                for (ii = 0; ii < size; ii++)
                   printf("%ld) %f \n", ii, gval[ii]);
                 */

            }
        }
        else
        {
            printf("sthing wrong\n");
        }
        fits_close_file(fptr, &status);
    }
    printf("\nStatistics of %ld x %ld x %ld x %ld image. Number of axis = %d. Bits per pixel=%d.\n", naxes[0], naxes[1], naxes[2], naxes[3], naxis, bitpix);

    /**
    bitpix variable can assume:
    BYTE_IMG   =   8   ( 8-bit byte pixels, 0 - 255)
    SHORT_IMG  =  16   (16 bit integer pixels)
    LONG_IMG   =  32   (32-bit integer pixels)
    FLOAT_IMG  = -32   (32-bit floating point pixels)
    DOUBLE_IMG = -64   (64-bit floating point pixels)
    **/
    return 1;
}


int WriteFITS3DCUBE(char *filename, char *inputimagefilename, greyval_t *out)
{
    fitsfile *outfptr, *infptr;  /* FITS file pointers */
    int status = 0;  /* CFITSIO status value MUST be initialized to zero! */

    char str1[100];
    strcpy(str1, "!"); // '!' symbol makes the output image, if existing, to be overwritten

    // create the new empty output file
    if (!fits_create_file(&outfptr, strcat(str1, filename), &status) )
    {
        // copy all the header keywords from original image to new output file
        fits_open_file(&infptr, inputimagefilename, READONLY, &status); // open input images
        fits_copy_header(infptr, outfptr, &status);
        fits_close_file(infptr, &status);

        long fpixel[3];
        fpixel[0] = fpixel[1] = fpixel[2] = 1; // start to copy from the first element in the array gval (row=1, col 1)
        fits_write_pix(outfptr, TFLOAT, fpixel, 1024*1024*1082, out, &status); // write new values to output image
        /**
        bitpix variable can assume:
        BYTE_IMG   =   8   ( 8-bit byte pixels, 0 - 255)
        SHORT_IMG  =  16   (16 bit integer pixels)
        LONG_IMG   =  32   (32-bit integer pixels)
        FLOAT_IMG  = -32   (32-bit floating point pixels)
        DOUBLE_IMG = -64   (64-bit floating point pixels)
        **/
    }
    fits_close_file(outfptr, &status);
    return 0;
}

int WriteFITS3D(char *filename, char *inputimagefilename, greyval_t *out)
{
    fitsfile *outfptr, *infptr;  /* FITS file pointers */
    int status = 0;  /* CFITSIO status value MUST be initialized to zero! */

    char str1[100];
    strcpy(str1, "!"); // '!' symbol makes the output image, if existing, to be overwritten

    // create the new empty output file
    if (!fits_create_file(&outfptr, strcat(str1, filename), &status) )
    {
        // copy all the header keywords from original image to new output file
        fits_open_file(&infptr, inputimagefilename, READONLY, &status); // open input images
        fits_copy_header(infptr, outfptr, &status);
        fits_close_file(infptr, &status);

        //int fits_write_imghdr / ffphps(fitsfile *fptr, int bitpix, int naxis, long *naxes, > int *status)
        //long naxes[3] = {200, 200, 4};
        //fits_write_imghdr(outfptr, -32, 3, naxes, &status);

        long fpixel[3];
        fpixel[0] = fpixel[1] = fpixel[2] = 1; // start to copy from the first element in the array gval (row=1, col 1)
        fits_write_pix(outfptr, TFLOAT, fpixel, width*depth*height, out, &status); // write new values to output image
        /**
        bitpix variable can assume:
        BYTE_IMG   =   8   ( 8-bit byte pixels, 0 - 255)
        SHORT_IMG  =  16   (16 bit integer pixels)
        LONG_IMG   =  32   (32-bit integer pixels)
        FLOAT_IMG  = -32   (32-bit floating point pixels)
        DOUBLE_IMG = -64   (64-bit floating point pixels)
        **/
    }
    fits_close_file(outfptr, &status);
    return 0;
}
