/* pngnq.c - quantize the colors in an alphamap down to 256 using
**  the Neuquant algorithm.
**
** Based on Greg Roelf's pngquant which was itself based on Jef Poskanzer's ppmquant.
** Uses Anthony Dekker's Neuquant algorithm extended to handle the alpha channel.
** Rewritten by Kornel Lesiński (2009)

**
** Copyright (C) 1989, 1991 by Jef Poskanzer.
** Copyright (C) 1997, 2000, 2002 by Greg Roelofs; based on an idea by
**                                Stefan Schneider.
** Copyright (C) 2004-2009 by Stuart Coyle
** Copyright (C) Kornel Lesiński (2009)

** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appear in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty.
*/

/* NeuQuant Neural-Net Quantization Algorithm
 * ------------------------------------------
 *
 * Copyright (c) 1994 Anthony Dekker
 *
 * NEUQUANT Neural-Net quantization algorithm by Anthony Dekker, 1994.
 * See "Kohonen neural networks for optimal colour quantization"
 * in "Network: Computation in Neural Systems" Vol. 5 (1994) pp 351-367.
 * for a discussion of the algorithm.
 * See also  http://members.ozemail.com.au/~dekker/NEUQUANT.HTML
 *
 * Any party obtaining a copy of these files from the author, directly or
 * indirectly, is granted, free of charge, a full and unrestricted irrevocable,
 * world-wide, paid up, royalty-free, nonexclusive right and license to deal
 * in this software and documentation files (the "Software"), including without
 * limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons who receive
 * copies from any such party to do so, with the only requirement being
 * that this copyright notice remain intact.
 *
 */


#include "png.h"
#include "neuquant32.h"
#include "rwpng.h"
#include "errors.h"
#include "Python.h"

#ifdef  bool
#undef  bool
#endif
#define bool        int

#ifdef  false
#undef  false
#endif
#define false       0

#ifdef  true
#undef  true
#endif
#define true        1

typedef struct {
    PyObject_HEAD
    char *filename;
    PyObject *config;
} PngNQObject;

static PyTypeObject PngNQType;


/* remap_floyd creates the 8 bit indexed pixel data for the output image by examining the 32 bit RGBA pixel data from the input
 * image and selecting the best matching colour for each pixel from the remap[] palette.  remap_floyd keeps track of the
 * accumulated colour error, and performs Floyd-Steinberg style dithering.
 */
static void remap_floyd(mainprog_info *rwpng_info,
                        int cols,
                        int rows,
                        unsigned char map[MAXNETSIZE][4],
                        unsigned int* remap,
                        uch **row_pointers,
                        int quantization_method,
                        int use_alpha_importance_heuristic)
{

    uch *outrow = NULL; /* Output image pixels */

    int i,row;
    #define CLAMP(a) ((a)>=0 ? ((a)<=255 ? (a) : 255)  : 0)

    /* Do each image row */
    for ( row = 0; (ulg)row < rows; ++row ) {
        int offset, nextoffset;
        outrow = rwpng_info->interlaced? row_pointers[row] :
        rwpng_info->indexed_data;

        int rederr=0;
        int blueerr=0;
        int greenerr=0;
        int alphaerr=0;

        offset = row*cols*4;
        nextoffset = offset; if (row+1<rows) nextoffset += cols*4;
        int increment = 4;

        if (0)//row&1)
        {
            offset += cols*4 - 4;
            nextoffset += cols*4 - 4;
            increment = -4;
        }

        for( i=0;i<cols;i++, offset+=increment, nextoffset+=increment)
        {
            int idx;
            unsigned int floyderr = rederr*rederr + greenerr*greenerr + blueerr*blueerr + alphaerr*alphaerr;

            idx = inxsearch(CLAMP(rwpng_info->rgba_data[offset+3] - alphaerr),
                            CLAMP(rwpng_info->rgba_data[offset+2] - blueerr),
                            CLAMP(rwpng_info->rgba_data[offset+1] - greenerr),
                            CLAMP(rwpng_info->rgba_data[offset]   - rederr  ));

            outrow[increment > 0 ? i : cols-i-1] = remap[idx];

            int alpha = MAX(map[idx][3],rwpng_info->rgba_data[offset+3]);
            int colorimp = 255;
            if(use_alpha_importance_heuristic) {
                colorimp = 255 - ((255-alpha) * (255-alpha) / 255);
            }

            int thisrederr=(map[idx][0] -   rwpng_info->rgba_data[offset]) * colorimp   / 255;
            int thisblueerr=(map[idx][1] - rwpng_info->rgba_data[offset+1]) * colorimp  / 255;
            int thisgreenerr=(map[idx][2] -  rwpng_info->rgba_data[offset+2]) * colorimp  / 255;
            int thisalphaerr=map[idx][3] - rwpng_info->rgba_data[offset+3];

            rederr += thisrederr;
            greenerr += thisblueerr;
            blueerr +=  thisgreenerr;
            alphaerr += thisalphaerr;

            unsigned int thiserr = (thisrederr*thisrederr + thisblueerr*thisblueerr + thisgreenerr*thisgreenerr + thisalphaerr*thisalphaerr)*2;
             floyderr = rederr*rederr + greenerr*greenerr + blueerr*blueerr + alphaerr*alphaerr;

            int L = quantization_method;
            while (rederr*rederr > L*L || greenerr*greenerr > L*L || blueerr*blueerr > L*L || alphaerr*alphaerr > L*L ||
                   floyderr > thiserr || floyderr > L*L*2)
            {
                rederr /=2;greenerr /=2;blueerr /=2;alphaerr /=2;
                floyderr = rederr*rederr + greenerr*greenerr + blueerr*blueerr + alphaerr*alphaerr;
            }

            if (i>0)
            {
                rwpng_info->rgba_data[nextoffset-increment+3]=CLAMP(rwpng_info->rgba_data[nextoffset-increment+3] - alphaerr*3/16);
                rwpng_info->rgba_data[nextoffset-increment+2]=CLAMP(rwpng_info->rgba_data[nextoffset-increment+2] - blueerr*3/16 );
                rwpng_info->rgba_data[nextoffset-increment+1]=CLAMP(rwpng_info->rgba_data[nextoffset-increment+1] - greenerr*3/16);
                rwpng_info->rgba_data[nextoffset-increment]  =CLAMP(rwpng_info->rgba_data[nextoffset-increment]   - rederr*3/16  );
            }
            if (i+1<cols)
            {
                rwpng_info->rgba_data[nextoffset+increment+3]=CLAMP(rwpng_info->rgba_data[nextoffset+increment+3] - alphaerr/16);
                rwpng_info->rgba_data[nextoffset+increment+2]=CLAMP(rwpng_info->rgba_data[nextoffset+increment+2] - blueerr/16 );
                rwpng_info->rgba_data[nextoffset+increment+1]=CLAMP(rwpng_info->rgba_data[nextoffset+increment+1] - greenerr/16);
                rwpng_info->rgba_data[nextoffset+increment]  =CLAMP(rwpng_info->rgba_data[nextoffset+increment]   - rederr/16  );
            }
            rwpng_info->rgba_data[nextoffset+3]=CLAMP(rwpng_info->rgba_data[nextoffset+3] - alphaerr*5/16);
            rwpng_info->rgba_data[nextoffset+2]=CLAMP(rwpng_info->rgba_data[nextoffset+2] - blueerr*5/16 );
            rwpng_info->rgba_data[nextoffset+1]=CLAMP(rwpng_info->rgba_data[nextoffset+1] - greenerr*5/16);
            rwpng_info->rgba_data[nextoffset]  =CLAMP(rwpng_info->rgba_data[nextoffset]   - rederr*5/16  );

            rederr = rederr*7/16; greenerr =greenerr*7/16; blueerr =blueerr*7/16; alphaerr =alphaerr*7/16;
        }


        /* if non-interlaced PNG, write row now */
        if (!rwpng_info->interlaced)
            rwpng_write_image_row(rwpng_info);
    }

}

/* remap_simple creates the 8 bit indexed pixel data for the output image by examining the 32 bit RGBA pixel data from the input
 * image and selecting the best matching colour for each pixel from the remap[] palette.  remap_simple does not attempt to
 * dither.
 */
static void remap_simple(mainprog_info *rwpng_info,
                        unsigned int cols,
                        unsigned int rows,
                        unsigned char map[MAXNETSIZE][4],
                        unsigned int* remap,
                        uch **row_pointers)
{
    uch *outrow = NULL; /* Output image pixels */

    unsigned int i,row;
    /* Do each image row */
    for ( row = 0; (ulg)row < rows; ++row )
    {
        unsigned int offset;
        outrow = rwpng_info->interlaced? row_pointers[row] : rwpng_info->indexed_data;
        /* Assign the new colors */
        offset = row*cols*4;
        for( i=0;i<cols;i++){
            outrow[i] = remap[inxsearch(rwpng_info->rgba_data[i*4+offset+3],
                                        rwpng_info->rgba_data[i*4+offset+2],
                                        rwpng_info->rgba_data[i*4+offset+1],
                                        rwpng_info->rgba_data[i*4+offset])];
        }

        /* if non-interlaced PNG, write row now */
        if (!rwpng_info->interlaced){
            rwpng_write_image_row(rwpng_info);
        }
    }
}

/* Methods */

void
pngnq__free(PngNQObject* self)
{
    if (self->filename) {
        free(self->filename);
        self->filename = NULL;
    }

    if (self->config) {
        PyDict_Clear(self->config);
        Py_XDECREF(self->config);
        self->config = NULL;
    }
}

static void
pngnq_dealloc(PngNQObject* self)
{
    if(self != NULL ){
        pngnq__free(self);
        PyObject_Del(self);
    }
}

static PyObject *
pngnq__close(register PngNQObject *self, PyObject *args)
{
    if (!PyArg_ParseTuple(args, ":close")){
        return NULL;
    }
    pngnq__free(self);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pngnq__load(register PngNQObject *self, PyObject *args)
{

    if(self->filename){
        free(self->filename);
    }

    if (!PyArg_ParseTuple(args, "s:load", &self->filename)){
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pngnq__save(register PngNQObject *self, PyObject *args)
{
    FILE                    *fp;
    char                    *filename;

    static mainprog_info    rwpng_info;
    static mainprog_info    rwpngpal_info;

    int                     bot_idx;
    int                     top_idx; /* for remapping of indices */
    unsigned int            remap[MAXNETSIZE];

    ulg                     cols;
    ulg                     rows;
    ulg                     row;

    unsigned char           map[MAXNETSIZE][4];
    int                     x;
    uch                     **row_pointers = NULL; /* Pointers to rows of pixels */

    double                  newcolors  = PyFloat_AsDouble(PyDict_GetItemString(self->config, "n_colours"));
    double                  colour_space  = PyFloat_AsDouble(PyDict_GetItemString(self->config, "colour_space"));
    PyObject                *use_floyd = PyDict_GetItemString(self->config, "use_floyd");
    double                  force_gamma = PyFloat_AsDouble(PyDict_GetItemString(self->config, "force_gamma"));
    double                  sample_factor = PyFloat_AsDouble(PyDict_GetItemString(self->config, "sample_factor"));
    PyObject                *strict_pal_rgba = PyDict_GetItemString(self->config, "strict_pal_rgba");
    double                  write_gamma = PyFloat_AsDouble(PyDict_GetItemString(self->config, "write_gamma"));

    double                  file_gamma;
    double                  quantization_gamma;

    int                     inpal_colours = 0;
    double                  pal_gamma = 1.8;

    char                    *input_palette_name = PyString_AsString(PyDict_GetItemString(self->config, "input_palette_name"));

    double                  unisolate = PyFloat_AsDouble(PyDict_GetItemString(self->config, "unisolate"));
    double                  exclusion_threshold = PyFloat_AsDouble(PyDict_GetItemString(self->config, "exclusion_threshold"));

    double                  r_sens = PyFloat_AsDouble(PyDict_GetItemString(self->config, "r_sens"));
    double                  g_sens = PyFloat_AsDouble(PyDict_GetItemString(self->config, "g_sens"));
    double                  b_sens = PyFloat_AsDouble(PyDict_GetItemString(self->config, "b_sens"));
    double                  a_sens = PyFloat_AsDouble(PyDict_GetItemString(self->config, "a_sens"));

    double                  remap_r_sens = PyFloat_AsDouble(PyDict_GetItemString(self->config, "remap_r_sens"));
    double                  remap_g_sens = PyFloat_AsDouble(PyDict_GetItemString(self->config, "remap_g_sens"));
    double                  remap_b_sens = PyFloat_AsDouble(PyDict_GetItemString(self->config, "remap_b_sens"));
    double                  remap_a_sens = PyFloat_AsDouble(PyDict_GetItemString(self->config, "remap_a_sens"));

    double                  alpha_class_correction = PyFloat_AsDouble(PyDict_GetItemString(self->config, "alpha_class_correction"));
    int                     force_alpha_class_correctness = PyFloat_AsDouble(PyDict_GetItemString(self->config, "force_alpha_class_correctness"));
    int                     use_alpha_importance_heuristic = PyFloat_AsDouble(PyDict_GetItemString(self->config, "use_alpha_importance_heuristic"));

    int                     quantization_method = use_floyd ? 10: 0;


    if (!PyArg_ParseTuple(args, "s:save", &filename)){
        return NULL;
    }

    fp = fopen(self->filename, "r");

    if( fp == NULL ){
        PyErr_SetString(PyExc_ValueError, "file is not exists.");
        return NULL;
    }

    rwpng_read_image(fp, &rwpng_info);
    fclose(fp);

    if (rwpng_info.retval) {
        PyErr_SetString(PyExc_ValueError, "rwpng_read_image() error");
        return NULL;
    }

    fp = fopen(filename, "w");

    if(input_palette_name) {
        inpal_colours = rwpngpal_info.width * rwpngpal_info.height;
    }

    cols = rwpng_info.width;
    rows = rwpng_info.height;

    if(!rwpng_info.rgba_data) {
        PyErr_SetString(PyExc_ValueError, "no pixel data found.");
        return NULL;
    }

    file_gamma = rwpng_info.gamma;
    if (force_gamma > 0) {
         quantization_gamma = force_gamma;
         file_gamma=0;
    } else if (file_gamma > 0) {
        quantization_gamma = file_gamma;
    } else {
        quantization_gamma = 1.8;
        file_gamma = 0;
    }

    rwpng_info.have_bg = false;

    if(write_gamma) {
       rwpng_info.gamma = quantization_gamma;
    } else {
       rwpng_info.gamma = 0;
    }

    if (sample_factor<1) {
        sample_factor = 1 + rows*cols / (512*512);
        if (sample_factor > 10) {
            sample_factor = 10;
        }
    }

    if(input_palette_name){
        if(strict_pal_rgba) {
            pal_gamma = quantization_gamma;
        } else if (force_gamma > 0) {
            pal_gamma = force_gamma;
        } else if (rwpngpal_info.gamma > 0) {
            pal_gamma = rwpngpal_info.gamma;
        } else {
            pal_gamma = 1.8;
        }

        if(pal_gamma != quantization_gamma)
        {
            PyErr_SetString(PyExc_ValueError, "rwpng_read_image() error");
            return NULL;
        }
    }

    /* Start neuquant */
    if(!input_palette_name) {
        palinitnet(NULL, 0, 1.0, (unsigned char*)rwpng_info.rgba_data,rows*cols*4,newcolors,
            colour_space, quantization_gamma, alpha_class_correction,
            force_alpha_class_correctness, r_sens, g_sens, b_sens, a_sens,
            remap_r_sens, remap_g_sens, remap_b_sens, remap_a_sens,
            exclusion_threshold, use_alpha_importance_heuristic);
    } else {
        palinitnet((unsigned char*)rwpngpal_info.rgba_data,inpal_colours,pal_gamma,
               (unsigned char*)rwpng_info.rgba_data,rows*cols*4,newcolors,
            colour_space, quantization_gamma, alpha_class_correction,
            force_alpha_class_correctness, r_sens, g_sens, b_sens, a_sens,
            remap_r_sens, remap_g_sens, remap_b_sens, remap_a_sens,
            exclusion_threshold, use_alpha_importance_heuristic);
    }

    learn(sample_factor, unisolate, false);
    getcolormap((unsigned char*)map, strict_pal_rgba ? 1: 0);

    for (top_idx = newcolors-1, bot_idx = x = 0;  x < newcolors;  ++x) {
        if (map[x][3] == 255) { /* maxval */
            remap[x] = top_idx--;
        } else {
            remap[x] = bot_idx++;
        }
    }

    if (bot_idx != top_idx + 1) {

        if (rwpng_info.row_pointers) {
            free(rwpng_info.row_pointers);
        }

        if (rwpng_info.rgba_data) {
            free(rwpng_info.rgba_data);
        }

        fclose(fp);
        PyErr_SetString(PyExc_ValueError, "Internal logic error.");
        return NULL;
    }

    /* Fill in the palette info in the pngrw structure. */
    rwpng_info.sample_depth = 8;
    rwpng_info.num_palette = newcolors;
    rwpng_info.num_trans = bot_idx;

    /* GRR TO DO:  if bot_idx == 0, check whether all RGB samples are gray
     and if so, whether grayscale sample_depth would be same
     => skip following palette section and go grayscale */

    /* Remap and make palette entries */
    for (x = 0; x < newcolors; ++x) {
        rwpng_info.palette[remap[x]].red  = map[x][0];
        rwpng_info.palette[remap[x]].green = map[x][1];
        rwpng_info.palette[remap[x]].blue = map[x][2];
        rwpng_info.trans[remap[x]] = map[x][3];
    }

    /* Allocate memory*/
    if (rwpng_info.interlaced) {
        if ((rwpng_info.indexed_data = (uch *)malloc(rows * cols)) != NULL) {
            if ((row_pointers = (uch **)malloc(rows * sizeof(uch *))) != NULL) {
                for (row = 0;  (ulg)row < rows;  ++row) {
                    row_pointers[row] = rwpng_info.indexed_data + row*cols;
                }
            }
        }
    } else {
        rwpng_info.indexed_data = (uch *)malloc(cols);
    }

    if (rwpng_info.indexed_data == NULL || (rwpng_info.interlaced && row_pointers == NULL)) {
        if (rwpng_info.row_pointers) {
            free(rwpng_info.row_pointers);
        }

        if (rwpng_info.rgba_data) {
            free(rwpng_info.rgba_data);
        }

        if (rwpng_info.indexed_data) {
            free(rwpng_info.indexed_data);
        }

        fclose(fp);

        PyErr_SetString(PyExc_ValueError, "Insufficient memory for indexed data and/or row pointers.");
        return NULL;
    }

    /* Write headers and such. */
    if (rwpng_write_image_init(fp, &rwpng_info) != 0) {
        if (rwpng_info.rgba_data) {
            free(rwpng_info.rgba_data);
        }

        if (rwpng_info.row_pointers) {
            free(rwpng_info.row_pointers);
        }

        if (rwpng_info.indexed_data) {
            free(rwpng_info.indexed_data);
        }

        if (row_pointers) {
            free(row_pointers);
        }

        fclose(fp);

        PyErr_SetString(PyExc_ValueError, "rwpng_write_image_init() error.");
        return NULL;
    }

    /* Actually build the quantized output image using the colour palette Neuquant has already selected for us. */
    if (quantization_method > 0) {
        remap_floyd(&rwpng_info, cols, rows, map, remap,
                    row_pointers, quantization_method,
                    use_alpha_importance_heuristic);
    } else {
        remap_simple(&rwpng_info, cols,rows,map,remap,row_pointers);
    }

    /* now we're done with the INPUT data and row_pointers, so free 'em */
    if (rwpng_info.rgba_data) {
        free(rwpng_info.rgba_data);
        rwpng_info.rgba_data = NULL;
    }

    if (rwpng_info.row_pointers) {
        free(rwpng_info.row_pointers);
        rwpng_info.row_pointers = NULL;
    }

    /* write entire interlaced palette PNG, or finish/flush noninterlaced one */
    if (rwpng_info.interlaced) {
        rwpng_info.row_pointers = row_pointers;   /* now for OUTPUT data */
        rwpng_write_image_whole(&rwpng_info);
    } else {
        rwpng_write_image_finish(&rwpng_info);
    }

    /* Have finished writing file */
    fclose(fp);

    /* now we're done with the OUTPUT data and row_pointers, too */
    if (rwpng_info.indexed_data) {
        free(rwpng_info.indexed_data);
        rwpng_info.indexed_data = NULL;
    }

    if (row_pointers) {
        free(row_pointers);
        row_pointers = rwpng_info.row_pointers = NULL;
    }

    Py_INCREF(Py_True);
    return Py_True;
}

static PyMethodDef pngnq_methods[] = {
    {"close", (PyCFunction)pngnq__close, METH_VARARGS,
     "close()\nClose the handler."},
    {"load", (PyCFunction)pngnq__load, METH_VARARGS,
     "load(filename)\nload new png image."},
    {"save", (PyCFunction)pngnq__save, METH_VARARGS,
     "save(filename)\n Return bool."},
    {NULL, NULL} /* sentinel */
};

static PyObject *
pngnq_getattr(PngNQObject *self, char *name)
{
    if(strcmp("config", name)==0){
        return self->config;
    }

    if (self->config != NULL) {
        PyObject * v = PyDict_GetItemString(self->config, name);
        if (v != NULL) {
            Py_INCREF(v);
            return v;
        }
    }

    return Py_FindMethod(pngnq_methods, (PyObject *)self, name);
}

static int
pngnq_setattr(PngNQObject *self, char *name, PyObject *val)
{
    if (self->config == NULL) {
        self->config = PyDict_New();
        if (self->config == NULL){
            return -1;
        }
    }

    if(PyDict_Check(val) && strcmp("config", name) == 0){
        return PyDict_Update(self->config, val);
    }

    if (val == NULL) {
        int rv = PyDict_DelItemString(self->config, name);
        if (rv < 0){
            PyErr_SetString(PyExc_AttributeError, "delete non-existing attribute");
        }
        return rv;
    }

    if(strcmp("colour_space", name) == 0 && PyInt_AsLong(val) != RGB
        && PyInt_AsLong(val) != YUV){
        PyErr_Format(PyExc_AttributeError, "colour_space must be %i or %i.", RGB, YUV);
        return -1;
    }

    return PyDict_SetItemString(self->config, name, val);
}

static PyTypeObject PngNQType = {
    PyObject_HEAD_INIT(NULL)
    0,
    "pngnq.PngNQ",
    sizeof(PngNQObject),
    0,
    (destructor)pngnq_dealloc,          /*tp_dealloc*/
    0,                                  /*tp_print*/
    (getattrfunc)pngnq_getattr,         /*tp_getattr*/
    (setattrfunc)pngnq_setattr,         /*tp_setattr*/
    0,                                  /*tp_compare*/
    0,                                  /*tp_repr*/
    0,                                  /*tp_as_number*/
    0,                                  /*tp_as_sequence*/
    0,                                  /*tp_as_mapping*/
    0,                                  /*tp_hash  */
    0,                                  /*tp_call */
    0,                                  /*tp_str */
    0,                                  /*tp_getattro */
    0,                                  /*tp_setattro */
    0,                                  /*tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                 /* tp_flags */
    "pngnq.PngNQ",                    /* tp_doc */
    0,                                  /* tp_traverse */
    0,                                  /* tp_clear */
    0,                                  /* tp_richcompare */
    0,                                  /* tp_weaklistoffset */
    0,                                  /* tp_iter */
    0,                                  /* tp_iternext */
    0,                                  /* tp_methods */
    0,                                  /* tp_members */
    0,                                  /* tp_getset */
    0,                                  /* tp_base */
    0,                                  /* tp_dict */
    0,                                  /* tp_descr_get */
    0,                                  /* tp_descr_set */
    0,                                  /* tp_dictoffset */
    0,                                  /* tp_init */
    0,                                  /* tp_alloc */
    0,                                  /* tp_new */
    0,                                  /* tp_free */
    0,                                  /* tp_is_gc */
};

/* ----------------------------------------------------------------- */
/* pngnq module                                                      */
/* ----------------------------------------------------------------- */

static PyObject *
pngnq_new(PyObject* self, PyObject* args)
{
    PngNQObject *dp;

    dp = PyObject_New(PngNQObject, &PngNQType);

    if (dp == NULL){
        return NULL;
    }

    dp->filename = NULL;
    dp->config = PyDict_New();
    Py_INCREF(dp->config);

    PyDict_SetItemString(dp->config, "input_palette_name", Py_BuildValue("s", ""));
    PyDict_SetItemString(dp->config, "n_colours", PyInt_FromLong(256));
    PyDict_SetItemString(dp->config, "use_floyd", PyBool_FromLong(0));
    PyDict_SetItemString(dp->config, "sample_factor", PyFloat_FromDouble(10));
    PyDict_SetItemString(dp->config, "colour_space", PyInt_FromLong(RGB));
    PyDict_SetItemString(dp->config, "force_gamma", PyFloat_FromDouble(1.0));
    PyDict_SetItemString(dp->config, "write_gamma", PyFloat_FromDouble(0));
    PyDict_SetItemString(dp->config, "strict_pal_rgba", PyBool_FromLong(0));

    PyDict_SetItemString(dp->config, "unisolate", PyFloat_FromDouble(0));
    PyDict_SetItemString(dp->config, "exclusion_threshold", PyFloat_FromDouble(1.0));

    PyDict_SetItemString(dp->config, "r_sens", PyFloat_FromDouble(1.0));
    PyDict_SetItemString(dp->config, "g_sens", PyFloat_FromDouble(1.0));
    PyDict_SetItemString(dp->config, "b_sens", PyFloat_FromDouble(1.0));
    PyDict_SetItemString(dp->config, "a_sens", PyFloat_FromDouble(1.0));

    PyDict_SetItemString(dp->config, "remap_r_sens", PyFloat_FromDouble(1.0));
    PyDict_SetItemString(dp->config, "remap_g_sens", PyFloat_FromDouble(1.0));
    PyDict_SetItemString(dp->config, "remap_b_sens", PyFloat_FromDouble(1.0));
    PyDict_SetItemString(dp->config, "remap_a_sens", PyFloat_FromDouble(1.0));

    PyDict_SetItemString(dp->config, "alpha_class_correction", PyFloat_FromDouble(0.0));
    PyDict_SetItemString(dp->config, "force_alpha_class_correctness", PyFloat_FromDouble(0));
    PyDict_SetItemString(dp->config, "use_alpha_importance_heuristic", PyFloat_FromDouble(0));

    if (!PyArg_ParseTuple(args, "s", &dp->filename)){
        return NULL;
    }

    if(!dp->filename){
        PyErr_SetString(PyExc_ValueError, "ValueError: filename.");
        return NULL;
    }

    return (PyObject *)dp;
}

static PyMethodDef pngnqmodule_methods[] = {
    { "PngNQ", (PyCFunction)pngnq_new, METH_VARARGS,
      "PngNQ() -> mapping\n"
      "Return a pngnq object."},
    { 0, 0 },
};

PyMODINIT_FUNC
initpngnq(void) {
    PyObject *m, *d, *s;

    PngNQType.ob_type = &PyType_Type;

    m = Py_InitModule("pngnq", pngnqmodule_methods);
    if (m == NULL){
        return;
    }

    d = PyModule_GetDict(m);
    s = PyString_FromString("PngNQ");

    if (s != NULL) {
        PyDict_SetItemString(d, "__doc__", s);
        Py_DECREF(s);
    }

    s = PyString_FromString("0.1");
    if (s != NULL) {
        PyDict_SetItemString(d, "__version__", s);
        Py_DECREF(s);
    }

}

