/*
 *
 * Copyright (c) 2012, Neurasmus B.V., The Netherlands,
 * web: www.neurasmus.com email: info@neurasmus.com
 *
 * Any use reproduction in whole or in parts is prohibited
 * without the written consent of the copyright owner.
 *
 * All Rights Reserved.
 *
 *
 * Author: Sebastian Isaza
 * Created: 10-04-2012
 * Modified: 06-06-2012
 *
 * Description : Top header file of the Inferior Olive model. It contains the
 * constant model conductances, the data structures that hold the cell state and
 * the function prototypes.
 *
 */


#ifndef MAIN_H_
#define MAIN_H_
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "variables.h" 

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif


#ifdef DEBUG
# define DEBUG_PRINT(x) printf x
#else
# define DEBUG_PRINT(x) do {} while (0)
#endif


/*** TYPEDEFS AND STRUCTS***/

typedef cl_double cl_mod_prec;
//typedef cl_float cl_mod_prec

typedef struct StepData{
	cl_mod_prec iApp;  // External input of the dendrite
	cl_int x;       // current position in dimension 1 of the IO network
	cl_int y;       // current position in dimension 2 of the IO network
} StepData;



#endif /* MAIN_H_ */
