#include "variables.h"
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

typedef double mod_prec;
//typedef float mod_prec;

// Private struct for each thread so that get_global_id does not need to be called everywhere
typedef struct StepData
{
    mod_prec iApp; // External input of the dendrite
    int i;
    int x;         // current position in dimension 1 of the IO network
    int y;         // current position in dimension 2 of the IO network
    int prevCellIdx;
    int newCellIdx;
} StepData;

// Return struct for 2 variables
typedef struct CompRet
{
    mod_prec ret1;
    mod_prec ret2;
} CompRet;

void ComputeOneCell(write_only image2d_t t_cellVDendPtr, global mod_prec *, StepData step);
void CompDend(write_only image2d_t t_cellVDendPtr, global mod_prec *cellCompParamsPtr, StepData step);
mod_prec DendHCurr(mod_prec v, mod_prec prevComp1);
mod_prec DendCaCurr(mod_prec v, mod_prec prevComp1);
mod_prec DendKCurr(mod_prec prevComp1, mod_prec prevComp2);
mod_prec DendCal(mod_prec prevComp1, mod_prec prevComp2);
CompRet DendCurrVolt(mod_prec I_c,
                     mod_prec I_app,
                     mod_prec prevV_dend,
                     mod_prec prevV_soma,
                     mod_prec q,
                     mod_prec r,
                     mod_prec s);
mod_prec IcNeighbors(global mod_prec *cellCompParamsPtr, mod_prec prevV_dend, StepData step);

void CompSoma(global mod_prec *cellCompParamsPtr, StepData step);
CompRet SomaCalcium(mod_prec vSoma, mod_prec prevComp1, mod_prec prevComp2);
CompRet SomaSodium(mod_prec vSoma, mod_prec prevComp1, mod_prec prevComp2);
CompRet SomaPotassium(mod_prec vSoma, mod_prec prevComp1, mod_prec prevComp2);
mod_prec SomaPotassiumX(mod_prec vSoma, mod_prec prevComp1);
mod_prec SomaCurrVolt(mod_prec g_CaL,
                      mod_prec prevV_dend,
                      mod_prec prevV_soma,
                      mod_prec prevV_axon,
                      mod_prec k,
                      mod_prec l,
                      mod_prec m,
                      mod_prec h,
                      mod_prec n,
                      mod_prec x_s);

void CompAxon(global mod_prec *cellCompParamsPtr, StepData step);
CompRet AxonSodium(mod_prec vAxon, mod_prec prevComp1);
mod_prec AxonPotassium(mod_prec vAxon, mod_prec prevComp1);
mod_prec AxonCurrVolt(mod_prec prevV_soma,
                      mod_prec prevV_axon,
                      mod_prec m_a,
                      mod_prec h_a,
                      mod_prec x_a);
