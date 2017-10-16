#include "variables.h"

typedef double mod_prec;
//typedef float mod_prec;

// Private struct for each thread so that get_global_id does not need to be called everywhere
typedef struct StepData
{
    mod_prec iApp; // External input of the dendrite
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

void ComputeOneCell(global mod_prec *, StepData step);
void CompDend(global mod_prec *cellCompParamsPtr, StepData step);
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
void SomaSodium(mod_prec *chPrms_v, mod_prec *chPrms_prevComp1, mod_prec *chPrms_prevComp2, mod_prec *chPrms_newComp1, mod_prec *chPrms_newComp2);
void SomaPotassium(mod_prec *chPrms_v, mod_prec *chPrms_prevComp1, mod_prec *chPrms_prevComp2, mod_prec *chPrms_newComp1, mod_prec *chPrms_newComp2);
void SomaPotassiumX(mod_prec *chPrms_v, mod_prec *chPrms_prevComp1, mod_prec *chPrms_newComp1);
void SomaCurrVolt(mod_prec *chComps_g_CaL, mod_prec *chComps_vDend, mod_prec *chComps_vSoma, mod_prec *chComps_newVSoma, mod_prec *chComps_vAxon, mod_prec *chComps_k, mod_prec *chComps_l, mod_prec *chComps_m, mod_prec *chComps_h, mod_prec *chComps_n, mod_prec *chComps_x_s);

void CompAxon(global mod_prec *cellCompParamsPtr, StepData step);
void AxonSodium(mod_prec *chPrms_v, mod_prec *chPrms_prevComp1, mod_prec *chPrms_newComp1, mod_prec *chPrms_newComp2);
void AxonPotassium(mod_prec *chPrms_v, mod_prec *chPrms_prevComp1, mod_prec *chPrms_newComp1);
void AxonCurrVolt(mod_prec *chComps_vSoma, mod_prec *chComps_vAxon, mod_prec *chComps_newVAxon, mod_prec *chComps_m_a, mod_prec *chComps_h_a, mod_prec *chComps_x_a);
