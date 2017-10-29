#include "variables.h"
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

typedef double mod_prec;
//typedef float mod_prec;

// private struct for each thread so that get_global_id does not need to be called everywhere
typedef struct StepData
{
    mod_prec iApp; // External input of the dendrite
    int i;         // Only for debug
    int x;         // Only for debug
    int y;         // Only for debug
    //int prevCellIdx;
    //int newCellIdx;
} StepData;

void ComputeOneCell(private mod_prec *cellCompParamsPtr);
void CompDend(private mod_prec *cellCompParamsPtr);
void DendHCurr(private mod_prec *chPrms_v, private mod_prec *chPrms_prevComp1, private mod_prec *chPrms_newComp1);
void DendCaCurr(private mod_prec *chPrms_v, private mod_prec *chPrms_prevComp1, private mod_prec *chPrms_newComp1);
void DendKCurr(private mod_prec *chPrms_prevComp1, private mod_prec *chPrms_prevComp2, private mod_prec *chPrms_newComp1);
void DendCal(private mod_prec *chPrms_prevComp1, private mod_prec *chPrms_prevComp2, private mod_prec *chPrms_newComp1);
void DendCurrVolt(mod_prec chComps_iC,
                  private mod_prec *chComps_iApp,
                  private mod_prec *chComps_vDend,
                  private mod_prec *chComps_newVDend,
                  private mod_prec *chComps_vSoma,
                  private mod_prec *chComps_q,
                  private mod_prec *chComps_r,
                  private mod_prec *chComps_s,
                  private mod_prec *chComps_newI_CaH);
mod_prec IcNeighbors(private mod_prec *cellCompParamsPtr, mod_prec prevV_dend);

void CompSoma(private mod_prec *cellCompParamsPtr);
void SomaCalcium(private mod_prec *chPrms_v, private mod_prec *chPrms_prevComp1, private mod_prec *chPrms_prevComp2, private mod_prec *chPrms_newComp1, private mod_prec *chPrms_newComp2);
void SomaSodium(private mod_prec *chPrms_v, private mod_prec *chPrms_prevComp1, private mod_prec *chPrms_prevComp2, private mod_prec *chPrms_newComp1, private mod_prec *chPrms_newComp2);
void SomaPotassium(private mod_prec *chPrms_v, private mod_prec *chPrms_prevComp1, private mod_prec *chPrms_prevComp2, private mod_prec *chPrms_newComp1, private mod_prec *chPrms_newComp2);
void SomaPotassiumX(private mod_prec *chPrms_v, private mod_prec *chPrms_prevComp1, private mod_prec *chPrms_newComp1);
void SomaCurrVolt(private mod_prec *chComps_g_CaL,
                  private mod_prec *chComps_vDend,
                  private mod_prec *chComps_vSoma,
                  private mod_prec *chComps_newVSoma,
                  private mod_prec *chComps_vAxon,
                  private mod_prec *chComps_k,
                  private mod_prec *chComps_l,
                  private mod_prec *chComps_m,
                  private mod_prec *chComps_h,
                  private mod_prec *chComps_n,
                  private mod_prec *chComps_x_s);

void CompAxon(private mod_prec *cellCompParamsPtr);
void AxonSodium(private mod_prec *chPrms_v, private mod_prec *chPrms_prevComp1, private mod_prec *chPrms_newComp1, private mod_prec *chPrms_newComp2);
void AxonPotassium(private mod_prec *chPrms_v, private mod_prec *chPrms_prevComp1, private mod_prec *chPrms_newComp1);
void AxonCurrVolt(private mod_prec *chComps_vSoma,
                  private mod_prec *chComps_vAxon,
                  private mod_prec *chComps_newVAxon,
                  private mod_prec *chComps_m_a,
                  private mod_prec *chComps_h_a,
                  private mod_prec *chComps_x_a);