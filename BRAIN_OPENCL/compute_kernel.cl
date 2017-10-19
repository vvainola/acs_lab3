#include "kernel.h"
//#include <math.h>


void ComputeOneCell(private mod_prec *cellCompParamsPtr, mod_prec iApp, StepData step)
{
    //The three compartments can be computed concurrently but only across a single sim step
    CompDend(cellCompParamsPtr, iApp, step);
    CompSoma(cellCompParamsPtr, step);
    CompAxon(cellCompParamsPtr, step);
    return;
}

void CompDend(private mod_prec *cellCompParamsPtr, mod_prec iApp, StepData step)
{
    private mod_prec *chPrms_v;
    private mod_prec *chPrms_prevComp1, *chPrms_prevComp2;
    private mod_prec *chPrms_newComp1; // *chPrms_newComp2;
    private mod_prec *chComps_iApp;
    mod_prec chComps_iC;
    private mod_prec *chComps_vDend;
    private mod_prec *chComps_vSoma;
    private mod_prec *chComps_q, *chComps_r, *chComps_s;
    private mod_prec *chComps_newVDend;
    private mod_prec *chComps_newI_CaH;

    //printf("Dendrite ");

    //Prepare pointers to inputs/outputs
    chPrms_v = &(cellCompParamsPtr[DEND_V]);                     //&cellCompParamsPtr->prevCellState->dend.V_dend;
    chPrms_prevComp1 = &(cellCompParamsPtr[DEND_H]);             //&cellCompParamsPtr->prevCellState->dend.Hcurrent_q;
    chPrms_newComp1 = &(cellCompParamsPtr[PARAM_SIZE + DEND_H]); //&cellCompParamsPtr->newCellState->dend.Hcurrent_q;
    //Compute
    DendHCurr(chPrms_v, chPrms_prevComp1, chPrms_newComp1);

    //Prepare pointers to inputs/outputs
    chPrms_v = &(cellCompParamsPtr[DEND_V]);                       //&cellCompParamsPtr->prevCellState->dend.V_dend;
    chPrms_prevComp1 = &(cellCompParamsPtr[DEND_CAL]);             //&cellCompParamsPtr->prevCellState->dend.Calcium_r;
    chPrms_newComp1 = &(cellCompParamsPtr[PARAM_SIZE + DEND_CAL]); //&cellCompParamsPtr->newCellState->dend.Calcium_r;
    //Compute
    DendCaCurr(chPrms_v, chPrms_prevComp1, chPrms_newComp1);

    //Prepare pointers to inputs/outputs
    chPrms_prevComp1 = &(cellCompParamsPtr[DEND_P]);             //&cellCompParamsPtr->prevCellState->dend.Potassium_s;
    chPrms_prevComp2 = &(cellCompParamsPtr[DEND_CA2]);           //&cellCompParamsPtr->prevCellState->dend.Ca2Plus;
    chPrms_newComp1 = &(cellCompParamsPtr[PARAM_SIZE + DEND_P]); //&cellCompParamsPtr->newCellState->dend.Potassium_s;
    //Compute
    DendKCurr(chPrms_prevComp1, chPrms_prevComp2, chPrms_newComp1);

    //Prepare pointers to inputs/outputs
    chPrms_prevComp1 = &(cellCompParamsPtr[DEND_CA2]);             //&cellCompParamsPtr->prevCellState->dend.Ca2Plus;
    chPrms_prevComp2 = &(cellCompParamsPtr[DEND_I]);               //&cellCompParamsPtr->prevCellState->dend.I_CaH;
    chPrms_newComp1 = &(cellCompParamsPtr[PARAM_SIZE + DEND_CA2]); //&cellCompParamsPtr->newCellState->dend.Ca2Plus;
    //Compute
    DendCal(chPrms_prevComp1, chPrms_prevComp2, chPrms_newComp1);

    chComps_iC = IcNeighbors(cellCompParamsPtr, cellCompParamsPtr[DEND_V]);
    //IcNeighbors(cellCompParamsPtr->neighVdend, cellCompParamsPtr->prevCellState->dend.V_dend);
    //chComps_iApp = &iApp;                                       //&cellCompParamsPtr->iAppIn;
    chComps_vDend = &(cellCompParamsPtr[DEND_V]);                 //&cellCompParamsPtr->prevCellState->dend.V_dend;
    chComps_newVDend = &(cellCompParamsPtr[PARAM_SIZE + DEND_V]); //&cellCompParamsPtr->newCellState->dend.V_dend;
    chComps_vSoma = &(cellCompParamsPtr[SOMA_V]);                 //&cellCompParamsPtr->prevCellState->soma.V_soma
    chComps_q = &(cellCompParamsPtr[PARAM_SIZE + DEND_H]);        //&cellCompParamsPtr->newCellState->dend.Hcurrent_q;
    chComps_r = &(cellCompParamsPtr[PARAM_SIZE + DEND_CAL]);      //&cellCompParamsPtr->newCellState->dend.Calcium_r;
    chComps_s = &(cellCompParamsPtr[PARAM_SIZE + DEND_P]);        //&cellCompParamsPtr->newCellState->dend.Potassium_s;
    chComps_newI_CaH = &(cellCompParamsPtr[PARAM_SIZE + DEND_I]); //&cellCompParamsPtr->newCellState->dend.I_CaH;
    DendCurrVolt(chComps_iC, iApp, chComps_vDend, chComps_newVDend, chComps_vSoma, chComps_q, chComps_r, chComps_s, chComps_newI_CaH);

    return;
}

void DendHCurr(private mod_prec *chPrms_v, private mod_prec *chPrms_prevComp1, private mod_prec *chPrms_newComp1)
{

    mod_prec q_inf, tau_q, dq_dt, q_private;

    //Get inputs
    mod_prec prevV_dend = *chPrms_v;             // *chPrms->v;
    mod_prec prevHcurrent_q = *chPrms_prevComp1; //*chPrms->prevComp1;

    // Update dendritic H current component
    q_inf = 1 / (1 + exp((prevV_dend + 80) / 4));
    tau_q = 1 / (exp(-0.086 * prevV_dend - 14.6) + exp(0.070 * prevV_dend - 1.87));
    dq_dt = (q_inf - prevHcurrent_q) / tau_q;
    q_private = DELTA * dq_dt + prevHcurrent_q;
    //Put result
    *chPrms_newComp1 = q_private;

    return;
}

void DendCaCurr(private mod_prec *chPrms_v, private mod_prec *chPrms_prevComp1, private mod_prec *chPrms_newComp1)
{

    mod_prec alpha_r, beta_r, r_inf, tau_r, dr_dt, r_private;

    //Get inputs
    mod_prec prevV_dend = *chPrms_v;            //*chPrms->v;
    mod_prec prevCalcium_r = *chPrms_prevComp1; //*chPrms->prevComp1;

    // Update dendritic high-threshold Ca current component
    alpha_r = 1.7 / (1 + exp(-(prevV_dend - 5) / 13.9));
    beta_r = 0.02 * (prevV_dend + 8.5) / (exp((prevV_dend + 8.5) / 5) - 1);
    r_inf = alpha_r / (alpha_r + beta_r);
    tau_r = 5 / (alpha_r + beta_r);
    dr_dt = (r_inf - prevCalcium_r) / tau_r;
    r_private = DELTA * dr_dt + prevCalcium_r;
    //Put result
    *chPrms_newComp1 = r_private; // *chPrms->newComp1

    return;
}

void DendKCurr(private mod_prec *chPrms_prevComp1, private mod_prec *chPrms_prevComp2, private mod_prec *chPrms_newComp1)
{

    mod_prec alpha_s = 0.01, beta_s, s_inf, tau_s, ds_dt, s_private;

    //Get inputs
    mod_prec prevPotassium_s = *chPrms_prevComp1; //*chPrms->prevComp1;
    mod_prec prevCa2Plus = *chPrms_prevComp2;     //*chPrms->prevComp2;

    // Update dendritic Ca-dependent K current component
    if ((0.00002 * prevCa2Plus) < 0.01)
        alpha_s = (0.00002 * prevCa2Plus);
    beta_s = 0.015;
    s_inf = alpha_s / (alpha_s + beta_s);
    tau_s = 1 / (alpha_s + beta_s);
    ds_dt = (s_inf - prevPotassium_s) / tau_s;
    s_private = DELTA * ds_dt + prevPotassium_s;
    //Put result
    *chPrms_newComp1 = s_private; //*chPrms->newComp1

    return;
}

//Consider merging DendCal into DendKCurr since DendCal's output doesn't go to DendCurrVolt but to DendKCurr
void DendCal(private mod_prec *chPrms_prevComp1, private mod_prec *chPrms_prevComp2, private mod_prec *chPrms_newComp1)
{
    mod_prec dCa_dt, Ca2Plus_private;

    //Get inputs
    mod_prec prevCa2Plus = *chPrms_prevComp1; //*chPrms->prevComp1;
    mod_prec prevI_CaH = *chPrms_prevComp2;   //*chPrms->prevComp2;

    // update Calcium concentration
    dCa_dt = -3 * prevI_CaH - 0.075 * prevCa2Plus;
    Ca2Plus_private = DELTA * dCa_dt + prevCa2Plus;
    //Put result
    *chPrms_newComp1 = Ca2Plus_private; //*chPrms->newComp1 //This state value is read in DendKCurr

    return;
}

void DendCurrVolt(mod_prec chComps_iC, mod_prec iApp, private mod_prec *chComps_vDend, private mod_prec *chComps_newVDend, private mod_prec *chComps_vSoma, private mod_prec *chComps_q, private mod_prec *chComps_r, private mod_prec *chComps_s, private mod_prec *chComps_newI_CaH)
{

    //Loca variables
    mod_prec I_sd, I_CaH, I_K_Ca, I_ld, I_h, dVd_dt;

    //Get inputs
    mod_prec I_c = chComps_iC;            //chComps->iC;
    mod_prec I_app = iApp;                //*chComps->iApp;
    mod_prec prevV_dend = *chComps_vDend; //*chComps->vDend;
    mod_prec prevV_soma = *chComps_vSoma; //*chComps->vSoma;
    mod_prec q = *chComps_q;              //*chComps->q;
    mod_prec r = *chComps_r;              //*chComps->r;
    mod_prec s = *chComps_s;              //*chComps->s;

    // DENDRITIC CURRENTS

    // Soma-dendrite interaction current I_sd
    I_sd = (G_INT / (1 - P1)) * (prevV_dend - prevV_soma);
    // Inward high-threshold Ca current I_CaH
    I_CaH = G_CAH * r * r * (prevV_dend - V_CA);
    // Outward Ca-dependent K current I_K_Ca
    I_K_Ca = G_K_CA * s * (prevV_dend - V_K);
    // Leakage current I_ld
    I_ld = G_LD * (prevV_dend - V_L);
    // Inward anomalous rectifier I_h
    I_h = G_H * q * (prevV_dend - V_H);

    dVd_dt = (-(I_CaH + I_sd + I_ld + I_K_Ca + I_c + I_h) + I_app) / C_M;

    //Put result (update V_dend)
    *chComps_newVDend = DELTA * dVd_dt + prevV_dend; //*chComps->newVDend
    *chComps_newI_CaH = I_CaH;                       //*chComps->newI_CaH //This is a state value read in DendCal
    return;
}

mod_prec IcNeighbors(private mod_prec *cellCompParamsPtr, mod_prec prevV_dend)
{

    int i;
    mod_prec f, V, I_c;
    //printf("Ic[0]= %f\n", neighVdend[0]);

    I_c = 0;
    for (i = 0; i < 8; i++)
    {
        //printf("%d prevdend: %0.10lf, neighVdend: %0.10lf\n",i, prevV_dend, *neighVdend );
        V = prevV_dend - cellCompParamsPtr[STATE_SIZE + i];
        f = 0.8 * exp(-1 * pow(V, 2) / 100) + 0.2; // SCHWEIGHOFER 2004 VERSION
        I_c = I_c + (CONDUCTANCE * f * V);
    }
    //printf("ja hallo hier is IC, met wie spreek ik: %0.10lf\n", I_c);
    return I_c;
}

void CompSoma(private mod_prec *cellCompParamsPtr, StepData step)
{
    private mod_prec *chPrms_v;
    private mod_prec *chPrms_prevComp1, *chPrms_prevComp2;
    private mod_prec *chPrms_newComp1, *chPrms_newComp2;
    private mod_prec *chComps_g_CaL;
    private mod_prec *chComps_vSoma;
    private mod_prec *chComps_vDend;
    private mod_prec *chComps_vAxon;
    private mod_prec *chComps_k, *chComps_l, *chComps_m, *chComps_h, *chComps_n, *chComps_x_s;
    private mod_prec *chComps_newVSoma;

    // update somatic components
    // SCHWEIGHOFER:

    //Prepare pointers to inputs/outputs
    chPrms_v = &(cellCompParamsPtr[SOMA_V]);                      //&cellCompParamsPtr->prevCellState->soma.V_soma;
    chPrms_prevComp1 = &(cellCompParamsPtr[SOMA_CK]);             //&cellCompParamsPtr->prevCellState->soma.Calcium_k;
    chPrms_prevComp2 = &(cellCompParamsPtr[SOMA_CL]);             //&cellCompParamsPtr->prevCellState->soma.Calcium_l;
    chPrms_newComp1 = &(cellCompParamsPtr[PARAM_SIZE + SOMA_CK]); //&cellCompParamsPtr->newCellState->soma.Calcium_k;
    chPrms_newComp2 = &(cellCompParamsPtr[PARAM_SIZE + SOMA_CL]); //&cellCompParamsPtr->newCellState->soma.Calcium_l;
    //Compute
    SomaCalcium(chPrms_v, chPrms_prevComp1, chPrms_prevComp2, chPrms_newComp1, chPrms_newComp2);

    //Prepare pointers to inputs/outputs
    chPrms_v = &(cellCompParamsPtr[SOMA_V]);                      //&cellCompParamsPtr->prevCellState->soma.V_soma;
    chPrms_prevComp1 = &(cellCompParamsPtr[SOMA_SM]);             //&cellCompParamsPtr->prevCellState->soma.Sodium_m;
    chPrms_prevComp2 = &(cellCompParamsPtr[SOMA_SH]);             //&cellCompParamsPtr->prevCellState->soma.Sodium_h;
    chPrms_newComp1 = &(cellCompParamsPtr[PARAM_SIZE + SOMA_SM]); //&cellCompParamsPtr->newCellState->soma.Sodium_m;
    chPrms_newComp2 = &(cellCompParamsPtr[PARAM_SIZE + SOMA_SH]); //&cellCompParamsPtr->newCellState->soma.Sodium_h;
    //Compute
    SomaSodium(chPrms_v, chPrms_prevComp1, chPrms_prevComp2, chPrms_newComp1, chPrms_newComp2);

    //Prepare pointers to inputs/outputs
    chPrms_v = &(cellCompParamsPtr[SOMA_V]);                      //&cellCompParamsPtr->prevCellState->soma.V_soma;
    chPrms_prevComp1 = &(cellCompParamsPtr[SOMA_PN]);             //&cellCompParamsPtr->prevCellState->soma.Potassium_n;
    chPrms_prevComp2 = &(cellCompParamsPtr[SOMA_PP]);             //&cellCompParamsPtr->prevCellState->soma.Potassium_p;
    chPrms_newComp1 = &(cellCompParamsPtr[PARAM_SIZE + SOMA_PN]); //&cellCompParamsPtr->newCellState->soma.Potassium_n;
    chPrms_newComp2 = &(cellCompParamsPtr[PARAM_SIZE + SOMA_PP]); //&cellCompParamsPtr->newCellState->soma.Potassium_p;
    //Compute
    SomaPotassium(chPrms_v, chPrms_prevComp1, chPrms_prevComp2, chPrms_newComp1, chPrms_newComp2);

    //Prepare pointers to inputs/outputs
    chPrms_v = &(cellCompParamsPtr[SOMA_V]);                       //&cellCompParamsPtr->prevCellState->soma.V_soma;
    chPrms_prevComp1 = &(cellCompParamsPtr[SOMA_PXS]);             //&cellCompParamsPtr->prevCellState->soma.Potassium_x_s;
    chPrms_newComp1 = &(cellCompParamsPtr[PARAM_SIZE + SOMA_PXS]); //&cellCompParamsPtr->newCellState->soma.Potassium_x_s;
    //Compute
    SomaPotassiumX(chPrms_v, chPrms_prevComp1, chPrms_newComp1);

    chComps_g_CaL = &(cellCompParamsPtr[SOMA_G]);                 //&cellCompParamsPtr->prevCellState->soma.g_CaL;
    chComps_vDend = &(cellCompParamsPtr[DEND_V]);                 //&cellCompParamsPtr->prevCellState->dend.V_dend;
    chComps_vSoma = &(cellCompParamsPtr[SOMA_V]);                 //&cellCompParamsPtr->prevCellState->soma.V_soma;
    chComps_newVSoma = &(cellCompParamsPtr[PARAM_SIZE + SOMA_V]); //&cellCompParamsPtr->newCellState->soma.V_soma;
    chComps_vAxon = &(cellCompParamsPtr[AXON_V]);                 //&cellCompParamsPtr->prevCellState->axon.V_axon;
    chComps_k = &(cellCompParamsPtr[PARAM_SIZE + SOMA_CK]);       //&cellCompParamsPtr->newCellState->soma.Calcium_k;
    chComps_l = &(cellCompParamsPtr[PARAM_SIZE + SOMA_CL]);       //&cellCompParamsPtr->newCellState->soma.Calcium_l;
    chComps_m = &(cellCompParamsPtr[PARAM_SIZE + SOMA_SM]);       //&cellCompParamsPtr->newCellState->soma.Sodium_m;
    chComps_h = &(cellCompParamsPtr[PARAM_SIZE + SOMA_SH]);       //&cellCompParamsPtr->newCellState->soma.Sodium_h;
    chComps_n = &(cellCompParamsPtr[PARAM_SIZE + SOMA_PN]);       //&cellCompParamsPtr->newCellState->soma.Potassium_n;
    chComps_x_s = &(cellCompParamsPtr[PARAM_SIZE + SOMA_PXS]);    // &cellCompParamsPtr->newCellState->soma.Potassium_x_s;
    SomaCurrVolt(chComps_g_CaL, chComps_vDend, chComps_vSoma, chComps_newVSoma, chComps_vAxon, chComps_k, chComps_l, chComps_m, chComps_h, chComps_n, chComps_x_s);

    return;
}

void SomaCalcium(private mod_prec *chPrms_v, private mod_prec *chPrms_prevComp1, private mod_prec *chPrms_prevComp2, private mod_prec *chPrms_newComp1, private mod_prec *chPrms_newComp2)
{
    mod_prec k_inf, l_inf, tau_k, tau_l, dk_dt, dl_dt, k_private, l_private;

    //Get inputs
    mod_prec prevV_soma = *chPrms_v;            //*chPrms->v;
    mod_prec prevCalcium_k = *chPrms_prevComp1; //*chPrms->prevComp1;
    mod_prec prevCalcium_l = *chPrms_prevComp2; //*chPrms->prevComp2;

    k_inf = (1 / (1 + exp(-1 * (prevV_soma + 61) / 4.2)));
    l_inf = (1 / (1 + exp((prevV_soma + 85.5) / 8.5)));
    tau_k = 1;
    tau_l = ((20 * exp((prevV_soma + 160) / 30) / (1 + exp((prevV_soma + 84) / 7.3))) + 35);
    dk_dt = (k_inf - prevCalcium_k) / tau_k;
    dl_dt = (l_inf - prevCalcium_l) / tau_l;
    k_private = DELTA * dk_dt + prevCalcium_k;
    l_private = DELTA * dl_dt + prevCalcium_l;
    //Put result
    *chPrms_newComp1 = k_private; //*chPrms->newComp1
    *chPrms_newComp2 = l_private; //*chPrms->newComp2

    return;
}

void SomaSodium(private mod_prec *chPrms_v, private mod_prec *chPrms_prevComp1, private mod_prec *chPrms_prevComp2, private mod_prec *chPrms_newComp1, private mod_prec *chPrms_newComp2)
{
    mod_prec m_inf, h_inf, tau_h, dh_dt, m_private, h_private;

    //Get inputs
    mod_prec prevV_soma = *chPrms_v; //*chPrms->v;
    //mod_prec prevSodium_m = *chPrms->prevComp1;
    mod_prec prevSodium_h = *chPrms_prevComp2; //*chPrms->prevComp2;

    // RAT THALAMOCORTICAL SODIUM:
    m_inf = 1 / (1 + (exp((-30 - prevV_soma) / 5.5)));
    h_inf = 1 / (1 + (exp((-70 - prevV_soma) / -5.8)));
    tau_h = 3 * exp((-40 - prevV_soma) / 33);
    dh_dt = (h_inf - prevSodium_h) / tau_h;
    m_private = m_inf;
    h_private = prevSodium_h + DELTA * dh_dt;
    //Put result
    *chPrms_newComp1 = m_private; //*chPrms->newComp1
    *chPrms_newComp2 = h_private; //*chPrms->newComp2

    return;
}

void SomaPotassium(private mod_prec *chPrms_v, private mod_prec *chPrms_prevComp1, private mod_prec *chPrms_prevComp2, private mod_prec *chPrms_newComp1, private mod_prec *chPrms_newComp2)
{
    mod_prec n_inf, p_inf, tau_n, tau_p, dn_dt, dp_dt, n_private, p_private;

    //Get inputs
    mod_prec prevV_soma = *chPrms_v;              //*chPrms->v;
    mod_prec prevPotassium_n = *chPrms_prevComp1; //*chPrms->prevComp1;
    mod_prec prevPotassium_p = *chPrms_prevComp2; //*chPrms->prevComp2;

    // NEOCORTICAL
    n_inf = 1 / (1 + exp((-3 - prevV_soma) / 10));
    p_inf = 1 / (1 + exp((-51 - prevV_soma) / -12));
    tau_n = 5 + (47 * exp(-(-50 - prevV_soma) / 900));
    tau_p = tau_n;
    dn_dt = (n_inf - prevPotassium_n) / tau_n;
    dp_dt = (p_inf - prevPotassium_p) / tau_p;
    n_private = DELTA * dn_dt + prevPotassium_n;
    p_private = DELTA * dp_dt + prevPotassium_p;
    //Put result
    *chPrms_newComp1 = n_private; //*chPrms->newComp1
    *chPrms_newComp2 = p_private; //*chPrms->newComp2

    return;
}

void SomaPotassiumX(private mod_prec *chPrms_v, private mod_prec *chPrms_prevComp1, private mod_prec *chPrms_newComp1)
{

    mod_prec alpha_x_s, beta_x_s, x_inf_s, tau_x_s, dx_dt_s, x_s_private;

    //Get inputs
    mod_prec prevV_soma = *chPrms_v;                //*chPrms->v;
    mod_prec prevPotassium_x_s = *chPrms_prevComp1; //*chPrms->prevComp1;

    // Voltage-dependent (fast) potassium
    alpha_x_s = 0.13 * (prevV_soma + 25) / (1 - exp(-(prevV_soma + 25) / 10));
    beta_x_s = 1.69 * exp(-0.0125 * (prevV_soma + 35));
    x_inf_s = alpha_x_s / (alpha_x_s + beta_x_s);
    tau_x_s = 1 / (alpha_x_s + beta_x_s);
    dx_dt_s = (x_inf_s - prevPotassium_x_s) / tau_x_s;
    x_s_private = 0.05 * dx_dt_s + prevPotassium_x_s;
    //Put result
    *chPrms_newComp1 = x_s_private; //*chPrms->newComp1

    return;
}

void SomaCurrVolt(private mod_prec *chComps_g_CaL, private mod_prec *chComps_vDend, private mod_prec *chComps_vSoma, private mod_prec *chComps_newVSoma, private mod_prec *chComps_vAxon, private mod_prec *chComps_k, private mod_prec *chComps_l, private mod_prec *chComps_m, private mod_prec *chComps_h, private mod_prec *chComps_n, private mod_prec *chComps_x_s)
{
    //private variables
    mod_prec I_ds, I_CaL, I_Na_s, I_ls, I_Kdr_s, I_K_s, I_as, dVs_dt;

    //Get inputs
    mod_prec g_CaL = *chComps_g_CaL;      //*chComps->g_CaL;
    mod_prec prevV_dend = *chComps_vDend; //*chComps->vDend;
    mod_prec prevV_soma = *chComps_vSoma; //*chComps->vSoma;
    mod_prec prevV_axon = *chComps_vAxon; //*chComps->vAxon;
    mod_prec k = *chComps_k;              //*chComps->k;
    mod_prec l = *chComps_l;              //*chComps->l;
    mod_prec m = *chComps_m;              //*chComps->m;
    mod_prec h = *chComps_h;              //*chComps->h;
    mod_prec n = *chComps_n;              //*chComps->n;
    mod_prec x_s = *chComps_x_s;          //*chComps->x_s;

    // SOMATIC CURRENTS

    // Dendrite-soma interaction current I_ds
    I_ds = (G_INT / P1) * (prevV_soma - prevV_dend);
    // Inward low-threshold Ca current I_CaL
    I_CaL = g_CaL * k * k * k * l * (prevV_soma - V_CA); //k^3
    // Inward Na current I_Na_s
    I_Na_s = G_NA_S * m * m * m * h * (prevV_soma - V_NA);
    // Leakage current I_ls
    I_ls = G_LS * (prevV_soma - V_L);
    // Outward delayed potassium current I_Kdr
    I_Kdr_s = G_KDR_S * n * n * n * n * (prevV_soma - V_K); // SCHWEIGHOFER
    // I_K_s
    I_K_s = G_K_S * pow(x_s, 4) * (prevV_soma - V_K);
    // Axon-soma interaction current I_as
    I_as = (G_INT / (1 - P2)) * (prevV_soma - prevV_axon);

    dVs_dt = (-(I_CaL + I_ds + I_as + I_Na_s + I_ls + I_Kdr_s + I_K_s)) / C_M;
    *chComps_newVSoma = DELTA * dVs_dt + prevV_soma; // *chComps->newVSoma

    return;
}

void CompAxon(private mod_prec *cellCompParamsPtr, StepData step)
{

    private mod_prec *chPrms_v;
    private mod_prec *chPrms_prevComp1; // *chPrms_prevComp2;
    private mod_prec *chPrms_newComp1, *chPrms_newComp2;
    private mod_prec *chComps_vSoma;
    private mod_prec *chComps_vAxon;
    private mod_prec *chComps_m_a, *chComps_h_a, *chComps_x_a;
    private mod_prec *chComps_newVAxon;

    // update somatic components
    // SCHWEIGHOFER:

    //Prepare pointers to inputs/outputs
    chPrms_v = &(cellCompParamsPtr[AXON_V]);                      //prevCellState->axon.V_axon;
    chPrms_prevComp1 = &(cellCompParamsPtr[AXON_SH]);             //prevCellState->axon.Sodium_h_a;
    chPrms_newComp1 = &(cellCompParamsPtr[PARAM_SIZE + AXON_SH]); //&cellCompParamsPtr->newCellState->axon.Sodium_h_a;
    chPrms_newComp2 = &(cellCompParamsPtr[PARAM_SIZE + AXON_SM]); //&cellCompParamsPtr->newCellState->axon.Sodium_m_a;
    //Compute
    AxonSodium(chPrms_v, chPrms_prevComp1, chPrms_newComp1, chPrms_newComp2);

    //Prepare pointers to inputs/outputs
    chPrms_v = &(cellCompParamsPtr[AXON_V]);                     //&cellCompParamsPtr->prevCellState->axon.V_axon;
    chPrms_prevComp1 = &(cellCompParamsPtr[AXON_P]);             //&cellCompParamsPtr->prevCellState->axon.Potassium_x_a;
    chPrms_newComp1 = &(cellCompParamsPtr[PARAM_SIZE + AXON_P]); //&cellCompParamsPtr->newCellState->axon.Potassium_x_a;
    //Compute
    AxonPotassium(chPrms_v, chPrms_prevComp1, chPrms_newComp1);

    //Get inputs
    chComps_vSoma = &(cellCompParamsPtr[SOMA_V]);                 //&cellCompParamsPtr->prevCellState->soma.V_soma;
    chComps_vAxon = &(cellCompParamsPtr[AXON_V]);                 //&cellCompParamsPtr->prevCellState->axon.V_axon;
    chComps_newVAxon = &(cellCompParamsPtr[PARAM_SIZE + AXON_V]); //&cellCompParamsPtr->newCellState->axon.V_axon;
    chComps_m_a = &(cellCompParamsPtr[PARAM_SIZE + AXON_SM]);     //&cellCompParamsPtr->newCellState->axon.Sodium_m_a;
    chComps_h_a = &(cellCompParamsPtr[PARAM_SIZE + AXON_SH]);     //&cellCompParamsPtr->newCellState->axon.Sodium_h_a;
    chComps_x_a = &(cellCompParamsPtr[PARAM_SIZE + AXON_P]);      //&cellCompParamsPtr->newCellState->axon.Potassium_x_a;
    AxonCurrVolt(chComps_vSoma, chComps_vAxon, chComps_newVAxon, chComps_m_a, chComps_h_a, chComps_x_a);

    return;
}

void AxonSodium(private mod_prec *chPrms_v, private mod_prec *chPrms_prevComp1, private mod_prec *chPrms_newComp1, private mod_prec *chPrms_newComp2)
{

    mod_prec m_inf_a, h_inf_a, tau_h_a, dh_dt_a, m_a_private, h_a_private;

    //Get inputs
    mod_prec prevV_axon = *chPrms_v;             //*chPrms->v;
    mod_prec prevSodium_h_a = *chPrms_prevComp1; //*chPrms->prevComp1;

    // Update axonal Na components
    // NOTE: current has shortened inactivation to account for high
    // firing frequencies in axon hillock
    m_inf_a = 1 / (1 + (exp((-30 - prevV_axon) / 5.5)));
    h_inf_a = 1 / (1 + (exp((-60 - prevV_axon) / (-5.8))));
    tau_h_a = 1.5 * exp((-40 - prevV_axon) / 33);
    dh_dt_a = (h_inf_a - prevSodium_h_a) / tau_h_a;
    m_a_private = m_inf_a;
    h_a_private = prevSodium_h_a + DELTA * dh_dt_a;
    //Put result
    *chPrms_newComp1 = h_a_private; //*chPrms->newComp1
    *chPrms_newComp2 = m_a_private; //*chPrms->newComp2

    return;
}

void AxonPotassium(private mod_prec *chPrms_v, private mod_prec *chPrms_prevComp1, private mod_prec *chPrms_newComp1)
{

    mod_prec alpha_x_a, beta_x_a, x_inf_a, tau_x_a, dx_dt_a, x_a_private;

    //Get inputs
    mod_prec prevV_axon = *chPrms_v;                //*chPrms->v;
    mod_prec prevPotassium_x_a = *chPrms_prevComp1; //*chPrms->prevComp1;

    // D'ANGELO 2001 -- Voltage-dependent potassium
    alpha_x_a = 0.13 * (prevV_axon + 25) / (1 - exp(-(prevV_axon + 25) / 10));
    beta_x_a = 1.69 * exp(-0.0125 * (prevV_axon + 35));
    x_inf_a = alpha_x_a / (alpha_x_a + beta_x_a);
    tau_x_a = 1 / (alpha_x_a + beta_x_a);
    dx_dt_a = (x_inf_a - prevPotassium_x_a) / tau_x_a;
    x_a_private = 0.05 * dx_dt_a + prevPotassium_x_a;
    //Put result
    *chPrms_newComp1 = x_a_private; //*chPrms->newComp1

    return;
}

void AxonCurrVolt(private mod_prec *chComps_vSoma, private mod_prec *chComps_vAxon, private mod_prec *chComps_newVAxon, private mod_prec *chComps_m_a, private mod_prec *chComps_h_a, private mod_prec *chComps_x_a)
{

    //private variable
    mod_prec I_Na_a, I_la, I_sa, I_K_a, dVa_dt;

    //Get inputs
    mod_prec prevV_soma = *chComps_vSoma; //*chComps->vSoma;
    mod_prec prevV_axon = *chComps_vAxon; //*chComps->vAxon;
    mod_prec m_a = *chComps_m_a;          //*chComps->m_a;
    mod_prec h_a = *chComps_h_a;          //*chComps->h_a;
    mod_prec x_a = *chComps_x_a;          //*chComps->x_a;

    // AXONAL CURRENTS
    // Sodium
    I_Na_a = G_NA_A * m_a * m_a * m_a * h_a * (prevV_axon - V_NA);
    // Leak
    I_la = G_LA * (prevV_axon - V_L);
    // Soma-axon interaction current I_sa
    I_sa = (G_INT / P2) * (prevV_axon - prevV_soma);
    // Potassium (transient)
    //I_K_a   = G_K_A * pow(x_a, 4) * (prevV_axon - V_K);
    I_K_a = G_K_A * x_a * x_a * x_a * x_a * (prevV_axon - V_K);
    dVa_dt = (-(I_K_a + I_sa + I_la + I_Na_a)) / C_M;
    *chComps_newVAxon = DELTA * dVa_dt + prevV_axon; //*chComps->newVAxon
    return;
}

/**
Input: cellCompParamsPtr, cellStatePtr, iApp ,i
cellCompParamsPtr: Array of struct which stores values of neighbours for each cell.
cellStatePtr: Array with values for each cell.
iApp: Extenal input of the dendrite
i: Current simulation step 

Retreive the external input of the dedrite 
and update the previous and new state of the current cell.
Then Compute the new variables of the current cell with ComputeOneCell.
**/
__kernel void compute_kernel(global mod_prec *cellStatePtr, global mod_prec *cellCompParamsPtr, mod_prec iApp, uint i)
{

    private mod_prec cellCompParamsPtrPrivate[LOCAL_PARAM_SIZE];

    int y = get_global_id(0);
    int x = get_global_id(1); 

    // Previous cell state in indices 0-26
    // Next cell state in indices 27-53
    for (int idx = 0; idx < STATE_SIZE; idx++)
    {
        // Previous cell state
        cellCompParamsPtrPrivate[idx] = cellStatePtr[(i % 2) * IO_NETWORK_SIZE * STATE_SIZE + (y * IO_NETWORK_DIM2 + x) * STATE_SIZE + idx];
        // Next cell state
        cellCompParamsPtrPrivate[PARAM_SIZE + idx] = cellStatePtr[((i % 2) ^ 1) * IO_NETWORK_SIZE * STATE_SIZE + (y * IO_NETWORK_DIM2 + x) * STATE_SIZE + idx];
    }

    // Copy neighbour data to private pointer
    for (int idx = STATE_SIZE; idx < PARAM_SIZE; idx++)
    {
        cellCompParamsPtrPrivate[idx] = cellCompParamsPtr[(y * IO_NETWORK_DIM2 + x) * LOCAL_PARAM_SIZE + idx];
    }

    StepData step;
    step.i = i;
    step.x = x;
    step.y = y;

    ComputeOneCell(cellCompParamsPtrPrivate, iApp, step);

    // Copy data back
    for (int idx = 0; idx < STATE_SIZE; idx++)
    {
        cellStatePtr[((i % 2) ^ 1) * IO_NETWORK_SIZE * STATE_SIZE + (y * IO_NETWORK_DIM2 + x) * STATE_SIZE + idx] = cellCompParamsPtrPrivate[PARAM_SIZE + idx];
    }
}
