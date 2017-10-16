#include "kernel.h"

//#include <math.h>
void ComputeOneCell(global mod_prec *cellCompParamsPtr, StepData step)
{

    //The three compartments can be computed concurrently but only across a single sim step
    CompDend(cellCompParamsPtr, step);
    CompSoma(cellCompParamsPtr, step);
    CompAxon(cellCompParamsPtr, step);
    return;
}

void CompDend(global mod_prec *cellCompParamsPtr, StepData step)
{

    //define variables
    mod_prec vDend;
    mod_prec prevComp1;
    mod_prec prevComp2;

    //printf("Dendrite ");

    //Prepare pointers to inputs/outputs
    vDend = cellCompParamsPtr[step.prevCellIdx + DEND_V];
    prevComp1 = cellCompParamsPtr[step.prevCellIdx + DEND_H];
    //Compute
    cellCompParamsPtr[step.newCellIdx + DEND_H] = DendHCurr(vDend, prevComp1);

    //Prepare pointers to inputs/outputs
    //vDend = cellCompParamsPtr[step.prevCellIdx + DEND_V];
    prevComp1 = cellCompParamsPtr[step.prevCellIdx + DEND_CAL];
    //Compute
    cellCompParamsPtr[step.newCellIdx + DEND_CAL] = DendCaCurr(vDend, prevComp1);

    //Prepare pointers to inputs/outputs
    prevComp1 = cellCompParamsPtr[step.prevCellIdx + DEND_P];
    prevComp2 = cellCompParamsPtr[step.prevCellIdx + DEND_CA2];
    //Compute
    cellCompParamsPtr[step.newCellIdx + DEND_P] = DendKCurr(prevComp1, prevComp2);

    //Prepare pointers to inputs/outputs
    prevComp1 = cellCompParamsPtr[step.prevCellIdx + DEND_CA2];
    prevComp2 = cellCompParamsPtr[step.prevCellIdx + DEND_I];
    //Compute
    cellCompParamsPtr[step.newCellIdx + DEND_P] = DendCal(prevComp1, prevComp2);

    //Prepare pointers to inputs/outputs
    mod_prec Ic = IcNeighbors(cellCompParamsPtr, vDend, step);
    //vDend = cellCompParamsPtr[step.prevCellIdx + DEND_V];
    mod_prec vSoma = cellCompParamsPtr[step.prevCellIdx + SOMA_V];
    mod_prec q = cellCompParamsPtr[step.newCellIdx + DEND_H];
    mod_prec r = cellCompParamsPtr[step.newCellIdx + DEND_CAL];
    mod_prec s = cellCompParamsPtr[step.newCellIdx + DEND_P];
    //Compute
    CompRet dendCurrRet = DendCurrVolt(Ic,
                                       step.iApp,
                                       vDend,
                                       vSoma,
                                       q,
                                       r,
                                       s);
    cellCompParamsPtr[step.newCellIdx + DEND_V] = dendCurrRet.ret1;
    cellCompParamsPtr[step.newCellIdx + DEND_I] = dendCurrRet.ret2;

    return;
}

mod_prec DendHCurr(mod_prec v, mod_prec prevComp1)
{

    mod_prec q_inf, tau_q, dq_dt, q_local;

    //Get inputs
    mod_prec prevV_dend = v;             // *chPrms->v;
    mod_prec prevHcurrent_q = prevComp1; //*chPrms->prevComp1;

    // Update dendritic H current component
    q_inf = 1 / (1 + exp((prevV_dend + 80) / 4));
    tau_q = 1 / (exp(-0.086 * prevV_dend - 14.6) + exp(0.070 * prevV_dend - 1.87));
    dq_dt = (q_inf - prevHcurrent_q) / tau_q;
    q_local = DELTA * dq_dt + prevHcurrent_q;

    //Return result
    return q_local;
}

mod_prec DendCaCurr(mod_prec v, mod_prec prevComp1)
{

    mod_prec alpha_r, beta_r, r_inf, tau_r, dr_dt, r_local;

    //Get inputs
    mod_prec prevV_dend = v;            //*chPrms->v;
    mod_prec prevCalcium_r = prevComp1; //*chPrms->prevComp1;

    // Update dendritic high-threshold Ca current component
    alpha_r = 1.7 / (1 + exp(-(prevV_dend - 5) / 13.9));
    beta_r = 0.02 * (prevV_dend + 8.5) / (exp((prevV_dend + 8.5) / 5) - 1);
    r_inf = alpha_r / (alpha_r + beta_r);
    tau_r = 5 / (alpha_r + beta_r);
    dr_dt = (r_inf - prevCalcium_r) / tau_r;
    r_local = DELTA * dr_dt + prevCalcium_r;
    //Put result
    return r_local; // *chPrms->newComp1
}

mod_prec DendKCurr(mod_prec prevComp1, mod_prec prevComp2)
{

    mod_prec alpha_s = 0.01, beta_s, s_inf, tau_s, ds_dt, s_local;

    //Get inputs
    mod_prec prevPotassium_s = prevComp1; //*chPrms->prevComp1;
    mod_prec prevCa2Plus = prevComp2;     //*chPrms->prevComp2;

    // Update dendritic Ca-dependent K current component
    if ((0.00002 * prevCa2Plus) < 0.01)
        alpha_s = (0.00002 * prevCa2Plus);
    beta_s = 0.015;
    s_inf = alpha_s / (alpha_s + beta_s);
    tau_s = 1 / (alpha_s + beta_s);
    ds_dt = (s_inf - prevPotassium_s) / tau_s;
    s_local = DELTA * ds_dt + prevPotassium_s;
    //Put result
    return s_local; //*chPrms->newComp1
}

//Consider merging DendCal into DendKCurr since DendCal's output doesn't go to DendCurrVolt but to DendKCurr
mod_prec DendCal(mod_prec prevComp1, mod_prec prevComp2)
{
    mod_prec dCa_dt, Ca2Plus_local;

    //Get inputs
    mod_prec prevCa2Plus = prevComp1; //*chPrms->prevComp1;
    mod_prec prevI_CaH = prevComp2;   //*chPrms->prevComp2;

    // update Calcium concentration
    dCa_dt = -3 * prevI_CaH - 0.075 * prevCa2Plus;
    Ca2Plus_local = DELTA * dCa_dt + prevCa2Plus;
    //Put result
    return Ca2Plus_local; //*chPrms->newComp1 //This state value is read in DendKCurr
}

CompRet DendCurrVolt(mod_prec I_c,
                     mod_prec I_app,
                     mod_prec prevV_dend,
                     mod_prec prevV_soma,
                     mod_prec q,
                     mod_prec r,
                     mod_prec s)
{

    //Local variables
    mod_prec I_sd, I_CaH, I_K_Ca, I_ld, I_h, dVd_dt;

    //Get inputs
    /* mod_prec I_c = chComps_iC;            //chComps->iC;
    mod_prec I_app = *chComps_iApp;       //*chComps->iApp;
    mod_prec prevV_dend = *chComps_vDend; //*chComps->vDend;
    mod_prec prevV_soma = *chComps_vSoma; //*chComps->vSoma;
    mod_prec q = *chComps_q;              //*chComps->q;
    mod_prec r = *chComps_r;              //*chComps->r;
    mod_prec s = *chComps_s;              //*chComps->s; */

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
    CompRet ret;
    ret.ret1 = DELTA * dVd_dt + prevV_dend; //*chComps->newVDend
    ret.ret2 = I_CaH;                       //*chComps->newI_CaH //This is a state value read in DendCal
    return ret;
}
mod_prec IcNeighbors(global mod_prec *cellCompParamsPtr, mod_prec prevV_dend, StepData step)
{

    int i;
    mod_prec f, V, I_c;
    //printf("Ic[0]= %f\n", neighVdend[0]);

    I_c = 0;
    for (i = 0; i < 8; i++)
    {
        //printf("%d prevdend: %0.10lf, neighVdend: %0.10lf\n",i, prevV_dend, *neighVdend );
        V = prevV_dend - cellCompParamsPtr[step.prevCellIdx + STATE_SIZE + i];
        f = 0.8 * exp(-1 * pow(V, 2) / 100) + 0.2; // SCHWEIGHOFER 2004 VERSION
        I_c = I_c + (CONDUCTANCE * f * V);
    }
    //printf("ja hallo hier is IC, met wie spreek ik: %0.10lf\n", I_c);
    return I_c;
}

void CompSoma(global mod_prec *cellCompParamsPtr, StepData step)
{
    //define variables
    mod_prec vSoma; //*chComps->vSoma;
    mod_prec prevComp1;
    mod_prec prevComp2;
    mod_prec g_CaL; //*chComps->g_CaL;
    mod_prec vDend; //*chComps->vDend;
    mod_prec vAxon; //*chComps->vAxon;
    mod_prec k;     //*chComps->k;
    mod_prec l;     //*chComps->l;
    mod_prec m;     //*chComps->m;
    mod_prec h;     //*chComps->h;
    mod_prec n;     //*chComps->n;
    mod_prec x_s;   //*chComps->x_s;
    CompRet retVals;

    // update somatic components
    // SCHWEIGHOFER:

    // SomaCalcium
    vSoma = cellCompParamsPtr[step.prevCellIdx + SOMA_V];
    prevComp1 = cellCompParamsPtr[step.prevCellIdx + SOMA_CK];
    prevComp2 = cellCompParamsPtr[step.prevCellIdx + SOMA_CL];
    retVals = SomaCalcium(vSoma, prevComp1, prevComp2);
    cellCompParamsPtr[step.newCellIdx + SOMA_CK] = retVals.ret1;
    cellCompParamsPtr[step.newCellIdx + SOMA_CL] = retVals.ret2;

    // SomaSodium
    // vSoma = cellCompParamsPtr[step.prevCellIdx + SOMA_V];
    prevComp1 = cellCompParamsPtr[step.prevCellIdx + SOMA_SM];
    prevComp2 = cellCompParamsPtr[step.prevCellIdx + SOMA_SH];
    retVals = SomaSodium(vSoma, prevComp1, prevComp2);
    cellCompParamsPtr[step.newCellIdx + SOMA_SM] = retVals.ret1;
    cellCompParamsPtr[step.newCellIdx + SOMA_SH] = retVals.ret2;

    // SomaPotassium
    // vSoma = cellCompParamsPtr[step.prevCellIdx + SOMA_V];
    prevComp1 = cellCompParamsPtr[step.prevCellIdx + SOMA_PN];
    prevComp2 = cellCompParamsPtr[step.prevCellIdx + SOMA_PP];
    retVals = SomaPotassium(vSoma, prevComp1, prevComp2);
    cellCompParamsPtr[step.newCellIdx + SOMA_PN] = retVals.ret1;
    cellCompParamsPtr[step.newCellIdx + SOMA_PP] = retVals.ret2;

    // SomaPotassiumX
    // vSoma = cellCompParamsPtr[step.prevCellIdx + SOMA_V];
    prevComp1 = cellCompParamsPtr[step.prevCellIdx + SOMA_PXS];
    cellCompParamsPtr[step.newCellIdx + SOMA_PXS] = SomaPotassiumX(vSoma, prevComp1);

    // SomaCurrVolt
    g_CaL = cellCompParamsPtr[step.prevCellIdx + SOMA_G]; //*chComps->g_CaL;
    vDend = cellCompParamsPtr[step.prevCellIdx + DEND_V]; //*chComps->vDend;
    // vSoma = cellCompParamsPtr[step.prevCellIdx + SOMA_V]; //*chComps->vSoma;
    vAxon = cellCompParamsPtr[step.prevCellIdx + AXON_V]; //*chComps->vAxon;
    k = cellCompParamsPtr[step.newCellIdx + SOMA_CK];     //*chComps->k;
    l = cellCompParamsPtr[step.newCellIdx + SOMA_CL];     //*chComps->l;
    m = cellCompParamsPtr[step.newCellIdx + SOMA_SM];     //*chComps->m;
    h = cellCompParamsPtr[step.newCellIdx + SOMA_SH];     //*chComps->h;
    n = cellCompParamsPtr[step.newCellIdx + SOMA_PN];     //*chComps->n;
    x_s = cellCompParamsPtr[step.newCellIdx + SOMA_PXS];  //*chComps->x_s;
    cellCompParamsPtr[step.newCellIdx + SOMA_V] = SomaCurrVolt(g_CaL,
                                                               vDend,
                                                               vSoma,
                                                               vAxon,
                                                               k,
                                                               l,
                                                               m,
                                                               h,
                                                               n,
                                                               x_s);

    return;
}

CompRet SomaCalcium(mod_prec vSoma, mod_prec prevComp1, mod_prec prevComp2)
{

    mod_prec k_inf, l_inf, tau_k, tau_l, dk_dt, dl_dt, k_local, l_local;

    //Get inputs
    mod_prec prevV_soma = vSoma;        //*chPrms->v;
    mod_prec prevCalcium_k = prevComp1; //*chPrms->prevComp1;
    mod_prec prevCalcium_l = prevComp2; //*chPrms->prevComp2;

    k_inf = (1 / (1 + exp(-1 * (prevV_soma + 61) / 4.2)));
    l_inf = (1 / (1 + exp((prevV_soma + 85.5) / 8.5)));
    tau_k = 1;
    tau_l = ((20 * exp((prevV_soma + 160) / 30) / (1 + exp((prevV_soma + 84) / 7.3))) + 35);
    dk_dt = (k_inf - prevCalcium_k) / tau_k;
    dl_dt = (l_inf - prevCalcium_l) / tau_l;
    k_local = DELTA * dk_dt + prevCalcium_k;
    l_local = DELTA * dl_dt + prevCalcium_l;
    //Put result
    CompRet retVals;
    retVals.ret1 = k_local; //*chPrms->newComp1
    retVals.ret2 = l_local; //*chPrms->newComp2
    return retVals;
}

CompRet SomaSodium(mod_prec vSoma, mod_prec prevComp1, mod_prec prevComp2)
{

    mod_prec m_inf, h_inf, tau_h, dh_dt, m_local, h_local;

    //Get inputs
    mod_prec prevV_soma = vSoma; //*chPrms->v;
    //mod_prec prevSodium_m = *chPrms->prevComp1;
    mod_prec prevSodium_h = prevComp2; //*chPrms->prevComp2;

    // RAT THALAMOCORTICAL SODIUM:
    m_inf = 1 / (1 + (exp((-30 - prevV_soma) / 5.5)));
    h_inf = 1 / (1 + (exp((-70 - prevV_soma) / -5.8)));
    tau_h = 3 * exp((-40 - prevV_soma) / 33);
    dh_dt = (h_inf - prevSodium_h) / tau_h;
    m_local = m_inf;
    h_local = prevSodium_h + DELTA * dh_dt;
    //Put result
    CompRet retVals;
    retVals.ret1 = m_local; //*chPrms->newComp1
    retVals.ret2 = h_local; //*chPrms->newComp2
    return retVals;
}

CompRet SomaPotassium(mod_prec vSoma, mod_prec prevComp1, mod_prec prevComp2)
{

    mod_prec n_inf, p_inf, tau_n, tau_p, dn_dt, dp_dt, n_local, p_local;

    //Get inputs
    mod_prec prevV_soma = vSoma;          //*chPrms->v;
    mod_prec prevPotassium_n = prevComp1; //*chPrms->prevComp1;
    mod_prec prevPotassium_p = prevComp2; //*chPrms->prevComp2;

    // NEOCORTICAL
    n_inf = 1 / (1 + exp((-3 - prevV_soma) / 10));
    p_inf = 1 / (1 + exp((-51 - prevV_soma) / -12));
    tau_n = 5 + (47 * exp(-(-50 - prevV_soma) / 900));
    tau_p = tau_n;
    dn_dt = (n_inf - prevPotassium_n) / tau_n;
    dp_dt = (p_inf - prevPotassium_p) / tau_p;
    n_local = DELTA * dn_dt + prevPotassium_n;
    p_local = DELTA * dp_dt + prevPotassium_p;
    //Put result
    CompRet retVals;
    retVals.ret1 = n_local; //*chPrms->newComp1
    retVals.ret2 = p_local; //*chPrms->newComp2
    return retVals;
}

mod_prec SomaPotassiumX(mod_prec vSoma, mod_prec prevComp1)
{

    mod_prec alpha_x_s, beta_x_s, x_inf_s, tau_x_s, dx_dt_s, x_s_local;

    //Get inputs
    mod_prec prevV_soma = vSoma;            //*chPrms->v;
    mod_prec prevPotassium_x_s = prevComp1; //*chPrms->prevComp1;

    // Voltage-dependent (fast) potassium
    alpha_x_s = 0.13 * (prevV_soma + 25) / (1 - exp(-(prevV_soma + 25) / 10));
    beta_x_s = 1.69 * exp(-0.0125 * (prevV_soma + 35));
    x_inf_s = alpha_x_s / (alpha_x_s + beta_x_s);
    tau_x_s = 1 / (alpha_x_s + beta_x_s);
    dx_dt_s = (x_inf_s - prevPotassium_x_s) / tau_x_s;
    x_s_local = 0.05 * dx_dt_s + prevPotassium_x_s;
    //Put result
    return x_s_local; //*chPrms->newComp1
}
mod_prec SomaCurrVolt(mod_prec g_CaL,
                      mod_prec prevV_dend,
                      mod_prec prevV_soma,
                      mod_prec prevV_axon,
                      mod_prec k,
                      mod_prec l,
                      mod_prec m,
                      mod_prec h,
                      mod_prec n,
                      mod_prec x_s)
{
    //Local variables
    mod_prec I_ds, I_CaL, I_Na_s, I_ls, I_Kdr_s, I_K_s, I_as, dVs_dt;

    //Get inputs
    /* mod_prec g_CaL = *chComps_g_CaL;      //*chComps->g_CaL;
    mod_prec prevV_dend = *chComps_vDend; //*chComps->vDend;
    mod_prec prevV_soma = *chComps_vSoma; //*chComps->vSoma;
    mod_prec prevV_axon = *chComps_vAxon; //*chComps->vAxon;
    mod_prec k = *chComps_k;              //*chComps->k;
    mod_prec l = *chComps_l;              //*chComps->l;
    mod_prec m = *chComps_m;              //*chComps->m;
    mod_prec h = *chComps_h;              //*chComps->h;
    mod_prec n = *chComps_n;              //*chComps->n;
    mod_prec x_s = *chComps_x_s;          //*chComps->x_s; */

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
    return (DELTA * dVs_dt + prevV_soma); // *chComps->newVSoma
}

void CompAxon(global mod_prec *cellCompParamsPtr, StepData step)
{

    // update somatic components
    // SCHWEIGHOFER:
    mod_prec vAxon;
    mod_prec prevComp1;
    mod_prec prevComp2;
    CompRet retVals;
    mod_prec vSoma;
    mod_prec m_a;
    mod_prec h_a;
    mod_prec x_a;

    // AxonSodium
    vAxon = cellCompParamsPtr[step.prevCellIdx + AXON_V];
    prevComp1 = cellCompParamsPtr[step.prevCellIdx + AXON_SH];
    retVals = AxonSodium(vAxon, prevComp1);
    cellCompParamsPtr[step.newCellIdx + AXON_SH] = retVals.ret1;
    cellCompParamsPtr[step.newCellIdx + AXON_SM] = retVals.ret2;

    // AxonPotassium
    // vAxon = cellCompParamsPtr[step.prevCellIdx + AXON_V];
    prevComp1 = cellCompParamsPtr[step.prevCellIdx + AXON_P];
    cellCompParamsPtr[step.newCellIdx + AXON_P] = AxonPotassium(vAxon, prevComp1);

    // AxonCurrVolt();
    vSoma = cellCompParamsPtr[step.prevCellIdx + SOMA_V]; //&cellCompParamsPtr->prevCellState->soma.V_soma;
    vAxon = cellCompParamsPtr[step.prevCellIdx + AXON_V]; //&cellCompParamsPtr->prevCellState->axon.V_axon;
    m_a = cellCompParamsPtr[step.newCellIdx + AXON_SM];   //&cellCompParamsPtr->newCellState->axon.Sodium_m_a;
    h_a = cellCompParamsPtr[step.newCellIdx + AXON_SH];   //&cellCompParamsPtr->newCellState->axon.Sodium_h_a;
    x_a = cellCompParamsPtr[step.newCellIdx + AXON_P];    //&cellCompParamsPtr->newCellState->axon.Potassium_x_a;
    cellCompParamsPtr[step.newCellIdx + AXON_V] = AxonCurrVolt(vSoma, vAxon, m_a, h_a, x_a);

    return;
}

CompRet AxonSodium(mod_prec vAxon, mod_prec prevComp1)
{

    mod_prec m_inf_a, h_inf_a, tau_h_a, dh_dt_a, m_a_local, h_a_local;

    //Get inputs
    mod_prec prevV_axon = vAxon;         //*chPrms->v;
    mod_prec prevSodium_h_a = prevComp1; //*chPrms->prevComp1;

    // Update axonal Na components
    // NOTE: current has shortened inactivation to account for high
    // firing frequencies in axon hillock
    m_inf_a = 1 / (1 + (exp((-30 - prevV_axon) / 5.5)));
    h_inf_a = 1 / (1 + (exp((-60 - prevV_axon) / (-5.8))));
    tau_h_a = 1.5 * exp((-40 - prevV_axon) / 33);
    dh_dt_a = (h_inf_a - prevSodium_h_a) / tau_h_a;
    m_a_local = m_inf_a;
    h_a_local = prevSodium_h_a + DELTA * dh_dt_a;
    //Put result
    CompRet retVals;
    retVals.ret1 = h_a_local; //*chPrms->newComp1
    retVals.ret2 = m_a_local; //*chPrms->newComp2
    return retVals;
}

mod_prec AxonPotassium(mod_prec vAxon, mod_prec prevComp1)
{

    mod_prec alpha_x_a, beta_x_a, x_inf_a, tau_x_a, dx_dt_a, x_a_local;

    //Get inputs
    mod_prec prevV_axon = vAxon;            //*chPrms->v;
    mod_prec prevPotassium_x_a = prevComp1; //*chPrms->prevComp1;

    // D'ANGELO 2001 -- Voltage-dependent potassium
    alpha_x_a = 0.13 * (prevV_axon + 25) / (1 - exp(-(prevV_axon + 25) / 10));
    beta_x_a = 1.69 * exp(-0.0125 * (prevV_axon + 35));
    x_inf_a = alpha_x_a / (alpha_x_a + beta_x_a);
    tau_x_a = 1 / (alpha_x_a + beta_x_a);
    dx_dt_a = (x_inf_a - prevPotassium_x_a) / tau_x_a;
    x_a_local = 0.05 * dx_dt_a + prevPotassium_x_a;
    //Put result
    return x_a_local; //*chPrms->newComp1
}

mod_prec AxonCurrVolt(mod_prec prevV_soma,
                      mod_prec prevV_axon,
                      mod_prec m_a,
                      mod_prec h_a,
                      mod_prec x_a)
{

    //Local variable
    mod_prec I_Na_a, I_la, I_sa, I_K_a, dVa_dt;

    //Get inputs
    /* mod_prec prevV_soma = *chComps_vSoma; //*chComps->vSoma;
    mod_prec prevV_axon = *chComps_vAxon; //*chComps->vAxon;
    mod_prec m_a = *chComps_m_a;          //*chComps->m_a;
    mod_prec h_a = *chComps_h_a;          //*chComps->h_a;
    mod_prec x_a = *chComps_x_a;          //*chComps->x_a; */

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
    return (DELTA * dVa_dt + prevV_axon); //*chComps->newVAxon
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
    int x = get_global_id(0);
    int y = get_global_id(1);

    // TODO smarter way of copying data? memcpy not supported in opencl
    for (int idx = 0; idx < STATE_SIZE; idx++)
    {
        // Previous cell state
        cellCompParamsPtr[(y * IO_NETWORK_DIM1 + x) * LOCAL_PARAM_SIZE + idx] = cellStatePtr[(i % 2) * IO_NETWORK_SIZE * STATE_SIZE + (y * IO_NETWORK_DIM1 + x) * STATE_SIZE + idx];
        // Next cell state
        cellCompParamsPtr[(y * IO_NETWORK_DIM1 + x) * LOCAL_PARAM_SIZE + PARAM_SIZE + idx] = cellStatePtr[((i % 2) ^ 1) * IO_NETWORK_SIZE * STATE_SIZE + (y * IO_NETWORK_DIM1 + x) * STATE_SIZE + idx];
    }
    // Previous cell state in indices 0-27
    // Next cell state in indices 28-54

    //Step data that needs to go with the ptr to access correct indices
    StepData step;
    step.iApp = iApp;
    step.x = x;
    step.y = y;
    step.prevCellIdx = (y * IO_NETWORK_DIM1 + x) * LOCAL_PARAM_SIZE;
    step.newCellIdx = step.prevCellIdx + PARAM_SIZE;

    ComputeOneCell(cellCompParamsPtr, step);
}