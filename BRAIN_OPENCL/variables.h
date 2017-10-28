/*** MACROS ***/
#define RAND_INIT 0  // make it zero to facilitate debugging
#define SIMTIME 10 // in ms, for when no input file is provided
//IO network size is IO_NETWORK_DIM1*IO_NETWORK_DIM2
#define IO_NETWORK_DIM1 32
#define IO_NETWORK_DIM2 32
#define IO_NETWORK_SIZE IO_NETWORK_DIM1 *IO_NETWORK_DIM2

#define IAPP_MAX_CHARS 6 //2 integer, the dot, 2 decimals and the delimiter

// Cell properties
#define DELTA 0.05
//Conductance for neighbors' coupling
#define CONDUCTANCE 0.04
// Capacitance
#define C_M 1
// Somatic conductances (mS/cm2)
#define G_NA_S 150  // Na gate conductance (=90 in Schweighofer code, 70 in paper) 120 too little
#define G_KDR_S 9.0 // K delayed rectifier gate conductance (alternative value: 18)
#define G_K_S 5     // Voltage-dependent (fast) potassium
#define G_LS 0.016  // Leak conductance (0.015)
// Dendritic conductances (mS/cm2)
#define G_K_CA 35  // Potassium gate conductance (35)
#define G_CAH 4.5  // High-threshold Ca gate conductance (4.5)
#define G_LD 0.016 // Dendrite leak conductance (0.015)
#define G_H 0.125  // H current gate conductance (1.5) (0.15 in SCHWEIGHOFER 2004)
// Axon hillock conductances (mS/cm2)
#define G_NA_A 240 // Na gate conductance (according to literature: 100 to 200 times as big as somatic conductance)
#define G_NA_R 0   // Na (resurgent) gate conductance
#define G_K_A 20   // K voltage-dependent
#define G_LA 0.016 // Leak conductance
// Cell morphology
#define P1 0.25    // Cell surface ratio soma/dendrite (0.2)
#define P2 0.15    // Cell surface ratio axon(hillock)/soma (0.1)
#define G_INT 0.13 // Cell internal conductance (0.13)
// Reversal potentials
#define V_NA 55  // Na reversal potential (55)
#define V_K -75  // K reversal potential
#define V_CA 120 // Ca reversal potential (120)
#define V_H -43  // H current reversal potential
#define V_L 10   // leak current

#define DEBUG 1
#define EXTRA_TIMING 1

#define WRITE_OUTPUT 1

#define BLOCKSIZEX 32
#define BLOCKSIZEY 1

#define STATE_SIZE 19
//#define CELL_STATE_SIZE 27
#define PARAM_SIZE 27
#define LOCAL_PARAM_SIZE 54
#define STATEADD 8
#define VNEIGHSTARTADD 1
#define PREVSTATESTARTADD 16
#define NEXTSTATESTARTADD 35
#define DEND_V 0
#define DEND_H 1
#define DEND_CAL 2
#define DEND_P 3
#define DEND_I 4
#define DEND_CA2 5
#define SOMA_G 6
#define SOMA_V 7
#define SOMA_SM 8
#define SOMA_SH 9
#define SOMA_CK 10
#define SOMA_CL 11
#define SOMA_PN 12
#define SOMA_PP 13
#define SOMA_PXS 14
#define AXON_V 15
#define AXON_SM 16
#define AXON_SH 17
#define AXON_P 18