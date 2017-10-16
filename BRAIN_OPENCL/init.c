#include "init.h" 

void mallocCells(cl_mod_prec **cellCompParamsPtr, cl_mod_prec **cellStatePtr){
    int k;
    DEBUG_PRINT(("cellStatePtr: %luB\n", 2*IO_NETWORK_SIZE*STATE_SIZE*sizeof(cl_mod_prec)));
    //Two cell state structs are needed so as to avoid having to synchronize all consumers before they start rewriting the cell state.
    (*cellStatePtr) = malloc(2*IO_NETWORK_SIZE*STATE_SIZE*sizeof(cl_mod_prec));//current and next state
    if((*cellStatePtr)==NULL){
        printf("Error: Couldn't malloc for cellStatePtr\n");
        exit(EXIT_FAILURE);
    }

    DEBUG_PRINT(("cellCompParamsPtr: %luB\n", IO_NETWORK_SIZE*LOCAL_PARAM_SIZE*sizeof(cl_mod_prec)));
    (*cellCompParamsPtr) = malloc(IO_NETWORK_SIZE*LOCAL_PARAM_SIZE*sizeof(cl_mod_prec));
    if((*cellCompParamsPtr) ==NULL){
        printf("Error: Couldn't malloc for cellCompParamsPtr\n");
        exit(EXIT_FAILURE);
    }
}

void InitState(cl_mod_prec *cellStatePtr){
    int i, b;
    cl_mod_prec  cellStateInit[STATE_SIZE];
    //Initial dendritic parameters
    cellStateInit[DEND_V]   = -60;  
    cellStateInit[DEND_H]   = 0.0337836;
    cellStateInit[DEND_CAL] = 0.0112788;
    cellStateInit[DEND_P]   = 0.0049291;
    cellStateInit[DEND_I]   = 0.5;
    cellStateInit[DEND_CA2] = 3.7152;

    cellStateInit[SOMA_G]   = 0.68;
    cellStateInit[SOMA_V]   = -60;
    cellStateInit[SOMA_SM]  = 1.0127807;
    cellStateInit[SOMA_SH]  = 0.3596066;
    cellStateInit[SOMA_CK]  = 0.7423159;
    cellStateInit[SOMA_CL]  = 0.0321349;
    cellStateInit[SOMA_PN]  = 0.2369847;
    cellStateInit[SOMA_PP]  = 0.2369847;
    cellStateInit[SOMA_PXS] = 0.1;

    cellStateInit[AXON_V]   = -60;
    cellStateInit[AXON_SM]  = 0.003596066;
    cellStateInit[AXON_SH]  = 0.9;
    cellStateInit[AXON_P]   = 0.2369847;

    //Copy init sate to all cell states
    for(i=0;i<IO_NETWORK_SIZE;i++){
        for(b=0;b<STATE_SIZE;b++){
            cellStatePtr[i*STATE_SIZE + b] = cellStateInit[b];
        }
    }

    return;
}

void init_g_CaL(cl_mod_prec *cellStatePtr){
    int seedvar, i;
    seedvar = 1;
    for(i=0;i<IO_NETWORK_SIZE;i++){
            srand(seedvar++);   // use this for debugging, now there is difference
            cellStatePtr[(IO_NETWORK_SIZE + i)*STATE_SIZE + SOMA_G] = cellStatePtr[i*STATE_SIZE + SOMA_G] = 0.68;
    }
}
