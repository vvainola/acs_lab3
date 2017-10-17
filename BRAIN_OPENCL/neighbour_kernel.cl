#include "kernel.h"

/**
Input: cellCompParamsPtr, cellStatePtr, i
cellCompParamsPtr: Array of struct which stores values of neighbours for each cell.
cellStatePtr: Array with values for each cell.
i: current simulation step


Retreive the voltage of the dendrite (V_dend) from each neighbour
**/
__kernel void neighbour_kernel(global mod_prec *cellStatePtr, global mod_prec *cellCompParamsPtr,  uint i){
    int n, p, q;
    n = 0;
    int y = get_global_id(0);
    int x = get_global_id(1);
    for (p = x - 1; p <= x + 1; p++)
    {
        for (q = y - 1; q <= y + 1; q++)
        {
            if(((p!=x)||(q!=y)) && ((p>=0)&&(q>=0)) && ((p<IO_NETWORK_DIM1)&&(q<IO_NETWORK_DIM2))){
                cellCompParamsPtr[(y*IO_NETWORK_DIM1 + x)*LOCAL_PARAM_SIZE + STATE_SIZE + n++] = cellStatePtr[(i%2)*IO_NETWORK_SIZE*STATE_SIZE + (q*IO_NETWORK_DIM1 + p)*STATE_SIZE + DEND_V];
            }else if(p==x && q==y){
                ;   // do nothing, this is the cell itself
            }else{
                //store same V_dend so that Ic becomes zero by the subtraction
                cellCompParamsPtr[(y*IO_NETWORK_DIM1 + x)*LOCAL_PARAM_SIZE + STATE_SIZE + n++] = cellStatePtr[(i%2)*IO_NETWORK_SIZE*STATE_SIZE + (y*IO_NETWORK_DIM1 + x)*STATE_SIZE + DEND_V];
            }
        }
    }

    // DEBUG
    
    //if (i < 3)
    //{
        printf("neighbour\n");
        int u;
        for (u = 0; u < 8; u++)
        {
            printf("%f\n", cellCompParamsPtr[(y*IO_NETWORK_DIM1 + x)*LOCAL_PARAM_SIZE + STATE_SIZE + u]);
        }
    //} 

}