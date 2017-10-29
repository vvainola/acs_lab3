#include "kernel.h"

/**
Input: cellCompParamsPtr, cellStatePtr, i
cellCompParamsPtr: Array of struct which stores values of neighbours for each
cell.
cellStatePtr: Array with values for each cell.
i: current simulation step


Retreive the voltage of the dendrite (V_dend) from each neighbour
**/

constant sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                             CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

__kernel void neighbour_kernel(read_only image2d_t t_cellVDendPtr,
                               global mod_prec *cellCompParamsPtr, uint i )
{
    int n, p, q;

    n = 0;
    int y = get_global_id(0);
    int x = get_global_id(1);
    // int offset = (y * IO_NETWORK_DIM1 + x) * LOCAL_PARAM_SIZE + STATE_SIZE;

    for (p = x - 1; p <= x + 1; p++)
    {
        for (q = y - 1; q <= y + 1; q++)
        {

            double2 temp =
                as_double2(read_imageui(t_cellVDendPtr, sampler, (int2)(p, q)));
            cellCompParamsPtr[(y * IO_NETWORK_DIM1 + x) * LOCAL_PARAM_SIZE +
                              STATE_SIZE + n++] = temp.lo;

            if (p == x && q == y)
                n = n - 1;
        }
    }
}
