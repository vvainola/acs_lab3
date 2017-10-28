#include "kernel.h"

/**
Input: cellCompParamsPtr, cellStatePtr, i 
cellStatePtr: Array with values for each cell.
i: current simulation step  


Retreive the voltage of the dendrite (V_dend) from each neighbour
**/

__kernel void neighbour_kernel(global mod_prec *cellStatePtr, global mod_prec *cellVDendPtr, int i)
{
	int n, p, q;

	n = 0;
	int y = get_global_id(0);
	int x = get_global_id(1);
	int offset = (y * IO_NETWORK_DIM2 + x) * PARAM_SIZE + STATE_SIZE;

	for (p = x - 1; p <= x + 1; p++)
	{
		for (q = y - 1; q <= y + 1; q++)
		{
			if (((p != x) || (q != y)) && ((p >= 0) && (q >= 0)) && ((p < IO_NETWORK_DIM1) && (q < IO_NETWORK_DIM2)))
			{
				cellStatePtr[offset + (n++)] = cellVDendPtr[q * IO_NETWORK_DIM2 + p];
			}
			else if (p == x && q == y)
			{
				; // do nothing, this is the cell itself
			}
			else
			{
				//store same V_dend so that Ic becomes zero by the subtraction
				cellStatePtr[offset + (n++)] = cellVDendPtr[y * IO_NETWORK_DIM2 + x];
			}
		}
	}

	/*
	//Debug print
	if (i < 5 && x == 0 && y == 0)
	{
		for (int idx = 0; idx < PARAM_SIZE; idx++)
		{
			printf("neighbour i=%d idx=%d val=%f\n", i, idx, cellStatePtr[offset - STATE_SIZE + idx]);
		}
	} */
}
