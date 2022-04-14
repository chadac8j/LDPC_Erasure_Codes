/******************************************************************************
*  @file    ldpc_erasure_decoder.cl
*  @author  Chad Cole
*  @date    11/2/2021
*  @version 1.0
*
*  @brief Device side OpenCL implemenation of LDPC Decoder.
*           
*
*  @section DESCRIPTION
*
*  This program sets up the device side interface for a FPGA implemetation 
*  of a LDPC erasure decoder. This kernel is a single
*  work item, single workgroup kernel. It is intended to be used for a 
*  all supported code rates with either short or normal frame sizes
*
* ASSUMPTIONS: 1) The symbol_type.is_erasure field is properly set to 1 for erased packets
*   and 0 for received packets.  
* 			   2) The symbol_type.symbol array is populated with all 0's for erased symbols since 
*	this field is XORed with the other symbols in a parity check regardless of erasure state
*******************************************************************************/
#include "n2000_k1000_no6cycle_ldpc_Vlist_device.h"

__kernel
void ldpc_erasure_decoder(const short num_iter)
{
	while(1)
	{
	
	// I guess have k_ldpc and n_ldpc as variables declared in n2000_k1000_no6cycle_ldpc_Vlist.h?
	// Might want to make them kernel input arguments so that they can be changed in the host code on the fly
	
	// bring in data from read kernel and store in codeword
	symbol_type codeword[n_ldpc];
	unsigned int num_parity_checks = n_ldpc - k_ldpc;
	for (int ii=0; ii<n_ldpc; ii++)
	{  
		symbol_type dataIn;
		dataIn = read_channel_intel(LDPC_DEC_DIN);
		codeword[ii] = dataIn;
	   	// for (int l=0; l<SYM_LEN; l++)
		// {
		//   	 codeword[ii].symbol[l] = dataIn.symbol[l];
		// }
		// codeword[ii].is_erasure = dataIn.is_erasure;
	}
	
// set up iteration loop
	for (int ii=0; ii<num_iter; ii++)
	{
		
		
		#pragma unroll 1
		for (int k=0; k<num_parity_checks; k++)
		{
// since erasures have been stored as '000..00' they can be XORed to the cumulative XOR value without changing it				
			symbol_type parity_accumulator;
			unsigned char num_erasures = 0; // assume we won't have more than 255 erasures for one parity check and wrap a uchar
			unsigned int erasure_ind = 0;
			#pragma unroll
			for (int l=0; l<SYM_LEN; l++)
			{
			 	parity_accumulator.symbol[l] = 0;
			}
				
			// must accumulate parity symbols according to parity_check_mat_Vlist
									
			for(int ii=0; ii<(parity_check_mat_Vlist[k][0]); ii++)
			{
				//remember Matlab indexing from parity_check_mat_Vlist is 1 based and C is 0 based
				#pragma unroll
				for (int l=0; l<SYM_LEN; l++)
				{
				 	parity_accumulator.symbol[l] = (parity_accumulator.symbol[l]^codeword[parity_check_mat_Vlist[k][ii+1]-1].symbol[l]);
				}
				if (codeword[parity_check_mat_Vlist[k][ii+1]-1].is_erasure == 1)
				{
					num_erasures = num_erasures + codeword[parity_check_mat_Vlist[k][ii+1]-1].is_erasure;
					erasure_ind = parity_check_mat_Vlist[k][ii+1]-1;  // store C-based index, not Matlab
				}
			}
			if (num_erasures == 1)  // we can correct this erasure
			{
				codeword[erasure_ind].is_erasure = 0;
				#pragma unroll
				for (int l=0; l<SYM_LEN; l++)
				{
				 	codeword[erasure_ind].symbol[l] = parity_accumulator.symbol[l];
				}
			}
//	    			printf("%d %lu\n", k, codeword[k].symbol[l]);
		}
	}
	
	// since this is a systematic code, the first k symbols are the source symbols
	// Here we don't have an early stopping criteria and don't verify if all check nodes are satisfied, just return systematic symbols
	for (int ii=0; ii<k_ldpc; ii++)
	{
		symbol_type dataOut;
		dataOut = codeword[ii];
		write_channel_intel(LDPC_DEC_DOUT, dataOut);
	}
	
	} // outer while(1)
}



	
		
	


