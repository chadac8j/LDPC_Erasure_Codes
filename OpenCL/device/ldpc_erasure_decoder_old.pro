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
	int frame_num = 0;  // use this when sending multiple frames through for performance testing
	int num_frame_errors = 0;
	int num_RS_frame_errors = 0;
	int n_RS = 250;
	int k_RS = 125;
	//int RS_multiplier = n_ldpc/n_RS;
	
	while(1)
	{
	
	// I guess have k_ldpc and n_ldpc as variables declared in n2000_k1000_no6cycle_ldpc_Vlist.h?
	// Might want to make them kernel input arguments so that they can be changed in the host code on the fly
	
	// bring in data from read kernel and store in codeword
	symbol_type codeword[n_ldpc];
	unsigned int num_parity_checks = n_ldpc - k_ldpc;
	int num_initial_erasures = 0;
	int num_RS_erasures = 0;
	for (int ii=0; ii<n_ldpc; ii++)
	{  
		symbol_type dataIn;
		dataIn = read_channel_intel(LDPC_DEC_DIN);
		codeword[ii] = dataIn;
		if (codeword[ii].is_erasure == 1) {
			num_initial_erasures += 1;
			num_RS_erasures += 1;
		}
		// For RS comparison performance
		if ((ii+1) % n_RS == 0) {
			if (num_RS_erasures > (n_RS - k_RS)) {
				num_RS_frame_errors += 1;
			}
			num_RS_erasures = 0;
		}
		
	   	// for (int l=0; l<SYM_LEN; l++)
		// {
		//   	 codeword[ii].symbol[l] = dataIn.symbol[l];
		// }
		// codeword[ii].is_erasure = dataIn.is_erasure;
	}
	//printf("In Frame number %d, there are still %d erasures BEFORE decoding.\n", frame_num, num_initial_erasures);
	
// set up iteration loop
//	for (int ii=0; ii<num_iter; ii++)
	int iter_ind = 0;
	int stop_sig = 0;
	while ((iter_ind<num_iter) && (stop_sig == 0))
	{
		
		
		#pragma unroll 4
		#pragma ivdep
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
				num_initial_erasures -= 1;  // decrement this each time we correct an erasure
			}
//	    			printf("%d %lu\n", k, codeword[k].symbol[l]);
		}
		// early stopping criteria
		if (num_initial_erasures == 0) {
			stop_sig = 1;
		}
		iter_ind += 1;
	}
	
	// since this is a systematic code, the first k symbols are the source symbols
	// Here we don't have an early stopping criteria and don't verify if all check nodes are satisfied, just return systematic symbols
	int num_final_erasures = 0;
	for (int ii=0; ii<k_ldpc; ii++)
	{
		if (codeword[ii].is_erasure == 1) {
			num_final_erasures += 1;
		}
		symbol_type dataOut;
		dataOut = codeword[ii];
		//write_channel_intel(LDPC_DEC_DOUT, dataOut);
	}
//	printf("In Frame number %d, there are still %d erasures after decoding.\n", frame_num, num_final_erasures);
	frame_num += 1;
	if (num_final_erasures > 0) {
		num_frame_errors += 1; }
		
	if (frame_num % 10000 == 0) {
		//printf("In Frame# %d, iteration count=%d, frame error rate is: %f, RS FER=%f\n", frame_num, iter_ind, (float)num_frame_errors/(float)frame_num, (float)num_RS_frame_errors/(RS_MULTIPLIER*(float)frame_num));
	}
	error_type error_out;
	error_out.num_LDPC_errors = num_frame_errors;
	error_out.num_RS_errors = num_RS_frame_errors;
	write_channel_intel(ERROR_STAT, error_out);
	
	} // outer while(1)
}



	
		
	


