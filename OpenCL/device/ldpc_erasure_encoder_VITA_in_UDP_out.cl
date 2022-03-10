/******************************************************************************
*  @file    ldpc_erasure_encoder.cl
*  @author  Chad Cole
*  @date    10/25/2021
*  @version 1.0
*
*  @brief Device side OpenCL implemenation of LDPC Encoder.
*           
*
*  @section DESCRIPTION
*
*  This program sets up the device side interface for a FPGA implemetation 
*  of a LDPC encoder used in the DVB-S2 waveform. This kernel is a single
*  work item, single workgroup kernel. It is intended to be used for a 
*  all supported code rates with either short or normal frame sizes
*
*******************************************************************************/
#include "n2000_k1000_no6cycle_ldpc_Vlist_device.h"

// #define SYM_LEN 128
// typedef struct {
// 	unsigned long symbol[SYM_LEN];
// 	unsigned char is_erasure;
// } symbol_type;

// typedef struct {
// 	unsigned long data_8bytes;
// 	unsigned char is_SOF;
// } VITA_data_type;

__kernel
void ldpc_erasure_encoder()
{
	
//while(1)
//{
	// I guess have k_ldpc and n_ldpc as variables declared in n2000_k1000_no6cycle_ldpc_Vlist.h?
	
	// FEC Header parameters
	unsigned char code_type = 0x01; // we'll say our LDPC code is type 1
	unsigned char block_counter = 0;  // this modulo 256 counter value keeps track of codeblock indexing
	unsigned short symbol_counter = 0;  // this modulo 2^16 counter value keeps track of symbol indexing within a codeblock

	while(1)
	{
		// How to signal that a new VITA packet is ready for encoder?  Have a separate channel for FEC header that is called at start of new VITA packet?
		unsigned long dataIn_HDR = read_channel_intel(LDPC_ENC_HDR_DIN);
		// since FEC header is 4 bytes, make it a uint?
		unsigned int dataOut_HDR = ((code_type & 0xFF) << 24) | ((block_counter & 0xFF) << 16) | ((symbol_counter & 0xFFFF) );
		write_channel_intel(LDPC_ENC_HDR_DOUT, dataOut_HDR);
		
		symbol_type codeword[n_ldpc];
		for (int ii=0; ii<n_ldpc; ii++)
		{  
	   	// codeword[ii].symbol = 0; 
	   	// codeword[ii].is_erasure = 0; 
			for (int l=0; l<SYM_LEN; l++)
			{
				codeword[ii].symbol[l] = 0; 
			}
			codeword[ii].is_erasure = 0; 
		}
	
		symbol_counter = 0;
		
		#pragma unroll 1
		for (int k=0; k<n_ldpc; k++)
		{
				
			symbol_type parity_accumulator;
			for (int l=0; l<SYM_LEN; l++)
			{
			 	parity_accumulator.symbol[l] = 0;
			}
				
			//unsigned char __attribute((register)) din[360];
						
			// since this is a systematic code, first k symbols are part of final codeword
			if(k<k_ldpc){
				symbol_type dataIn;
				dataIn = read_channel_intel(LDPC_ENC_DIN);
				write_channel_intel(LDPC_ENC_DOUT, dataIn);
				// codeword[k].symbol = dataIn.symbol;
				for (int l=0; l<SYM_LEN; l++)
				{
				 	codeword[k].symbol[l] = dataIn.symbol[l];
	    			printf("%d %lu\n", k, codeword[k].symbol[l]);
				}
			} else {  // must accumulate parity symbols according to parity_check_mat_Vlist
									
				#pragma unroll
				for(int ii=0; ii<(parity_check_mat_Vlist[k-k_ldpc][0]-1); ii++){  // Don't use last non-zero col in H because it is part of the triangle
				//remember Matlab indexing from parity_check_mat_Vlist is 1 based and C is 0 based
					// parity_accumulator = (parity_accumulator^codeword[parity_check_mat_Vlist[k-k_ldpc][ii+1]-1]);

					for (int l=0; l<SYM_LEN; l++)
					{
					 	parity_accumulator.symbol[l] = (parity_accumulator.symbol[l]^codeword[parity_check_mat_Vlist[k-k_ldpc][ii+1]-1].symbol[l]);
					}
				}
				write_channel_intel(LDPC_ENC_DOUT, parity_accumulator);
				// codeword[k] = parity_accumulator;
				for (int l=0; l<SYM_LEN; l++)
				{
				 	codeword[k].symbol[l] = parity_accumulator.symbol[l];
	    			printf("%d %lu\n", k, codeword[k].symbol[l]);
				}
			}
			symbol_counter += 1;
		}
		block_counter += 1;
	}
//} // end outer main while
}



	
		
	


