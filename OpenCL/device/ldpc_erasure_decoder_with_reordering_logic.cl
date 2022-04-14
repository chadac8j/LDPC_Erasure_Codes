/******************************************************************************
*  @file    ldpc_erasure_decoder_with_reordering_logic.cl
*  @author  Chad Cole
*  @date    4/14/2022
*  @version 1.0
*
*  @brief Device side OpenCL implemenation of LDPC Decoder.  Includes logic to read data from UDP kernel and reorder packets
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
	// I guess have k_ldpc and n_ldpc as variables declared in n2000_k1000_no6cycle_ldpc_Vlist.h?
	// Might want to make them kernel input arguments so that they can be changed in the host code on the fly

__kernel
void ldpc_erasure_decoder(const short num_iter)
{
	// FEC Header parameters
	unsigned long FECClassCode=0x0001;	// First byte of FEC header contains code info; we'll say our LDPC code is type 1
	unsigned long blockNum = 0;  // this modulo 256 counter value keeps track of codeblock indexing
	unsigned long symbolNum = 0; // this modulo 2^16 counter value keeps track of symbol indexing within a codeblock

	symbol_type codeword_first[n_ldpc];
	symbol_type codeword_second[n_ldpc];

	int is_first_codeword = 1;  //since we have 2 codeword buffers, need a flag to designate current one
	int cur_block_num = -1;  // store two codewords during assembly process, where to initialize these?
	int next_block_num = -1;
	int cur_block_num_cnt = 0;  // keep track of how many symbols have arrived for current codewords
	int next_block_num_cnt = 0;
	
	int desired_parity_rx = round((n_ldpc - k_ldpc)*0.8);  // if next block begins to appear, can try to decode with this many received symbols
	int min_parity_rx = round((n_ldpc - k_ldpc)*0.2);  // don't even try to decode without this many received symbols
	
	int num_longs_used = 0;  // set this based on symbol size, such as VITA 49 packet size
	unsigned int num_parity_checks = n_ldpc - k_ldpc;
	int ready_to_decode = 0;
	
	for (int ii=0; ii<n_ldpc; ii++)
	{  
		for (int l=0; l<num_longs_used; l++)
		{
			codeword_first[ii].symbol[l] = 0; 
			codeword_second[ii].symbol[l] = 0; 
		}
		codeword_first[ii].is_erasure = 1; // for decoder, assume all symbols are erased until they are produced in a packet by UDP
		codeword_second[ii].is_erasure = 1; 
	}

	while(1)
	{
	// non-blocking read for header
	bool udprdy = false;
	udpTxHeader = read_channel_nb_intel(USERRXDATAGRAMS_HDR, &udprdy);  
	//strip off UDP header and do something with it
	//mem_fence(CLK_CHANNEL_MEM_FENCE); //not sure what this does?
	if (udprdy == true)
	{
	// Now get FEC header
		udpFECHeaderDataIn = read_channel_intel(USERRXDATAGRAMS);  // get FEC header, should be first ulong in UDP packet
		FECClassCode = ((udpFECHeaderDataIn >> 24) & 0xff);
		blockNum = ((udpFECHeaderDataIn >> 16) & 0xff);
		symbolNum = (udpFECHeaderDataIn & 0xffff);
	}
	// initialize cur_block_num/next_block_num here?  this code should only execute once
	if (cur_block_num == -1) && (next_block_num == -1) {
		cur_block_num = blockNum;
		next_block_num = cur_block_num + 1;
	}
	
	// Grab one FEC symbol - should be the rest of the UDP packet
	for (int ss=0; ss<num_longs_used; ss++)
	{  
		udpDataIn = read_channel_intel(USERRXDATAGRAMS);  //read a ulong at a time
	// need to package these longs into a symbol_type
		if (blockNum == cur_block_num) 
		{
			if (is_first_codeword == 1) {
				codeword_first[symbolNum].symbol[ss] = udpDataIn;
			}
			else {
				codeword_second[symbolNum].symbol[ss] = udpDataIn;				
			}
		}
		elseif (blockNum == next_block_num) // ignore (i.e. drop) all blockNum's that are not current or next
		{
			if (is_first_codeword == 0) {
				codeword_first[symbolNum].symbol[ss] = udpDataIn;
			}
			else {
				codeword_second[symbolNum].symbol[ss] = udpDataIn;				
			}
		}
	}
	if (blockNum == cur_block_num) 
	{
		cur_block_num_cnt += 1;
		if (is_first_codeword == 1) {
			codeword_first[symbolNum].is_erasure = 0; //since we received this packet, it is no longer erased
		else {
			codeword_second[symbolNum].is_erasure = 0;
		}
	}
	elseif (blockNum == next_block_num) // ignore (i.e. drop) all blockNum's that are not current or next
	{
		next_block_num_cnt += 1;
		if (is_first_codeword == 0) {
			codeword_first[symbolNum].is_erasure = 0; //since we received this packet, it is no longer erased
		else {
			codeword_second[symbolNum].is_erasure = 0;
		}
	}
		
// Once enough packets have arrived, begin decoding
// Need logic to determine when to hand off codeword_first or codeword_second to decoder
// This logic contains a number of parameters that will need to be tuned
	if ((cur_block_num_cnt == n_ldpc) || ((cur_block_num_cnt > (k_ldpc + desired_parity_rx)) && (next_block_num_cnt > 10)) || ((cur_block_num_cnt > (k_ldpc + min_parity_rx)) && (next_block_num_cnt > 100)) ){
		ready_to_decode = 1;
	}

// set up iteration loop
	if (ready_to_decode == 1)
		ready_to_decode = 0;
	{
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
			for (int l=0; l<num_longs_used; l++)
			{
			 	parity_accumulator.symbol[l] = 0;
			}
				
			// must accumulate parity symbols according to parity_check_mat_Vlist
									
			for(int ii=0; ii<(parity_check_mat_Vlist[k][0]); ii++)
			{
				//remember Matlab indexing from parity_check_mat_Vlist is 1 based and C is 0 based
				#pragma unroll
				for (int l=0; l<num_longs_used; l++)
				{
					if (is_first_codeword == 1) {
						parity_accumulator.symbol[l] = (parity_accumulator.symbol[l]^codeword_first[parity_check_mat_Vlist[k][ii+1]-1].symbol[l]);
					} else {
						parity_accumulator.symbol[l] = (parity_accumulator.symbol[l]^codeword_second[parity_check_mat_Vlist[k][ii+1]-1].symbol[l]);					
					}
				}
				if ( ((is_first_codeword == 1) && (codeword_first[parity_check_mat_Vlist[k][ii+1]-1].is_erasure == 1))  ||
					((is_first_codeword == 0) && (codeword_second[parity_check_mat_Vlist[k][ii+1]-1].is_erasure == 1)) )
				{
					num_erasures = num_erasures + 1;
					erasure_ind = parity_check_mat_Vlist[k][ii+1]-1;  // store C-based index, not Matlab
				}
			}
			if (num_erasures == 1)  // we can correct this erasure
			{
				if (is_first_codeword == 1) {
					codeword_first[erasure_ind].is_erasure = 0;
				} else {
					codeword_second[erasure_ind].is_erasure = 0;
				}

				#pragma unroll
				for (int l=0; l<num_longs_used; l++)
				{
					if (is_first_codeword == 1) {
						codeword_first[erasure_ind].symbol[l] = parity_accumulator.symbol[l];
					} else {
						codeword_second[erasure_ind].symbol[l] = parity_accumulator.symbol[l];
					}
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
	
	// prepare codeword memory for next incoming codeword
		cur_block_num = next_block_num;  // move up code block
		next_block_num += 1;  // assume sequential block numbers?
		cur_block_num_cnt = next_block_num_cnt;  // keep track of how many symbols have arrived for current codewords
		next_block_num_cnt = 0;
		for (int ii=0; ii<n_ldpc; ii++)
		{  
			for (int l=0; l<num_longs_used; l++)
			{
				if (is_first_codeword == 1)
				{
					codeword_first[ii].symbol[l] = 0; 
				} else {
					codeword_second[ii].symbol[l] = 0; 
				}
			}
			if (is_first_codeword == 1)
			{
				codeword_first[ii].is_erasure = 1; 
			} else {
				codeword_second[ii].is_erasure = 1; 
			}
		}
		is_first_codeword = is_first_codeword ? 0:1;  // switch input buffers
	} // end if decode
	
	
	} // outer while(1)
}



	
		
	


