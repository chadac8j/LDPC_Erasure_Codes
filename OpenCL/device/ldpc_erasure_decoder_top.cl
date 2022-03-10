/******************************************************************************
*  @file    ldpc_erasure_decoder_top.cl
*  @author  Chad Cole
*  @date    11/2/2021
*  @version 1.0
*
*  @brief Device side OpenCL implemenation of LDPC Decoder.
*                                                            
*                       -------------         
*   ----------------    |           |         -------------- 
*   |              |    |           |         |            |  
*   |              |    |  LDPC     |  parity | wr to mem  | 
*   |    data_in   |--->|  Decoder  |-------->|            |      
*   |  hard code   |    |           |         |            |
*   |   all 0's cw |    |           |         --------------
*   | also send rand   
*   ----------------    |           |
*                       |           |         
*                       -------------         
*                                            
*         
*
*  @section DESCRIPTION
*
*  This is the top level device side interface for a FPGA implementation
*  of a LDPC erasure decoder. It uses channels to 
*  stitch together three OpenCL kernels.  One for reading data from memory, 
*  one for encoding the message, the third for writing the parity bits backIt
*  to memory.  
*
*******************************************************************************/

#include "threefry.h"
#include "LDPC_Vlist_data.h"

#pragma OPENCL EXTENSION cl_intel_channels : enable

#define SYM_LEN 128
//#define RS_MULTIPLIER 8

typedef struct {
	unsigned long symbol[SYM_LEN];
	unsigned char is_erasure;
} symbol_type;

typedef struct {
	int num_LDPC_errors;
	int num_RS_errors;
} error_type;

// Channel declarations
channel symbol_type LDPC_DEC_DIN          __attribute__((depth(1024)));
channel symbol_type LDPC_DEC_DOUT         __attribute__((depth(1024)));

channel error_type ERROR_STAT         __attribute__((depth(1024)));

__kernel 
void data_in(
					__global  symbol_type* data_in, 
					const unsigned short nldpc,
					const int seed,
					const int PER_numerator_div_64,
					const int code_ind,
					const long numFrames  // number of frames for throughput testing
					) 
{
//	unsigned useed = 9999; //maybe make this a function of system time so that it changes each run
	unsigned useed = seed; //get this from host where it is randomly created
	unsigned tid = 1;
    int k_ldpc = ldpc_params[code_ind][1];  // (n_ldpc - k_ldpc) here must match the number of rows in parity_check_mat_Vlist below
    int n_ldpc = ldpc_params[code_ind][0];
		
	printf("Enter data_in,  nldpc=%d\n", nldpc);
	threefry4x32_key_t key = {{tid, useed}};
    threefry4x32_ctr_t c = {{}}; // start counter from 0
	
	symbol_type temp_dIn;
//	#pragma unroll
	for (int l=0; l<SYM_LEN; l++)
	{
		temp_dIn.symbol[l] = 0x0;
	}
	
	for(long itr=0; itr<numFrames; itr++)  //this outer itr loop is for throughput testing
	{	
		symbol_type dIn;
	
		#pragma unroll 1
		for (int k=0; k<n_ldpc; k++)
		{
//			dIn = temp_dIn;
			union {
				threefry4x32_ctr_t c;
				int4 i;
			} u;
			c.v[0]++;
			u.c = threefry4x32(c, key);
			long rv = u.i.x; //, y1 = u.i.y;
			//long x2 = u.i.z, y2 = u.i.w;
			//printf("RNG produces =%lu\n", rv);

//			if ((rv & 0x000000000000003F) < 16) {  // 25% packet error rate
//			if ((rv & 0x000000000000003F) < 20) {  // 31.25% packet error rate
//			if ((rv & 0x000000000000003F) < 24) {  // 37.5% packet error rate
			if ((rv & 0x000000000000003F) < PER_numerator_div_64) {  // (PER_numerator_div_64/64)*100% packet error rate
				dIn.is_erasure = 1; 
			}
			else {
				dIn.is_erasure = 0;
			}
		
			 // need RNG in kernel?
			write_channel_intel(LDPC_DEC_DIN, dIn);	
			//printf("%d\n", k);
				
		}
	}
		

}


__kernel 
void data_out(	global symbol_type *data_out,
				const int code_ind,
				const long numFrames  // number of frames for throughput testing
			)
{
    int k_ldpc = ldpc_params[code_ind][1];  // (n_ldpc - k_ldpc) here must match the number of rows in parity_check_mat_Vlist below
    int n_ldpc = ldpc_params[code_ind][0];
	int RS_multiplier = n_ldpc/ldpc_params[code_ind][4];
	symbol_type dout;
	error_type error_out;
	int num_frame_errors = 0;
	int num_RS_frame_errors = 0;

	for(long itr=0; itr<numFrames; itr++)  //this outer itr loop is for throughput testing
	{
		//printf("itr %d : \n", itr);
		/*for(unsigned short i=0; i<kldpc; i++)
		{
			dout = read_channel_intel(LDPC_DEC_DOUT);
			data_out[i] = dout;
			//printf("row %d : ", i);
			for (int l=0; l<SYM_LEN; l++)
			{
				//printf(" %lu", dout.symbol[l]);
			}
			//printf("\n");
		} */
		error_out = read_channel_intel(ERROR_STAT);
		
	}
	num_frame_errors = error_out.num_LDPC_errors;
	num_RS_frame_errors = error_out.num_RS_errors;
	printf("In data_out, frame error rate is: %f, RS FER=%f\n", (float)num_frame_errors/(float)numFrames, (float)num_RS_frame_errors/(RS_multiplier*(float)numFrames));

}

// Include the datapath kernels
#include "ldpc_erasure_decoder.cl"


