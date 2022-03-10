/******************************************************************************
*  @file    ldpc_erasure_encoder_top.cl
*  @author  Chad Cole
*  @date    10/25/2021
*  @version 1.0
*
*  @brief Device side OpenCL implemenation of LDPC Encoder.
*                                                            
*                       -------------         
*   ----------------    |           |         -------------- 
*   |              |    |           |         |            |  
*   |              |    |  LDPC     |  parity | wr to mem  | 
*   |    data_in   |--->|  Encoder  |-------->|            |      
*   |  rd from mem |    |           |         |            |
*   |              |    |           |         --------------
*   ----------------    |           |
*                       |           |         
*                       -------------         
*                                            
*         
*
*  @section DESCRIPTION
*
*  This is the top level device side interface for a FPGA implemetation 
*  of a LDPC erasure codec. This kernel generates random data for the encoder.
*
*******************************************************************************/

#pragma OPENCL EXTENSION cl_intel_channels : enable

#define SYM_LEN 128

typedef struct {
	unsigned long symbol[SYM_LEN];
	unsigned char is_erasure;
} symbol_type;

// Channel declarations
channel symbol_type LDPC_ENC_DIN          __attribute__((depth(1024)));
channel symbol_type LDPC_ENC_DOUT         __attribute__((depth(1024)));


__kernel 
void data_in(
					__global  symbol_type* data_in, 
					const unsigned int numItr,
					const unsigned short kldpc,
					const unsigned short nldpc
					) 
// 					__constant unsigned short * restrict parity_addr_table,

{
	
	printf("k=%d, n=%d\n", kldpc, nldpc);
	for(uint itr=0; itr<numItr; itr++)
	{	
		symbol_type dIn;
	
		//char data[360];
		#pragma unroll 1
		for (int k=0; k<kldpc; k++)
		{
			dIn = data_in[k];
			write_channel_intel(LDPC_ENC_DIN, dIn);	
			//printf("%d\n", k);
				
		}
	}
		

}


__kernel 
void data_out(	global symbol_type *data_out,
				const unsigned int   numItr,
				const unsigned short nldpc
			)
//				const unsigned short kldpc,

{
	symbol_type dout;
	for(uint itr=0; itr<numItr; itr++)  //not sure what this outer itr loop is for
	{
		for(unsigned short i=0; i<nldpc; i++)
		{
			dout = read_channel_intel(LDPC_ENC_DOUT); // & 0x1;
			data_out[i] = dout;
			//printf("%d %d\n", i, dout);
		}
	}	
}

// Include the datapath kernels
#include "ldpc_erasure_encoder.cl"


