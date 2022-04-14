/******************************************************************************
*  @file    ldpc_erasure_encoder.cl
*  @author  Chad Cole
*  @date    1/5/2022
*  @version 1.0
*
*  @brief Device side OpenCL implemenation of LDPC Encoder.  Gets data from VITA and gives it to UDP
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

#define RTP_DEBUG 1

#define SYM_LEN 256 // overprovision memory and then only index into the data that is actually populated
typedef struct {
 	unsigned long symbol[SYM_LEN];
 	unsigned char is_erasure;
} symbol_type;

// typedef struct {
// 	unsigned long data_8bytes;
// 	unsigned char is_SOF;
// } VITA_data_type;

__kernel
void ldpc_erasure_encoder
(
				uint streamID,             // stream identifier
				ulong dataPerContext,      // number of packets between a context packet
				unsigned char testModeEn,  // Packet type selector
				uint disableContextPackets // 1: disable context packets
)
{
	
	unsigned long pktClassCode=0x000a;	
	unsigned long dataPerContextCnt = 0;
	
	uint din;
	ulong dout;
		
	udpIpTxHeader udpTxHeader;	
	
	udpTxHeader.destPort=4991;
	udpTxHeader.srcPort=4991;
	
	ulong dataOut;
	ulong packetLen;
	
	// FEC Header parameters
	unsigned long FECClassCode=0x0001;	// First byte of FEC header contains code info; we'll say our LDPC code is type 1
	unsigned long blockNum = 0;  // this modulo 256 counter value keeps track of codeblock indexing
	unsigned long symbolNum = 0; // this modulo 2^16 counter value keeps track of symbol indexing within a codeblock
	
	symbol_type codeword_first[n_ldpc];
	symbol_type codeword_second[n_ldpc];
	int is_first_codeword = 1;
	int num_longs_used = 0;
	int FEC_header_length_in_words = 2;
	
	for (int ii=0; ii<n_ldpc; ii++)
	{  
		for (int l=0; l<num_longs_used; l++)
		{
			codeword_first[ii].symbol[l] = 0; 
			codeword_second[ii].symbol[l] = 0; 
		}
		codeword_first[ii].is_erasure = 0; 
		codeword_second[ii].is_erasure = 0; 
	}

	int word_cnt = 0; // keeps track of the index within symbol_type symbol structure

	while(1){
		
	// k_ldpc and n_ldpc are variables declared in n2000_k1000_no6cycle_ldpc_Vlist.h
		if (symbolNum == (k_ldpc)) {  // Time to do encoding procedure
			for (int repair_sym_ind=symbolNum; repair_sym_ind<n_ldpc; repair_sym_ind++)
			{
				symbol_type parity_accumulator;
				for (int l=0; l<num_longs_used; l++)
				{
					parity_accumulator.symbol[l] = 0;
				}
				#pragma unroll
				for(int ii=0; ii<(parity_check_mat_Vlist[repair_sym_ind-k_ldpc][0]-1); ii++){  // Don't use last non-zero col in H because it is part of the triangle
				//remember Matlab indexing from parity_check_mat_Vlist is 1 based and C is 0 based
					for (int l=0; l<num_longs_used; l++)
					{
						
						if (is_first_codeword == 1) {
							parity_accumulator.symbol[l] = (parity_accumulator.symbol[l]^codeword_first[parity_check_mat_Vlist[repair_sym_ind-k_ldpc][ii+1]-1].symbol[l]);
						}
						else {
							parity_accumulator.symbol[l] = (parity_accumulator.symbol[l]^codeword_second[parity_check_mat_Vlist[repair_sym_ind-k_ldpc][ii+1]-1].symbol[l]);
						}
					}
				}
				
			// Each repair packet will be sent as a separate UDP packet, so need to create a new UDP Header.  Can reuse same header info
				write_channel_intel(USERTXDATAGRAMS_HDR, udpTxHeader);
			// ------------------
			//    FEC Header - Make it one ulong (8 bytes) for ease of implementation.  Can repeat the 4 byte header
			// ------------------
				dout = 0x00000000ffffffff & (((FECClassCode & 0xff)<<24) | ((blockNum & 0xff)<<16) | (repair_sym_ind & 0xffff));
				dataOut = ((dout<<32) & 0xffffffff00000000) | (dout & 0x00000000ffffffff);
				write_channel_intel(USERTXDATAGRAMS, dataOut);
				for (int l=0; l<num_longs_used; l++)
				{
					write_channel_intel(USERTXDATAGRAMS, parity_accumulator.symbol[l]);
#if (RTP_DEBUG == 1)
					//printf("In repair code, repair_sym_ind:codeword segment:data = %d:%d:%lu \n", repair_sym_ind, l, parity_accumulator.symbol[l]);
					if(l%10 == 0) printf("In repair code, repair_sym_ind:codeword segment:data = %d:%d:%lu \n", repair_sym_ind, l, parity_accumulator.symbol[l]);
#endif
				// Need to store repair packet in FEC structure
					if (is_first_codeword == 1) {
						codeword_first[repair_sym_ind].symbol[l] = parity_accumulator.symbol[l];
					}
					else {
						codeword_second[repair_sym_ind].symbol[l] = parity_accumulator.symbol[l];
					}
				}

				
			}
			symbolNum = 0;
			blockNum += 1;
			is_first_codeword = is_first_codeword ? 0:1;  // switch input buffers
		}
		
		// Shouldn't I read the VITA packet length from the VITA packet header?  Instead of hard coding it.
		// The following read should be the first word of VITA packet, which contains header
		din = read_channel_intel(VRT49_EXT_DATA_OUT);
		packetLen  =  (din & 0xffff) + 2;   // UDP payload in words = vrt packetLen + FEC Header (FEC_header_length_in_words = 2 words)
		if (dataPerContextCnt==dataPerContext){
		    dataPerContextCnt = 0;
		    if (testModeEn==1){
				pktClassCode  = 0x0008;
//				packetLen  = 28;
		    }else{
				pktClassCode  = 0x000B;  
//				packetLen  = 367;
		    }
		}else{
			dataPerContextCnt++;
			if (testModeEn==1){
				pktClassCode  = 0x0006;
//				packetLen  = 367;  // UDP payload in words = vrt packetLen + FEC Header (FEC_header_length_in_words = 2 words)
			}else{
				pktClassCode  = 0x000A;  
//				packetLen  = 367;  // UDP payload in words = vrt packetLen + FEC Header (FEC_header_length_in_words = 2 words)
			}
		}
		udpTxHeader.payloadSize=packetLen*4;
		num_longs_used = ceil(((float)packetLen - 2)/2); // this is how long the FEC payload is in 8 byte units
		//printf("num_longs_used %d ", num_longs_used);
		//printf("symbolNum %lu ", symbolNum);
		//printf("packetLen %lu ", packetLen);
		
		//if context packets are turned off, overwrite context packet counter
		if (disableContextPackets==1) dataPerContextCnt=1;

		write_channel_intel(USERTXDATAGRAMS_HDR, udpTxHeader);

			// ------------------
			//    FEC Header - Make it one ulong (8 bytes) for ease of implementation.  Can repeat the 4 byte header to align at 8 bytes
			// ------------------
		dout = 0x00000000ffffffff & (((FECClassCode & 0xff)<<24) | ((blockNum & 0xff)<<16) | (symbolNum & 0xffff));
		dataOut = ((dout<<32) & 0xffffffff00000000) | (dout & 0x00000000ffffffff);
		write_channel_intel(USERTXDATAGRAMS, dataOut);
		
		word_cnt = 0;
//		for(uint i=0; i < packetLen; i++) {	
		for(uint i=1; i < packetLen-2; i++) {	 //index from 1 since we already pulled out VITA header first and subtract 2 for FEC header
			din = read_channel_intel(VRT49_EXT_DATA_OUT);
#if (RTP_DEBUG == 1)
				if(i%(packetLen-3) == 0) printf("In encoder kernel, symbol num=%lu, i=%d   %08X \n", symbolNum, i, din);

				//printf("In encoder kernel, symbol num=%lu, i=%d   %08X ", symbolNum, i, din);
            	//if(i%10 == 0) printf("\n");
#endif
			dout = ((ulong)(din)&0x00000000ffffffff);

			if( (i & 0x1)==0x0 ){
				dataOut = (dout<<32) & 0xffffffff00000000;
			}else{
				dataOut |= dout & 0x00000000ffffffff;
				write_channel_intel(USERTXDATAGRAMS, dataOut);
				// Need to store packet in FEC structure
				if (is_first_codeword == 1) {
					codeword_first[symbolNum].symbol[word_cnt] = dataOut;
				}
				else {
					codeword_second[symbolNum].symbol[word_cnt] = dataOut;				
				}
				word_cnt += 1;
				
#if (RTP_DEBUG == 2)
				printf("%d:%08lX ", i, dataOut);
            	if(i%10 == 0) printf("\n");
#endif
			}
		}
		symbolNum += 1;
	} // END WHILE


}  // end kernel