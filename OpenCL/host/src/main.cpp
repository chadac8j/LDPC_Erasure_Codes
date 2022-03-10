/******************************************************************************
*  @file    main.cpp
*  @author  Chad Cole
*  @date    11/02/2021
*  @version 1.0
*
*  @brief Host side OpenCL implementation of a LDPC Erasure decoder 
*
*  @section DESCRIPTION
*
*  This program sets up the host side interface for a FPGA implementation 
*  of a LDPC Erasure decoder that will decode packet level data for transmission
*  over lossy internet channels.  
*
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include "CL/opencl.h"
#include "AOCLUtils/aocl_utils.h"
#include <malloc.h>
#include "AOCLUtils/options.h"

// These include files contain the parity check matrix in Vlist format
//#include "ldpc_parity_matrix_n1024_k768_.h"
#include "../inc/Main_LDPC_header.h"

#define AOCL_ALIGNMENT 64
#define VERBOSE

using namespace aocl_utils;
using namespace std;


#define EMULATION_PLAT	0
#define HARDWARE_PLAT	1
const int SYM_LEN = 128;

typedef struct {
	unsigned long symbol[SYM_LEN];
	unsigned char is_erasure;
} symbol_type;

enum KERNELS {
	K_DATA_IN,
	K_LDPC_ERASURE_DECODER,
	K_DATA_OUT,
	K_NUM_KERNELS
};

static const char* kernel_names[K_NUM_KERNELS] =
{
	"data_in",
	"ldpc_erasure_decoder",
	"data_out"
};


int k_LEN = 2000; // n,k should be defined in a .h file
int n_LEN = 4000; // initialize these to a maximum value so that memory reserved is enough, but then update these values with command line input

const char *device_kernel;
const char *input_data_file="LDPC_ErasureDecoder_IN_n2000_Shorts_35PercentPER.txt";
const char *output_data_file="LDPC_ErasureDecoder_OUT_k1000_Shorts.txt";

// OpenCL runtime configuration
cl_platform_id platform = NULL;
unsigned num_devices = 0;
cl_device_id device; 
cl_context context = NULL;
cl_command_queue queue[K_NUM_KERNELS];
cl_program program = NULL;
cl_kernel kernel[K_NUM_KERNELS];

cl_mem input_a_buf; 
cl_mem output_buf; 

// Problem data.
void *din_array_ptr = alignedMalloc(n_LEN*sizeof(symbol_type));
//void *din_erasure_ind_array_ptr = alignedMalloc(n_LEN*sizeof(char));
void *dout_vector_ptr = alignedMalloc(k_LEN*sizeof(symbol_type));
symbol_type *din_array = (symbol_type *)din_array_ptr;    //matlab generated input data
//char *din_erasure_ind_array = (char *)din_erasure_ind_array_ptr;    //matlab generated input data
symbol_type *dout_vector = (symbol_type *)dout_vector_ptr;  //matlab generated results used for comparision


void *dout_ptr = alignedMalloc(k_LEN*sizeof(symbol_type));
symbol_type *dout = (symbol_type *)dout_ptr;

char ptype = EMULATION_PLAT;
//char frameType = NORMAL_FRAME;

//short numItr = 10;  // number of iterations in Message Passing Algorithm
short numItr = 50;  // number of iterations in Message Passing Algorithm
long numFrames = 1000000;  // number of iterations for throughput speed tests
int seed = 0;
int PER_numerator_div_64 = 0;
int code_ind = 0;

// Function prototypes
float rand_float();
int read_test_vector_file_noisy_packets(const char *filename, symbol_type *din_array, const int length);
//int read_test_vector_file_noisy_packets(const char *filename, symbol_type *din_array, char *din_erasure_ind_array , const int length);
int read_test_vector_file_packets(const char *filename, symbol_type *din_array, const int length);
// int read_test_vector_file_short(const char *filename, unsigned short *din_array);
bool init_opencl();
void run();
void cleanup();
int verify_output();


//****************************************
// Option Parser structures and routines
//****************************************
struct Arg : public option::Arg
{
	static void printError(const char* msg1, const option::Option& opt, const char* msg2)
	{
		fprintf(stderr, "%s", msg1);
		fwrite(opt.name, opt.namelen, 1, stderr);
		fprintf(stderr, "%s", msg2);
	}
		static option::ArgStatus Unknown(const option::Option& option, bool msg)
	{
		if (msg) printError("Unknown option '", option, "'\n");
		return option::ARG_ILLEGAL;
	}
		static option::ArgStatus Required(const option::Option& option, bool msg)
	{
		if (option.arg != 0)
			return option::ARG_OK;
			if (msg) printError("Option '", option, "' requires an argument\n");
		return option::ARG_ILLEGAL;
	}
		static option::ArgStatus NonEmpty(const option::Option& option, bool msg)
	{
		if (option.arg != 0 && option.arg[0] != 0)
			return option::ARG_OK;
			if (msg) printError("Option '", option, "' requires a non-empty argument\n");
		return option::ARG_ILLEGAL;
	}
		static option::ArgStatus Numeric(const option::Option& option, bool msg)
	{
		char* endptr = 0;
		if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
		if (endptr != option.arg && *endptr == 0)
			return option::ARG_OK;
			if (msg) printError("Option '", option, "' requires a numeric argument\n");
		return option::ARG_ILLEGAL;
	}
};
enum  optionIndex { UNKNOWN, HELP, N_FRAMES, PER, NFRAME, EMODE, HMODE, ITERATIONS, CODE };
const option::Descriptor usage[] = {
	{ UNKNOWN, 0, "", "", Arg::Unknown, "USAGE: example_arg [options]\n\n"
	"Options:" },
	{ HELP, 0, "", "help", Arg::None, "  \t--help  \tPrint usage and exit." },
	{ PER, 0, "p", "PER cases are: ", Arg::Required, " -p <arg>, \t--required=<arg>, (arg/64)*100% packet error rate " },
	{ NFRAME, 0, "n", "specify how many frames to send through simulation", Arg::Required, " -n <arg>, \t--required=<arg>  " },
	{ EMODE, 0, "e", "run in emulation mode", Arg::None, "  -e\t\tRun in emulation mode" },
	{ HMODE, 0, "h", "run on hardware", Arg::None, "  -h\t\tRun on hardware" },
	{ ITERATIONS, 0, "i", "number of iterations", Arg::Required, "  -i <arg>, \t--required=<arg>  \tNumber of DVB-S2 frames to encode" },
	{ CODE, 0, "c", "code type ", Arg::Required, "  -c <arg>, \t--required=<arg>  \tcode rate\n \t\t0 = (2000, 1000),\t\t1 = (2040, 1530),\n \t\t2 = 2 / 5,"
	"\t\t3 = 1 / 2,\n \t\t4 = 3 / 5,\t\t5 = 2 / 3,\n \t\t6 = 3 / 4,\t\t7 = 4 / 5,\n \t\t8 = 5 / 6,\t\t9 = 8 / 9,\n \t\t10 = 9 / 10,"
	"\t\t35 = 22 / 30,\n \t\t36 = 135 / 180,\t\t37 = 140 / 180,\n \t\t38 = 7 / 9,\t\t39 = 154 / 180\n"},
	{ 0, 0, 0, 0, 0, 0 } };
	
/*************************************************************************

	@brief The main function serves as the launch point for initializing
	         an OpenCL application.  The function finds a device, sets up 
			 the context and read/write buffers.  It then launches the
			 kernel and verifes the results against Matlab generated 
			 test vectors.

	@param argc argument count
	@param argv argument varaibles
	@return int if a negative value is returned the function failed

**************************************************************************/
int main(int argc, char **argv) 
{
	int stat = 0;
	bool emulation=0;

	argc -= (argc>0); argv += (argc>0); // skip program name argv[0] if present
	option::Stats stats(usage, argc, argv);

#ifdef __GNUC__
	// GCC supports C99 VLAs for C++ with proper constructor calls.
	option::Option options[stats.options_max], buffer[stats.buffer_max];
#else
	// use calloc() to allocate 0-initialized memory. It's not the same
	// as properly constructed elements, but good enough. Obviously in an
	// ordinary C++ program you'd use new[], but this file demonstrates that
	// TLMC++OP can be used without any dependency on the C++ standard library.
	option::Option* options = (option::Option*)calloc(stats.options_max, sizeof(option::Option));
	option::Option* buffer = (option::Option*)calloc(stats.buffer_max, sizeof(option::Option));
#endif

	option::Parser parse(usage, argc, argv, options, buffer);

	if (parse.error())
		return 1;

	if (options[HELP] || argc == 0)
	{
		int columns = getenv("COLUMNS") ? atoi(getenv("COLUMNS")) : 80;
		option::printUsage(fwrite, stdout, usage, columns);
		return 0;
	}

	for (int i = 0; i < parse.optionsCount(); ++i)
	{
		option::Option& opt = buffer[i];
		switch (opt.index())
		{
		case HELP:
		case PER:
		 	PER_numerator_div_64 = atoi(opt.arg);
		 	break;
		case NFRAME:
			numFrames = atoi(opt.arg);
			break;
		case CODE:
			code_ind = atoi(opt.arg);
			break;
		case ITERATIONS:
			numItr = atoi(opt.arg);
			break;
		case EMODE:
			ptype = EMULATION_PLAT;
			break;
		case HMODE:
			ptype = HARDWARE_PLAT;
			break;
		case UNKNOWN:
			// not possible because Arg::Unknown returns ARG_ILLEGAL
			// which aborts the parse with an error
			break;
		}
	}

	for (int i = 0; i < parse.nonOptionsCount(); ++i)
		fprintf(stdout, "Non-option argument #%d is %s\n", i, parse.nonOption(i));


	if (ptype == EMULATION_PLAT)
		device_kernel = "ldpc_erasure_decoder_em";
	else
		device_kernel = "ldpc_erasure_decoder";

// Use data input from command line:
    k_LEN = ldpc_params[code_ind][1];  // (n_ldpc - k_ldpc) here must match the number of rows in parity_check_mat_Vlist below
    n_LEN = ldpc_params[code_ind][0];

	float codeRate = float(k_LEN)/float(n_LEN);

  //**********************
  // Initialize OpenCL.
  //**********************
  if(!init_opencl()) {
    return -1;
  }

  //*********************************************
  //get input and output test vectors from files
  //*********************************************
  // if (read_test_vector_file_short(input_data_file, din_array) < 0)
  if (read_test_vector_file_noisy_packets(input_data_file, din_array, SYM_LEN) < 0)
  {
	  printf("Error opening input vector file\n");
	  return -1;
  }
  // if (read_test_vector_file_short(output_data_file, dout_vector) < 0)
  if (read_test_vector_file_packets(output_data_file, dout_vector, SYM_LEN) < 0)
  {
	  printf("Error opening output vector file\n");
	  return -1;
  }

  // for(int i=0; i<k_LEN; i++){
	  // printf("%d\n", din_array[i]);
  // }

  //**********************
  // Run the kernel.
  //**********************
  run();

  //**********************
  // Verify results.
  //**********************
  if (verify_output() >= 0)
  {
	  printf("PASSED\n");
  }
  else
  {
	  printf("FAILED!\n");
  }

  //******************************
  // Free the resources allocated
  //******************************
  cleanup();

  return 0;
}



/************************************************************************
	
	@brief The read_test_vector_file function reads a text file of input 
		parameters intended to be used as inputs to an OpenCL module
	
	@param filename is a const char pointer containg the name of the 
				file to be parsed
	@param din_array is a char array that holds the fiels parsed data
	@return int if a negative value is returned the function failed

*************************************************************************/

int read_test_vector_file_packets(const char *filename, symbol_type *din_array, const int SYM_LEN)
{
	FILE* file = fopen(filename, "rt");
	char line[256];
	const char *s = ", ";
	char *token;
	unsigned long tmp_short;
	unsigned long tmp_long;
	int cnt = 0;

// Since I don't want to regenerate my Matlab input file to be SYM_LEN long's wide, just take
// existing single shorts and duplicate
	while (fgets(line, sizeof(line), file))
	{
		token = strtok(line, s);

		tmp_short = atoi(token);
		tmp_long = ((tmp_short & 0xFFFF) << 48) | ((tmp_short & 0xFFFF) << 32) | ((tmp_short & 0xFFFF) << 16) | ((tmp_short & 0xFFFF) );

		for (int i = 0; i < SYM_LEN; i++)
		{
			din_array[cnt].symbol[i] = tmp_long;
	  }
	  // din_array++;
		cnt++;

	}

	fclose(file);
	return cnt;
}

//int read_test_vector_file_noisy_packets(const char *filename, symbol_type *din_array, char *din_erasure_ind_array, const int SYM_LEN)
int read_test_vector_file_noisy_packets(const char *filename, symbol_type *din_array, const int SYM_LEN)
{
	FILE* file = fopen(filename, "rt");
	char line[256];
	const char *s = ", ";
	char *token;
	unsigned long tmp_short;
	unsigned long tmp_long;
	int cnt = 0;

// Since I don't want to regenerate my Matlab input file to be SYM_LEN long's wide, just take
// existing single shorts and duplicate
	while (fgets(line, sizeof(line), file))
	{
		token = strtok(line, s);

		tmp_short = atoi(token);
// Assume an all zeros input is an erasure? Make sure to reflect this in Matlab generated input file
		if (tmp_short == 0)
		{
			din_array[cnt].is_erasure = 1; 
		}
		else 
		{
			din_array[cnt].is_erasure = 0; 
		}
		tmp_long = ((tmp_short & 0xFFFF) << 48) | ((tmp_short & 0xFFFF) << 32) | ((tmp_short & 0xFFFF) << 16) | ((tmp_short & 0xFFFF) );

		for (int i = 0; i < SYM_LEN; i++)
		{
			din_array[cnt].symbol[i] = tmp_long;
	  }
	  // din_array++;
		cnt++;

	}

	fclose(file);
	return cnt;
}

/**************************************************************

@brief The verify_output function verifies the kernels outputs 
	against a text file containing outputs generated in Matlab.
	This function only compares the parity bits generated by the 
	LDPC Encoder.

@return int if less than 0 an error has occurred

**************************************************************/
int verify_output()
{
	for (int i = 0; i < k_LEN; i++)
	{
		//if (i < kldpc){
	  for (int j = 0; j < SYM_LEN; j++)
	  {
			if (dout_vector[i].symbol[j] != dout[i].symbol[j])
				return -1;
	  }
	}
	return 0;
}

/*************************************************************************

@brief The init_opencl function intializes the OpenCL objects.
	Essentially it looks for an Altera OpenCl device, creates a 
	context for that device.  The AOC compiled kernel is pointed to 
	using the global variable "device_kernel".  Provided a valid .aocx
	compiled kernel is found, the FPGA is programmed and command and data
	queues are generated.

@return bool True if successful, otherwise false

**************************************************************************/
bool init_opencl() {
  cl_int status;
#ifdef VERBOSE
  printf("Initializing OpenCL\n");
#endif
  if(!setCwdToExeDir()) {
    return false;
  }

  //****************************
  // Get the OpenCL platform.
  //****************************

  // platform = findPlatform("Intel(R) FPGA SDK for OpenCL(TM)");
  // printf("INFO: The platform is %s\n", getPlatformName(platform).c_str());

  if (ptype == EMULATION_PLAT) {
             // new 'fast' emulator
        platform = findPlatform("Intel(R) FPGA Emulation Platform for OpenCL(TM)");
        if (platform == NULL) {
            printf("ERROR: Unable to find Intel(R) FPGA Emulation Platform for OpenCL(TM).\n");
         // For legacy emulator
            platform = findPlatform("Intel(R) FPGA SDK for OpenCL(TM)");
            if (platform == NULL) {
                printf("ERROR: Unable to find Intel(R) FPGA Legacy Emulation Platform for OpenCL(TM).\n");
                return false;
            }
        }
  }else{    //use for hardware and simulation
        platform = findPlatform("Intel(R) FPGA SDK for OpenCL(TM)");
        if (platform == NULL) {
            printf("ERROR: Unable to find Intel(R) FPGA SDK for OpenCL(TM).\n");
            return false;
        }
  }

  //*************************************
  // Query the available OpenCL devices.
  //*************************************
  scoped_array<cl_device_id> devices;
  devices.reset(getDevices(platform, CL_DEVICE_TYPE_ALL, &num_devices));
  device = devices[0];
#ifdef VERBOSE
  printf("Platform: %s\n", getPlatformName(platform).c_str());
  printf("Using %d device(s)\n", num_devices);
  for(unsigned i = 0; i < num_devices; ++i) {
    printf("  %s\n", getDeviceName(device).c_str());
  }
#endif

  //*********************
  // Create the context.
  //*********************
  context = clCreateContext(NULL, num_devices, &device, &oclContextCallback, NULL, &status);
  checkError(status, "Failed to create context");

  //************************************
  // Create the program for all device. 
  //************************************
  std::string binary_file = getBoardBinaryFile(device_kernel, device);
#ifdef VERBOSE
  printf("Using AOCX: %s\n", binary_file.c_str());
#endif
  program = createProgramFromBinary(context, binary_file.c_str(), &device, num_devices);

  //*******************************************
  // Build the program that was just created.
  //*******************************************
  status = clBuildProgram(program, 0, NULL, "", NULL, NULL);
  checkError(status, "Failed to build program");

  //**********************
  // Create Command queue.
  //**********************
  for (int i = 0; i < K_NUM_KERNELS; ++i)
  {
	  queue[i] = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
	  checkError(status, "Failed to create command queue %d", i);
  }

  //*****************
  // Create Kernel.
  //*****************
  for (int i = 0; i < K_NUM_KERNELS; ++i)
  {
	  kernel[i] = clCreateKernel(program, kernel_names[i], &status);
	  checkError(status, "Failed to create kernel %s", kernel_names[i]);
  }


  //**********************
  // Create Input buffers.
  //**********************
  input_a_buf = clCreateBuffer(context, CL_MEM_READ_ONLY, 
	n_LEN * sizeof(symbol_type), NULL, &status);
  checkError(status, "Failed to create buffer for input A");

  //**********************  
  // Create Output buffer.
  //**********************
  output_buf = clCreateBuffer(context, CL_MEM_WRITE_ONLY, 
	  k_LEN * sizeof(symbol_type), NULL, &status);
  checkError(status, "Failed to create buffer for output");
  
  return true;
}


/*************************************************************************

@brief The run function sets kernel arguments and then launches 
	the kernels.  Data is then read from the kernels write buffer.

@return void

**************************************************************************/
void run() {
  cl_int status;
  cl_event kernel_event;
  

  const double start_time = getCurrentTimestamp();
  seed = ((int)round(start_time)) % 1000000;
  printf("The seed used in this run is: %d\n", seed);
  
  //***********************************
  // Copy data from host to device
  //***********************************
  status = clEnqueueWriteBuffer(queue[0], input_a_buf, CL_TRUE,
	  0, n_LEN * sizeof(symbol_type), din_array, 0, NULL, NULL);
  checkError(status, "Failed to transfer input A");



  //***********************************
  // Set kernel arguments.
  //***********************************   
    
  //Read 
  status = clSetKernelArg(kernel[K_DATA_IN], 0, sizeof(cl_mem), &input_a_buf);
  checkError(status, "Failed to set kernel_rd arg 0");
  status = clSetKernelArg(kernel[K_DATA_IN], 1, sizeof(unsigned short), &n_LEN);
  checkError(status, "Failed to set kernel_rd arg 1");
  status = clSetKernelArg(kernel[K_DATA_IN], 2, sizeof(int), &seed);
  checkError(status, "Failed to set kernel_rd arg 2");
  status = clSetKernelArg(kernel[K_DATA_IN], 3, sizeof(int), &PER_numerator_div_64);
  checkError(status, "Failed to set kernel_rd arg 3");
  status = clSetKernelArg(kernel[K_DATA_IN], 4, sizeof(int), &code_ind);
  checkError(status, "Failed to set kernel_rd arg 4");
  status = clSetKernelArg(kernel[K_DATA_IN], 5, sizeof(long), &numFrames);
  checkError(status, "Failed to set kernel_rd arg 5");
  
  
  //decoder
  status = clSetKernelArg(kernel[K_LDPC_ERASURE_DECODER], 0, sizeof(short), &numItr);
  checkError(status, "Failed to set kernel_decoder arg 0");
  status = clSetKernelArg(kernel[K_LDPC_ERASURE_DECODER], 1, sizeof(int), &code_ind);
  checkError(status, "Failed to set kernel_decoder arg 1");
  
  //Write  
  status = clSetKernelArg(kernel[K_DATA_OUT], 0, sizeof(cl_mem), &output_buf);
  checkError(status, "Failed to set kernel_wr arg 0");
  status = clSetKernelArg(kernel[K_DATA_OUT], 1, sizeof(int), &code_ind);
  checkError(status, "Failed to set K_WRITER arg 1");
  status = clSetKernelArg(kernel[K_DATA_OUT], 2, sizeof(long), &numFrames);
  checkError(status, "Failed to set K_WRITER arg 2");

  

  //***********************************
  // Enqueue kernel.
  //***********************************
  size_t global_work_size = n_LEN;
#ifdef VERBOSE
  printf("Launching for device %d (%d elements)\n", 1, (int)global_work_size);
#endif

  //Read 
  status = clEnqueueTask(queue[K_DATA_IN], kernel[K_DATA_IN], 0, NULL, NULL);
  checkError(status, "Failed to launch ldpc_erasure_encoder");
  
  //LDPC Decoder
  status = clEnqueueTask(queue[K_LDPC_ERASURE_DECODER], kernel[K_LDPC_ERASURE_DECODER], 0, NULL, NULL);
  checkError(status, "Failed to launch K_LDPC_ERASURE_DECODER");

  //Write
  status = clEnqueueTask(queue[K_DATA_OUT], kernel[K_DATA_OUT], 0, NULL, &kernel_event);
  checkError(status, "Failed to launch kernel_write");


  //***************************************************
  // Wait for command queue to complete pending events
  //***************************************************
  status = clFinish(queue[K_DATA_OUT]);
    checkError(status, "Failed to finish (%d: %s)", K_DATA_OUT, kernel_names[K_DATA_OUT]);

  //********************************************
  // Read the result. This the final operation.
  //********************************************
  status = clEnqueueReadBuffer(queue[0], output_buf, CL_TRUE,
	  0, (k_LEN)*sizeof(symbol_type), dout, 0, NULL, NULL);

  
  const double end_time = getCurrentTimestamp();

  for(int i=0; i<n_LEN; i++){
	  //printf("%d\n", dout[i]);
  }

  // Wall-clock time taken.
  //printf("\nTime: %0.3f ms\n", (end_time - start_time) * 1e3);

  // Get kernel times using the OpenCL event profiling API.
  cl_ulong time_ns = getStartEndTime(kernel_event);
#ifdef VERBOSE
	printf("Kernel time: %0.3f ms\n", double(time_ns) * 1e-6);
	printf("The throughput in information bits/sec: %f\t", (float(SYM_LEN)*8.0*8.0*float(numFrames)*float(k_LEN))/(double(time_ns)* 1e-9));
#else
	printf("%f\t", (numFrames*n_LEN)/(double(time_ns)* 1e-9));
#endif
}


/*********************************************************

@brief Free the resources allocated during initialization

@return void
*********************************************************/
void cleanup() {

  for (int i = 0; i<K_NUM_KERNELS; ++i) {
	if (kernel[i])
	  clReleaseKernel(kernel[i]);
  }
  for (int i = 0; i<K_NUM_KERNELS; ++i) {
	if (queue[i])
	  clReleaseCommandQueue(queue[i]);
  }

  if(input_a_buf && input_a_buf) {
    clReleaseMemObject(input_a_buf);
  }
  if(output_buf && output_buf) {
    clReleaseMemObject(output_buf);
  }
  if(program) {
    clReleaseProgram(program);
  }
  if(context) {
    clReleaseContext(context);
  }
}


