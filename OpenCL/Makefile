#############################################################################################
# Copy this file to a location that is on the same level as dir: bin, device, host, etc.
# 
# STEP 1 - Configurations
#   Change K_SRC, K_DIR, K_AOCX_EMU_NAME, K_AOCX_HW_NAME, and K_BOARD below as appropriate.
#
# STEP 2 - Usage
#
#  COMPILATION:
#   To compile kernel code for legacy emulator, run:
#   > make emul
#   To compile kernel code for emulator, run:
#   > make emu
#   To compile kernel code for hardware, run:
#   > make hw_fast
#   OR:
#   > make hw
#
#  CLEAN:
#   To remove host executable in 'bin' dir, run:
#   > make clean
#   To remove all emulation files in 'bin' dir, run:
#   > make clean_emu
#   To remove compiled files in 'bin' dir, run:
#   > make clean_all
#  
#  PRINT SUGGESTED RUN COMMANDS
#   > make run
#############################################################################################

#--------------------------------------------
# USER CONFIGURATIONS FOR SPECIFIC KERNEL
# 
# File name for Kernel source code
K_SRC := ldpc_erasure_decoder_top.cl
# Relative path to folder where the kernel is at, from this Makefile, example: device
K_DIR := device
# File name of generated aocx file for simulation.
K_AOCX_SIM_NAME := mldpc_decoder_sim
# File name of generated aocx file for emulation. Same for legacy & normal emulator
K_AOCX_EMU_NAME := ldpc_erasure_decoder_em
# File name of generated aocx file for hardware
K_AOCX_HW_NAME := ldpc_erasure_decoder
# Target board name
K_BOARD := evc_p520_hpc_sg280h
# Relative path to inside 'common' dir from this Makefile, example: ../common
PATH_TO_COMMON := ../common
#--------------------------------------------
DEBUG=1

ifeq ($(VERBOSE),1)
ECHO := 
else
ECHO := @
endif

# Where is the Intel(R) FPGA SDK for OpenCL(TM) software?
#ifeq ($(wildcard $(INTELFPGAOCLSDKROOT)),)
#$(error Set INTELFPGAOCLSDKROOT to the root directory of the Intel(R) FPGA SDK for OpenCL(TM) software installation)
#endif
#ifeq ($(wildcard $(INTELFPGAOCLSDKROOT)/host/include/CL/opencl.h),)
#$(error Set INTELFPGAOCLSDKROOT to the root directory of the Intel(R) FPGA SDK for OpenCL(TM) software installation.)
#endif

# OpenCL compile and link flags.
AOCL_COMPILE_CONFIG := $(shell aocl compile-config )
AOCL_LINK_LIBS := $(shell aocl ldlibs )
AOCL_LINK_FLAGS := $(shell aocl ldflags )
# Linking with defences enabled
AOCL_LINK_FLAGS += -z noexecstack
AOCL_LINK_FLAGS += -Wl,-z,relro,-z,now
AOCL_LINK_FLAGS += -Wl,-Bsymbolic
AOCL_LINK_FLAGS += -pie
AOCL_LINK_CONFIG := $(AOCL_LINK_FLAGS) $(AOCL_LINK_LIBS)

# Compilation flags
ifeq ($(DEBUG),1)
CXXFLAGS += -g
else
CXXFLAGS += -O2
endif

# Compiling with defences enabled
CXXFLAGS += -fstack-protector
CXXFLAGS += -D_FORTIFY_SOURCE=2
CXXFLAGS += -Wformat -Wformat-security
CXXFLAGS += -fPIE

# We must force GCC to never assume that it can shove in its own
# sse2/sse3 versions of strlen and strcmp because they will CRASH.
# Very hard to debug!
CXXFLAGS += -fPIC -std=c++11 

# Compiler
CXX := g++ -g -O0

# Target
TARGET := host
TARGET_DIR := bin

# Directories
INC_DIRS := $(PATH_TO_COMMON)/inc ./host/inc 
LIB_DIRS := 

# Files
INCS := $(wildcard )
SRCS := $(wildcard host/src/*.cpp $(PATH_TO_COMMON)/src/AOCLUtils/*.cpp)
LIBS := rt pthread

#----------------------------------------
# Compilation
#----------------------------------------
# File name of generated aocx file for simulation.
K_AOCX_SIM := $(K_AOCX_SIM_NAME).aocx
# File name of generated aocx file for emulation. Same for legacy & normal emulator
K_AOCX_EMU := $(K_AOCX_EMU_NAME).aocx
# File name of generated aocx file for hardware
K_AOCX_HW := $(K_AOCX_HW_NAME).aocx

# simulator (-ghdl: enable collecting the waveform during the simulation)
sim : $(K_DIR)/$(K_SRC)
	${INTELFPGAOCLSDKROOT}/bin/aoc -march=simulator -ghdl -v $(K_DIR)/$(K_SRC) -o $(TARGET_DIR)/$(K_AOCX_SIM)

# legacy emulator
emul : $(K_DIR)/$(K_SRC)
	${INTELFPGAOCLSDKROOT}/bin/aoc -march=emulator -legacy-emulator -v $(K_DIR)/$(K_SRC) -o $(TARGET_DIR)/$(K_AOCX_EMU)

# emulator
emu : $(K_DIR)/$(K_SRC)
	${INTELFPGAOCLSDKROOT}/bin/aoc -march=emulator -v $(K_DIR)/$(K_SRC) -o $(TARGET_DIR)/$(K_AOCX_EMU)

# compile for hardware using fast-compile option
hw_fast : $(K_DIR)/$(K_SRC)
	${INTELFPGAOCLSDKROOT}/bin/aoc -v -board=$(K_BOARD) $(K_DIR)/$(K_SRC) -fast-compile -report -o $(TARGET_DIR)/$(K_AOCX_HW)

# compile for hardware
hw : $(K_DIR)/$(K_SRC)
	${INTELFPGAOCLSDKROOT}/bin/aoc -v -board=$(K_BOARD) -board-package=${AOCL_BOARD_PACKAGE_ROOT} $(K_DIR)/$(K_SRC) -seed=54321 -report -o $(TARGET_DIR)/$(K_AOCX_HW)

# Make it all!
#all : $(TARGET_DIR)/$(TARGET)
# compile host and kernel for emulation
all : host emul

# Host executable target.
#$(TARGET_DIR)/$(TARGET) : Makefile $(SRCS) $(INCS) $(TARGET_DIR)
host : Makefile $(SRCS) $(INCS) $(TARGET_DIR)
	$(ECHO)$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(EXTRACXXFLAGS) -fPIC $(foreach D,$(INC_DIRS),-I$D) \
			$(AOCL_COMPILE_CONFIG) $(SRCS) $(AOCL_LINK_CONFIG) \
			$(foreach D,$(LIB_DIRS),-L$D) \
			$(foreach L,$(LIBS),-l$L) \
			-o $(TARGET_DIR)/$(TARGET)

$(TARGET_DIR) :
	$(ECHO)mkdir $(TARGET_DIR)

#----------------------------------------
# PRINTOUT RUN COMMAND
#----------------------------------------
run :
	$(ECHO)echo 
	$(ECHO)echo Go into bin dir, and execute command below
	$(ECHO)echo FOR run on SIMULATOR:
	$(ECHO)echo unset CL_CONTEXT_EMULATOR_DEVICE_INTELFPGA
	$(ECHO)echo env CL_CONTEXT_MPSIM_DEVICE_INTELFPGA=1 ./$(TARGET) -e
	$(ECHO)echo set CL_CONTEXT_COMPILER_MODE_INTELFPGA=3 if the host program cannot find the simulator device
	$(ECHO)echo 
	$(ECHO)echo FOR run on LEGACY EMULATOR:
	$(ECHO)echo env CL_CONTEXT_EMULATOR_DEVICE_INTELFPGA=1 ./$(TARGET) -e
	$(ECHO)echo 
	$(ECHO)echo FOR run on EMULATOR:
	$(ECHO)echo env CL_CONFIG_CPU_EMULATE_DEVICES=1 ./$(TARGET) -e
	$(ECHO)echo 
	$(ECHO)echo FOR run on HARDWARE:
	$(ECHO)echo unset CL_CONFIG_CPU_EMULATE_DEVICES
	$(ECHO)echo unset CL_CONTEXT_EMULATOR_DEVICE_INTELFPGA
	$(ECHO)echo ./$(TARGET) -h
	$(ECHO)echo 

#----------------------------------------
# CLEAN
#----------------------------------------	
# Standard make targets
# delete host executable
clean :
	$(ECHO)rm -f $(TARGET_DIR)/$(TARGET)

# delete the generated comiled files/dir for emulation
clean_emu :
	$(ECHO)rm -rf $(TARGET_DIR)/$(K_AOCX_EMU_NAME)
	$(ECHO)rm -f $(TARGET_DIR)/$(K_AOCX_EMU_NAME).aoc*

# delete the generated compiled files/dir for hw
clean_hw :
	$(ECHO)rm -rf $(TARGET_DIR)/$(K_AOCX_HW_NAME)
	$(ECHO)rm -f $(TARGET_DIR)/$(K_AOCX_HW_NAME).aoc*

clean_all : clean clean_emu clean_hw

.PHONY : all clean
