###########################################################################
## Makefile generated for component 'windEmulatorStep4_WECSim'. 
## 
## Makefile     : windEmulatorStep4_WECSim.mk
## Generated on : Wed Jul 22 11:35:00 2026
## Final product: $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim
## Product type : executable
## 
###########################################################################

###########################################################################
## MACROS
###########################################################################

# Macro Descriptions:
# PRODUCT_NAME            Name of the system to build
# MAKEFILE                Name of this makefile

PRODUCT_NAME              = windEmulatorStep4_WECSim
MAKEFILE                  = windEmulatorStep4_WECSim.mk
MATLAB_ROOT               = C:/PROGRA~1/MATLAB/R2025b
MATLAB_BIN                = C:/PROGRA~1/MATLAB/R2025b/bin
MATLAB_ARCH_BIN           = $(MATLAB_BIN)/win64
START_DIR                 = C:/Users/jleonqu/Documents/GitHub/WEC-Sim_Applications/WEC-Sim_HIL/WindEmulatorFiles_edited
SOLVER                    = 
SOLVER_OBJ                = 
CLASSIC_INTERFACE         = 0
TGT_FCN_LIB               = ISO_C++
MODEL_HAS_DYNAMICALLY_LOADED_SFCNS = 0
RELATIVE_PATH_TO_ANCHOR   = ../..
C_STANDARD_OPTS           = 
CPP_STANDARD_OPTS         = 
LIBSSC_SLI_SLRT_X64_OBJS  = 
LIBSM_SSCI_SLRT_X64_OBJS  = 
LIBSSC_CORE_SLRT_X64_OBJS = 
LIBSM_SLRT_X64_OBJS       = 
LIBPM_MATH_SLRT_X64_OBJS  = 
LIBPM_ST_SLRT_X64_OBJS    = 
LIBMC_SLRT_X64_OBJS       = 
LIBEX_SLRT_X64_OBJS       = 
LIBPM_SLRT_X64_OBJS       = 

###########################################################################
## TOOLCHAIN SPECIFICATIONS
###########################################################################

# Toolchain Name:          Simulink Real-Time Toolchain
# Supported Version(s):    
# ToolchainInfo Version:   2025b
# Specification Revision:  1.0
# 
#-------------------------------------------
# Macros assumed to be defined elsewhere
#-------------------------------------------

# SLREALTIME_QNX_SP_ROOT
# SLREALTIME_QNX_VERSION

#-----------
# MACROS
#-----------

QCC_TARGET             = gcc_ntox86_64

TOOLCHAIN_SRCS = 
TOOLCHAIN_INCS = 
TOOLCHAIN_LIBS = -L$(MATLAB_ROOT)/toolbox/slrealtime/target/win64/target/lib -ltraceparser -lpps -lslrealtime_kernel -lslrealtime_platform -lslrealtime_rtps -lsocket -lboost_system -lboost_log -lpci -lopenblas -lpcap

#------------------------
# BUILD TOOL COMMANDS
#------------------------

# C Compiler: QNX C Compiler
CC = qcc

# Linker: QCC Linker
LD = q++

# C++ Compiler: QNX C++ Compiler
CPP = q++

# C++ Linker: QCC C++ Linker
CPP_LD = q++

# Archiver: QNX Archiver
AR = ntox86_64-gcc-ar

# Builder: GMAKE Utility
MAKE = make


#-------------------------
# Directives/Utilities
#-------------------------

CDEBUG              = -g -O0 -finstrument-functions
C_OUTPUT_FLAG       = -o
LDDEBUG             = -g
OUTPUT_FLAG         = -o
CPPDEBUG            = -g -O0 -finstrument-functions
CPP_OUTPUT_FLAG     = -o
CPPLDDEBUG          = -g
OUTPUT_FLAG         = -o
ARDEBUG             =
STATICLIB_OUTPUT_FLAG =
RM                  = @del /F
ECHO                = @echo
MV                  = @move
RUN                 =

#--------------------------------------
# "Faster Runs" Build Configuration
#--------------------------------------

ARFLAGS              = ruvs
CFLAGS               = -c -V$(QCC_TARGET) -g \
                       -O2 -fwrapv
CPPFLAGS             = -c -V$(QCC_TARGET) -g -std=gnu++14 -stdlib=libstdc++ \
                       -O2 -fwrapv
CPP_LDFLAGS          = -V$(QCC_TARGET) -g -std=gnu++14 -stdlib=libstdc++
CPP_SHAREDLIB_LDFLAGS  = -V$(QCC_TARGET) -shared -Wl,--no-undefined -g
LDFLAGS              = -V$(QCC_TARGET) -g -std=gnu++14 -stdlib=libstdc++
MAKE_FLAGS           = -f $(MAKEFILE)
SHAREDLIB_LDFLAGS    = -V$(QCC_TARGET) -shared -Wl,--no-undefined -g



###########################################################################
## OUTPUT INFO
###########################################################################

PRODUCT = $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim
PRODUCT_TYPE = "executable"
BUILD_TYPE = "Top-Level Standalone Executable"

###########################################################################
## INCLUDE PATHS
###########################################################################

INCLUDES_BUILDINFO = -I$(START_DIR) -I$(START_DIR)/windEmulatorStep4_WECSim_sg_rtw -I$(MATLAB_ROOT)/toolbox/slrealtime/simulink/blocks/dist/include -I$(MATLAB_ROOT)/toolbox/slrealtime/target/kernel/dist/include -I$(MATLAB_ROOT)/extern/include -I$(MATLAB_ROOT)/simulink/include -I$(MATLAB_ROOT)/rtw/c/src -I$(MATLAB_ROOT)/rtw/c/src/ext_mode/common -I$(MATLAB_ROOT)/extern/physmod/win64/ex/include -I$(MATLAB_ROOT)/extern/physmod/win64/mc/include -I$(MATLAB_ROOT)/extern/physmod/win64/pd/include -I$(MATLAB_ROOT)/extern/physmod/win64/pm/include -I$(MATLAB_ROOT)/extern/physmod/win64/pm_log/include -I$(MATLAB_ROOT)/extern/physmod/win64/pm_math/include -I$(MATLAB_ROOT)/extern/physmod/win64/pm_st/include -I$(MATLAB_ROOT)/extern/physmod/win64/sm/include -I$(MATLAB_ROOT)/extern/physmod/win64/sm_ssci/include -I$(MATLAB_ROOT)/extern/physmod/win64/ssc_core/include -I$(MATLAB_ROOT)/extern/physmod/win64/ssc_dae/include -I$(MATLAB_ROOT)/extern/physmod/win64/ssc_ds/include -I$(MATLAB_ROOT)/extern/physmod/win64/ssc_sli/include -IC:/ProgramData/Speedgoat/speedgoatlib/R2025b/10.0.1/sg_blocks/common/libsg -I$(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/instrumented

INCLUDES = $(INCLUDES_BUILDINFO)

###########################################################################
## DEFINES
###########################################################################

DEFINES_ = -DSIMULINK_REAL_TIME -D_QNX_SOURCE -DFMU_CG_TARGET=30
DEFINES_BUILD_ARGS = -DCLASSIC_INTERFACE=0 -DALLOCATIONFCN=0 -DEXT_MODE=1 -DMAT_FILE=0 -DONESTEPFCN=1 -DTERMFCN=1 -DMULTI_INSTANCE_CODE=0 -DINTEGER_CODE=0 -DMT=0
DEFINES_CUSTOM = 
DEFINES_OPTS = -DTID01EQ=1
DEFINES_STANDARD = -DMODEL=windEmulatorStep4_WECSim -DNUMST=2 -DNCSTATES=44 -DHAVESTDIO -DRT -DUSE_RTMODEL

DEFINES = $(DEFINES_) $(DEFINES_BUILD_ARGS) $(DEFINES_CUSTOM) $(DEFINES_OPTS) $(DEFINES_STANDARD)

###########################################################################
## SOURCE FILES
###########################################################################

SRCS = $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_dxf_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_tdxf_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_tdxy_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_dxy_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_vdf.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_duf.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_obs_exp.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_dxf.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_zc.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_vmf.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_obs_act.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_dnf_v_x.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_assert.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_dnf_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_imin.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_log.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_obs_all.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_obs_il.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_mcon_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_f.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_acon_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_acon.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_imax.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_mode.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_dxicr_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_tduf_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_gateway.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_create.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_setParameters.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_asserts.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_deriv.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_checkDynamics.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_compOutputsKin.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_compOutputsDyn.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_asm_delegate.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_sim_delegate.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_mode_zero_crossings.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_logging.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_geometries.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_gateway.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_y.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_ic.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_f.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dxy_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_slv.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_mcon_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_slf.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dxicr_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dxicrm_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_obs_exp.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_mode.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_zc.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_eq_tol.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dxm_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_duf.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_lv.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dxm.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_acon_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_imin.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_mnl_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_vdf.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_qx_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_mdxy_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_apar_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_vpf.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_imax.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_vsf.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_nldv.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_sclv.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dxf.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_tdxy_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_obs_act.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dxy.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_var_tol.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dnf_v_x.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_mnl.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_tdxf_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_log.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_icrm_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_assert.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_obs_all.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_mpar_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_obs_il.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_slf0.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_vmf.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dxf_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_acon.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dnf_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dxcer_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dxdelt_p.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_gateway.c $(MATLAB_ROOT)/rtw/c/src/rt_matrx.c $(MATLAB_ROOT)/rtw/c/src/rt_printf.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/rt_backsubrr_dbl.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/rt_forwardsubrr_dbl.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/rt_lu_real.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/rt_matrixlib_dbl.c $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/rtGetInf.cpp $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/rtGetNaN.cpp $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/rt_nonfinite.cpp $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/slrealtime_datatype_ground.cpp $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim.cpp $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_cal.cpp $(START_DIR)/ecat_config_xml_0.cpp $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/rte_windEmulatorStep4_WECSim_parameters.cpp $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/main.cpp $(MATLAB_ROOT)/toolbox/slrealtime/simulink/blocks/mex/slrealtimeenablelogging.cpp $(MATLAB_ROOT)/toolbox/slrealtime/simulink/blocks/mex/slrtEcatUtils.cpp host_timer_x86.c slrealtime_code_profiling_utility_functions.cpp

ALL_SRCS = $(SRCS)

###########################################################################
## OBJECTS
###########################################################################

OBJS = windEmulatorStep4_WECSim_1e9c788f_1_ds.o windEmulatorStep4_WECSim_1e9c788f_1_ds_dxf_p.o windEmulatorStep4_WECSim_1e9c788f_1_ds_tdxf_p.o windEmulatorStep4_WECSim_1e9c788f_1_ds_tdxy_p.o windEmulatorStep4_WECSim_1e9c788f_1_ds_dxy_p.o windEmulatorStep4_WECSim_1e9c788f_1_ds_vdf.o windEmulatorStep4_WECSim_1e9c788f_1_ds_duf.o windEmulatorStep4_WECSim_1e9c788f_1_ds_obs_exp.o windEmulatorStep4_WECSim_1e9c788f_1_ds_dxf.o windEmulatorStep4_WECSim_1e9c788f_1_ds_zc.o windEmulatorStep4_WECSim_1e9c788f_1_ds_vmf.o windEmulatorStep4_WECSim_1e9c788f_1_ds_obs_act.o windEmulatorStep4_WECSim_1e9c788f_1_ds_dnf_v_x.o windEmulatorStep4_WECSim_1e9c788f_1_ds_assert.o windEmulatorStep4_WECSim_1e9c788f_1_ds_dnf_p.o windEmulatorStep4_WECSim_1e9c788f_1_ds_imin.o windEmulatorStep4_WECSim_1e9c788f_1_ds_log.o windEmulatorStep4_WECSim_1e9c788f_1_ds_obs_all.o windEmulatorStep4_WECSim_1e9c788f_1_ds_obs_il.o windEmulatorStep4_WECSim_1e9c788f_1_ds_mcon_p.o windEmulatorStep4_WECSim_1e9c788f_1_ds_f.o windEmulatorStep4_WECSim_1e9c788f_1_ds_acon_p.o windEmulatorStep4_WECSim_1e9c788f_1_ds_acon.o windEmulatorStep4_WECSim_1e9c788f_1_ds_imax.o windEmulatorStep4_WECSim_1e9c788f_1_ds_mode.o windEmulatorStep4_WECSim_1e9c788f_1_ds_dxicr_p.o windEmulatorStep4_WECSim_1e9c788f_1_ds_tduf_p.o windEmulatorStep4_WECSim_1e9c788f_1.o windEmulatorStep4_WECSim_1e9c788f_1_gateway.o windEmulatorStep4_WECSim_dfbb7ac7_1.o windEmulatorStep4_WECSim_dfbb7ac7_1_create.o windEmulatorStep4_WECSim_dfbb7ac7_1_setParameters.o windEmulatorStep4_WECSim_dfbb7ac7_1_asserts.o windEmulatorStep4_WECSim_dfbb7ac7_1_deriv.o windEmulatorStep4_WECSim_dfbb7ac7_1_checkDynamics.o windEmulatorStep4_WECSim_dfbb7ac7_1_compOutputsKin.o windEmulatorStep4_WECSim_dfbb7ac7_1_compOutputsDyn.o windEmulatorStep4_WECSim_dfbb7ac7_1_asm_delegate.o windEmulatorStep4_WECSim_dfbb7ac7_1_sim_delegate.o windEmulatorStep4_WECSim_dfbb7ac7_1_mode_zero_crossings.o windEmulatorStep4_WECSim_dfbb7ac7_1_logging.o windEmulatorStep4_WECSim_dfbb7ac7_1_geometries.o windEmulatorStep4_WECSim_dfbb7ac7_1_gateway.o windEmulatorStep4_WECSim_5bdcd402_1_ds.o windEmulatorStep4_WECSim_5bdcd402_1_ds_y.o windEmulatorStep4_WECSim_5bdcd402_1_ds_ic.o windEmulatorStep4_WECSim_5bdcd402_1_ds_f.o windEmulatorStep4_WECSim_5bdcd402_1_ds_dxy_p.o windEmulatorStep4_WECSim_5bdcd402_1_ds_slv.o windEmulatorStep4_WECSim_5bdcd402_1_ds_mcon_p.o windEmulatorStep4_WECSim_5bdcd402_1_ds_slf.o windEmulatorStep4_WECSim_5bdcd402_1_ds_dxicr_p.o windEmulatorStep4_WECSim_5bdcd402_1_ds_dxicrm_p.o windEmulatorStep4_WECSim_5bdcd402_1_ds_obs_exp.o windEmulatorStep4_WECSim_5bdcd402_1_ds_mode.o windEmulatorStep4_WECSim_5bdcd402_1_ds_zc.o windEmulatorStep4_WECSim_5bdcd402_1_ds_eq_tol.o windEmulatorStep4_WECSim_5bdcd402_1_ds_dxm_p.o windEmulatorStep4_WECSim_5bdcd402_1_ds_duf.o windEmulatorStep4_WECSim_5bdcd402_1_ds_lv.o windEmulatorStep4_WECSim_5bdcd402_1_ds_dxm.o windEmulatorStep4_WECSim_5bdcd402_1_ds_acon_p.o windEmulatorStep4_WECSim_5bdcd402_1_ds_imin.o windEmulatorStep4_WECSim_5bdcd402_1_ds_mnl_p.o windEmulatorStep4_WECSim_5bdcd402_1_ds_vdf.o windEmulatorStep4_WECSim_5bdcd402_1_ds_qx_p.o windEmulatorStep4_WECSim_5bdcd402_1_ds_mdxy_p.o windEmulatorStep4_WECSim_5bdcd402_1_ds_apar_p.o windEmulatorStep4_WECSim_5bdcd402_1_ds_vpf.o windEmulatorStep4_WECSim_5bdcd402_1_ds_imax.o windEmulatorStep4_WECSim_5bdcd402_1_ds_vsf.o windEmulatorStep4_WECSim_5bdcd402_1_ds_nldv.o windEmulatorStep4_WECSim_5bdcd402_1_ds_sclv.o windEmulatorStep4_WECSim_5bdcd402_1_ds_dxf.o windEmulatorStep4_WECSim_5bdcd402_1_ds_tdxy_p.o windEmulatorStep4_WECSim_5bdcd402_1_ds_obs_act.o windEmulatorStep4_WECSim_5bdcd402_1_ds_dxy.o windEmulatorStep4_WECSim_5bdcd402_1_ds_var_tol.o windEmulatorStep4_WECSim_5bdcd402_1_ds_dnf_v_x.o windEmulatorStep4_WECSim_5bdcd402_1_ds_mnl.o windEmulatorStep4_WECSim_5bdcd402_1_ds_tdxf_p.o windEmulatorStep4_WECSim_5bdcd402_1_ds_log.o windEmulatorStep4_WECSim_5bdcd402_1_ds_icrm_p.o windEmulatorStep4_WECSim_5bdcd402_1_ds_assert.o windEmulatorStep4_WECSim_5bdcd402_1_ds_obs_all.o windEmulatorStep4_WECSim_5bdcd402_1_ds_mpar_p.o windEmulatorStep4_WECSim_5bdcd402_1_ds_obs_il.o windEmulatorStep4_WECSim_5bdcd402_1_ds_slf0.o windEmulatorStep4_WECSim_5bdcd402_1_ds_vmf.o windEmulatorStep4_WECSim_5bdcd402_1_ds_dxf_p.o windEmulatorStep4_WECSim_5bdcd402_1_ds_acon.o windEmulatorStep4_WECSim_5bdcd402_1_ds_dnf_p.o windEmulatorStep4_WECSim_5bdcd402_1_ds_dxcer_p.o windEmulatorStep4_WECSim_5bdcd402_1_ds_dxdelt_p.o windEmulatorStep4_WECSim_5bdcd402_1.o windEmulatorStep4_WECSim_5bdcd402_1_gateway.o rt_matrx.o rt_printf.o rt_backsubrr_dbl.o rt_forwardsubrr_dbl.o rt_lu_real.o rt_matrixlib_dbl.o rtGetInf.o rtGetNaN.o rt_nonfinite.o slrealtime_datatype_ground.o windEmulatorStep4_WECSim.o windEmulatorStep4_WECSim_cal.o ecat_config_xml_0.o rte_windEmulatorStep4_WECSim_parameters.o main.o slrealtimeenablelogging.o slrtEcatUtils.o host_timer_x86.o slrealtime_code_profiling_utility_functions.o

ALL_OBJS = $(OBJS)

###########################################################################
## PREBUILT OBJECT FILES
###########################################################################

PREBUILT_OBJS = 

###########################################################################
## LIBRARIES
###########################################################################

LIBS = $(MATLAB_ROOT)/extern/physmod/win64/ssc_sli/lib/ssc_sli_slrt_x64.a $(MATLAB_ROOT)/extern/physmod/win64/sm_ssci/lib/sm_ssci_slrt_x64.a $(MATLAB_ROOT)/extern/physmod/win64/ssc_core/lib/ssc_core_slrt_x64.a $(MATLAB_ROOT)/extern/physmod/win64/sm/lib/sm_slrt_x64.a $(MATLAB_ROOT)/extern/physmod/win64/pm_math/lib/pm_math_slrt_x64.a $(MATLAB_ROOT)/extern/physmod/win64/pm_st/lib/pm_st_slrt_x64.a $(MATLAB_ROOT)/extern/physmod/win64/mc/lib/mc_slrt_x64.a $(MATLAB_ROOT)/extern/physmod/win64/ex/lib/ex_slrt_x64.a $(MATLAB_ROOT)/extern/physmod/win64/pm/lib/pm_slrt_x64.a $(MATLAB_ROOT)/toolbox/slrealtime/simulink/blocks/dist/win64/lib/libslrealtime_libsrc_ecatapi_slrt_x64.a $(MATLAB_ROOT)/toolbox/slrealtime/simulink/blocks/dist/win64/lib/libslrealtime_libsrc_ecatstack_slrt_x64.a $(MATLAB_ROOT)/toolbox/slrealtime/simulink/blocks/dist/win64/lib/libslrealtime_libsrc_pervar_slrt_x64.a

###########################################################################
## SYSTEM LIBRARIES
###########################################################################

SYSTEM_LIBS = 

###########################################################################
## ADDITIONAL TOOLCHAIN FLAGS
###########################################################################

#---------------
# C Compiler
#---------------

CFLAGS_BASIC = $(DEFINES) $(INCLUDES)

CFLAGS += $(CFLAGS_BASIC)

#-----------------
# C++ Compiler
#-----------------

CPPFLAGS_BASIC = $(DEFINES) $(INCLUDES)

CPPFLAGS += $(CPPFLAGS_BASIC)

#---------------
# C++ Linker
#---------------

CPP_LDFLAGS_ = -lsg_qnx710_x86_64 -LC:/ProgramData/Speedgoat/speedgoatlib/R2025b/10.0.1/sg_blocks/common/libsg

CPP_LDFLAGS += $(CPP_LDFLAGS_)

#------------------------------
# C++ Shared Library Linker
#------------------------------

CPP_SHAREDLIB_LDFLAGS_ = -lsg_qnx710_x86_64 -LC:/ProgramData/Speedgoat/speedgoatlib/R2025b/10.0.1/sg_blocks/common/libsg

CPP_SHAREDLIB_LDFLAGS += $(CPP_SHAREDLIB_LDFLAGS_)

#-----------
# Linker
#-----------

LDFLAGS_ = -lsg_qnx710_x86_64 -LC:/ProgramData/Speedgoat/speedgoatlib/R2025b/10.0.1/sg_blocks/common/libsg

LDFLAGS += $(LDFLAGS_)

#--------------------------
# Shared Library Linker
#--------------------------

SHAREDLIB_LDFLAGS_ = -lsg_qnx710_x86_64 -LC:/ProgramData/Speedgoat/speedgoatlib/R2025b/10.0.1/sg_blocks/common/libsg

SHAREDLIB_LDFLAGS += $(SHAREDLIB_LDFLAGS_)

###########################################################################
## INLINED COMMANDS
###########################################################################

###########################################################################
## PHONY TARGETS
###########################################################################

.PHONY : all build buildobj clean info prebuild


all : build
	@echo "### Successfully generated all binary outputs."


build : prebuild $(PRODUCT)


buildobj : prebuild $(OBJS) $(PREBUILT_OBJS) $(LIBS)
	@echo "### Successfully generated all binary outputs."


prebuild : 


###########################################################################
## FINAL TARGET
###########################################################################

#-------------------------------------------
# Create a standalone executable            
#-------------------------------------------

$(PRODUCT) : $(OBJS) $(PREBUILT_OBJS) $(LIBS)
	@echo "### Creating standalone executable "$(PRODUCT)" ..."
	$(CPP_LD) $(CPP_LDFLAGS) -o $(PRODUCT) $(OBJS) -Wl,--start-group $(LIBS) -Wl,--end-group $(SYSTEM_LIBS) $(TOOLCHAIN_LIBS)
	@echo "### Created: $(PRODUCT)"


###########################################################################
## INTERMEDIATE TARGETS
###########################################################################

#---------------------
# SOURCE-TO-OBJECT
#---------------------

%.o : $(RELATIVE_PATH_TO_ANCHOR)/%.c
	$(CC) $(CFLAGS) -o $@ $<


%.o : $(RELATIVE_PATH_TO_ANCHOR)/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/extern/physmod/win64/ex/src/%.c
	$(CC) $(CFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/extern/physmod/win64/ex/src/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/extern/physmod/win64/mc/src/%.c
	$(CC) $(CFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/extern/physmod/win64/mc/src/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/extern/physmod/win64/pm/src/%.c
	$(CC) $(CFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/extern/physmod/win64/pm/src/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/extern/physmod/win64/pm_math/src/%.c
	$(CC) $(CFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/extern/physmod/win64/pm_math/src/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/extern/physmod/win64/pm_st/src/%.c
	$(CC) $(CFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/extern/physmod/win64/pm_st/src/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/extern/physmod/win64/sm/src/%.c
	$(CC) $(CFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/extern/physmod/win64/sm/src/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/extern/physmod/win64/sm_ssci/src/%.c
	$(CC) $(CFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/extern/physmod/win64/sm_ssci/src/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/extern/physmod/win64/ssc_core/src/%.c
	$(CC) $(CFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/extern/physmod/win64/ssc_core/src/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/extern/physmod/win64/ssc_sli/src/%.c
	$(CC) $(CFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/extern/physmod/win64/ssc_sli/src/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/toolbox/slrealtime/simulink/blocks/mex/%.c
	$(CC) $(CFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/toolbox/slrealtime/simulink/blocks/mex/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


%.o : $(START_DIR)/%.c
	$(CC) $(CFLAGS) -o $@ $<


%.o : $(START_DIR)/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


%.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/%.c
	$(CC) $(CFLAGS) -o $@ $<


%.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/rtw/c/src/%.c
	$(CC) $(CFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/rtw/c/src/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/simulink/src/%.c
	$(CC) $(CFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/simulink/src/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/toolbox/simulink/blocks/src/%.c
	$(CC) $(CFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/toolbox/simulink/blocks/src/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


%.o : ../%.c
	$(CC) $(CFLAGS) -o $@ $<


%.o : ../%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/toolbox/coder/profile/src/%.c
	$(CC) $(CFLAGS) -o $@ $<


%.o : $(MATLAB_ROOT)/toolbox/coder/profile/src/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


%.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/instrumented/%.c
	$(CC) $(CFLAGS) -o $@ $<


%.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/instrumented/%.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_dxf_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_dxf_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_tdxf_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_tdxf_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_tdxy_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_tdxy_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_dxy_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_dxy_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_vdf.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_vdf.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_duf.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_duf.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_obs_exp.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_obs_exp.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_dxf.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_dxf.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_zc.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_zc.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_vmf.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_vmf.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_obs_act.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_obs_act.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_dnf_v_x.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_dnf_v_x.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_assert.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_assert.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_dnf_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_dnf_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_imin.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_imin.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_log.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_log.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_obs_all.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_obs_all.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_obs_il.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_obs_il.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_mcon_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_mcon_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_f.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_f.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_acon_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_acon_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_acon.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_acon.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_imax.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_imax.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_mode.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_mode.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_dxicr_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_dxicr_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_ds_tduf_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_ds_tduf_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_1e9c788f_1_gateway.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_1e9c788f_1_gateway.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_dfbb7ac7_1.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_dfbb7ac7_1_create.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_create.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_dfbb7ac7_1_setParameters.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_setParameters.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_dfbb7ac7_1_asserts.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_asserts.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_dfbb7ac7_1_deriv.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_deriv.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_dfbb7ac7_1_checkDynamics.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_checkDynamics.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_dfbb7ac7_1_compOutputsKin.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_compOutputsKin.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_dfbb7ac7_1_compOutputsDyn.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_compOutputsDyn.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_dfbb7ac7_1_asm_delegate.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_asm_delegate.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_dfbb7ac7_1_sim_delegate.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_sim_delegate.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_dfbb7ac7_1_mode_zero_crossings.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_mode_zero_crossings.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_dfbb7ac7_1_logging.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_logging.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_dfbb7ac7_1_geometries.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_geometries.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_dfbb7ac7_1_gateway.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_dfbb7ac7_1_gateway.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_y.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_y.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_ic.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_ic.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_f.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_f.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_dxy_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dxy_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_slv.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_slv.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_mcon_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_mcon_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_slf.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_slf.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_dxicr_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dxicr_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_dxicrm_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dxicrm_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_obs_exp.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_obs_exp.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_mode.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_mode.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_zc.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_zc.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_eq_tol.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_eq_tol.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_dxm_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dxm_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_duf.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_duf.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_lv.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_lv.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_dxm.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dxm.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_acon_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_acon_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_imin.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_imin.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_mnl_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_mnl_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_vdf.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_vdf.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_qx_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_qx_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_mdxy_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_mdxy_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_apar_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_apar_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_vpf.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_vpf.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_imax.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_imax.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_vsf.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_vsf.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_nldv.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_nldv.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_sclv.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_sclv.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_dxf.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dxf.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_tdxy_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_tdxy_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_obs_act.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_obs_act.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_dxy.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dxy.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_var_tol.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_var_tol.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_dnf_v_x.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dnf_v_x.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_mnl.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_mnl.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_tdxf_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_tdxf_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_log.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_log.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_icrm_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_icrm_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_assert.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_assert.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_obs_all.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_obs_all.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_mpar_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_mpar_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_obs_il.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_obs_il.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_slf0.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_slf0.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_vmf.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_vmf.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_dxf_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dxf_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_acon.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_acon.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_dnf_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dnf_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_dxcer_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dxcer_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_ds_dxdelt_p.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_ds_dxdelt_p.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1.c
	$(CC) $(CFLAGS) -o $@ $<


windEmulatorStep4_WECSim_5bdcd402_1_gateway.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_5bdcd402_1_gateway.c
	$(CC) $(CFLAGS) -o $@ $<


rt_matrx.o : $(MATLAB_ROOT)/rtw/c/src/rt_matrx.c
	$(CC) $(CFLAGS) -o $@ $<


rt_printf.o : $(MATLAB_ROOT)/rtw/c/src/rt_printf.c
	$(CC) $(CFLAGS) -o $@ $<


rt_backsubrr_dbl.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/rt_backsubrr_dbl.c
	$(CC) $(CFLAGS) -o $@ $<


rt_forwardsubrr_dbl.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/rt_forwardsubrr_dbl.c
	$(CC) $(CFLAGS) -o $@ $<


rt_lu_real.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/rt_lu_real.c
	$(CC) $(CFLAGS) -o $@ $<


rt_matrixlib_dbl.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/rt_matrixlib_dbl.c
	$(CC) $(CFLAGS) -o $@ $<


rtGetInf.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/rtGetInf.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


rtGetNaN.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/rtGetNaN.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


rt_nonfinite.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/rt_nonfinite.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


slrealtime_datatype_ground.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/slrealtime_datatype_ground.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


windEmulatorStep4_WECSim.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


windEmulatorStep4_WECSim_cal.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/windEmulatorStep4_WECSim_cal.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


ecat_config_xml_0.o : $(START_DIR)/ecat_config_xml_0.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


rte_windEmulatorStep4_WECSim_parameters.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/rte_windEmulatorStep4_WECSim_parameters.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


main.o : $(START_DIR)/windEmulatorStep4_WECSim_sg_rtw/main.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


slrealtimeenablelogging.o : $(MATLAB_ROOT)/toolbox/slrealtime/simulink/blocks/mex/slrealtimeenablelogging.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


slrtEcatUtils.o : $(MATLAB_ROOT)/toolbox/slrealtime/simulink/blocks/mex/slrtEcatUtils.cpp
	$(CPP) $(CPPFLAGS) -o $@ $<


#------------------------
# BUILDABLE LIBRARIES
#------------------------

$(MATLAB_ROOT)/extern/physmod/win64/ssc_sli/lib/ssc_sli_slrt_x64.a : $(LIBSSC_SLI_SLRT_X64_OBJS)
	@echo "### Creating static library $@ ..."
	$(AR) $(ARFLAGS)  $@ $(LIBSSC_SLI_SLRT_X64_OBJS)


$(MATLAB_ROOT)/extern/physmod/win64/sm_ssci/lib/sm_ssci_slrt_x64.a : $(LIBSM_SSCI_SLRT_X64_OBJS)
	@echo "### Creating static library $@ ..."
	$(AR) $(ARFLAGS)  $@ $(LIBSM_SSCI_SLRT_X64_OBJS)


$(MATLAB_ROOT)/extern/physmod/win64/ssc_core/lib/ssc_core_slrt_x64.a : $(LIBSSC_CORE_SLRT_X64_OBJS)
	@echo "### Creating static library $@ ..."
	$(AR) $(ARFLAGS)  $@ $(LIBSSC_CORE_SLRT_X64_OBJS)


$(MATLAB_ROOT)/extern/physmod/win64/sm/lib/sm_slrt_x64.a : $(LIBSM_SLRT_X64_OBJS)
	@echo "### Creating static library $@ ..."
	$(AR) $(ARFLAGS)  $@ $(LIBSM_SLRT_X64_OBJS)


$(MATLAB_ROOT)/extern/physmod/win64/pm_math/lib/pm_math_slrt_x64.a : $(LIBPM_MATH_SLRT_X64_OBJS)
	@echo "### Creating static library $@ ..."
	$(AR) $(ARFLAGS)  $@ $(LIBPM_MATH_SLRT_X64_OBJS)


$(MATLAB_ROOT)/extern/physmod/win64/pm_st/lib/pm_st_slrt_x64.a : $(LIBPM_ST_SLRT_X64_OBJS)
	@echo "### Creating static library $@ ..."
	$(AR) $(ARFLAGS)  $@ $(LIBPM_ST_SLRT_X64_OBJS)


$(MATLAB_ROOT)/extern/physmod/win64/mc/lib/mc_slrt_x64.a : $(LIBMC_SLRT_X64_OBJS)
	@echo "### Creating static library $@ ..."
	$(AR) $(ARFLAGS)  $@ $(LIBMC_SLRT_X64_OBJS)


$(MATLAB_ROOT)/extern/physmod/win64/ex/lib/ex_slrt_x64.a : $(LIBEX_SLRT_X64_OBJS)
	@echo "### Creating static library $@ ..."
	$(AR) $(ARFLAGS)  $@ $(LIBEX_SLRT_X64_OBJS)


$(MATLAB_ROOT)/extern/physmod/win64/pm/lib/pm_slrt_x64.a : $(LIBPM_SLRT_X64_OBJS)
	@echo "### Creating static library $@ ..."
	$(AR) $(ARFLAGS)  $@ $(LIBPM_SLRT_X64_OBJS)


###########################################################################
## DEPENDENCIES
###########################################################################

$(ALL_OBJS) : rtw_proj.tmw $(MAKEFILE)


###########################################################################
## MISCELLANEOUS TARGETS
###########################################################################

info : 
	@echo "### PRODUCT = $(PRODUCT)"
	@echo "### PRODUCT_TYPE = $(PRODUCT_TYPE)"
	@echo "### BUILD_TYPE = $(BUILD_TYPE)"
	@echo "### INCLUDES = $(INCLUDES)"
	@echo "### DEFINES = $(DEFINES)"
	@echo "### ALL_SRCS = $(ALL_SRCS)"
	@echo "### ALL_OBJS = $(ALL_OBJS)"
	@echo "### LIBS = $(LIBS)"
	@echo "### MODELREF_LIBS = $(MODELREF_LIBS)"
	@echo "### SYSTEM_LIBS = $(SYSTEM_LIBS)"
	@echo "### TOOLCHAIN_LIBS = $(TOOLCHAIN_LIBS)"
	@echo "### CFLAGS = $(CFLAGS)"
	@echo "### LDFLAGS = $(LDFLAGS)"
	@echo "### SHAREDLIB_LDFLAGS = $(SHAREDLIB_LDFLAGS)"
	@echo "### CPPFLAGS = $(CPPFLAGS)"
	@echo "### CPP_LDFLAGS = $(CPP_LDFLAGS)"
	@echo "### CPP_SHAREDLIB_LDFLAGS = $(CPP_SHAREDLIB_LDFLAGS)"
	@echo "### ARFLAGS = $(ARFLAGS)"
	@echo "### MAKE_FLAGS = $(MAKE_FLAGS)"


clean : 
	$(ECHO) "### Deleting all derived files ..."
	$(RM) $(subst /,\,$(PRODUCT))
	$(RM) $(subst /,\,$(ALL_OBJS))
	$(ECHO) "### Deleted all derived files."


