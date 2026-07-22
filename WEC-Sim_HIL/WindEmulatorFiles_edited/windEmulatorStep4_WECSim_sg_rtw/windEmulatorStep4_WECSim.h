/*
 * windEmulatorStep4_WECSim.h
 *
 * Code generation for model "windEmulatorStep4_WECSim".
 *
 * Model version              : 10.8
 * Simulink Coder version : 25.2 (R2025b) 28-Jul-2025
 * C++ source code generated on : Wed Jul 22 11:34:26 2026
 *
 * Target selection: speedgoat.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef windEmulatorStep4_WECSim_h_
#define windEmulatorStep4_WECSim_h_
#include <logsrv.h>
#include "rtwtypes.h"
#include "simstruc.h"
#include "fixedpoint.h"
#include "rtw_extmode.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "stddef.h"
#include "stdlib.h"
#include "string"
#include "Logger.hpp"
#include "slecatutils.h"
#include "sl_types_def.h"
#include "verify/verifyIntrf.h"
#include "StartCallbackAPI.h"
#include "sysran_types.h"
#include "nesl_rtw_rtp.h"
#include "windEmulatorStep4_WECSim_1e9c788f_1_gateway.h"
#include "nesl_rtw.h"
#include "windEmulatorStep4_WECSim_dfbb7ac7_1_gateway.h"
#include "windEmulatorStep4_WECSim_5bdcd402_1_gateway.h"
#include "windEmulatorStep4_WECSim_types.h"
#include "abbState.h"
#include "expType.h"

extern "C"
{

#include "rt_nonfinite.h"

}

#include <cstring>
#include <stddef.h>
#include <cfloat>
#include <cmath>

extern "C"
{

#include "rtGetInf.h"

}

extern "C"
{

#include "rtGetNaN.h"

}

#include "windEmulatorStep4_WECSim_cal.h"
#include "rt_matrixlib.h"
#include "zero_crossing_types.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContStateDisabled
#define rtmGetContStateDisabled(rtm)   ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
#define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
#define rtmGetContStates(rtm)          ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
#define rtmSetContStates(rtm, val)     ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
#define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
#define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetFinalTime
#define rtmGetFinalTime(rtm)           ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetIntgData
#define rtmGetIntgData(rtm)            ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
#define rtmSetIntgData(rtm, val)       ((rtm)->intgData = (val))
#endif

#ifndef rtmGetMassMatrixIr
#define rtmGetMassMatrixIr(rtm)        ((rtm)->massMatrixIr)
#endif

#ifndef rtmSetMassMatrixIr
#define rtmSetMassMatrixIr(rtm, val)   ((rtm)->massMatrixIr = (val))
#endif

#ifndef rtmGetMassMatrixJc
#define rtmGetMassMatrixJc(rtm)        ((rtm)->massMatrixJc)
#endif

#ifndef rtmSetMassMatrixJc
#define rtmSetMassMatrixJc(rtm, val)   ((rtm)->massMatrixJc = (val))
#endif

#ifndef rtmGetMassMatrixNzMax
#define rtmGetMassMatrixNzMax(rtm)     ((rtm)->massMatrixNzMax)
#endif

#ifndef rtmSetMassMatrixNzMax
#define rtmSetMassMatrixNzMax(rtm, val) ((rtm)->massMatrixNzMax = (val))
#endif

#ifndef rtmGetMassMatrixPr
#define rtmGetMassMatrixPr(rtm)        ((rtm)->massMatrixPr)
#endif

#ifndef rtmSetMassMatrixPr
#define rtmSetMassMatrixPr(rtm, val)   ((rtm)->massMatrixPr = (val))
#endif

#ifndef rtmGetMassMatrixType
#define rtmGetMassMatrixType(rtm)      ((rtm)->massMatrixType)
#endif

#ifndef rtmSetMassMatrixType
#define rtmSetMassMatrixType(rtm, val) ((rtm)->massMatrixType = (val))
#endif

#ifndef rtmGetOdeDELTA
#define rtmGetOdeDELTA(rtm)            ((rtm)->odeDELTA)
#endif

#ifndef rtmSetOdeDELTA
#define rtmSetOdeDELTA(rtm, val)       ((rtm)->odeDELTA = (val))
#endif

#ifndef rtmGetOdeDFDX
#define rtmGetOdeDFDX(rtm)             ((rtm)->odeDFDX)
#endif

#ifndef rtmSetOdeDFDX
#define rtmSetOdeDFDX(rtm, val)        ((rtm)->odeDFDX = (val))
#endif

#ifndef rtmGetOdeE
#define rtmGetOdeE(rtm)                ((rtm)->odeE)
#endif

#ifndef rtmSetOdeE
#define rtmSetOdeE(rtm, val)           ((rtm)->odeE = (val))
#endif

#ifndef rtmGetOdeEDOT
#define rtmGetOdeEDOT(rtm)             ((rtm)->odeEDOT)
#endif

#ifndef rtmSetOdeEDOT
#define rtmSetOdeEDOT(rtm, val)        ((rtm)->odeEDOT = (val))
#endif

#ifndef rtmGetOdeF0
#define rtmGetOdeF0(rtm)               ((rtm)->odeF0)
#endif

#ifndef rtmSetOdeF0
#define rtmSetOdeF0(rtm, val)          ((rtm)->odeF0 = (val))
#endif

#ifndef rtmGetOdeF1
#define rtmGetOdeF1(rtm)               ((rtm)->odeF1)
#endif

#ifndef rtmSetOdeF1
#define rtmSetOdeF1(rtm, val)          ((rtm)->odeF1 = (val))
#endif

#ifndef rtmGetOdeFAC
#define rtmGetOdeFAC(rtm)              ((rtm)->odeFAC)
#endif

#ifndef rtmSetOdeFAC
#define rtmSetOdeFAC(rtm, val)         ((rtm)->odeFAC = (val))
#endif

#ifndef rtmGetOdeFMXDOT
#define rtmGetOdeFMXDOT(rtm)           ((rtm)->odeFMXDOT)
#endif

#ifndef rtmSetOdeFMXDOT
#define rtmSetOdeFMXDOT(rtm, val)      ((rtm)->odeFMXDOT = (val))
#endif

#ifndef rtmGetOdeMASSMATRIX_M
#define rtmGetOdeMASSMATRIX_M(rtm)     ((rtm)->odeMASSMATRIX_M)
#endif

#ifndef rtmSetOdeMASSMATRIX_M
#define rtmSetOdeMASSMATRIX_M(rtm, val) ((rtm)->odeMASSMATRIX_M = (val))
#endif

#ifndef rtmGetOdeMASSMATRIX_M1
#define rtmGetOdeMASSMATRIX_M1(rtm)    ((rtm)->odeMASSMATRIX_M1)
#endif

#ifndef rtmSetOdeMASSMATRIX_M1
#define rtmSetOdeMASSMATRIX_M1(rtm, val) ((rtm)->odeMASSMATRIX_M1 = (val))
#endif

#ifndef rtmGetOdePIVOTS
#define rtmGetOdePIVOTS(rtm)           ((rtm)->odePIVOTS)
#endif

#ifndef rtmSetOdePIVOTS
#define rtmSetOdePIVOTS(rtm, val)      ((rtm)->odePIVOTS = (val))
#endif

#ifndef rtmGetOdeW
#define rtmGetOdeW(rtm)                ((rtm)->odeW)
#endif

#ifndef rtmSetOdeW
#define rtmSetOdeW(rtm, val)           ((rtm)->odeW = (val))
#endif

#ifndef rtmGetOdeX0
#define rtmGetOdeX0(rtm)               ((rtm)->odeX0)
#endif

#ifndef rtmSetOdeX0
#define rtmSetOdeX0(rtm, val)          ((rtm)->odeX0 = (val))
#endif

#ifndef rtmGetOdeX1START
#define rtmGetOdeX1START(rtm)          ((rtm)->odeX1START)
#endif

#ifndef rtmSetOdeX1START
#define rtmSetOdeX1START(rtm, val)     ((rtm)->odeX1START = (val))
#endif

#ifndef rtmGetOdeXDOT
#define rtmGetOdeXDOT(rtm)             ((rtm)->odeXDOT)
#endif

#ifndef rtmSetOdeXDOT
#define rtmSetOdeXDOT(rtm, val)        ((rtm)->odeXDOT = (val))
#endif

#ifndef rtmGetOdeXTMP
#define rtmGetOdeXTMP(rtm)             ((rtm)->odeXTMP)
#endif

#ifndef rtmSetOdeXTMP
#define rtmSetOdeXTMP(rtm, val)        ((rtm)->odeXTMP = (val))
#endif

#ifndef rtmGetOdeZTMP
#define rtmGetOdeZTMP(rtm)             ((rtm)->odeZTMP)
#endif

#ifndef rtmSetOdeZTMP
#define rtmSetOdeZTMP(rtm, val)        ((rtm)->odeZTMP = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
#define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
#define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
#define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
#define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetSampleHitArray
#define rtmGetSampleHitArray(rtm)      ((rtm)->Timing.sampleHitArray)
#endif

#ifndef rtmGetStepSize
#define rtmGetStepSize(rtm)            ((rtm)->Timing.stepSize)
#endif

#ifndef rtmGetZCCacheNeedsReset
#define rtmGetZCCacheNeedsReset(rtm)   ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
#define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGet_TimeOfLastOutput
#define rtmGet_TimeOfLastOutput(rtm)   ((rtm)->Timing.timeOfLastOutput)
#endif

#ifndef rtmGetdX
#define rtmGetdX(rtm)                  ((rtm)->derivs)
#endif

#ifndef rtmSetdX
#define rtmSetdX(rtm, val)             ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
#define rtmGetStopRequested(rtm)       ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
#define rtmSetStopRequested(rtm, val)  ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
#define rtmGetStopRequestedPtr(rtm)    (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTFinal
#define rtmGetTFinal(rtm)              ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

#ifndef rtmGetTStart
#define rtmGetTStart(rtm)              ((rtm)->Timing.tStart)
#endif

#ifndef rtmGetTimeOfLastOutput
#define rtmGetTimeOfLastOutput(rtm)    ((rtm)->Timing.timeOfLastOutput)
#endif

/* Block signals for system '<S41>/parseCtrlWord' */
struct B_parseCtrlWord_windEmulatorS_T {
  uint16_T off1Ctrl;                   /* '<S41>/parseCtrlWord' */
  uint16_T off2Ctrl;                   /* '<S41>/parseCtrlWord' */
  uint16_T off3Ctrl;                   /* '<S41>/parseCtrlWord' */
  uint16_T enableOperation;            /* '<S41>/parseCtrlWord' */
  uint16_T rampOutZero;                /* '<S41>/parseCtrlWord' */
  uint16_T rampHold;                   /* '<S41>/parseCtrlWord' */
  uint16_T rampInZero;                 /* '<S41>/parseCtrlWord' */
  uint16_T reset;                      /* '<S41>/parseCtrlWord' */
  uint16_T inching1;                   /* '<S41>/parseCtrlWord' */
  uint16_T inching2;                   /* '<S41>/parseCtrlWord' */
  uint16_T remoteCmd;                  /* '<S41>/parseCtrlWord' */
  uint16_T extCtrlLoc;                 /* '<S41>/parseCtrlWord' */
};

/* Block states (default storage) for system '<S41>/parseCtrlWord' */
struct DW_parseCtrlWord_windEmulator_T {
  int32_T sfEvent;                     /* '<S41>/parseCtrlWord' */
  boolean_T doneDoubleBufferReInit;    /* '<S41>/parseCtrlWord' */
};

/* Block signals for system '<S43>/Moving Average' */
struct B_MovingAverage_windEmulatorS_T {
  real_T csumrev[2499];
  real_T MovingAverage;                /* '<S43>/Moving Average' */
};

/* Block states (default storage) for system '<S43>/Moving Average' */
struct DW_MovingAverage_windEmulator_T {
  dsp_simulink_MovingAverage_wi_T obj; /* '<S43>/Moving Average' */
  boolean_T objisempty;                /* '<S43>/Moving Average' */
};

/* Block signals for system '<S44>/Parse Status Word' */
struct B_ParseStatusWord_windEmulato_T {
  uint16_T rdy_on;                     /* '<S44>/Parse Status Word' */
  uint16_T rdy_run;                    /* '<S44>/Parse Status Word' */
  uint16_T rdy_ref;                    /* '<S44>/Parse Status Word' */
  uint16_T tripped;                    /* '<S44>/Parse Status Word' */
  uint16_T off2;                       /* '<S44>/Parse Status Word' */
  uint16_T off3;                       /* '<S44>/Parse Status Word' */
  uint16_T swc_on_inhib;               /* '<S44>/Parse Status Word' */
  uint16_T alarm;                      /* '<S44>/Parse Status Word' */
  uint16_T at_setpoint;                /* '<S44>/Parse Status Word' */
  uint16_T remote;                     /* '<S44>/Parse Status Word' */
  uint16_T above_limit;                /* '<S44>/Parse Status Word' */
  uint16_T ext_ctrl_loc;               /* '<S44>/Parse Status Word' */
  uint16_T ext_run_enable;             /* '<S44>/Parse Status Word' */
  uint16_T msw_b13;                    /* '<S44>/Parse Status Word' */
  uint16_T msw_b14;                    /* '<S44>/Parse Status Word' */
  uint16_T comm_err;                   /* '<S44>/Parse Status Word' */
};

/* Block states (default storage) for system '<S44>/Parse Status Word' */
struct DW_ParseStatusWord_windEmulat_T {
  int32_T sfEvent;                     /* '<S44>/Parse Status Word' */
  boolean_T doneDoubleBufferReInit;    /* '<S44>/Parse Status Word' */
};

/* Block signals for system '<S51>/Moving Average' */
struct B_MovingAverage_windEmulato_c_T {
  real_T csumrev[2499];
  real_T MovingAverage;                /* '<S51>/Moving Average' */
};

/* Block states (default storage) for system '<S51>/Moving Average' */
struct DW_MovingAverage_windEmulat_f_T {
  dsp_simulink_MovingAverage_wi_T obj; /* '<S51>/Moving Average' */
  boolean_T objisempty;                /* '<S51>/Moving Average' */
};

/* Block signals for system '<S60>/Nonlinear Wave Elevation' */
struct B_NonlinearWaveElevation_wind_T {
  real_T x_cg[6];                      /* '<S65>/Add' */
  real_T zero;                         /* '<S80>/zero' */
};

/* Block signals for system '<S81>/quaternion2EulXYZ' */
struct B_quaternion2EulXYZ_windEmula_T {
  real_T E[3];                         /* '<S81>/quaternion2EulXYZ' */
};

/* Block states (default storage) for system '<S81>/quaternion2EulXYZ' */
struct DW_quaternion2EulXYZ_windEmul_T {
  int32_T sfEvent;                     /* '<S81>/quaternion2EulXYZ' */
  boolean_T doneDoubleBufferReInit;    /* '<S81>/quaternion2EulXYZ' */
};

/* Block signals for system '<S126>/MATLAB Function1' */
struct B_MATLABFunction1_windEmulato_T {
  real_T dispPhase;                    /* '<S126>/MATLAB Function1' */
};

/* Block states (default storage) for system '<S126>/MATLAB Function1' */
struct DW_MATLABFunction1_windEmulat_T {
  int32_T sfEvent;                     /* '<S126>/MATLAB Function1' */
  boolean_T doneDoubleBufferReInit;    /* '<S126>/MATLAB Function1' */
};

/* Block signals for system '<S70>/Yaw Force Transforms' */
struct B_YawForceTransforms_windEmul_T {
  real_T TmpSignalConversionAtSFunctionI[6];/* '<S70>/Yaw Force Transforms' */
  real_T F_Excitation[6];              /* '<S70>/Yaw Force Transforms' */
  real_T F_RadiationDamping[6];        /* '<S70>/Yaw Force Transforms' */
  real_T F_AddedMass[6];               /* '<S70>/Yaw Force Transforms' */
  real_T F_Restoring[6];               /* '<S70>/Yaw Force Transforms' */
};

/* Block states (default storage) for system '<S70>/Yaw Force Transforms' */
struct DW_YawForceTransforms_windEmu_T {
  int32_T sfEvent;                     /* '<S70>/Yaw Force Transforms' */
  boolean_T doneDoubleBufferReInit;    /* '<S70>/Yaw Force Transforms' */
};

/* Block signals for system '<S133>/Yaw Kinematic Transforms' */
struct B_YawKinematicTransforms_wind_T {
  real_T TmpSignalConversionAtSFunctionI[6];/* '<S133>/Yaw Kinematic Transforms' */
  real_T TmpSignalConversionAtSFunctio_e[6];/* '<S133>/Yaw Kinematic Transforms' */
  real_T dispLoc[6];                   /* '<S133>/Yaw Kinematic Transforms' */
  real_T velLoc[6];                    /* '<S133>/Yaw Kinematic Transforms' */
  real_T accLoc[6];                    /* '<S133>/Yaw Kinematic Transforms' */
};

/* Block states (default storage) for system '<S133>/Yaw Kinematic Transforms' */
struct DW_YawKinematicTransforms_win_T {
  int32_T sfEvent;                     /* '<S133>/Yaw Kinematic Transforms' */
  boolean_T doneDoubleBufferReInit;    /* '<S133>/Yaw Kinematic Transforms' */
};

/* Block signals (default storage) */
struct B_windEmulatorStep4_WECSim_T {
  invPowerBus BusAssignment;           /* '<S9>/Bus Assignment' */
  invPowerBus BusAssignment_l;         /* '<S11>/Bus Assignment' */
  hptoSignalBus BusAssignment_n;       /* '<S6>/Bus Assignment' */
  acs880SignalBus BusAssignment_a;     /* '<S10>/Bus Assignment' */
  acs800SignalBus BusAssignment_h;     /* '<S8>/Bus Assignment' */
  shaftSignalBus BusAssignment_g;      /* '<S12>/Bus Assignment' */
  sidInfoBus BusAssignment_j;          /* '<S555>/Bus Assignment' */
  hptoCtrlBus BusAssignment_c;         /* '<S7>/Bus Assignment' */
  expCtrlBus BusAssignment_b;          /* '<S4>/Bus Assignment' */
  acs880CtrlBus BusAssignment_k;       /* '<S2>/Bus Assignment' */
  acs800CtrlBus BusAssignment_kc;      /* '<S1>/Bus Assignment' */
  sidCtrlBus BusAssignment_f;          /* '<S13>/Bus Assignment' */
  real_T CastToDouble;                 /* '<S10>/Cast To Double' */
  real_T Gain;                         /* '<S10>/Gain' */
  real_T CastToDouble1;                /* '<S10>/Cast To Double1' */
  real_T Gain2;                        /* '<S10>/Gain2' */
  real_T CastToDouble2;                /* '<S10>/Cast To Double2' */
  real_T Gain3;                        /* '<S10>/Gain3' */
  real_T CastToDouble3;                /* '<S10>/Cast To Double3' */
  real_T Gain4;                        /* '<S10>/Gain4' */
  real_T DataTypeConversion1;          /* '<S10>/Data Type Conversion1' */
  real_T Gain1;                        /* '<S10>/Gain1' */
  real_T DataTypeConversion2;          /* '<S10>/Data Type Conversion2' */
  real_T Gain5;                        /* '<S10>/Gain5' */
  real_T rpmrads;                      /* '<S48>/rpm -> rad//s' */
  real_T Product;                      /* '<S48>/Product' */
  real_T shaftPowerAverage_W;          /* '<S48>/shaftPowerAverage_W' */
  real_T shaftPower_W;                 /* '<S48>/shaftPower_W' */
  real_T frequency_Hz;                 /* '<S28>/frequency_Hz' */
  real_T motorCurrent_A;               /* '<S28>/motorCurrent_A' */
  real_T motorSpeed_rpm;               /* '<S28>/motorSpeed_rpm' */
  real_T motorTorque_Nm;               /* '<S28>/motorTorque_Nm' */
  real_T motorVoltage_V;               /* '<S28>/motorVoltage_V' */
  real_T shaftPower_W_a;               /* '<S28>/shaftPower_W' */
  real_T CastToDouble_o;               /* '<S2>/Cast To Double' */
  real_T expRunTime;                   /* '<S4>/expRunTime' */
  real_T rampValue;                    /* '<S555>/Switch' */
  real_T Switch2;                      /* '<S609>/Switch2' */
  real_T torqueSlewRate;               /* '<S609>/torqueSlewRate' */
  real_T MultiportSwitch;              /* '<S555>/Multiport Switch' */
  real_T Product_g;                    /* '<S555>/Product' */
  real_T Switch1;                      /* '<S609>/Switch1' */
  real_T speedSlewRate;                /* '<S609>/speedSlewRate' */
  real_T MultiportSwitch1;             /* '<S555>/Multiport Switch1' */
  real_T Product1;                     /* '<S555>/Product1' */
  real_T Sum;                          /* '<S13>/Sum' */
  real_T PProdOut;                     /* '<S595>/PProd Out' */
  real_T Integrator;                   /* '<S590>/Integrator' */
  real_T Sum_k;                        /* '<S600>/Sum' */
  real_T Switch;                       /* '<S598>/Switch' */
  real_T Switch2_h;                    /* '<S598>/Switch2' */
  real_T RateLimiter;                  /* '<S498>/Rate Limiter' */
  real_T INPUT_1_1_1[4];               /* '<S541>/INPUT_1_1_1' */
  real_T Internal;                     /* '<S531>/Internal' */
  real_T INPUT_2_1_1[4];               /* '<S541>/INPUT_2_1_1' */
  real_T Internal_j;                   /* '<S545>/Internal' */
  real_T Gain_l;                       /* '<S510>/Gain' */
  real_T INPUT_4_1_1[4];               /* '<S541>/INPUT_4_1_1' */
  real_T RateLimiter_j;                /* '<S7>/Rate Limiter' */
  real_T fw;                           /* '<S7>/f->w' */
  real_T Product4;                     /* '<S7>/Product4' */
  real_T Sin2;                         /* '<S7>/Sin2' */
  real_T excForceAmpNow_N;             /* '<S7>/excForceAmpNow_N' */
  real_T waveBotExcitationForce_N;     /* '<S7>/ExcitationForce_N' */
  real_T rampValue_a;                  /* '<S7>/Switch' */
  real_T Product_m;                    /* '<S7>/Product' */
  real_T DataTypeConversion4;          /* '<S8>/Data Type Conversion4' */
  real_T Gain3_g;                      /* '<S8>/Gain3' */
  real_T DataTypeConversion3;          /* '<S8>/Data Type Conversion3' */
  real_T Gain4_p;                      /* '<S8>/Gain4' */
  real_T DataTypeConversion5;          /* '<S8>/Data Type Conversion5' */
  real_T Gain5_g;                      /* '<S8>/Gain5' */
  real_T DataTypeConversion;           /* '<S8>/Data Type Conversion' */
  real_T Gain_c;                       /* '<S8>/Gain' */
  real_T DataTypeConversion1_e;        /* '<S8>/Data Type Conversion1' */
  real_T Gain1_f;                      /* '<S8>/Gain1' */
  real_T DataTypeConversion2_i;        /* '<S8>/Data Type Conversion2' */
  real_T Gain2_d;                      /* '<S8>/Gain2' */
  real_T Gain_d;                       /* '<S366>/Gain' */
  real_T INPUT_5_1_1[4];               /* '<S541>/INPUT_5_1_1' */
  real_T Internal_h;                   /* '<S542>/Internal' */
  real_T Gain_g;                       /* '<S509>/Gain' */
  real_T INPUT_3_1_1[4];               /* '<S541>/INPUT_3_1_1' */
  real_T RTP_1;                        /* '<S508>/RTP_1' */
  real_T STATE_1[37];                  /* '<S541>/STATE_1' */
  real_T OUTPUT_1_0[11];               /* '<S541>/OUTPUT_1_0' */
  real_T Pressure;                     /* '<S504>/Gain' */
  real_T psibar;                       /* '<S6>/psi -> bar' */
  real_T ShaftSpeedPump;               /* '<S511>/Gain' */
  real_T Step;                         /* '<S370>/Step' */
  real_T Clock;                        /* '<S370>/Clock' */
  real_T Sum_m;                        /* '<S370>/Sum' */
  real_T Product_a;                    /* '<S370>/Product' */
  real_T Output;                       /* '<S370>/Output' */
  real_T Saturation;                   /* '<S365>/Saturation' */
  real_T kDampingNow;                  /* '<S371>/kDampingNow' */
  real_T Product_mo;                   /* '<S371>/Product' */
  real_T kSpringNow;                   /* '<S371>/kSpringNow' */
  real_T Product1_g;                   /* '<S371>/Product1' */
  real_T ForceD;                       /* '<S371>/Add' */
  real_T Product2;                     /* '<S367>/Product2' */
  real_T Gain_n;                       /* '<S367>/Gain' */
  real_T TorqueInputRef;               /* '<S365>/Product' */
  real_T Abs;                          /* '<S365>/Abs' */
  real_T Sum_f;                        /* '<S437>/Sum' */
  real_T DiscreteTimeIntegrator;       /* '<S437>/Discrete-Time Integrator' */
  real_T Gain1_b;                      /* '<S437>/Gain1' */
  real_T ControlSignal1;               /* '<S365>/Switch' */
  real_T Abs2;                         /* '<S365>/Abs2' */
  real_T PressureRef;                  /* '<S365>/Switch2' */
  real_T Sum_me;                       /* '<S372>/Sum' */
  real_T DiscreteTimeIntegrator_i;     /* '<S372>/Discrete-Time Integrator' */
  real_T Gain1_j;                      /* '<S372>/Gain1' */
  real_T Sum_d;                        /* '<S435>/Sum' */
  real_T DiscreteTimeIntegrator_e;     /* '<S435>/Discrete-Time Integrator' */
  real_T Gain1_h;                      /* '<S435>/Gain1' */
  real_T ControlSignal2;               /* '<S365>/Switch' */
  real_T FlowMotor1;                   /* '<S501>/Gain' */
  real_T MultiportSwitch_p;            /* '<S3>/Multiport Switch' */
  real_T acs880RateLim;                /* '<S2>/acs880RateLim' */
  real_T ACS880Setpoint;               /* '<S2>/Saturation' */
  real_T Nm;                           /* '<S2>/Nm -> %' */
  real_T torqueSetpoint_Nm;            /* '<S26>/torqueSetpoint_Nm' */
  real_T torqueSetpoint_percent;       /* '<S26>/torqueSetpoint_percent' */
  real_T rpmrads_o;                    /* '<S43>/rpm -> rad//s' */
  real_T Product_o;                    /* '<S43>/Product' */
  real_T shaftPowerAverage_W_a;        /* '<S43>/shaftPowerAverage_W' */
  real_T shaftPower_W_e;               /* '<S43>/shaftPower_W' */
  real_T dcBusVoltage_V;               /* '<S24>/dcBusVoltage_V' */
  real_T frequency_Hz_j;               /* '<S24>/frequency_Hz' */
  real_T motorSpeed_rpm_j;             /* '<S24>/motorSpeed_rpm' */
  real_T motorTorque_Nm_p;             /* '<S24>/motorTorque_Nm' */
  real_T shaftPower_W_p;               /* '<S24>/shaftPower_W' */
  real_T temperature;                  /* '<S24>/temperature' */
  real_T CastToDouble_j;               /* '<S1>/Cast To Double' */
  real_T MultiportSwitch1_m;           /* '<S3>/Multiport Switch1' */
  real_T Nm_j;                         /* '<S1>/Nm -> %' */
  real_T torqueSetpoint_Nm_j;          /* '<S22>/torqueSetpoint_Nm' */
  real_T torqueSetpoint_percent_c;     /* '<S22>/torqueSetpoint_percent' */
  real_T barPa;                        /* '<S51>/bar->Pa' */
  real_T lmm3s;                        /* '<S51>/l//m -> m3//s' */
  real_T Product_gm;                   /* '<S51>/Product' */
  real_T rpmrads_oh;                   /* '<S51>/rpm -> rad//s' */
  real_T Product1_i;                   /* '<S51>/Product1' */
  real_T excShaftPowerAverage_W;       /* '<S51>/excShaftPowerAverage_W' */
  real_T excShaftPower_W;              /* '<S51>/excShaftPower_W' */
  real_T hydrPowerAverage_W;           /* '<S51>/hydrPowerAverage_W' */
  real_T hydrPower_W;                  /* '<S51>/hydrPower_W' */
  real_T ctrlSignal1;                  /* '<S35>/ctrlSignal1' */
  real_T ctrlSignal2;                  /* '<S35>/ctrlSignal2' */
  real_T excShaftSpeed_rpm;            /* '<S35>/excShaftSpeed_rpm' */
  real_T excShaftTorque_Nm;            /* '<S35>/excShaftTorque_Nm' */
  real_T genPumpFlow_lpm;              /* '<S35>/genPumpFlow_lpm' */
  real_T genShaftSpeed_rpm;            /* '<S35>/genShaftSpeed_rpm' */
  real_T genTorqueCmd_Nm;              /* '<S35>/genTorqueCmd_Nm' */
  real_T hmOutputShaftTorque_Nm;       /* '<S35>/hmOutputShaftTorque_Nm' */
  real_T pressure_bar;                 /* '<S35>/pressure_bar' */
  real_T excForce_N;                   /* '<S33>/excForce_N' */
  real_T genSpeedActual;               /* '<S33>/genSpeedActual' */
  real_T speedRef_rpm;                 /* '<S33>/speedRef_rpm' */
  real_T ramp;                         /* '<S31>/ramp' */
  real_T time;                         /* '<S31>/time' */
  real_T CastToDouble_a;               /* '<S12>/Cast To Double' */
  real_T Gain_nm;                      /* '<S12>/Gain' */
  real_T CastToDouble3_b;              /* '<S552>/Cast To Double3' */
  real_T encoderCountsToRad;           /* '<S552>/encoderCountsToRad' */
  real_T CastToDouble1_m;              /* '<S552>/Cast To Double1' */
  real_T Gain_b;                       /* '<S552>/Gain' */
  real_T Add2;                         /* '<S552>/Add2' */
  real_T TSamp;                        /* '<S553>/TSamp' */
  real_T Uk1;                          /* '<S553>/UD' */
  real_T Diff;                         /* '<S553>/Diff' */
  real_T radsrpm;                      /* '<S12>/rad//s->rpm' */
  real_T absEncoderPosition_rad;       /* '<S39>/absEncoderPosition_rad' */
  real_T absEncoderSpeed_rpm;          /* '<S39>/absEncoderSpeed_rpm' */
  real_T torqueActual_Nm;              /* '<S39>/torqueActual_Nm' */
  real_T L1Voltage;                    /* '<S9>/L1Voltage' */
  real_T L1Current;                    /* '<S9>/L1Current' */
  real_T L1PowFactor;                  /* '<S9>/L1PowFactor' */
  real_T L1ActivePow;                  /* '<S9>/L1ActivePow' */
  real_T L1THDu;                       /* '<S9>/L1THDu' */
  real_T L1THDi;                       /* '<S9>/L1THDi' */
  real_T L2Voltage;                    /* '<S9>/L2Voltage' */
  real_T L2Current;                    /* '<S9>/L2Current' */
  real_T L2PowFactor;                  /* '<S9>/L2PowFactor' */
  real_T L2ActivePow;                  /* '<S9>/L2ActivePow' */
  real_T L2THDu;                       /* '<S9>/L2THDu' */
  real_T L2THDi;                       /* '<S9>/L2THDi' */
  real_T L3Voltage;                    /* '<S9>/L3Voltage' */
  real_T L3Current;                    /* '<S9>/L3Current' */
  real_T L3PowFactor;                  /* '<S9>/L3PowFactor' */
  real_T L3ActivePow;                  /* '<S9>/L3ActivePow' */
  real_T L3THDu;                       /* '<S9>/L3THDu' */
  real_T L3THDi;                       /* '<S9>/L3THDi' */
  real_T totalFrequency;               /* '<S9>/totalFrequency' */
  real_T totalPowFactor;               /* '<S9>/totalPowFactor' */
  real_T totalActivePow;               /* '<S9>/totalActivePow' */
  real_T L1L2Voltage;                  /* '<S9>/L1L2Voltage' */
  real_T L2L3Voltage;                  /* '<S9>/L2L3Voltage' */
  real_T L3L1Voltage;                  /* '<S9>/L3L1Voltage' */
  real_T L1Voltage_f;                  /* '<S11>/L1Voltage' */
  real_T L1Current_a;                  /* '<S11>/L1Current' */
  real_T L1PowFactor_l;                /* '<S11>/L1PowFactor' */
  real_T L1ActivePow_m;                /* '<S11>/L1ActivePow' */
  real_T L1THDu_l;                     /* '<S11>/L1THDu' */
  real_T L1THDi_e;                     /* '<S11>/L1THDi' */
  real_T L2Voltage_j;                  /* '<S11>/L2Voltage' */
  real_T L2Current_h;                  /* '<S11>/L2Current' */
  real_T L2PowFactor_p;                /* '<S11>/L2PowFactor' */
  real_T L2ActivePow_b;                /* '<S11>/L2ActivePow' */
  real_T L2THDu_i;                     /* '<S11>/L2THDu' */
  real_T L2THDi_m;                     /* '<S11>/L2THDi' */
  real_T L3Voltage_n;                  /* '<S11>/L3Voltage' */
  real_T L3Current_h;                  /* '<S11>/L3Current' */
  real_T L3PowFactor_o;                /* '<S11>/L3PowFactor' */
  real_T L3ActivePow_l;                /* '<S11>/L3ActivePow' */
  real_T L3THDu_i;                     /* '<S11>/L3THDu' */
  real_T L3THDi_j;                     /* '<S11>/L3THDi' */
  real_T totalFrequency_g;             /* '<S11>/totalFrequency' */
  real_T totalPowFactor_e;             /* '<S11>/totalPowFactor' */
  real_T totalActivePow_g;             /* '<S11>/totalActivePow' */
  real_T L1L2Voltage_n;                /* '<S11>/L1L2Voltage' */
  real_T L2L3Voltage_h;                /* '<S11>/L2L3Voltage' */
  real_T L3L1Voltage_m;                /* '<S11>/L3L1Voltage' */
  real_T Switch3;                      /* '<S609>/Switch3' */
  real_T caseCounterSignalsNow;        /* '<S609>/caseCounterSignalsNow' */
  real_T DataTypeConversion_f;         /* '<S29>/Data Type Conversion' */
  real_T time_s;                       /* '<S29>/time_s' */
  real_T Product_k;                    /* '<S611>/Product' */
  real_T Product_d;                    /* '<S610>/Product' */
  real_T Product_l;                    /* '<S612>/Product' */
  real_T STATE_1_p[2];                 /* '<S216>/STATE_1' */
  real_T OUTPUT_1_1[28];               /* '<S216>/OUTPUT_1_1' */
  real_T velocity[6];                  /* '<S59>/Assignment1' */
  real_T Switch_l;                     /* '<S58>/Switch' */
  real_T Step_g;                       /* '<S58>/Step' */
  real_T INPUT_2_1_1_c[4];             /* '<S332>/INPUT_2_1_1' */
  real_T INPUT_3_1_1_d[4];             /* '<S332>/INPUT_3_1_1' */
  real_T DelayOneStep;                 /* '<S58>/Delay One Step' */
  real_T Sum1;                         /* '<S58>/Sum1' */
  real_T ProportionalGain;             /* '<S287>/Proportional Gain' */
  real_T Integrator_f;                 /* '<S282>/Integrator' */
  real_T DerivativeGain;               /* '<S273>/Derivative Gain' */
  real_T Tsamp;                        /* '<S277>/Tsamp' */
  real_T UD;                           /* '<S275>/UD' */
  real_T Diff_o;                       /* '<S275>/Diff' */
  real_T Sum_g;                        /* '<S291>/Sum' */
  real_T INPUT_1_1_1_n[4];             /* '<S332>/INPUT_1_1_1' */
  real_T STATE_1_d[62];                /* '<S332>/STATE_1' */
  real_T OUTPUT_1_0_m[24];             /* '<S332>/OUTPUT_1_0' */
  real_T Abs_p;                        /* '<S58>/Abs' */
  real_T Gain_f;                       /* '<S58>/Gain' */
  real_T ptoTorqueHydraulic;           /* '<S58>/Product5' */
  real_T ptoPowerMech;                 /* '<S58>/Product4' */
  real_T Switch1_k;                    /* '<S58>/Switch1' */
  real_T Abs1;                         /* '<S58>/Abs1' */
  real_T pistonPowerMech;              /* '<S58>/Product1' */
  real_T shaftSpeed;                   /* '<S58>/Product' */
  real_T shaftPower;                   /* '<S58>/Product2' */
  real_T powerHM;                      /* '<S58>/Product3' */
  real_T IntegralGain;                 /* '<S279>/Integral Gain' */
  real_T flowRateHCB;                  /* '<S219>/Gain' */
  real_T flowRateCV4;                  /* '<S220>/Gain' */
  real_T flowRateCV1;                  /* '<S221>/Gain' */
  real_T flowRateCV3;                  /* '<S222>/Gain' */
  real_T flowRateHCA;                  /* '<S223>/Gain' */
  real_T flowRateC;                    /* '<S224>/Gain' */
  real_T flowRateAccHP;                /* '<S225>/Gain' */
  real_T flowRateHMin;                 /* '<S226>/Gain' */
  real_T flowRateHMout;                /* '<S227>/Gain' */
  real_T flowRateAccLP;                /* '<S228>/Gain' */
  real_T flowRateD;                    /* '<S229>/Gain' */
  real_T flowRateCV2;                  /* '<S230>/Gain' */
  real_T pressureA;                    /* '<S241>/Gain' */
  real_T pressureB;                    /* '<S242>/Gain' */
  real_T pressureC;                    /* '<S243>/Gain' */
  real_T pressureHM;                   /* '<S244>/Gain' */
  real_T pressureD;                    /* '<S245>/Gain' */
  real_T position[6];                  /* '<S59>/Assignment ' */
  real_T INPUT_5_1_1_c[4];             /* '<S216>/INPUT_5_1_1' */
  real_T TransportDelay[6];            /* '<S60>/Transport Delay' */
  real_T F_SingleFrequency[6];         /* '<S131>/Product' */
  real_T F_AddedMass[6];               /* '<S69>/Product1' */
  real_T Clock_n;                      /* '<S57>/Clock' */
  real_T frequency;                    /* '<S124>/Divide' */
  real_T R;                            /* '<S68>/Switch' */
  real_T Product3[6];                  /* '<S126>/Product3' */
  real_T Product_p;                    /* '<S126>/Product' */
  real_T x_cg[6];                      /* '<S126>/Add3' */
  real_T Add;                          /* '<S126>/Add' */
  real_T coswt;                        /* '<S126>/Sine Wave Function1' */
  real_T Product1_b[6];                /* '<S126>/Product1' */
  real_T Add2_i;                       /* '<S126>/Add2' */
  real_T sinwt;                        /* '<S126>/Sine Wave Function' */
  real_T Product2_o[6];                /* '<S126>/Product2' */
  real_T Add1[6];                      /* '<S126>/Add1' */
  real_T F_wave[6];                    /* '<S68>/Add' */
  real_T F_Excitation[6];              /* '<S68>/Product' */
  real_T F_Gravity;                    /* '<S75>/Product' */
  real_T F_Buoyancy;                   /* '<S75>/Product1' */
  real_T Add1_a;                       /* '<S75>/Add1' */
  real_T VerticalBuoyancyForce[6];
               /* '<S75>/Assignment (Add Net Bouyancy Force  to Z-Direction)' */
  real_T AssignmentAddNetBouyancyForceto[3];
              /* '<S75>/Assignment (Add Net Bouyancy Force  to Z-Direction)1' */
  real_T Add3[3];                      /* '<S75>/Add3' */
  real_T Elementproduct[6];            /* '<S76>/Element product' */
  real_T Add3_d[3];                    /* '<S76>/Add3' */
  real_T Rotationalbuoyancyforce[6];
              /* '<S75>/Assignment (Add Net Bouyancy Force  to Z-Direction)2' */
  real_T Netbuoyancyforce[6];          /* '<S75>/Add2' */
  real_T x_cg_b[6];                    /* '<S63>/Add' */
  real_T LinearRestoringForce[6];      /* '<S74>/Product2' */
  real_T F_Restoring[6];               /* '<S74>/Add' */
  real_T v[6];
  real_T Abs1_e[6];                    /* '<S64>/Abs1' */
  real_T vv[6];                        /* '<S64>/Product' */
  real_T F_quadraticViscous[6];        /* '<S64>/Product1' */
  real_T F_MorisonAndViscous[6];       /* '<S64>/VisSum ' */
  real_T F_LinearDamping[6];           /* '<S61>/Product1' */
  real_T F_Total[6];                   /* '<S60>/Sum' */
  real_T INPUT_1_1_1_p[4];             /* '<S216>/INPUT_1_1_1' */
  real_T INPUT_1_1_2[4];               /* '<S216>/INPUT_1_1_2' */
  real_T INPUT_1_1_3[4];               /* '<S216>/INPUT_1_1_3' */
  real_T INPUT_2_1_1_l[4];             /* '<S216>/INPUT_2_1_1' */
  real_T INPUT_2_1_2[4];               /* '<S216>/INPUT_2_1_2' */
  real_T INPUT_2_1_3[4];               /* '<S216>/INPUT_2_1_3' */
  real_T TransportDelay_b[6];          /* '<S139>/Transport Delay' */
  real_T F_SingleFrequency_b[6];       /* '<S210>/Product' */
  real_T F_AddedMass_k[6];             /* '<S148>/Product1' */
  real_T frequency_g;                  /* '<S203>/Divide' */
  real_T R_p;                          /* '<S147>/Switch' */
  real_T Product3_p[6];                /* '<S205>/Product3' */
  real_T Product_oc;                   /* '<S205>/Product' */
  real_T x_cg_j[6];                    /* '<S205>/Add3' */
  real_T Add_m;                        /* '<S205>/Add' */
  real_T coswt_k;                      /* '<S205>/Sine Wave Function1' */
  real_T Product1_a[6];                /* '<S205>/Product1' */
  real_T Add2_f;                       /* '<S205>/Add2' */
  real_T sinwt_g;                      /* '<S205>/Sine Wave Function' */
  real_T Product2_p[6];                /* '<S205>/Product2' */
  real_T Add1_e[6];                    /* '<S205>/Add1' */
  real_T F_wave_i[6];                  /* '<S147>/Add' */
  real_T F_Excitation_f[6];            /* '<S147>/Product' */
  real_T F_Gravity_h;                  /* '<S154>/Product' */
  real_T F_Buoyancy_k;                 /* '<S154>/Product1' */
  real_T Add1_l;                       /* '<S154>/Add1' */
  real_T VerticalBuoyancyForce_g[6];
              /* '<S154>/Assignment (Add Net Bouyancy Force  to Z-Direction)' */
  real_T AssignmentAddNetBouyancyForce_c[3];
             /* '<S154>/Assignment (Add Net Bouyancy Force  to Z-Direction)1' */
  real_T Add3_m[3];                    /* '<S154>/Add3' */
  real_T Elementproduct_p[6];          /* '<S155>/Element product' */
  real_T Add3_h[3];                    /* '<S155>/Add3' */
  real_T Rotationalbuoyancyforce_l[6];
             /* '<S154>/Assignment (Add Net Bouyancy Force  to Z-Direction)2' */
  real_T Netbuoyancyforce_c[6];        /* '<S154>/Add2' */
  real_T x_cg_o[6];                    /* '<S142>/Add' */
  real_T LinearRestoringForce_l[6];    /* '<S153>/Product2' */
  real_T F_Restoring_a[6];             /* '<S153>/Add' */
  real_T v_j[6];
  real_T Abs1_b[6];                    /* '<S143>/Abs1' */
  real_T vv_j[6];                      /* '<S143>/Product' */
  real_T F_quadraticViscous_f[6];      /* '<S143>/Product1' */
  real_T F_MorisonAndViscous_d[6];     /* '<S143>/VisSum ' */
  real_T F_LinearDamping_k[6];         /* '<S140>/Product1' */
  real_T F_Total_c[6];                 /* '<S139>/Sum' */
  real_T INPUT_3_1_1_k[4];             /* '<S216>/INPUT_3_1_1' */
  real_T INPUT_3_1_2[4];               /* '<S216>/INPUT_3_1_2' */
  real_T INPUT_3_1_3[4];               /* '<S216>/INPUT_3_1_3' */
  real_T INPUT_4_1_1_k[4];             /* '<S216>/INPUT_4_1_1' */
  real_T INPUT_4_1_2[4];               /* '<S216>/INPUT_4_1_2' */
  real_T INPUT_4_1_3[4];               /* '<S216>/INPUT_4_1_3' */
  real_T OUTPUT_1_0_k[32];             /* '<S216>/OUTPUT_1_0' */
  real_T acceleration[6];              /* '<S59>/Assignment2' */
  real_T TmpSignalConversionAtAssignment[3];
  real_T Gain6;                        /* '<S59>/Gain6' */
  real_T Assignment6[3];               /* '<S59>/Assignment6' */
  real_T TmpSignalConversionAtAssignme_e[3];
  real_T Gain7;                        /* '<S59>/Gain7' */
  real_T Assignment7[3];               /* '<S59>/Assignment7' */
  real_T forceActuation[6];            /* '<S59>/Assignment3' */
  real_T TmpSignalConversionAtAssignme_p[3];
  real_T Gain4_j;                      /* '<S59>/Gain4' */
  real_T Assignment4[3];               /* '<S59>/Assignment4' */
  real_T TmpSignalConversionAtAssignme_f[3];
  real_T Gain5_j;                      /* '<S59>/Gain5' */
  real_T Assignment5[3];               /* '<S59>/Assignment5' */
  real_T forceInternalMechanics[6];    /* '<S59>/Add' */
  real_T powerInternalMechanics[6];    /* '<S59>/Product' */
  real_T TmpSignalConversionAtToWorkspac[48];
  real_T TmpSignalConversionAtToWorksp_c[60];
  real_T TransportDelay_o[6];          /* '<S64>/Transport Delay' */
  real_T position_d[6];                /* '<S55>/Constant' */
  real_T velocity_g[6];                /* '<S55>/Constant1' */
  real_T acceleration_b[6];            /* '<S55>/Constant2' */
  real_T TmpSignalConversionAtToWorksp_m[24];
  real_T TmpSignalConversionAtToWorks_ms[60];
  real_T TransportDelay_e[6];          /* '<S143>/Transport Delay' */
  real_T Add_mw;                       /* '<S376>/Add' */
  real_T ControlSignal31;              /* '<S431>/Product' */
  real_T RateLimiter1;                 /* '<S373>/Rate Limiter1' */
  real_T Add1_p;                       /* '<S376>/Add1' */
  real_T ControlSignal31_o;            /* '<S432>/Product' */
  real_T Switch_g;                     /* '<S376>/Switch' */
  real_T Gain_lr;                      /* '<S376>/Gain' */
  real_T RateLimiter_b;                /* '<S373>/Rate Limiter' */
  real_T wError;                       /* '<S375>/Sum' */
  real_T ProportionalGain_p;           /* '<S417>/Proportional Gain' */
  real_T Integrator_b;                 /* '<S412>/Integrator' */
  real_T DerivativeGain_a;             /* '<S405>/Derivative Gain' */
  real_T Filter;                       /* '<S407>/Filter' */
  real_T SumD;                         /* '<S407>/SumD' */
  real_T FilterCoefficient;            /* '<S415>/Filter Coefficient' */
  real_T Sum_c;                        /* '<S422>/Sum' */
  real_T Switch_i;                     /* '<S420>/Switch' */
  real_T Switch2_k;                    /* '<S420>/Switch2' */
  real_T ContolTorque;                 /* '<S375>/Gain2' */
  real_T IntegralGain_h;               /* '<S409>/Integral Gain' */
  real_T Add_mo;                       /* '<S439>/Add' */
  real_T ControlSignal31_d;            /* '<S494>/Product' */
  real_T RateLimiter1_m;               /* '<S436>/Rate Limiter1' */
  real_T Add1_f;                       /* '<S439>/Add1' */
  real_T ControlSignal31_m;            /* '<S495>/Product' */
  real_T Switch_n;                     /* '<S439>/Switch' */
  real_T Gain_fp;                      /* '<S439>/Gain' */
  real_T RateLimiter_a;                /* '<S436>/Rate Limiter' */
  real_T wError_c;                     /* '<S438>/Sum' */
  real_T ProportionalGain_h;           /* '<S480>/Proportional Gain' */
  real_T Integrator_l;                 /* '<S475>/Integrator' */
  real_T DerivativeGain_g;             /* '<S468>/Derivative Gain' */
  real_T Filter_j;                     /* '<S470>/Filter' */
  real_T SumD_c;                       /* '<S470>/SumD' */
  real_T FilterCoefficient_g;          /* '<S478>/Filter Coefficient' */
  real_T Sum_n;                        /* '<S485>/Sum' */
  real_T Switch_c;                     /* '<S483>/Switch' */
  real_T Switch2_m;                    /* '<S483>/Switch2' */
  real_T ContolTorque_f;               /* '<S438>/Gain2' */
  real_T IntegralGain_hk;              /* '<S472>/Integral Gain' */
  real_T ControlTorqueLoad;            /* '<S365>/Switch' */
  real_T SwitchLogic;                  /* '<S365>/Switch1' */
  real_T FlowPump1;                    /* '<S499>/m3toL' */
  real_T Sum_a;                        /* '<S522>/Sum' */
  real_T FlowAccumulator;              /* '<S500>/Gain' */
  real_T IProdOut;                     /* '<S587>/IProd Out' */
  real_T CastToDouble_g;               /* '<S609>/Cast To Double' */
  real_T CastToDouble1_i;              /* '<S609>/Cast To Double1' */
  real_T Divide1;                      /* '<S609>/Divide1' */
  real_T vecPercent;                   /* '<S609>/vecPercent' */
  real_T fromFileTorqueNow_Nm;         /* '<S609>/fromFileTorqueNow_Nm' */
  real_T fromFileSpeedNow_rpm;         /* '<S609>/fromFileSpeedNow_rpm' */
  real_T Gain_k;                       /* '<S365>/Gain' */
  real_T Abs1_a;                       /* '<S365>/Abs1' */
  real_T Switch_f;                     /* '<S436>/Switch' */
  real_T Switch_cn;                    /* '<S373>/Switch' */
  real_T Gain_j;                       /* '<S435>/Gain' */
  real_T Add_k;                        /* '<S435>/Add' */
  real_T Saturation_f;                 /* '<S435>/Saturation' */
  real_T Gain_ge;                      /* '<S372>/Gain' */
  real_T Add_i;                        /* '<S372>/Add' */
  real_T Saturation_a;                 /* '<S372>/Saturation' */
  real_T Gain2_n;                      /* '<S437>/Gain2' */
  real_T Divide;                       /* '<S437>/Divide' */
  real_T Product_aq;                   /* '<S437>/Product' */
  real_T Gain_a;                       /* '<S437>/Gain' */
  real_T Add_iu;                       /* '<S437>/Add' */
  real_T Add1_lh;                      /* '<S437>/Add1' */
  real_T Saturation_j;                 /* '<S437>/Saturation' */
  real_T Switch_ne;                    /* '<S374>/Switch' */
  real_T Switch1_g;                    /* '<S439>/Switch1' */
  real_T ControlSignal3;               /* '<S495>/Switch' */
  real_T Saturation_p;                 /* '<S495>/Saturation' */
  real_T ControlSignal3_h;             /* '<S494>/Switch' */
  real_T Saturation_af;                /* '<S494>/Saturation' */
  real_T Switch1_o;                    /* '<S376>/Switch1' */
  real_T ControlSignal3_e;             /* '<S432>/Switch' */
  real_T Saturation_k;                 /* '<S432>/Saturation' */
  real_T ControlSignal3_f;             /* '<S431>/Switch' */
  real_T Saturation_a3;                /* '<S431>/Saturation' */
  real_T Product2_c;                   /* '<S203>/Product2' */
  real_T Add2_g;                       /* '<S203>/Add2' */
  real_T SineWaveFunction;             /* '<S203>/Sine Wave Function' */
  real_T Add1_lhe;                     /* '<S203>/Add1' */
  real_T Product1_d;                   /* '<S203>/Product1' */
  real_T Product2_cu;                  /* '<S124>/Product2' */
  real_T Add2_e;                       /* '<S124>/Add2' */
  real_T SineWaveFunction_g;           /* '<S124>/Sine Wave Function' */
  real_T Add1_h;                       /* '<S124>/Add1' */
  real_T Product1_k;                   /* '<S124>/Product1' */
  real_T time_c;                       /* '<S4>/FexcRamp' */
  real_T ramp_l;                       /* '<S4>/FexcRamp' */
  uint32_T MultiportSwitch_h;          /* '<S609>/Multiport Switch' */
  uint32_T Mod;                        /* '<S609>/Mod' */
  uint32_T runCounter;                 /* '<S31>/runCounter' */
  uint32_T stepCounter;                /* '<S31>/stepCounter' */
  uint32_T readEncoderCounter;         /* '<S12>/readEncoderCounter' */
  uint32_T absEncoderCounts;           /* '<S39>/absEncoderCounts' */
  uint32_T CastTouint32;               /* '<S555>/Cast To uint32' */
  uint32_T Memory;                     /* '<S29>/Memory' */
  uint32_T Sum_e;                      /* '<S29>/Sum' */
  uint32_T loopCounter;                /* '<S29>/loopCounter' */
  uint32_T vecIndex;                   /* '<S609>/vecIndex' */
  uint32_T fileSamples;                /* '<S609>/fileSamples' */
  uint32_T runCounter_b;               /* '<S4>/FexcRamp' */
  uint32_T stepCounter_n;              /* '<S4>/FexcRamp' */
  real32_T L1VoltageRead;              /* '<S9>/L1VoltageRead' */
  real32_T L1CurrentRead;              /* '<S9>/L1CurrentRead' */
  real32_T L1PowFactorRead;            /* '<S9>/L1PowFactorRead' */
  real32_T L1ActivePowRead;            /* '<S9>/L1ActivePowRead' */
  real32_T L1THDuRead;                 /* '<S9>/L1THDuRead' */
  real32_T L1THDiRead;                 /* '<S9>/L1THDiRead' */
  real32_T L2VoltageRead;              /* '<S9>/L2VoltageRead' */
  real32_T L2CurrentRead;              /* '<S9>/L2CurrentRead' */
  real32_T L2PowFactorRead;            /* '<S9>/L2PowFactorRead' */
  real32_T L2ActivePowRead;            /* '<S9>/L2ActivePowRead' */
  real32_T L2THDuRead;                 /* '<S9>/L2THDuRead' */
  real32_T L2THDiRead;                 /* '<S9>/L2THDiRead' */
  real32_T L3VoltageRead;              /* '<S9>/L3VoltageRead' */
  real32_T L3CurrentRead;              /* '<S9>/L3CurrentRead' */
  real32_T L3PowFactorRead;            /* '<S9>/L3PowFactorRead' */
  real32_T L3ActivePowRead;            /* '<S9>/L3ActivePowRead' */
  real32_T L3THDuRead;                 /* '<S9>/L3THDuRead' */
  real32_T L3THDiRead;                 /* '<S9>/L3THDiRead' */
  real32_T FrequencyRead;              /* '<S9>/FrequencyRead' */
  real32_T totalPowFactorRead;         /* '<S9>/totalPowFactorRead' */
  real32_T totalActivePowRead;         /* '<S9>/totalActivePowRead' */
  real32_T L1L2VoltageRead;            /* '<S9>/L1L2VoltageRead' */
  real32_T L2L3VoltageRead;            /* '<S9>/L2L3VoltageRead' */
  real32_T L3L1VoltageRead;            /* '<S9>/L3L1VoltageRead' */
  real32_T L1VoltageRead_i;            /* '<S11>/L1VoltageRead' */
  real32_T L1CurrentRead_c;            /* '<S11>/L1CurrentRead' */
  real32_T L1PowFactorRead_k;          /* '<S11>/L1PowFactorRead' */
  real32_T L1ActivePowRead_h;          /* '<S11>/L1ActivePowRead' */
  real32_T L1THDuRead_a;               /* '<S11>/L1THDuRead' */
  real32_T L1THDiRead_a;               /* '<S11>/L1THDiRead' */
  real32_T L2VoltageRead_p;            /* '<S11>/L2VoltageRead' */
  real32_T L2CurrentRead_a;            /* '<S11>/L2CurrentRead' */
  real32_T L2PowFactorRead_o;          /* '<S11>/L2PowFactorRead' */
  real32_T L2ActivePowRead_e;          /* '<S11>/L2ActivePowRead' */
  real32_T L2THDuRead_l;               /* '<S11>/L2THDuRead' */
  real32_T L2THDiRead_j;               /* '<S11>/L2THDiRead' */
  real32_T L3VoltageRead_h;            /* '<S11>/L3VoltageRead' */
  real32_T L3CurrentRead_h;            /* '<S11>/L3CurrentRead' */
  real32_T L3PowFactorRead_o;          /* '<S11>/L3PowFactorRead' */
  real32_T L3ActivePowRead_m;          /* '<S11>/L3ActivePowRead' */
  real32_T L3THDuRead_b;               /* '<S11>/L3THDuRead' */
  real32_T L3THDiRead_m;               /* '<S11>/L3THDiRead' */
  real32_T FrequencyRead_g;            /* '<S11>/FrequencyRead' */
  real32_T totalPowFactorRead_n;       /* '<S11>/totalPowFactorRead' */
  real32_T totalActivePowRead_o;       /* '<S11>/totalActivePowRead' */
  real32_T L1L2VoltageRead_m;          /* '<S11>/L1L2VoltageRead' */
  real32_T L2L3VoltageRead_e;          /* '<S11>/L2L3VoltageRead' */
  real32_T L3L1VoltageRead_n;          /* '<S11>/L3L1VoltageRead' */
  int32_T EtherCATInit[6];             /* '<Root>/EtherCAT Init' */
  int32_T CastToDouble1_a;             /* '<S2>/Cast To Double1' */
  int32_T state;                       /* '<S26>/state' */
  int32_T CastToDouble1_ia;            /* '<S1>/Cast To Double1' */
  int32_T state_i;                     /* '<S22>/state' */
  int32_T readTorqueInput;             /* '<S12>/readTorqueInput' */
  int32_T lastRawCounts;               /* '<S552>/lastRawCounts' */
  int32_T CastToDouble_p;              /* '<S552>/Cast To Double' */
  int32_T Add_d;                       /* '<S552>/Add' */
  int32_T Abs_a;                       /* '<S552>/Abs' */
  int32_T Switch_g0;                   /* '<S552>/Switch' */
  int32_T lastTurn;                    /* '<S552>/lastTurn' */
  int32_T Add1_p2;                     /* '<S552>/Add1' */
  int32_T absEncoderTurns;             /* '<S39>/absEncoderTurns' */
  int32_T Sign;                        /* '<S552>/Sign' */
  abbStateEnum state_e;                /* '<S18>/ABB Fieldbus Control' */
  abbStateEnum state_ed;               /* '<S16>/ABB Fieldbus Control' */
  uint16_T EtherCATPDOReceive;         /* '<S10>/EtherCAT PDO Receive' */
  uint16_T statusWord;                 /* '<S28>/statusWord' */
  uint16_T ToUint16;                   /* '<S4>/ToUint16' */
  uint16_T ACS800StatusWord;           /* '<S8>/EtherCAT PDO Receive7' */
  uint16_T ctrlWord;                   /* '<S26>/ctrlWord' */
  uint16_T statusWord_m;               /* '<S24>/statusWord' */
  uint16_T ctrlWord_g;                 /* '<S22>/ctrlWord' */
  uint16_T expType;                    /* '<S31>/expType' */
  uint16_T ControlWord;                /* '<S18>/ABB Fieldbus Control' */
  uint16_T ControlWord_l;              /* '<S16>/ABB Fieldbus Control' */
  int16_T ACS880MotorVoltage;          /* '<S10>/EtherCAT PDO Receive1' */
  int16_T ACS880MotorCurrent;          /* '<S10>/EtherCAT PDO Receive2' */
  int16_T ACS880OutputFreq;            /* '<S10>/EtherCAT PDO Receive3' */
  int16_T ACS880MotorSpeed;            /* '<S10>/EtherCAT PDO Receive4' */
  int16_T ACS880MotorTorque;           /* '<S10>/EtherCAT PDO Receive5' */
  int16_T ACS880MotorShaftPower;       /* '<S10>/EtherCAT PDO Receive6' */
  int16_T ACS800DcBusVoltage;          /* '<S8>/EtherCAT PDO Receive8' */
  int16_T ACS800Frequency;             /* '<S8>/EtherCAT PDO Receive9' */
  int16_T ACS800Temperature;           /* '<S8>/EtherCAT PDO Receive10' */
  int16_T ACS800ActualFeedback;        /* '<S8>/EtherCAT PDO Receive11' */
  int16_T ACS800Torque;                /* '<S8>/EtherCAT PDO Receive12' */
  int16_T ACS800Power;                 /* '<S8>/EtherCAT PDO Receive13' */
  int16_T ACS800Ref2Int16;             /* '<S14>/ACS800Ref2Int16' */
  int16_T ACS800Ref1Int16;             /* '<S14>/ACS800Ref1Int16' */
  int16_T ACS880Ref2Int16;             /* '<S15>/ACS880Ref2Int16' */
  int16_T ACS880Ref1Int16;             /* '<S15>/ACS880Ref1Int16' */
  expTypeEnum expType_a;               /* '<S4>/expType' */
  expTypeEnum toExpTypeEnum;           /* '<S3>/toExpTypeEnum' */
  uint8_T encoderStatus[2];            /* '<S12>/EtherCAT PDO Receive7' */
  uint8_T CastToDouble1_f;             /* '<S12>/Cast To Double1' */
  uint8_T CastToDouble2_d;             /* '<S12>/Cast To Double2' */
  uint8_T absEncoderStatus1;           /* '<S39>/absEncoderStatus1' */
  uint8_T absEncoderStatus2;           /* '<S39>/absEncoderStatus2' */
  uint8_T speedCtrlReset;              /* '<S33>/speedCtrlReset' */
  boolean_T aboveLimit;                /* '<S49>/aboveLimit' */
  boolean_T alarm;                     /* '<S49>/alarm' */
  boolean_T atSetpoint;                /* '<S49>/atSetpoint' */
  boolean_T commErr;                   /* '<S49>/commErr' */
  boolean_T extCtrlLoc;                /* '<S49>/extCtrlLoc' */
  boolean_T extRunEnable;              /* '<S49>/extRunEnable' */
  boolean_T mswB13;                    /* '<S49>/mswB13' */
  boolean_T mswB14;                    /* '<S49>/mswB14' */
  boolean_T off2;                      /* '<S49>/off2' */
  boolean_T off3;                      /* '<S49>/off3' */
  boolean_T rdyOn;                     /* '<S49>/rdyOn' */
  boolean_T rdyRef;                    /* '<S49>/rdyRef' */
  boolean_T rdyRun;                    /* '<S49>/rdyRun' */
  boolean_T remote;                    /* '<S49>/remote' */
  boolean_T switchOnInhibit;           /* '<S49>/switchOnInhibit' */
  boolean_T tripped;                   /* '<S49>/tripped' */
  boolean_T Memory_b;                  /* '<S2>/Memory' */
  boolean_T NotEqual;                  /* '<S2>/NotEqual' */
  boolean_T Memory1;                   /* '<S2>/Memory1' */
  boolean_T NotEqual1;                 /* '<S2>/NotEqual1' */
  boolean_T Memory2;                   /* '<S2>/Memory2' */
  boolean_T NotEqual2;                 /* '<S2>/NotEqual2' */
  boolean_T Memory_i;                  /* '<S4>/Memory' */
  boolean_T NotEqual_i;                /* '<S4>/NotEqual' */
  boolean_T ACS880CtrlMode;            /* '<S2>/ACS880CtrlMode' */
  boolean_T Equal;                     /* '<S4>/Equal' */
  boolean_T Equal1;                    /* '<S4>/Equal1' */
  boolean_T Memory1_n;                 /* '<S4>/Memory1' */
  boolean_T NotEqual1_a;               /* '<S4>/NotEqual1' */
  boolean_T Memory2_c;                 /* '<S4>/Memory2' */
  boolean_T NotEqual2_c;               /* '<S4>/NotEqual2' */
  boolean_T Outofbounds;               /* '<S609>/Out of bounds' */
  boolean_T LowerRelop1;               /* '<S598>/LowerRelop1' */
  boolean_T UpperRelop;                /* '<S598>/UpperRelop' */
  boolean_T enableOperation;           /* '<S46>/enableOperation' */
  boolean_T extCtrlLoc_f;              /* '<S46>/extCtrlLoc' */
  boolean_T inching1;                  /* '<S46>/inching1' */
  boolean_T inching2;                  /* '<S46>/inching2' */
  boolean_T off1Ctrl;                  /* '<S46>/off1Ctrl' */
  boolean_T off2Ctrl;                  /* '<S46>/off2Ctrl' */
  boolean_T off3Ctrl;                  /* '<S46>/off3Ctrl' */
  boolean_T rampHold;                  /* '<S46>/rampHold' */
  boolean_T rampInZero;                /* '<S46>/rampInZero' */
  boolean_T rampOutZero;               /* '<S46>/rampOutZero' */
  boolean_T remoteCmd;                 /* '<S46>/remoteCmd' */
  boolean_T reset;                     /* '<S46>/reset' */
  boolean_T aboveLimit_c;              /* '<S44>/aboveLimit' */
  boolean_T alarm_c;                   /* '<S44>/alarm' */
  boolean_T atSetpoint_a;              /* '<S44>/atSetpoint' */
  boolean_T commErr_i;                 /* '<S44>/commErr' */
  boolean_T extCtrlLoc_e;              /* '<S44>/extCtrlLoc' */
  boolean_T extRunEnable_k;            /* '<S44>/extRunEnable' */
  boolean_T mswB13_a;                  /* '<S44>/mswB13' */
  boolean_T mswB14_m;                  /* '<S44>/mswB14' */
  boolean_T off2_n;                    /* '<S44>/off2' */
  boolean_T off3_p;                    /* '<S44>/off3' */
  boolean_T rdyOn_b;                   /* '<S44>/rdyOn' */
  boolean_T rdyRef_f;                  /* '<S44>/rdyRef' */
  boolean_T rdyRun_n;                  /* '<S44>/rdyRun' */
  boolean_T remote_d;                  /* '<S44>/remote' */
  boolean_T switchOnInhibit_p;         /* '<S44>/switchOnInhibit' */
  boolean_T tripped_p;                 /* '<S44>/tripped' */
  boolean_T Memory_j;                  /* '<S1>/Memory' */
  boolean_T NotEqual_f;                /* '<S1>/NotEqual' */
  boolean_T Memory1_f;                 /* '<S1>/Memory1' */
  boolean_T NotEqual1_m;               /* '<S1>/NotEqual1' */
  boolean_T Memory2_h;                 /* '<S1>/Memory2' */
  boolean_T NotEqual2_h;               /* '<S1>/NotEqual2' */
  boolean_T ACS800CtrlMode;            /* '<S1>/ACS800CtrlMode' */
  boolean_T enableOperation_c;         /* '<S41>/enableOperation' */
  boolean_T extCtrlLoc_ec;             /* '<S41>/extCtrlLoc' */
  boolean_T inching1_o;                /* '<S41>/inching1' */
  boolean_T inching2_j;                /* '<S41>/inching2' */
  boolean_T off1Ctrl_a;                /* '<S41>/off1Ctrl' */
  boolean_T off2Ctrl_h;                /* '<S41>/off2Ctrl' */
  boolean_T off3Ctrl_a;                /* '<S41>/off3Ctrl' */
  boolean_T rampHold_p;                /* '<S41>/rampHold' */
  boolean_T rampInZero_e;              /* '<S41>/rampInZero' */
  boolean_T rampOutZero_p;             /* '<S41>/rampOutZero' */
  boolean_T remoteCmd_p;               /* '<S41>/remoteCmd' */
  boolean_T reset_b;                   /* '<S41>/reset' */
  boolean_T resetHilIntegrator;        /* '<S31>/resetHilIntegrator' */
  boolean_T resetSidIntegrator;        /* '<S31>/resetSidIntegrator' */
  boolean_T runHil;                    /* '<S31>/runHil' */
  boolean_T runSid;                    /* '<S31>/runSid' */
  boolean_T RelationalOperator;        /* '<S552>/Relational Operator' */
  boolean_T L1InaccurateURead;         /* '<S9>/L1InaccurateURead' */
  boolean_T L1InaccurateU;             /* '<S9>/L1InaccurateU' */
  boolean_T L1InaccurateIRead;         /* '<S9>/L1InaccurateIRead' */
  boolean_T L1InaccurateI;             /* '<S9>/L1InaccurateI' */
  boolean_T L2InaccurateURead;         /* '<S9>/L2InaccurateURead' */
  boolean_T L2InaccurateU;             /* '<S9>/L2InaccurateU' */
  boolean_T L2InaccurateIRead;         /* '<S9>/L2InaccurateIRead' */
  boolean_T L2InaccurateI;             /* '<S9>/L2InaccurateI' */
  boolean_T L3InaccurateURead;         /* '<S9>/L3InaccurateURead' */
  boolean_T L3InaccurateU;             /* '<S9>/L3InaccurateU' */
  boolean_T L3InaccurateIRead;         /* '<S9>/L3InaccurateIRead' */
  boolean_T L3InaccurateI;             /* '<S9>/L3InaccurateI' */
  boolean_T L1InaccurateURead_l;       /* '<S11>/L1InaccurateURead' */
  boolean_T L1InaccurateU_f;           /* '<S11>/L1InaccurateU' */
  boolean_T L1InaccurateIRead_o;       /* '<S11>/L1InaccurateIRead' */
  boolean_T L1InaccurateI_f;           /* '<S11>/L1InaccurateI' */
  boolean_T L2InaccurateURead_h;       /* '<S11>/L2InaccurateURead' */
  boolean_T L2InaccurateU_p;           /* '<S11>/L2InaccurateU' */
  boolean_T L2InaccurateIRead_f;       /* '<S11>/L2InaccurateIRead' */
  boolean_T L2InaccurateI_p;           /* '<S11>/L2InaccurateI' */
  boolean_T L3InaccurateURead_o;       /* '<S11>/L3InaccurateURead' */
  boolean_T L3InaccurateU_n;           /* '<S11>/L3InaccurateU' */
  boolean_T L3InaccurateIRead_p;       /* '<S11>/L3InaccurateIRead' */
  boolean_T L3InaccurateI_m;           /* '<S11>/L3InaccurateI' */
  boolean_T Constant;                  /* '<S5>/Constant' */
  boolean_T LessThan;                  /* '<S68>/Less Than' */
  boolean_T LessThan_k;                /* '<S147>/Less Than' */
  boolean_T RelationalOperator_l;      /* '<S431>/Relational Operator' */
  boolean_T Memory_f;                  /* '<S433>/Memory' */
  boolean_T Logic[2];                  /* '<S433>/Logic' */
  boolean_T RelationalOperator_k;      /* '<S432>/Relational Operator' */
  boolean_T Memory_in;                 /* '<S434>/Memory' */
  boolean_T Logic_g[2];                /* '<S434>/Logic' */
  boolean_T LowerRelop1_g;             /* '<S420>/LowerRelop1' */
  boolean_T UpperRelop_g;              /* '<S420>/UpperRelop' */
  boolean_T RelationalOperator_e;      /* '<S494>/Relational Operator' */
  boolean_T Memory_o;                  /* '<S496>/Memory' */
  boolean_T Logic_c[2];                /* '<S496>/Logic' */
  boolean_T RelationalOperator_g;      /* '<S495>/Relational Operator' */
  boolean_T Memory_a;                  /* '<S497>/Memory' */
  boolean_T Logic_p[2];                /* '<S497>/Logic' */
  boolean_T LowerRelop1_h;             /* '<S483>/LowerRelop1' */
  boolean_T UpperRelop_m;              /* '<S483>/UpperRelop' */
  boolean_T resetHilIntegrator_h;      /* '<S4>/FexcRamp' */
  boolean_T resetSidIntegrator_h;      /* '<S4>/FexcRamp' */
  B_YawKinematicTransforms_wind_T sf_YawKinematicTransforms_l;/* '<S212>/Yaw Kinematic Transforms' */
  B_YawForceTransforms_windEmul_T sf_YawForceTransforms_i;/* '<S149>/Yaw Force Transforms' */
  B_MATLABFunction1_windEmulato_T sf_MATLABFunction1_e;/* '<S205>/MATLAB Function1' */
  B_quaternion2EulXYZ_windEmula_T sf_quaternion2EulXYZ_c;/* '<S160>/quaternion2EulXYZ' */
  B_NonlinearWaveElevation_wind_T NonlinearWaveElevation_j;/* '<S139>/Nonlinear Wave Elevation' */
  B_YawKinematicTransforms_wind_T sf_YawKinematicTransforms;/* '<S133>/Yaw Kinematic Transforms' */
  B_YawForceTransforms_windEmul_T sf_YawForceTransforms;/* '<S70>/Yaw Force Transforms' */
  B_MATLABFunction1_windEmulato_T sf_MATLABFunction1;/* '<S126>/MATLAB Function1' */
  B_quaternion2EulXYZ_windEmula_T sf_quaternion2EulXYZ;/* '<S81>/quaternion2EulXYZ' */
  B_NonlinearWaveElevation_wind_T NonlinearWaveElevation;/* '<S60>/Nonlinear Wave Elevation' */
  B_MovingAverage_windEmulato_c_T MovingAverage1;/* '<S51>/Moving Average' */
  B_MovingAverage_windEmulato_c_T MovingAverage_pn;/* '<S51>/Moving Average' */
  B_ParseStatusWord_windEmulato_T sf_ParseStatusWord_h;/* '<S49>/Parse Status Word' */
  B_MovingAverage_windEmulatorS_T MovingAverage_p;/* '<S43>/Moving Average' */
  B_parseCtrlWord_windEmulatorS_T sf_parseCtrlWord_h;/* '<S46>/parseCtrlWord' */
  B_ParseStatusWord_windEmulato_T sf_ParseStatusWord;/* '<S44>/Parse Status Word' */
  B_MovingAverage_windEmulatorS_T MovingAverage;/* '<S43>/Moving Average' */
  B_parseCtrlWord_windEmulatorS_T sf_parseCtrlWord;/* '<S41>/parseCtrlWord' */
};

/* Block states (default storage) for system '<Root>' */
struct DW_windEmulatorStep4_WECSim_T {
  real_T INPUT_1_1_1_Discrete_142920204[2];/* '<S541>/INPUT_1_1_1' */
  real_T INPUT_2_1_1_Discrete_1327804636[2];/* '<S541>/INPUT_2_1_1' */
  real_T INPUT_4_1_1_Discrete_3227796860[2];/* '<S541>/INPUT_4_1_1' */
  real_T INPUT_5_1_1_Discrete_4244925644[2];/* '<S541>/INPUT_5_1_1' */
  real_T INPUT_3_1_1_Discrete_1917098348[2];/* '<S541>/INPUT_3_1_1' */
  real_T STATE_1_Discrete_1041191992[22];/* '<S541>/STATE_1' */
  real_T DiscreteTimeIntegrator_DSTATE;/* '<S437>/Discrete-Time Integrator' */
  real_T DiscreteTimeIntegrator_DSTATE_n;/* '<S372>/Discrete-Time Integrator' */
  real_T DiscreteTimeIntegrator_DSTATE_l;/* '<S435>/Discrete-Time Integrator' */
  real_T UD_DSTATE;                    /* '<S553>/UD' */
  real_T INPUT_2_1_1_Discrete_3353911369[2];/* '<S332>/INPUT_2_1_1' */
  real_T INPUT_3_1_1_Discrete_4203252217;/* '<S332>/INPUT_3_1_1' */
  real_T INPUT_3_1_1_FirstOutput_4203252;/* '<S332>/INPUT_3_1_1' */
  real_T DelayOneStep_DSTATE;          /* '<S58>/Delay One Step' */
  real_T Integrator_DSTATE;            /* '<S282>/Integrator' */
  real_T UD_DSTATE_j;                  /* '<S275>/UD' */
  real_T INPUT_1_1_1_Discrete_2152258201[2];/* '<S332>/INPUT_1_1_1' */
  real_T STATE_1_Discrete_208214823[6];/* '<S332>/STATE_1' */
  real_T INPUT_5_1_1_Discrete_4076664036[2];/* '<S216>/INPUT_5_1_1' */
  real_T INPUT_1_1_1_Discrete_125588004[2];/* '<S216>/INPUT_1_1_1' */
  real_T INPUT_1_1_2_Discrete_2658468766[2];/* '<S216>/INPUT_1_1_2' */
  real_T INPUT_1_1_3_Discrete_3916575496[2];/* '<S216>/INPUT_1_1_3' */
  real_T INPUT_2_1_1_Discrete_1088170228[2];/* '<S216>/INPUT_2_1_1' */
  real_T INPUT_2_1_2_Discrete_3654646094[2];/* '<S216>/INPUT_2_1_2' */
  real_T INPUT_2_1_3_Discrete_2933017048[2];/* '<S216>/INPUT_2_1_3' */
  real_T INPUT_3_1_1_Discrete_2109473092[2];/* '<S216>/INPUT_3_1_1' */
  real_T INPUT_3_1_2_Discrete_3837087998[2];/* '<S216>/INPUT_3_1_2' */
  real_T INPUT_3_1_3_Discrete_2477940840[2];/* '<S216>/INPUT_3_1_3' */
  real_T INPUT_4_1_1_Discrete_3483163988[2];/* '<S216>/INPUT_4_1_1' */
  real_T INPUT_4_1_2_Discrete_1452641518[2];/* '<S216>/INPUT_4_1_2' */
  real_T INPUT_4_1_3_Discrete_563264632[2];/* '<S216>/INPUT_4_1_3' */
  real_T Integrator_DSTATE_d;          /* '<S412>/Integrator' */
  real_T Filter_DSTATE;                /* '<S407>/Filter' */
  real_T Integrator_DSTATE_e;          /* '<S475>/Integrator' */
  real_T Filter_DSTATE_b;              /* '<S470>/Filter' */
  real_T PrevY;                        /* '<S609>/torqueSlewRate' */
  real_T LastMajorTime;                /* '<S609>/torqueSlewRate' */
  real_T PrevY_a;                      /* '<S609>/speedSlewRate' */
  real_T LastMajorTime_j;              /* '<S609>/speedSlewRate' */
  real_T PrevY_k;                      /* '<S498>/Rate Limiter' */
  real_T PrevY_f;                      /* '<S7>/Rate Limiter' */
  real_T STATE_1_ZcValueStore;         /* '<S541>/STATE_1' */
  real_T OUTPUT_1_0_Discrete;          /* '<S541>/OUTPUT_1_0' */
  real_T OUTPUT_1_0_ZcValueStore;      /* '<S541>/OUTPUT_1_0' */
  real_T PrevY_b;                      /* '<S2>/acs880RateLim' */
  real_T LastMajorTime_a;              /* '<S2>/acs880RateLim' */
  real_T STATE_1_Discrete;             /* '<S216>/STATE_1' */
  real_T STATE_1_ZcValueStore_d;       /* '<S216>/STATE_1' */
  real_T OUTPUT_1_1_Discrete;          /* '<S216>/OUTPUT_1_1' */
  real_T OUTPUT_1_1_ZcValueStore;      /* '<S216>/OUTPUT_1_1' */
  real_T STATE_1_ZcValueStore_k;       /* '<S332>/STATE_1' */
  real_T OUTPUT_1_0_Discrete_k;        /* '<S332>/OUTPUT_1_0' */
  real_T OUTPUT_1_0_ZcValueStore_l;    /* '<S332>/OUTPUT_1_0' */
  real_T OUTPUT_1_0_Discrete_p;        /* '<S216>/OUTPUT_1_0' */
  real_T OUTPUT_1_0_ZcValueStore_f;    /* '<S216>/OUTPUT_1_0' */
  real_T PrevY_m;                      /* '<S373>/Rate Limiter1' */
  real_T PrevY_l;                      /* '<S373>/Rate Limiter' */
  real_T LastMajorTime_d;              /* '<S373>/Rate Limiter' */
  real_T PrevY_fq;                     /* '<S436>/Rate Limiter1' */
  real_T PrevY_e;                      /* '<S436>/Rate Limiter' */
  real_T LastMajorTime_k;              /* '<S436>/Rate Limiter' */
  real_T rampLast;                     /* '<S4>/FexcRamp' */
  real_T rampUpTime;                   /* '<S4>/FexcRamp' */
  real_T runTime;                      /* '<S4>/FexcRamp' */
  real_T swRDY_ON;                     /* '<S18>/ABB Fieldbus Control' */
  real_T swRDY_RUN;                    /* '<S18>/ABB Fieldbus Control' */
  real_T swRDY_REF;                    /* '<S18>/ABB Fieldbus Control' */
  real_T swTRIPPED;                    /* '<S18>/ABB Fieldbus Control' */
  real_T swOFF_2_STA;                  /* '<S18>/ABB Fieldbus Control' */
  real_T swOFF_3_STA;                  /* '<S18>/ABB Fieldbus Control' */
  real_T swSWC_ON_INHIB;               /* '<S18>/ABB Fieldbus Control' */
  real_T swAT_SETPOINT;                /* '<S18>/ABB Fieldbus Control' */
  real_T swEXT_RUN_ENABLE;             /* '<S18>/ABB Fieldbus Control' */
  real_T cwOFF2_CONTROL;               /* '<S18>/ABB Fieldbus Control' */
  real_T cwOFF3_CONTROL;               /* '<S18>/ABB Fieldbus Control' */
  real_T cwENABLE_OPERATION;           /* '<S18>/ABB Fieldbus Control' */
  real_T cwRAMP_OUT_ZERO;              /* '<S18>/ABB Fieldbus Control' */
  real_T cwRAMP_HOLD;                  /* '<S18>/ABB Fieldbus Control' */
  real_T cwRAMP_IN_ZERO;               /* '<S18>/ABB Fieldbus Control' */
  real_T cwRESET;                      /* '<S18>/ABB Fieldbus Control' */
  real_T cwREMOTE_CMD;                 /* '<S18>/ABB Fieldbus Control' */
  real_T cwOFF1_CONTROL;               /* '<S18>/ABB Fieldbus Control' */
  real_T swEXT_CTRL_LOC;               /* '<S18>/ABB Fieldbus Control' */
  real_T swREMOTE;                     /* '<S18>/ABB Fieldbus Control' */
  real_T swWARNING;                    /* '<S18>/ABB Fieldbus Control' */
  real_T swABOVE_LIMIT;                /* '<S18>/ABB Fieldbus Control' */
  real_T swMSW_B13;                    /* '<S18>/ABB Fieldbus Control' */
  real_T swMSW_B14;                    /* '<S18>/ABB Fieldbus Control' */
  real_T swCOMM_ERR;                   /* '<S18>/ABB Fieldbus Control' */
  real_T swRDY_ON_f;                   /* '<S16>/ABB Fieldbus Control' */
  real_T swRDY_RUN_f;                  /* '<S16>/ABB Fieldbus Control' */
  real_T swRDY_REF_a;                  /* '<S16>/ABB Fieldbus Control' */
  real_T swTRIPPED_e;                  /* '<S16>/ABB Fieldbus Control' */
  real_T swOFF_2_STA_c;                /* '<S16>/ABB Fieldbus Control' */
  real_T swOFF_3_STA_i;                /* '<S16>/ABB Fieldbus Control' */
  real_T swSWC_ON_INHIB_l;             /* '<S16>/ABB Fieldbus Control' */
  real_T swAT_SETPOINT_b;              /* '<S16>/ABB Fieldbus Control' */
  real_T swEXT_RUN_ENABLE_h;           /* '<S16>/ABB Fieldbus Control' */
  real_T cwOFF2_CONTROL_b;             /* '<S16>/ABB Fieldbus Control' */
  real_T cwOFF3_CONTROL_n;             /* '<S16>/ABB Fieldbus Control' */
  real_T cwENABLE_OPERATION_l;         /* '<S16>/ABB Fieldbus Control' */
  real_T cwRAMP_OUT_ZERO_m;            /* '<S16>/ABB Fieldbus Control' */
  real_T cwRAMP_HOLD_j;                /* '<S16>/ABB Fieldbus Control' */
  real_T cwRAMP_IN_ZERO_m;             /* '<S16>/ABB Fieldbus Control' */
  real_T cwRESET_e;                    /* '<S16>/ABB Fieldbus Control' */
  real_T cwREMOTE_CMD_c;               /* '<S16>/ABB Fieldbus Control' */
  real_T cwOFF1_CONTROL_i;             /* '<S16>/ABB Fieldbus Control' */
  real_T swEXT_CTRL_LOC_c;             /* '<S16>/ABB Fieldbus Control' */
  real_T swREMOTE_j;                   /* '<S16>/ABB Fieldbus Control' */
  real_T swWARNING_k;                  /* '<S16>/ABB Fieldbus Control' */
  real_T swABOVE_LIMIT_a;              /* '<S16>/ABB Fieldbus Control' */
  real_T swMSW_B13_l;                  /* '<S16>/ABB Fieldbus Control' */
  real_T swMSW_B14_i;                  /* '<S16>/ABB Fieldbus Control' */
  real_T swCOMM_ERR_p;                 /* '<S16>/ABB Fieldbus Control' */
  struct {
    real_T EXECRATIO;
  } EtherCATInit_RWORK;                /* '<Root>/EtherCAT Init' */

  real_T TransportDelay_RWORK[12289];  /* '<S60>/Transport Delay' */
  real_T TransportDelay_RWORK_k[12289];/* '<S139>/Transport Delay' */
  real_T TransportDelay_RWORK_l[12289];/* '<S64>/Transport Delay' */
  real_T TransportDelay_RWORK_f[12289];/* '<S143>/Transport Delay' */
  struct {
    void *AQHandles;
    void *SLRTSigHandles;
  } TAQSigLogging_InsertedFor_acs88;   /* synthesized block */

  void* RTP_1_RtpManager;              /* '<S508>/RTP_1' */
  void* STATE_1_Simulator;             /* '<S541>/STATE_1' */
  void* STATE_1_SimData;               /* '<S541>/STATE_1' */
  void* STATE_1_DiagMgr;               /* '<S541>/STATE_1' */
  void* STATE_1_ZcLogger;              /* '<S541>/STATE_1' */
  void* STATE_1_TsInfo;                /* '<S541>/STATE_1' */
  void* OUTPUT_1_0_Simulator;          /* '<S541>/OUTPUT_1_0' */
  void* OUTPUT_1_0_SimData;            /* '<S541>/OUTPUT_1_0' */
  void* OUTPUT_1_0_DiagMgr;            /* '<S541>/OUTPUT_1_0' */
  void* OUTPUT_1_0_ZcLogger;           /* '<S541>/OUTPUT_1_0' */
  void* OUTPUT_1_0_TsInfo;             /* '<S541>/OUTPUT_1_0' */
  struct {
    void *AQHandles;
    void *SLRTSigHandles;
  } TAQSigLogging_InsertedFor_acs_a;   /* synthesized block */

  struct {
    void *AQHandles;
    void *SLRTSigHandles;
  } TAQSigLogging_InsertedFor_acs80;   /* synthesized block */

  struct {
    void *AQHandles;
    void *SLRTSigHandles;
  } TAQSigLogging_InsertedFor_acs_l;   /* synthesized block */

  struct {
    void *AQHandles;
    void *SLRTSigHandles;
  } TAQSigLogging_InsertedFor_hptoS;   /* synthesized block */

  struct {
    void *AQHandles;
    void *SLRTSigHandles;
  } TAQSigLogging_InsertedFor_hptoC;   /* synthesized block */

  struct {
    void *AQHandles;
    void *SLRTSigHandles;
  } TAQSigLogging_InsertedFor_expCt;   /* synthesized block */

  struct {
    void *AQHandles;
    void *SLRTSigHandles;
  } TAQSigLogging_InsertedFor_shaft;   /* synthesized block */

  struct {
    void *AQHandles;
    void *SLRTSigHandles;
  } TAQSigLogging_InsertedFor_invPo;   /* synthesized block */

  struct {
    void *AQHandles;
    void *SLRTSigHandles;
  } TAQSigLogging_InsertedFor_inv_p;   /* synthesized block */

  struct {
    void *AQHandles;
    void *SLRTSigHandles;
  } TAQSigLogging_InsertedFor_sidIn;   /* synthesized block */

  void* STATE_1_Simulator_f;           /* '<S216>/STATE_1' */
  void* STATE_1_SimData_h;             /* '<S216>/STATE_1' */
  void* STATE_1_DiagMgr_o;             /* '<S216>/STATE_1' */
  void* STATE_1_ZcLogger_b;            /* '<S216>/STATE_1' */
  void* STATE_1_TsInfo_d;              /* '<S216>/STATE_1' */
  void* OUTPUT_1_1_Simulator;          /* '<S216>/OUTPUT_1_1' */
  void* OUTPUT_1_1_SimData;            /* '<S216>/OUTPUT_1_1' */
  void* OUTPUT_1_1_DiagMgr;            /* '<S216>/OUTPUT_1_1' */
  void* OUTPUT_1_1_ZcLogger;           /* '<S216>/OUTPUT_1_1' */
  void* OUTPUT_1_1_TsInfo;             /* '<S216>/OUTPUT_1_1' */
  struct {
    void *LoggedData;
  } Scope2_PWORK;                      /* '<S58>/Scope2' */

  void* STATE_1_Simulator_i;           /* '<S332>/STATE_1' */
  void* STATE_1_SimData_a;             /* '<S332>/STATE_1' */
  void* STATE_1_DiagMgr_g;             /* '<S332>/STATE_1' */
  void* STATE_1_ZcLogger_i;            /* '<S332>/STATE_1' */
  void* STATE_1_TsInfo_g;              /* '<S332>/STATE_1' */
  void* OUTPUT_1_0_Simulator_d;        /* '<S332>/OUTPUT_1_0' */
  void* OUTPUT_1_0_SimData_l;          /* '<S332>/OUTPUT_1_0' */
  void* OUTPUT_1_0_DiagMgr_i;          /* '<S332>/OUTPUT_1_0' */
  void* OUTPUT_1_0_ZcLogger_m;         /* '<S332>/OUTPUT_1_0' */
  void* OUTPUT_1_0_TsInfo_c;           /* '<S332>/OUTPUT_1_0' */
  struct {
    void *LoggedData;
  } Scope5_PWORK;                      /* '<S58>/Scope5' */

  struct {
    void *LoggedData;
  } Scope7_PWORK;                      /* '<S58>/Scope7' */

  struct {
    void *AQHandles;
  } TAQSigLogging_InsertedFor_Produ;   /* synthesized block */

  struct {
    void *AQHandles;
  } TAQSigLogging_InsertedFor_Pro_j;   /* synthesized block */

  struct {
    void *AQHandles;
  } TAQSigLogging_InsertedFor_Selec;   /* synthesized block */

  struct {
    void *LoggedData;
  } Scope3_PWORK;                      /* '<S58>/Scope3' */

  struct {
    void *LoggedData;
  } Scope6_PWORK;                      /* '<S58>/Scope6' */

  struct {
    void *LoggedData;
  } Scope8_PWORK;                      /* '<S58>/Scope8' */

  struct {
    void *AQHandles;
  } TAQSigLogging_InsertedFor_PSSim;   /* synthesized block */

  struct {
    void *AQHandles;
  } TAQSigLogging_InsertedFor_Pro_h;   /* synthesized block */

  struct {
    void *LoggedData;
  } Scope_PWORK;                       /* '<S240>/Scope' */

  struct {
    void *AQHandles;
  } TAQSigLogging_InsertedFor_PSS_h;   /* synthesized block */

  struct {
    void *LoggedData;
  } Scope_PWORK_a;                     /* '<S58>/Scope' */

  struct {
    void *LoggedData;
  } Scope1_PWORK;                      /* '<S58>/Scope1' */

  struct {
    void *AQHandles;
  } TAQSigLogging_InsertedFor_Discr;   /* synthesized block */

  struct {
    void *AQHandles;
  } TAQSigLogging_InsertedFor_PS_h4;   /* synthesized block */

  struct {
    void *AQHandles;
  } TAQSigLogging_InsertedFor_Pro_p;   /* synthesized block */

  struct {
    void *AQHandles;
  } TAQSigLogging_InsertedFor_Pr_hl;   /* synthesized block */

  struct {
    void *AQHandles;
  } TAQSigLogging_InsertedFor_Pro_k;   /* synthesized block */

  struct {
    void *LoggedData;
  } Scope4_PWORK;                      /* '<S58>/Scope4' */

  struct {
    void *AQHandles;
  } TAQSigLogging_InsertedFor_Gain_;   /* synthesized block */

  struct {
    void *AQHandles;
  } TAQSigLogging_InsertedFor_Gai_a;   /* synthesized block */

  struct {
    void *AQHandles;
  } TAQSigLogging_InsertedFor_Ga_ad;   /* synthesized block */

  struct {
    void *AQHandles;
  } TAQSigLogging_InsertedFor_G_ado;   /* synthesized block */

  struct {
    void *AQHandles;
  } TAQSigLogging_InsertedFor__adom;   /* synthesized block */

  struct {
    void *AQHandles;
  } TAQSigLogging_InsertedFor_adom5;   /* synthesized block */

  struct {
    void *AQHandles;
  } TAQSigLogging_InsertedFo_adom5n;   /* synthesized block */

  struct {
    void *AQHandles;
  } TAQSigLogging_InsertedF_adom5n1;   /* synthesized block */

  struct {
    void *AQHandles;
  } TAQSigLogging_Inserted_adom5n1n;   /* synthesized block */

  struct {
    void *AQHandles;
  } TAQSigLogging_Inserte_adom5n1nb;   /* synthesized block */

  struct {
    void *AQHandles;
  } TAQSigLogging_Insert_adom5n1nbj;   /* synthesized block */

  struct {
    void *AQHandles;
  } TAQSigLogging_Inser_adom5n1nbj0;   /* synthesized block */

  struct {
    void *LoggedData;
  } Scope_PWORK_d;                     /* '<S241>/Scope' */

  struct {
    void *AQHandles;
  } TAQSigLogging_Inse_adom5n1nbj0i;   /* synthesized block */

  struct {
    void *LoggedData;
  } Scope_PWORK_dk;                    /* '<S242>/Scope' */

  struct {
    void *AQHandles;
  } TAQSigLogging_Ins_adom5n1nbj0ib;   /* synthesized block */

  struct {
    void *AQHandles;
  } TAQSigLogging_In_adom5n1nbj0ibk;   /* synthesized block */

  struct {
    void *LoggedData;
  } Scope_PWORK_k;                     /* '<S244>/Scope' */

  struct {
    void *AQHandles;
  } TAQSigLogging_I_adom5n1nbj0ibkg;   /* synthesized block */

  struct {
    void *LoggedData;
  } Scope_PWORK_c;                     /* '<S245>/Scope' */

  struct {
    void *AQHandles;
  } TAQSigLogging_InsertedFor_Gai_f;   /* synthesized block */

  void *TransportDelay_PWORK[12];      /* '<S60>/Transport Delay' */
  void *TransportDelay_PWORK_f[12];    /* '<S139>/Transport Delay' */
  void* OUTPUT_1_0_Simulator_l;        /* '<S216>/OUTPUT_1_0' */
  void* OUTPUT_1_0_SimData_i;          /* '<S216>/OUTPUT_1_0' */
  void* OUTPUT_1_0_DiagMgr_a;          /* '<S216>/OUTPUT_1_0' */
  void* OUTPUT_1_0_ZcLogger_i;         /* '<S216>/OUTPUT_1_0' */
  void* OUTPUT_1_0_TsInfo_n;           /* '<S216>/OUTPUT_1_0' */
  struct {
    void *LoggedData;
  } ToWorkspace_PWORK;                 /* '<S59>/To Workspace' */

  struct {
    void *LoggedData;
  } ToWorkspace_PWORK_m;               /* '<S60>/To Workspace' */

  void *TransportDelay_PWORK_j[12];    /* '<S64>/Transport Delay' */
  struct {
    void *LoggedData;
  } ToWorkspace_PWORK_p;               /* '<S55>/To Workspace' */

  struct {
    void *LoggedData;
  } ToWorkspace_PWORK_h;               /* '<S139>/To Workspace' */

  void *TransportDelay_PWORK_n[12];    /* '<S143>/Transport Delay' */
  struct {
    void *LoggedData;
  } ToWorkspace_PWORK_b;               /* '<S57>/To Workspace' */

  void* SINK_1_RtwLogger;              /* '<S216>/SINK_1' */
  void* SINK_1_RtwLogBuffer;           /* '<S216>/SINK_1' */
  void* SINK_1_RtwLogFcnManager;       /* '<S216>/SINK_1' */
  void* SINK_1_InstRtwLogger;          /* '<S216>/SINK_1' */
  void* SINK_1_InstRtwLogBuffer;       /* '<S216>/SINK_1' */
  struct {
    void *LoggedData;
  } Scope_PWORK_l;                     /* '<S499>/Scope' */

  struct {
    void *LoggedData;
  } Scope_PWORK_e;                     /* '<S504>/Scope' */

  struct {
    void *LoggedData;
  } Scope_PWORK_b;                     /* '<S500>/Scope' */

  int32_T lastRawCounts_PreviousInput; /* '<S552>/lastRawCounts' */
  int32_T lastTurn_PreviousInput;      /* '<S552>/lastTurn' */
  int32_T sfEvent;                     /* '<S4>/FexcRamp' */
  int32_T sfEvent_l;                   /* '<S18>/ABB Fieldbus Control' */
  int32_T sfEvent_d;                   /* '<S16>/ABB Fieldbus Control' */
  uint32_T Memory_PreviousInput;       /* '<S29>/Memory' */
  uint32_T is_c3_windEmulatorStep4_WECSim;/* '<S4>/FexcRamp' */
  uint32_T temporalCounter_i1;         /* '<S4>/FexcRamp' */
  uint32_T is_UpdateStateMachine;      /* '<S18>/ABB Fieldbus Control' */
  uint32_T is_UpdateStateMachine_g;    /* '<S16>/ABB Fieldbus Control' */
  int_T EtherCATPDOReceive1_IWORK[7];  /* '<S10>/EtherCAT PDO Receive1' */
  int_T EtherCATPDOReceive2_IWORK[7];  /* '<S10>/EtherCAT PDO Receive2' */
  int_T EtherCATPDOReceive3_IWORK[7];  /* '<S10>/EtherCAT PDO Receive3' */
  int_T EtherCATPDOReceive4_IWORK[7];  /* '<S10>/EtherCAT PDO Receive4' */
  int_T EtherCATPDOReceive5_IWORK[7];  /* '<S10>/EtherCAT PDO Receive5' */
  int_T EtherCATPDOReceive6_IWORK[7];  /* '<S10>/EtherCAT PDO Receive6' */
  int_T EtherCATPDOReceive_IWORK[7];   /* '<S10>/EtherCAT PDO Receive' */
  int_T EtherCATPDOReceive8_IWORK[7];  /* '<S8>/EtherCAT PDO Receive8' */
  int_T EtherCATPDOReceive9_IWORK[7];  /* '<S8>/EtherCAT PDO Receive9' */
  int_T EtherCATPDOReceive10_IWORK[7]; /* '<S8>/EtherCAT PDO Receive10' */
  int_T EtherCATPDOReceive11_IWORK[7]; /* '<S8>/EtherCAT PDO Receive11' */
  int_T EtherCATPDOReceive12_IWORK[7]; /* '<S8>/EtherCAT PDO Receive12' */
  int_T EtherCATPDOReceive13_IWORK[7]; /* '<S8>/EtherCAT PDO Receive13' */
  int_T EtherCATPDOReceive7_IWORK[7];  /* '<S8>/EtherCAT PDO Receive7' */
  int_T STATE_1_Modes[15];             /* '<S541>/STATE_1' */
  int_T OUTPUT_1_0_Modes;              /* '<S541>/OUTPUT_1_0' */
  int_T readTorqueInput_IWORK[7];      /* '<S12>/readTorqueInput' */
  int_T readEncoderCounter_IWORK[7];   /* '<S12>/readEncoderCounter' */
  int_T EtherCATPDOReceive7_IWORK_l[7];/* '<S12>/EtherCAT PDO Receive7' */
  int_T L1InaccurateURead_IWORK[7];    /* '<S9>/L1InaccurateURead' */
  int_T L1InaccurateIRead_IWORK[7];    /* '<S9>/L1InaccurateIRead' */
  int_T L1VoltageRead_IWORK[7];        /* '<S9>/L1VoltageRead' */
  int_T L1CurrentRead_IWORK[7];        /* '<S9>/L1CurrentRead' */
  int_T L1PowFactorRead_IWORK[7];      /* '<S9>/L1PowFactorRead' */
  int_T L1ActivePowRead_IWORK[7];      /* '<S9>/L1ActivePowRead' */
  int_T L1THDuRead_IWORK[7];           /* '<S9>/L1THDuRead' */
  int_T L1THDiRead_IWORK[7];           /* '<S9>/L1THDiRead' */
  int_T L2InaccurateURead_IWORK[7];    /* '<S9>/L2InaccurateURead' */
  int_T L2InaccurateIRead_IWORK[7];    /* '<S9>/L2InaccurateIRead' */
  int_T L2VoltageRead_IWORK[7];        /* '<S9>/L2VoltageRead' */
  int_T L2CurrentRead_IWORK[7];        /* '<S9>/L2CurrentRead' */
  int_T L2PowFactorRead_IWORK[7];      /* '<S9>/L2PowFactorRead' */
  int_T L2ActivePowRead_IWORK[7];      /* '<S9>/L2ActivePowRead' */
  int_T L2THDuRead_IWORK[7];           /* '<S9>/L2THDuRead' */
  int_T L2THDiRead_IWORK[7];           /* '<S9>/L2THDiRead' */
  int_T L3InaccurateURead_IWORK[7];    /* '<S9>/L3InaccurateURead' */
  int_T L3InaccurateIRead_IWORK[7];    /* '<S9>/L3InaccurateIRead' */
  int_T L3VoltageRead_IWORK[7];        /* '<S9>/L3VoltageRead' */
  int_T L3CurrentRead_IWORK[7];        /* '<S9>/L3CurrentRead' */
  int_T L3PowFactorRead_IWORK[7];      /* '<S9>/L3PowFactorRead' */
  int_T L3ActivePowRead_IWORK[7];      /* '<S9>/L3ActivePowRead' */
  int_T L3THDuRead_IWORK[7];           /* '<S9>/L3THDuRead' */
  int_T L3THDiRead_IWORK[7];           /* '<S9>/L3THDiRead' */
  int_T FrequencyRead_IWORK[7];        /* '<S9>/FrequencyRead' */
  int_T totalPowFactorRead_IWORK[7];   /* '<S9>/totalPowFactorRead' */
  int_T totalActivePowRead_IWORK[7];   /* '<S9>/totalActivePowRead' */
  int_T L1L2VoltageRead_IWORK[7];      /* '<S9>/L1L2VoltageRead' */
  int_T L2L3VoltageRead_IWORK[7];      /* '<S9>/L2L3VoltageRead' */
  int_T L3L1VoltageRead_IWORK[7];      /* '<S9>/L3L1VoltageRead' */
  int_T L1InaccurateURead_IWORK_o[7];  /* '<S11>/L1InaccurateURead' */
  int_T L1InaccurateIRead_IWORK_o[7];  /* '<S11>/L1InaccurateIRead' */
  int_T L1VoltageRead_IWORK_e[7];      /* '<S11>/L1VoltageRead' */
  int_T L1CurrentRead_IWORK_f[7];      /* '<S11>/L1CurrentRead' */
  int_T L1PowFactorRead_IWORK_l[7];    /* '<S11>/L1PowFactorRead' */
  int_T L1ActivePowRead_IWORK_c[7];    /* '<S11>/L1ActivePowRead' */
  int_T L1THDuRead_IWORK_i[7];         /* '<S11>/L1THDuRead' */
  int_T L1THDiRead_IWORK_m[7];         /* '<S11>/L1THDiRead' */
  int_T L2InaccurateURead_IWORK_i[7];  /* '<S11>/L2InaccurateURead' */
  int_T L2InaccurateIRead_IWORK_a[7];  /* '<S11>/L2InaccurateIRead' */
  int_T L2VoltageRead_IWORK_h[7];      /* '<S11>/L2VoltageRead' */
  int_T L2CurrentRead_IWORK_e[7];      /* '<S11>/L2CurrentRead' */
  int_T L2PowFactorRead_IWORK_c[7];    /* '<S11>/L2PowFactorRead' */
  int_T L2ActivePowRead_IWORK_a[7];    /* '<S11>/L2ActivePowRead' */
  int_T L2THDuRead_IWORK_p[7];         /* '<S11>/L2THDuRead' */
  int_T L2THDiRead_IWORK_e[7];         /* '<S11>/L2THDiRead' */
  int_T L3InaccurateURead_IWORK_b[7];  /* '<S11>/L3InaccurateURead' */
  int_T L3InaccurateIRead_IWORK_i[7];  /* '<S11>/L3InaccurateIRead' */
  int_T L3VoltageRead_IWORK_m[7];      /* '<S11>/L3VoltageRead' */
  int_T L3CurrentRead_IWORK_a[7];      /* '<S11>/L3CurrentRead' */
  int_T L3PowFactorRead_IWORK_m[7];    /* '<S11>/L3PowFactorRead' */
  int_T L3ActivePowRead_IWORK_a[7];    /* '<S11>/L3ActivePowRead' */
  int_T L3THDuRead_IWORK_h[7];         /* '<S11>/L3THDuRead' */
  int_T L3THDiRead_IWORK_n[7];         /* '<S11>/L3THDiRead' */
  int_T FrequencyRead_IWORK_k[7];      /* '<S11>/FrequencyRead' */
  int_T totalPowFactorRead_IWORK_n[7]; /* '<S11>/totalPowFactorRead' */
  int_T totalActivePowRead_IWORK_j[7]; /* '<S11>/totalActivePowRead' */
  int_T L1L2VoltageRead_IWORK_f[7];    /* '<S11>/L1L2VoltageRead' */
  int_T L2L3VoltageRead_IWORK_k[7];    /* '<S11>/L2L3VoltageRead' */
  int_T L3L1VoltageRead_IWORK_b[7];    /* '<S11>/L3L1VoltageRead' */
  int_T ACS800CtrlWord_IWORK[7];       /* '<S14>/ACS800CtrlWord' */
  int_T ACS800TorqueSetpoint_IWORK[7]; /* '<S14>/ACS800TorqueSetpoint' */
  int_T ACS800SpeedSetpoint_IWORK[7];  /* '<S14>/ACS800SpeedSetpoint' */
  int_T ACS880CtrlWord_IWORK[7];       /* '<S15>/ACS880CtrlWord' */
  int_T ACS880TorqueSetpoint_IWORK[7]; /* '<S15>/ACS880TorqueSetpoint' */
  int_T ACS880SpeedSetpoint_IWORK[7];  /* '<S15>/ACS880SpeedSetpoint' */
  int_T STATE_1_Modes_j;               /* '<S216>/STATE_1' */
  int_T OUTPUT_1_1_Modes;              /* '<S216>/OUTPUT_1_1' */
  int_T STATE_1_Modes_i[21];           /* '<S332>/STATE_1' */
  int_T OUTPUT_1_0_Modes_b;            /* '<S332>/OUTPUT_1_0' */
  int_T TransportDelay_IWORK[24];      /* '<S60>/Transport Delay' */
  int_T TransportDelay_IWORK_c[24];    /* '<S139>/Transport Delay' */
  int_T OUTPUT_1_0_Modes_j;            /* '<S216>/OUTPUT_1_0' */
  int_T TransportDelay_IWORK_o[24];    /* '<S64>/Transport Delay' */
  int_T TransportDelay_IWORK_c1[24];   /* '<S143>/Transport Delay' */
  int32_T STATE_1_MASS_MATRIX_PR;      /* '<S332>/STATE_1' */
  int8_T Integrator_PrevResetState;    /* '<S412>/Integrator' */
  int8_T Filter_PrevResetState;        /* '<S407>/Filter' */
  int8_T Integrator_PrevResetState_g;  /* '<S475>/Integrator' */
  int8_T Filter_PrevResetState_g;      /* '<S470>/Filter' */
  uint8_T STATE_1_ZcSignalDir;         /* '<S541>/STATE_1' */
  uint8_T STATE_1_ZcStateStore;        /* '<S541>/STATE_1' */
  uint8_T OUTPUT_1_0_ZcSignalDir;      /* '<S541>/OUTPUT_1_0' */
  uint8_T OUTPUT_1_0_ZcStateStore;     /* '<S541>/OUTPUT_1_0' */
  uint8_T STATE_1_ZcSignalDir_o;       /* '<S216>/STATE_1' */
  uint8_T STATE_1_ZcStateStore_g;      /* '<S216>/STATE_1' */
  uint8_T OUTPUT_1_1_ZcSignalDir;      /* '<S216>/OUTPUT_1_1' */
  uint8_T OUTPUT_1_1_ZcStateStore;     /* '<S216>/OUTPUT_1_1' */
  uint8_T STATE_1_ZcSignalDir_f;       /* '<S332>/STATE_1' */
  uint8_T STATE_1_ZcStateStore_f;      /* '<S332>/STATE_1' */
  uint8_T OUTPUT_1_0_ZcSignalDir_b;    /* '<S332>/OUTPUT_1_0' */
  uint8_T OUTPUT_1_0_ZcStateStore_o;   /* '<S332>/OUTPUT_1_0' */
  uint8_T OUTPUT_1_0_ZcSignalDir_g;    /* '<S216>/OUTPUT_1_0' */
  uint8_T OUTPUT_1_0_ZcStateStore_o1;  /* '<S216>/OUTPUT_1_0' */
  uint8_T is_active_c3_windEmulatorStep4_;/* '<S4>/FexcRamp' */
  uint8_T is_active_c7_windEmulatorStep4_;/* '<S18>/ABB Fieldbus Control' */
  uint8_T is_active_UpdateStateMachine;/* '<S18>/ABB Fieldbus Control' */
  uint8_T is_active_UpdateControlWord; /* '<S18>/ABB Fieldbus Control' */
  uint8_T temporalCounter_i1_l;        /* '<S18>/ABB Fieldbus Control' */
  uint8_T is_active_c9_windEmulatorStep4_;/* '<S16>/ABB Fieldbus Control' */
  uint8_T is_active_UpdateStateMachine_a;/* '<S16>/ABB Fieldbus Control' */
  uint8_T is_active_UpdateControlWord_p;/* '<S16>/ABB Fieldbus Control' */
  uint8_T temporalCounter_i1_g;        /* '<S16>/ABB Fieldbus Control' */
  boolean_T Memory_PreviousInput_i;    /* '<S2>/Memory' */
  boolean_T Memory1_PreviousInput;     /* '<S2>/Memory1' */
  boolean_T Memory2_PreviousInput;     /* '<S2>/Memory2' */
  boolean_T Memory_PreviousInput_k;    /* '<S4>/Memory' */
  boolean_T Memory1_PreviousInput_d;   /* '<S4>/Memory1' */
  boolean_T Memory2_PreviousInput_l;   /* '<S4>/Memory2' */
  boolean_T PrevLimited;               /* '<S609>/torqueSlewRate' */
  boolean_T PrevLimited_d;             /* '<S609>/speedSlewRate' */
  boolean_T RTP_1_SetParametersNeeded; /* '<S508>/RTP_1' */
  boolean_T STATE_1_FirstOutput;       /* '<S541>/STATE_1' */
  boolean_T OUTPUT_1_0_FirstOutput;    /* '<S541>/OUTPUT_1_0' */
  boolean_T PrevLimited_o;             /* '<S2>/acs880RateLim' */
  boolean_T Memory_PreviousInput_d;    /* '<S1>/Memory' */
  boolean_T Memory1_PreviousInput_p;   /* '<S1>/Memory1' */
  boolean_T Memory2_PreviousInput_h;   /* '<S1>/Memory2' */
  boolean_T STATE_1_FirstOutput_i;     /* '<S216>/STATE_1' */
  boolean_T OUTPUT_1_1_FirstOutput;    /* '<S216>/OUTPUT_1_1' */
  boolean_T STATE_1_FirstOutput_a;     /* '<S332>/STATE_1' */
  boolean_T OUTPUT_1_0_FirstOutput_f;  /* '<S332>/OUTPUT_1_0' */
  boolean_T OUTPUT_1_0_FirstOutput_e;  /* '<S216>/OUTPUT_1_0' */
  boolean_T Memory_PreviousInput_kk;   /* '<S433>/Memory' */
  boolean_T Memory_PreviousInput_h;    /* '<S434>/Memory' */
  boolean_T PrevLimited_g;             /* '<S373>/Rate Limiter' */
  boolean_T Memory_PreviousInput_g;    /* '<S496>/Memory' */
  boolean_T Memory_PreviousInput_n;    /* '<S497>/Memory' */
  boolean_T PrevLimited_dz;            /* '<S436>/Rate Limiter' */
  DW_YawKinematicTransforms_win_T sf_YawKinematicTransforms_l;/* '<S212>/Yaw Kinematic Transforms' */
  DW_YawForceTransforms_windEmu_T sf_YawForceTransforms_i;/* '<S149>/Yaw Force Transforms' */
  DW_MATLABFunction1_windEmulat_T sf_MATLABFunction1_e;/* '<S205>/MATLAB Function1' */
  DW_quaternion2EulXYZ_windEmul_T sf_quaternion2EulXYZ_c;/* '<S160>/quaternion2EulXYZ' */
  DW_YawKinematicTransforms_win_T sf_YawKinematicTransforms;/* '<S133>/Yaw Kinematic Transforms' */
  DW_YawForceTransforms_windEmu_T sf_YawForceTransforms;/* '<S70>/Yaw Force Transforms' */
  DW_MATLABFunction1_windEmulat_T sf_MATLABFunction1;/* '<S126>/MATLAB Function1' */
  DW_quaternion2EulXYZ_windEmul_T sf_quaternion2EulXYZ;/* '<S81>/quaternion2EulXYZ' */
  DW_MovingAverage_windEmulat_f_T MovingAverage1;/* '<S51>/Moving Average' */
  DW_MovingAverage_windEmulat_f_T MovingAverage_pn;/* '<S51>/Moving Average' */
  DW_ParseStatusWord_windEmulat_T sf_ParseStatusWord_h;/* '<S49>/Parse Status Word' */
  DW_MovingAverage_windEmulator_T MovingAverage_p;/* '<S43>/Moving Average' */
  DW_parseCtrlWord_windEmulator_T sf_parseCtrlWord_h;/* '<S46>/parseCtrlWord' */
  DW_ParseStatusWord_windEmulat_T sf_ParseStatusWord;/* '<S44>/Parse Status Word' */
  DW_MovingAverage_windEmulator_T MovingAverage;/* '<S43>/Moving Average' */
  DW_parseCtrlWord_windEmulator_T sf_parseCtrlWord;/* '<S41>/parseCtrlWord' */
};

/* Continuous states (default storage) */
struct X_windEmulatorStep4_WECSim_T {
  real_T Integrator_CSTATE;            /* '<S590>/Integrator' */
  real_T Internal_CSTATE[3];           /* '<S531>/Internal' */
  real_T Internal_CSTATE_j;            /* '<S545>/Internal' */
  real_T Internal_CSTATE_a;            /* '<S542>/Internal' */
  real_T windEmulatorStep4_WECSimhptoSim[2];/* '<S216>/STATE_1' */
  real_T windEmulatorStep4_WECSimhptoS_n;/* '<S332>/INPUT_3_1_1' */
  real_T windEmulatorStep4_WECSimhptoS_h[35];/* '<S332>/STATE_1' */
};

/* State derivatives (default storage) */
struct XDot_windEmulatorStep4_WECSim_T {
  real_T Integrator_CSTATE;            /* '<S590>/Integrator' */
  real_T Internal_CSTATE[3];           /* '<S531>/Internal' */
  real_T Internal_CSTATE_j;            /* '<S545>/Internal' */
  real_T Internal_CSTATE_a;            /* '<S542>/Internal' */
  real_T windEmulatorStep4_WECSimhptoSim[2];/* '<S216>/STATE_1' */
  real_T windEmulatorStep4_WECSimhptoS_n;/* '<S332>/INPUT_3_1_1' */
  real_T windEmulatorStep4_WECSimhptoS_h[35];/* '<S332>/STATE_1' */
};

/* State disabled  */
struct XDis_windEmulatorStep4_WECSim_T {
  boolean_T Integrator_CSTATE;         /* '<S590>/Integrator' */
  boolean_T Internal_CSTATE[3];        /* '<S531>/Internal' */
  boolean_T Internal_CSTATE_j;         /* '<S545>/Internal' */
  boolean_T Internal_CSTATE_a;         /* '<S542>/Internal' */
  boolean_T windEmulatorStep4_WECSimhptoSim[2];/* '<S216>/STATE_1' */
  boolean_T windEmulatorStep4_WECSimhptoS_n;/* '<S332>/INPUT_3_1_1' */
  boolean_T windEmulatorStep4_WECSimhptoS_h[35];/* '<S332>/STATE_1' */
};

/* Zero-crossing (trigger) state */
struct PrevZCX_windEmulatorStep4_WECSim_T {
  ZCSigState Integrator_Reset_ZCE;     /* '<S590>/Integrator' */
};

/* Mass Matrix (global) */
struct MassMatrix_windEmulatorStep4_WECSim_T {
  int_T ir[26];
  int_T jc[45];
  real_T pr[26];
};

#ifndef ODE14X_INTG
#define ODE14X_INTG

/* ODE14X Integration Data */
struct ODE14X_IntgData {
  real_T *x0;
  real_T *f0;
  real_T *x1start;
  real_T *f1;
  real_T *Delta;
  real_T *E;
  real_T *fac;
  real_T *DFDX;
  real_T *W;
  int_T *pivots;
  real_T *xtmp;
  real_T *ztmp;
  real_T *M;
  real_T *M1;
  real_T *Edot;
  real_T *xdot;
  real_T *fminusMxdot;
  boolean_T isFirstStep;
};

#endif

/* External inputs (root inport signals with default storage) */
struct ExtU_windEmulatorStep4_WECSim_T {
  real_T inportSpeed_rpm;              /* '<Root>/inportSpeed_rpm' */
  real_T inportTorque_Nm;              /* '<Root>/inportTorque_Nm' */
  real_T inportCaseCounter;            /* '<Root>/inportCaseCounter' */
};

/* Real-time Model Data Structure */
struct tag_RTM_windEmulatorStep4_WECSim_T {
  struct SimStruct_tag * *childSfunctions;
  const char_T *errorStatus;
  SS_SimMode simMode;
  RTWSolverInfo solverInfo;
  RTWSolverInfo *solverInfoPtr;
  void *sfcnInfo;

  /*
   * NonInlinedSFcns:
   * The following substructure contains information regarding
   * non-inlined s-functions used in the model.
   */
  struct {
    RTWSfcnInfo sfcnInfo;
    time_T *taskTimePtrs[2];
    SimStruct childSFunctions[1];
    SimStruct *childSFunctionPtrs[1];
    struct _ssBlkInfo2 blkInfo2[1];
    struct _ssSFcnModelMethods2 methods2[1];
    struct _ssSFcnModelMethods3 methods3[1];
    struct _ssSFcnModelMethods4 methods4[1];
    struct _ssStatesInfo2 statesInfo2[1];
    ssPeriodicStatesInfo periodicStatesInfo[1];
    struct _ssPortInfo2 inputOutputPortInfo2[1];
    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      struct _ssInPortUnit inputPortUnits[1];
      struct _ssInPortCoSimAttribute inputPortCoSimAttribute[1];
    } Sfcn0;
  } NonInlinedSFcns;

  X_windEmulatorStep4_WECSim_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  XDis_windEmulatorStep4_WECSim_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  int_T massMatrixType;
  int_T massMatrixNzMax;
  int_T *massMatrixIr;
  int_T *massMatrixJc;
  real_T *massMatrixPr;
  real_T odeX0[44];
  real_T odeF0[44];
  real_T odeX1START[44];
  real_T odeF1[44];
  real_T odeDELTA[44];
  real_T odeE[4*44];
  real_T odeFAC[44];
  real_T odeDFDX[44*44];
  real_T odeW[44*44];
  int_T odePIVOTS[44];
  real_T odeXTMP[44];
  real_T odeZTMP[44];
  real_T odeMASSMATRIX_M[26];
  real_T odeMASSMATRIX_M1[26];
  real_T odeEDOT[4*44];
  real_T odeXDOT[44];
  real_T odeFMXDOT[44];
  ODE14X_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    uint32_T options;
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numU;
    int_T numY;
    int_T numSampTimes;
    int_T numBlocks;
    int_T numBlockIO;
    int_T numBlockPrms;
    int_T numDwork;
    int_T numSFcnPrms;
    int_T numSFcns;
    int_T numIports;
    int_T numOports;
    int_T numNonSampZCs;
    int_T sysDirFeedThru;
    int_T rtwGenSfcn;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    time_T stepSize;
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    time_T stepSize1;
    time_T tStart;
    time_T tFinal;
    time_T timeOfLastOutput;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *sampleTimes;
    time_T *offsetTimes;
    int_T *sampleTimeTaskIDPtr;
    int_T *sampleHits;
    int_T *perTaskSampleHits;
    time_T *t;
    time_T sampleTimesArray[2];
    time_T offsetTimesArray[2];
    int_T sampleTimeTaskIDArray[2];
    int_T sampleHitArray[2];
    int_T perTaskSampleHitsArray[4];
    time_T tArray[2];
  } Timing;
};

/* Block signals (default storage) */
#ifdef __cplusplus

extern "C"
{

#endif

  extern struct B_windEmulatorStep4_WECSim_T windEmulatorStep4_WECSim_B;

#ifdef __cplusplus

}

#endif

/* Continuous states (default storage) */
extern X_windEmulatorStep4_WECSim_T windEmulatorStep4_WECSim_X;

/* Disabled states (default storage) */
extern XDis_windEmulatorStep4_WECSim_T windEmulatorStep4_WECSim_XDis;

/* Block states (default storage) */
extern struct DW_windEmulatorStep4_WECSim_T windEmulatorStep4_WECSim_DW;

/* Zero-crossing (trigger) state */
extern PrevZCX_windEmulatorStep4_WECSim_T windEmulatorStep4_WECSim_PrevZCX;

/* global MassMatrix */
extern MassMatrix_windEmulatorStep4_WECSim_T windEmulatorStep4_WECSim_MassMatrix;

#ifdef __cplusplus

extern "C"
{

#endif

  /* External inputs (root inport signals with default storage) */
  extern struct ExtU_windEmulatorStep4_WECSim_T windEmulatorStep4_WECSim_U;

#ifdef __cplusplus

}

#endif

#ifdef __cplusplus

extern "C"
{

#endif

  /* Model entry point functions */
  extern void windEmulatorStep4_WECSim_initialize(void);
  extern void windEmulatorStep4_WECSim_step(void);
  extern void windEmulatorStep4_WECSim_terminate(void);

#ifdef __cplusplus

}

#endif

/* Real-time Model object */
#ifdef __cplusplus

extern "C"
{

#endif

  extern RT_MODEL_windEmulatorStep4_WECSim_T *const windEmulatorStep4_WECSim_M;

#ifdef __cplusplus

}

#endif

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'windEmulatorStep4_WECSim'
 * '<S1>'   : 'windEmulatorStep4_WECSim/acs800Ctrl'
 * '<S2>'   : 'windEmulatorStep4_WECSim/acs880Ctrl'
 * '<S3>'   : 'windEmulatorStep4_WECSim/ctrlSignalSelector'
 * '<S4>'   : 'windEmulatorStep4_WECSim/expCtrl'
 * '<S5>'   : 'windEmulatorStep4_WECSim/fileAndUI'
 * '<S6>'   : 'windEmulatorStep4_WECSim/hptoSim'
 * '<S7>'   : 'windEmulatorStep4_WECSim/processHptoInputs'
 * '<S8>'   : 'windEmulatorStep4_WECSim/readAcs800Pdos'
 * '<S9>'   : 'windEmulatorStep4_WECSim/readAcs800Power'
 * '<S10>'  : 'windEmulatorStep4_WECSim/readAcs880Pdos'
 * '<S11>'  : 'windEmulatorStep4_WECSim/readAcs880Power'
 * '<S12>'  : 'windEmulatorStep4_WECSim/readShaftSignals'
 * '<S13>'  : 'windEmulatorStep4_WECSim/sidCtrl'
 * '<S14>'  : 'windEmulatorStep4_WECSim/writeAcs800Pdos'
 * '<S15>'  : 'windEmulatorStep4_WECSim/writeAcs880Pdos'
 * '<S16>'  : 'windEmulatorStep4_WECSim/acs800Ctrl/ACS880FieldbusControl'
 * '<S17>'  : 'windEmulatorStep4_WECSim/acs800Ctrl/ACS880FieldbusControl/ABB Fieldbus Control'
 * '<S18>'  : 'windEmulatorStep4_WECSim/acs880Ctrl/ACS880FieldbusControl'
 * '<S19>'  : 'windEmulatorStep4_WECSim/acs880Ctrl/ACS880FieldbusControl/ABB Fieldbus Control'
 * '<S20>'  : 'windEmulatorStep4_WECSim/expCtrl/FexcRamp'
 * '<S21>'  : 'windEmulatorStep4_WECSim/fileAndUI/acs800CtrlSignalsToFile'
 * '<S22>'  : 'windEmulatorStep4_WECSim/fileAndUI/acs800CtrlSignalsUI'
 * '<S23>'  : 'windEmulatorStep4_WECSim/fileAndUI/acs800SignalsToFile'
 * '<S24>'  : 'windEmulatorStep4_WECSim/fileAndUI/acs800SignalsUI'
 * '<S25>'  : 'windEmulatorStep4_WECSim/fileAndUI/acs880CtrlSignalsToFile'
 * '<S26>'  : 'windEmulatorStep4_WECSim/fileAndUI/acs880CtrlSignalsUI'
 * '<S27>'  : 'windEmulatorStep4_WECSim/fileAndUI/acs880SignalsToFile'
 * '<S28>'  : 'windEmulatorStep4_WECSim/fileAndUI/acs880SignalsUI'
 * '<S29>'  : 'windEmulatorStep4_WECSim/fileAndUI/diagnostics'
 * '<S30>'  : 'windEmulatorStep4_WECSim/fileAndUI/expCtrlSignalsToFile'
 * '<S31>'  : 'windEmulatorStep4_WECSim/fileAndUI/expCtrlSignalsUI'
 * '<S32>'  : 'windEmulatorStep4_WECSim/fileAndUI/hptoCtrlSignals'
 * '<S33>'  : 'windEmulatorStep4_WECSim/fileAndUI/hptoCtrlSignalsUI'
 * '<S34>'  : 'windEmulatorStep4_WECSim/fileAndUI/hptoSignalsToFile'
 * '<S35>'  : 'windEmulatorStep4_WECSim/fileAndUI/hptoSignalsUI'
 * '<S36>'  : 'windEmulatorStep4_WECSim/fileAndUI/invPowerAcs800ToFile'
 * '<S37>'  : 'windEmulatorStep4_WECSim/fileAndUI/invPowerAcs880ToFile'
 * '<S38>'  : 'windEmulatorStep4_WECSim/fileAndUI/shaftSignalsToFile'
 * '<S39>'  : 'windEmulatorStep4_WECSim/fileAndUI/shaftSignalsUI'
 * '<S40>'  : 'windEmulatorStep4_WECSim/fileAndUI/sidInfoSignalsToFile'
 * '<S41>'  : 'windEmulatorStep4_WECSim/fileAndUI/acs800CtrlSignalsUI/ctrlWordDetail'
 * '<S42>'  : 'windEmulatorStep4_WECSim/fileAndUI/acs800CtrlSignalsUI/ctrlWordDetail/parseCtrlWord'
 * '<S43>'  : 'windEmulatorStep4_WECSim/fileAndUI/acs800SignalsUI/powerCals'
 * '<S44>'  : 'windEmulatorStep4_WECSim/fileAndUI/acs800SignalsUI/statusWordDetail'
 * '<S45>'  : 'windEmulatorStep4_WECSim/fileAndUI/acs800SignalsUI/statusWordDetail/Parse Status Word'
 * '<S46>'  : 'windEmulatorStep4_WECSim/fileAndUI/acs880CtrlSignalsUI/ctrlWordDetail'
 * '<S47>'  : 'windEmulatorStep4_WECSim/fileAndUI/acs880CtrlSignalsUI/ctrlWordDetail/parseCtrlWord'
 * '<S48>'  : 'windEmulatorStep4_WECSim/fileAndUI/acs880SignalsUI/powerCals'
 * '<S49>'  : 'windEmulatorStep4_WECSim/fileAndUI/acs880SignalsUI/statusWordDetail'
 * '<S50>'  : 'windEmulatorStep4_WECSim/fileAndUI/acs880SignalsUI/statusWordDetail/Parse Status Word'
 * '<S51>'  : 'windEmulatorStep4_WECSim/fileAndUI/hptoSignalsUI/powerCalcs'
 * '<S52>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel'
 * '<S53>'  : 'windEmulatorStep4_WECSim/hptoSim/hptoModel'
 * '<S54>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base'
 * '<S55>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Fixed'
 * '<S56>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap'
 * '<S57>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Global Reference Frame'
 * '<S58>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim'
 * '<S59>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Rotational PTO Actuation Torque'
 * '<S60>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body'
 * '<S61>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Additional Linear Damping Force Calculation'
 * '<S62>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/B2B Subsystem'
 * '<S63>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Hydrostatic Restoring Force Calculation'
 * '<S64>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Morison Element and Viscous Damping Force Calculation'
 * '<S65>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Nonlinear Wave Elevation'
 * '<S66>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure'
 * '<S67>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Variable Hydrodynamics Control'
 * '<S68>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation'
 * '<S69>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Wave Radiation Forces Calculation'
 * '<S70>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Yaw Force Transforms'
 * '<S71>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Yaw Kinematic Transforms'
 * '<S72>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/B2B Subsystem/noB2B'
 * '<S73>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Hydrostatic Restoring Force Calculation/Linear and Nonlinear Restoring Force Variant Subsystem'
 * '<S74>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Hydrostatic Restoring Force Calculation/Linear and Nonlinear Restoring Force Variant Subsystem/Linear Hydrostatic Restoring Force'
 * '<S75>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Hydrostatic Restoring Force Calculation/Linear and Nonlinear Restoring Force Variant Subsystem/Linear Hydrostatic Restoring Force/Net Buoyancy Force'
 * '<S76>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Hydrostatic Restoring Force Calculation/Linear and Nonlinear Restoring Force Variant Subsystem/Linear Hydrostatic Restoring Force/Net Buoyancy Force/Cross Product'
 * '<S77>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Morison Element and Viscous Damping Force Calculation/Morison Element Variant Subsystem'
 * '<S78>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Morison Element and Viscous Damping Force Calculation/Morison Element Variant Subsystem/Morison Element Off'
 * '<S79>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Nonlinear Wave Elevation/Linear or Instantaneous Free Surface'
 * '<S80>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Nonlinear Wave Elevation/Linear or Instantaneous Free Surface/Mean Water Free Surface'
 * '<S81>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor'
 * '<S82>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Simulink-PS Converter1'
 * '<S83>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Simulink-PS Converter2'
 * '<S84>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter'
 * '<S85>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter1'
 * '<S86>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter10'
 * '<S87>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter11'
 * '<S88>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter12'
 * '<S89>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter13'
 * '<S90>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter14'
 * '<S91>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter15'
 * '<S92>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter2'
 * '<S93>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter3'
 * '<S94>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter4'
 * '<S95>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter5'
 * '<S96>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter6'
 * '<S97>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter7'
 * '<S98>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter8'
 * '<S99>'  : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter9'
 * '<S100>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/quaternion2EulXYZ'
 * '<S101>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter/EVAL_KEY'
 * '<S102>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter1/EVAL_KEY'
 * '<S103>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter10/EVAL_KEY'
 * '<S104>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter11/EVAL_KEY'
 * '<S105>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter12/EVAL_KEY'
 * '<S106>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter13/EVAL_KEY'
 * '<S107>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter14/EVAL_KEY'
 * '<S108>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter15/EVAL_KEY'
 * '<S109>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter2/EVAL_KEY'
 * '<S110>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter3/EVAL_KEY'
 * '<S111>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter4/EVAL_KEY'
 * '<S112>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter5/EVAL_KEY'
 * '<S113>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter6/EVAL_KEY'
 * '<S114>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter7/EVAL_KEY'
 * '<S115>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter8/EVAL_KEY'
 * '<S116>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter9/EVAL_KEY'
 * '<S117>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Simulink-PS Converter1/EVAL_KEY'
 * '<S118>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Structure/Simulink-PS Converter2/EVAL_KEY'
 * '<S119>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Variable Hydrodynamics Control/Variable Hydrodynamics Subsystem'
 * '<S120>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Variable Hydrodynamics Control/Variable Hydrodynamics Subsystem/noVariableHydro'
 * '<S121>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation/FexcInputs'
 * '<S122>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation/Linear Wave Excitation Force Variant Subsystem'
 * '<S123>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation/Nonlinear FK Force Variant Subsystem'
 * '<S124>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation/Ramp function'
 * '<S125>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation/Second Order Excitation Force Variant Subsystem'
 * '<S126>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation/Linear Wave Excitation Force Variant Subsystem/Regular Wave  Excitation Force'
 * '<S127>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation/Linear Wave Excitation Force Variant Subsystem/Regular Wave  Excitation Force/MATLAB Function1'
 * '<S128>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation/Nonlinear FK Force Variant Subsystem/No Nonlinear FroudeKrylov Force'
 * '<S129>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation/Second Order Excitation Force Variant Subsystem/No Second Order Excitation Force'
 * '<S130>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Wave Radiation Forces Calculation/SS CI and Constant-Damping-CoeVariant Subsystem'
 * '<S131>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Wave Radiation Forces Calculation/SS CI and Constant-Damping-CoeVariant Subsystem/Constant Coefficients'
 * '<S132>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Yaw Force Transforms/Yaw Force Transforms'
 * '<S133>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Yaw Kinematic Transforms/No B2B'
 * '<S134>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Base/Hydrodynamic Body/Yaw Kinematic Transforms/No B2B/Yaw Kinematic Transforms'
 * '<S135>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Fixed/PS-Simulink Converter'
 * '<S136>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Fixed/PS-Simulink Converter1'
 * '<S137>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Fixed/PS-Simulink Converter/EVAL_KEY'
 * '<S138>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Fixed/PS-Simulink Converter1/EVAL_KEY'
 * '<S139>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body'
 * '<S140>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Additional Linear Damping Force Calculation'
 * '<S141>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/B2B Subsystem'
 * '<S142>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Hydrostatic Restoring Force Calculation'
 * '<S143>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Morison Element and Viscous Damping Force Calculation'
 * '<S144>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Nonlinear Wave Elevation'
 * '<S145>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure'
 * '<S146>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Variable Hydrodynamics Control'
 * '<S147>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation'
 * '<S148>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Wave Radiation Forces Calculation'
 * '<S149>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Yaw Force Transforms'
 * '<S150>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Yaw Kinematic Transforms'
 * '<S151>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/B2B Subsystem/noB2B'
 * '<S152>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Hydrostatic Restoring Force Calculation/Linear and Nonlinear Restoring Force Variant Subsystem'
 * '<S153>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Hydrostatic Restoring Force Calculation/Linear and Nonlinear Restoring Force Variant Subsystem/Linear Hydrostatic Restoring Force'
 * '<S154>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Hydrostatic Restoring Force Calculation/Linear and Nonlinear Restoring Force Variant Subsystem/Linear Hydrostatic Restoring Force/Net Buoyancy Force'
 * '<S155>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Hydrostatic Restoring Force Calculation/Linear and Nonlinear Restoring Force Variant Subsystem/Linear Hydrostatic Restoring Force/Net Buoyancy Force/Cross Product'
 * '<S156>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Morison Element and Viscous Damping Force Calculation/Morison Element Variant Subsystem'
 * '<S157>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Morison Element and Viscous Damping Force Calculation/Morison Element Variant Subsystem/Morison Element Off'
 * '<S158>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Nonlinear Wave Elevation/Linear or Instantaneous Free Surface'
 * '<S159>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Nonlinear Wave Elevation/Linear or Instantaneous Free Surface/Mean Water Free Surface'
 * '<S160>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor'
 * '<S161>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Simulink-PS Converter1'
 * '<S162>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Simulink-PS Converter2'
 * '<S163>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter'
 * '<S164>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter1'
 * '<S165>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter10'
 * '<S166>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter11'
 * '<S167>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter12'
 * '<S168>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter13'
 * '<S169>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter14'
 * '<S170>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter15'
 * '<S171>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter2'
 * '<S172>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter3'
 * '<S173>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter4'
 * '<S174>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter5'
 * '<S175>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter6'
 * '<S176>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter7'
 * '<S177>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter8'
 * '<S178>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter9'
 * '<S179>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/quaternion2EulXYZ'
 * '<S180>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter/EVAL_KEY'
 * '<S181>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter1/EVAL_KEY'
 * '<S182>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter10/EVAL_KEY'
 * '<S183>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter11/EVAL_KEY'
 * '<S184>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter12/EVAL_KEY'
 * '<S185>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter13/EVAL_KEY'
 * '<S186>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter14/EVAL_KEY'
 * '<S187>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter15/EVAL_KEY'
 * '<S188>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter2/EVAL_KEY'
 * '<S189>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter3/EVAL_KEY'
 * '<S190>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter4/EVAL_KEY'
 * '<S191>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter5/EVAL_KEY'
 * '<S192>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter6/EVAL_KEY'
 * '<S193>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter7/EVAL_KEY'
 * '<S194>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter8/EVAL_KEY'
 * '<S195>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Motion Sensor/PS-Simulink Converter9/EVAL_KEY'
 * '<S196>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Simulink-PS Converter1/EVAL_KEY'
 * '<S197>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Structure/Simulink-PS Converter2/EVAL_KEY'
 * '<S198>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Variable Hydrodynamics Control/Variable Hydrodynamics Subsystem'
 * '<S199>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Variable Hydrodynamics Control/Variable Hydrodynamics Subsystem/noVariableHydro'
 * '<S200>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation/FexcInputs'
 * '<S201>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation/Linear Wave Excitation Force Variant Subsystem'
 * '<S202>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation/Nonlinear FK Force Variant Subsystem'
 * '<S203>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation/Ramp function'
 * '<S204>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation/Second Order Excitation Force Variant Subsystem'
 * '<S205>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation/Linear Wave Excitation Force Variant Subsystem/Regular Wave  Excitation Force'
 * '<S206>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation/Linear Wave Excitation Force Variant Subsystem/Regular Wave  Excitation Force/MATLAB Function1'
 * '<S207>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation/Nonlinear FK Force Variant Subsystem/No Nonlinear FroudeKrylov Force'
 * '<S208>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation/Second Order Excitation Force Variant Subsystem/No Second Order Excitation Force'
 * '<S209>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Wave Radiation Forces Calculation/SS CI and Constant-Damping-CoeVariant Subsystem'
 * '<S210>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Wave Radiation Forces Calculation/SS CI and Constant-Damping-CoeVariant Subsystem/Constant Coefficients'
 * '<S211>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Yaw Force Transforms/Yaw Force Transforms'
 * '<S212>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Yaw Kinematic Transforms/No B2B'
 * '<S213>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Flap/Hydrodynamic Body/Yaw Kinematic Transforms/No B2B/Yaw Kinematic Transforms'
 * '<S214>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Global Reference Frame/Solver Configuration'
 * '<S215>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Global Reference Frame/waveVis'
 * '<S216>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Global Reference Frame/Solver Configuration/EVAL_KEY'
 * '<S217>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Global Reference Frame/waveVis/waveVisOff'
 * '<S218>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller'
 * '<S219>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor1'
 * '<S220>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor10'
 * '<S221>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor13'
 * '<S222>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor16'
 * '<S223>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor2'
 * '<S224>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor3'
 * '<S225>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor4'
 * '<S226>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor5'
 * '<S227>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor6'
 * '<S228>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor7'
 * '<S229>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor8'
 * '<S230>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor9'
 * '<S231>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/PS-Simulink Converter'
 * '<S232>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/PS-Simulink Converter1'
 * '<S233>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/PS-Simulink Converter2'
 * '<S234>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/PS-Simulink Converter3'
 * '<S235>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/PS-Simulink Converter4'
 * '<S236>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Simulink-PS Converter'
 * '<S237>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Simulink-PS Converter2'
 * '<S238>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Simulink-PS Converter3'
 * '<S239>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Solver Configuration'
 * '<S240>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/pistonPosition'
 * '<S241>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/pressureSensorA'
 * '<S242>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/pressureSensorB'
 * '<S243>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/pressureSensorC'
 * '<S244>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/pressureSensorC1'
 * '<S245>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/pressureSensorD'
 * '<S246>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Anti-windup'
 * '<S247>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/D Gain'
 * '<S248>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/External Derivative'
 * '<S249>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Filter'
 * '<S250>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Filter ICs'
 * '<S251>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/I Gain'
 * '<S252>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Ideal P Gain'
 * '<S253>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Ideal P Gain Fdbk'
 * '<S254>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Integrator'
 * '<S255>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Integrator ICs'
 * '<S256>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/N Copy'
 * '<S257>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/N Gain'
 * '<S258>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/P Copy'
 * '<S259>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Parallel P Gain'
 * '<S260>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Reset Signal'
 * '<S261>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Saturation'
 * '<S262>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Saturation Fdbk'
 * '<S263>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Sum'
 * '<S264>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Sum Fdbk'
 * '<S265>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Tracking Mode'
 * '<S266>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Tracking Mode Sum'
 * '<S267>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Tsamp - Integral'
 * '<S268>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Tsamp - Ngain'
 * '<S269>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/postSat Signal'
 * '<S270>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/preInt Signal'
 * '<S271>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/preSat Signal'
 * '<S272>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Anti-windup/Passthrough'
 * '<S273>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/D Gain/Internal Parameters'
 * '<S274>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/External Derivative/Error'
 * '<S275>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Filter/Differentiator'
 * '<S276>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Filter/Differentiator/Tsamp'
 * '<S277>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Filter/Differentiator/Tsamp/Internal Ts'
 * '<S278>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Filter ICs/Internal IC - Differentiator'
 * '<S279>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/I Gain/Internal Parameters'
 * '<S280>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Ideal P Gain/Passthrough'
 * '<S281>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Ideal P Gain Fdbk/Disabled'
 * '<S282>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Integrator/Discrete'
 * '<S283>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Integrator ICs/Internal IC'
 * '<S284>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/N Copy/Disabled wSignal Specification'
 * '<S285>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/N Gain/Passthrough'
 * '<S286>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/P Copy/Disabled'
 * '<S287>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Parallel P Gain/Internal Parameters'
 * '<S288>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Reset Signal/Disabled'
 * '<S289>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Saturation/Passthrough'
 * '<S290>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Saturation Fdbk/Disabled'
 * '<S291>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Sum/Sum_PID'
 * '<S292>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Sum Fdbk/Disabled'
 * '<S293>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Tracking Mode/Disabled'
 * '<S294>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Tracking Mode Sum/Passthrough'
 * '<S295>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Tsamp - Integral/TsSignalSpecification'
 * '<S296>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/Tsamp - Ngain/Passthrough'
 * '<S297>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/postSat Signal/Forward_Path'
 * '<S298>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/preInt Signal/Internal PreInt'
 * '<S299>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Discrete PID Controller/preSat Signal/Forward_Path'
 * '<S300>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor1/PS-Simulink Converter2'
 * '<S301>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor1/PS-Simulink Converter2/EVAL_KEY'
 * '<S302>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor10/PS-Simulink Converter2'
 * '<S303>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor10/PS-Simulink Converter2/EVAL_KEY'
 * '<S304>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor13/PS-Simulink Converter2'
 * '<S305>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor13/PS-Simulink Converter2/EVAL_KEY'
 * '<S306>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor16/PS-Simulink Converter2'
 * '<S307>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor16/PS-Simulink Converter2/EVAL_KEY'
 * '<S308>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor2/PS-Simulink Converter2'
 * '<S309>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor2/PS-Simulink Converter2/EVAL_KEY'
 * '<S310>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor3/PS-Simulink Converter2'
 * '<S311>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor3/PS-Simulink Converter2/EVAL_KEY'
 * '<S312>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor4/PS-Simulink Converter2'
 * '<S313>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor4/PS-Simulink Converter2/EVAL_KEY'
 * '<S314>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor5/PS-Simulink Converter2'
 * '<S315>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor5/PS-Simulink Converter2/EVAL_KEY'
 * '<S316>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor6/PS-Simulink Converter2'
 * '<S317>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor6/PS-Simulink Converter2/EVAL_KEY'
 * '<S318>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor7/PS-Simulink Converter2'
 * '<S319>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor7/PS-Simulink Converter2/EVAL_KEY'
 * '<S320>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor8/PS-Simulink Converter2'
 * '<S321>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor8/PS-Simulink Converter2/EVAL_KEY'
 * '<S322>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor9/PS-Simulink Converter2'
 * '<S323>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/FlowSensor9/PS-Simulink Converter2/EVAL_KEY'
 * '<S324>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/PS-Simulink Converter/EVAL_KEY'
 * '<S325>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/PS-Simulink Converter1/EVAL_KEY'
 * '<S326>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/PS-Simulink Converter2/EVAL_KEY'
 * '<S327>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/PS-Simulink Converter3/EVAL_KEY'
 * '<S328>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/PS-Simulink Converter4/EVAL_KEY'
 * '<S329>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Simulink-PS Converter/EVAL_KEY'
 * '<S330>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Simulink-PS Converter2/EVAL_KEY'
 * '<S331>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Simulink-PS Converter3/EVAL_KEY'
 * '<S332>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Solver Configuration/EVAL_KEY'
 * '<S333>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/pistonPosition/PS-Simulink Converter'
 * '<S334>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/pistonPosition/PS-Simulink Converter2'
 * '<S335>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/pistonPosition/PS-Simulink Converter/EVAL_KEY'
 * '<S336>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/pistonPosition/PS-Simulink Converter2/EVAL_KEY'
 * '<S337>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/pressureSensorA/PS-Simulink Converter2'
 * '<S338>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/pressureSensorA/PS-Simulink Converter2/EVAL_KEY'
 * '<S339>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/pressureSensorB/PS-Simulink Converter2'
 * '<S340>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/pressureSensorB/PS-Simulink Converter2/EVAL_KEY'
 * '<S341>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/pressureSensorC/PS-Simulink Converter2'
 * '<S342>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/pressureSensorC/PS-Simulink Converter2/EVAL_KEY'
 * '<S343>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/pressureSensorC1/PS-Simulink Converter2'
 * '<S344>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/pressureSensorC1/PS-Simulink Converter2/EVAL_KEY'
 * '<S345>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/pressureSensorD/PS-Simulink Converter2'
 * '<S346>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/pressureSensorD/PS-Simulink Converter2/EVAL_KEY'
 * '<S347>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Rotational PTO Actuation Torque/PS-Simulink Converter'
 * '<S348>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Rotational PTO Actuation Torque/PS-Simulink Converter1'
 * '<S349>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Rotational PTO Actuation Torque/PS-Simulink Converter2'
 * '<S350>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Rotational PTO Actuation Torque/PS-Simulink Converter3'
 * '<S351>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Rotational PTO Actuation Torque/PS-Simulink Converter4'
 * '<S352>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Rotational PTO Actuation Torque/PS-Simulink Converter5'
 * '<S353>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Rotational PTO Actuation Torque/PS-Simulink Converter6'
 * '<S354>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Rotational PTO Actuation Torque/PS-Simulink Converter7'
 * '<S355>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Rotational PTO Actuation Torque/Simulink-PS Converter1'
 * '<S356>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Rotational PTO Actuation Torque/PS-Simulink Converter/EVAL_KEY'
 * '<S357>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Rotational PTO Actuation Torque/PS-Simulink Converter1/EVAL_KEY'
 * '<S358>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Rotational PTO Actuation Torque/PS-Simulink Converter2/EVAL_KEY'
 * '<S359>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Rotational PTO Actuation Torque/PS-Simulink Converter3/EVAL_KEY'
 * '<S360>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Rotational PTO Actuation Torque/PS-Simulink Converter4/EVAL_KEY'
 * '<S361>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Rotational PTO Actuation Torque/PS-Simulink Converter5/EVAL_KEY'
 * '<S362>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Rotational PTO Actuation Torque/PS-Simulink Converter6/EVAL_KEY'
 * '<S363>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Rotational PTO Actuation Torque/PS-Simulink Converter7/EVAL_KEY'
 * '<S364>' : 'windEmulatorStep4_WECSim/hptoSim/WECSimModel/Rotational PTO Actuation Torque/Simulink-PS Converter1/EVAL_KEY'
 * '<S365>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm'
 * '<S366>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO'
 * '<S367>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/ForceToTorque'
 * '<S368>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin'
 * '<S369>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin'
 * '<S370>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/Ramp'
 * '<S371>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/ForceToTorque/DesiredForce'
 * '<S372>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/PressureControl1'
 * '<S373>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1'
 * '<S374>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/TorqueControl1'
 * '<S375>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem'
 * '<S376>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem1'
 * '<S377>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller'
 * '<S378>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Anti-windup'
 * '<S379>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/D Gain'
 * '<S380>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/External Derivative'
 * '<S381>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Filter'
 * '<S382>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Filter ICs'
 * '<S383>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/I Gain'
 * '<S384>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Ideal P Gain'
 * '<S385>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Ideal P Gain Fdbk'
 * '<S386>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Integrator'
 * '<S387>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Integrator ICs'
 * '<S388>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/N Copy'
 * '<S389>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/N Gain'
 * '<S390>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/P Copy'
 * '<S391>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Parallel P Gain'
 * '<S392>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Reset Signal'
 * '<S393>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Saturation'
 * '<S394>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Saturation Fdbk'
 * '<S395>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Sum'
 * '<S396>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Sum Fdbk'
 * '<S397>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Tracking Mode'
 * '<S398>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Tracking Mode Sum'
 * '<S399>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Tsamp - Integral'
 * '<S400>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Tsamp - Ngain'
 * '<S401>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/postSat Signal'
 * '<S402>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/preInt Signal'
 * '<S403>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/preSat Signal'
 * '<S404>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Anti-windup/Passthrough'
 * '<S405>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/D Gain/Internal Parameters'
 * '<S406>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/External Derivative/Error'
 * '<S407>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Filter/Disc. Forward Euler Filter'
 * '<S408>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Filter ICs/Internal IC - Filter'
 * '<S409>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/I Gain/Internal Parameters'
 * '<S410>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Ideal P Gain/Passthrough'
 * '<S411>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Ideal P Gain Fdbk/Disabled'
 * '<S412>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Integrator/Discrete'
 * '<S413>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Integrator ICs/Internal IC'
 * '<S414>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/N Copy/Disabled'
 * '<S415>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/N Gain/Internal Parameters'
 * '<S416>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/P Copy/Disabled'
 * '<S417>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Parallel P Gain/Internal Parameters'
 * '<S418>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Reset Signal/External Reset'
 * '<S419>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Saturation/External'
 * '<S420>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Saturation/External/Saturation Dynamic'
 * '<S421>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Saturation Fdbk/Disabled'
 * '<S422>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Sum/Sum_PID'
 * '<S423>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Sum Fdbk/Disabled'
 * '<S424>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Tracking Mode/Disabled'
 * '<S425>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Tracking Mode Sum/Passthrough'
 * '<S426>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Tsamp - Integral/TsSignalSpecification'
 * '<S427>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Tsamp - Ngain/Passthrough'
 * '<S428>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/postSat Signal/Forward_Path'
 * '<S429>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/preInt Signal/Internal PreInt'
 * '<S430>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/preSat Signal/Forward_Path'
 * '<S431>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem1/Subsystem'
 * '<S432>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem1/Subsystem1'
 * '<S433>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem1/Subsystem/S-R Flip-Flop'
 * '<S434>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem1/Subsystem1/S-R Flip-Flop'
 * '<S435>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/PressureControl1'
 * '<S436>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1'
 * '<S437>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/TorqueControl1'
 * '<S438>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem'
 * '<S439>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem1'
 * '<S440>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller'
 * '<S441>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Anti-windup'
 * '<S442>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/D Gain'
 * '<S443>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/External Derivative'
 * '<S444>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Filter'
 * '<S445>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Filter ICs'
 * '<S446>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/I Gain'
 * '<S447>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Ideal P Gain'
 * '<S448>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Ideal P Gain Fdbk'
 * '<S449>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Integrator'
 * '<S450>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Integrator ICs'
 * '<S451>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/N Copy'
 * '<S452>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/N Gain'
 * '<S453>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/P Copy'
 * '<S454>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Parallel P Gain'
 * '<S455>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Reset Signal'
 * '<S456>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Saturation'
 * '<S457>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Saturation Fdbk'
 * '<S458>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Sum'
 * '<S459>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Sum Fdbk'
 * '<S460>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Tracking Mode'
 * '<S461>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Tracking Mode Sum'
 * '<S462>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Tsamp - Integral'
 * '<S463>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Tsamp - Ngain'
 * '<S464>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/postSat Signal'
 * '<S465>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/preInt Signal'
 * '<S466>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/preSat Signal'
 * '<S467>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Anti-windup/Passthrough'
 * '<S468>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/D Gain/Internal Parameters'
 * '<S469>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/External Derivative/Error'
 * '<S470>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Filter/Disc. Forward Euler Filter'
 * '<S471>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Filter ICs/Internal IC - Filter'
 * '<S472>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/I Gain/Internal Parameters'
 * '<S473>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Ideal P Gain/Passthrough'
 * '<S474>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Ideal P Gain Fdbk/Disabled'
 * '<S475>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Integrator/Discrete'
 * '<S476>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Integrator ICs/Internal IC'
 * '<S477>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/N Copy/Disabled'
 * '<S478>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/N Gain/Internal Parameters'
 * '<S479>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/P Copy/Disabled'
 * '<S480>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Parallel P Gain/Internal Parameters'
 * '<S481>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Reset Signal/External Reset'
 * '<S482>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Saturation/External'
 * '<S483>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Saturation/External/Saturation Dynamic'
 * '<S484>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Saturation Fdbk/Disabled'
 * '<S485>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Sum/Sum_PID'
 * '<S486>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Sum Fdbk/Disabled'
 * '<S487>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Tracking Mode/Disabled'
 * '<S488>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Tracking Mode Sum/Passthrough'
 * '<S489>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Tsamp - Integral/TsSignalSpecification'
 * '<S490>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Tsamp - Ngain/Passthrough'
 * '<S491>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/postSat Signal/Forward_Path'
 * '<S492>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/preInt Signal/Internal PreInt'
 * '<S493>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/preSat Signal/Forward_Path'
 * '<S494>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem1/Subsystem'
 * '<S495>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem1/Subsystem1'
 * '<S496>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem1/Subsystem/S-R Flip-Flop'
 * '<S497>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem1/Subsystem1/S-R Flip-Flop'
 * '<S498>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/ChargeCircuit'
 * '<S499>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/FlowSensor'
 * '<S500>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/FlowSensor1'
 * '<S501>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/FlowSensor2'
 * '<S502>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/LinearVelInput'
 * '<S503>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/PS-Simulink Converter'
 * '<S504>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/PressureSensor'
 * '<S505>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Simulink-PS Converter'
 * '<S506>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Simulink-PS Converter1'
 * '<S507>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Simulink-PS Converter2'
 * '<S508>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Solver Configuration'
 * '<S509>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Subsystem1'
 * '<S510>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Subsystem2'
 * '<S511>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Subsystem3'
 * '<S512>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/ChargeCircuit/Simulink-PS Converter1'
 * '<S513>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/ChargeCircuit/Simulink-PS Converter1/EVAL_KEY'
 * '<S514>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/FlowSensor/PS-Simulink Converter'
 * '<S515>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/FlowSensor/PS-Simulink Converter/EVAL_KEY'
 * '<S516>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/FlowSensor1/PS-Simulink Converter'
 * '<S517>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/FlowSensor1/PS-Simulink Converter/EVAL_KEY'
 * '<S518>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/FlowSensor2/PS-Simulink Converter'
 * '<S519>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/FlowSensor2/PS-Simulink Converter/EVAL_KEY'
 * '<S520>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/LinearVelInput/PS-Simulink Converter'
 * '<S521>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/LinearVelInput/Subsystem'
 * '<S522>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/LinearVelInput/WavebotModel'
 * '<S523>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/LinearVelInput/PS-Simulink Converter/EVAL_KEY'
 * '<S524>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/LinearVelInput/Subsystem/PS-Simulink Converter'
 * '<S525>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/LinearVelInput/Subsystem/PS-Simulink Converter1'
 * '<S526>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/LinearVelInput/Subsystem/PS-Simulink Converter2'
 * '<S527>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/LinearVelInput/Subsystem/PS-Simulink Converter/EVAL_KEY'
 * '<S528>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/LinearVelInput/Subsystem/PS-Simulink Converter1/EVAL_KEY'
 * '<S529>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/LinearVelInput/Subsystem/PS-Simulink Converter2/EVAL_KEY'
 * '<S530>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/LinearVelInput/WavebotModel/Simulink-PS Converter'
 * '<S531>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/LinearVelInput/WavebotModel/WaveBotTF'
 * '<S532>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/LinearVelInput/WavebotModel/Simulink-PS Converter/EVAL_KEY'
 * '<S533>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/LinearVelInput/WavebotModel/WaveBotTF/Input Delay'
 * '<S534>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/LinearVelInput/WavebotModel/WaveBotTF/Output Delay'
 * '<S535>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/PS-Simulink Converter/EVAL_KEY'
 * '<S536>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/PressureSensor/PS-Simulink Converter'
 * '<S537>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/PressureSensor/PS-Simulink Converter/EVAL_KEY'
 * '<S538>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Simulink-PS Converter/EVAL_KEY'
 * '<S539>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Simulink-PS Converter1/EVAL_KEY'
 * '<S540>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Simulink-PS Converter2/EVAL_KEY'
 * '<S541>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Solver Configuration/EVAL_KEY'
 * '<S542>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Subsystem1/LTI System2'
 * '<S543>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Subsystem1/LTI System2/Input Delay'
 * '<S544>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Subsystem1/LTI System2/Output Delay'
 * '<S545>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Subsystem2/LTI System2'
 * '<S546>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Subsystem2/LTI System2/Input Delay'
 * '<S547>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Subsystem2/LTI System2/Output Delay'
 * '<S548>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Subsystem3/PS-Simulink Converter'
 * '<S549>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Subsystem3/PS-Simulink Converter1'
 * '<S550>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Subsystem3/PS-Simulink Converter/EVAL_KEY'
 * '<S551>' : 'windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Subsystem3/PS-Simulink Converter1/EVAL_KEY'
 * '<S552>' : 'windEmulatorStep4_WECSim/readShaftSignals/countsToRads'
 * '<S553>' : 'windEmulatorStep4_WECSim/readShaftSignals/posToVel'
 * '<S554>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller'
 * '<S555>' : 'windEmulatorStep4_WECSim/sidCtrl/setpointGenerator'
 * '<S556>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Anti-windup'
 * '<S557>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/D Gain'
 * '<S558>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/External Derivative'
 * '<S559>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Filter'
 * '<S560>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Filter ICs'
 * '<S561>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/I Gain'
 * '<S562>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Ideal P Gain'
 * '<S563>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Ideal P Gain Fdbk'
 * '<S564>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Integrator'
 * '<S565>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Integrator ICs'
 * '<S566>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/N Copy'
 * '<S567>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/N Gain'
 * '<S568>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/P Copy'
 * '<S569>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Parallel P Gain'
 * '<S570>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Reset Signal'
 * '<S571>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Saturation'
 * '<S572>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Saturation Fdbk'
 * '<S573>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Sum'
 * '<S574>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Sum Fdbk'
 * '<S575>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Tracking Mode'
 * '<S576>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Tracking Mode Sum'
 * '<S577>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Tsamp - Integral'
 * '<S578>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Tsamp - Ngain'
 * '<S579>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/postSat Signal'
 * '<S580>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/preInt Signal'
 * '<S581>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/preSat Signal'
 * '<S582>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Anti-windup/Passthrough'
 * '<S583>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/D Gain/Disabled'
 * '<S584>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/External Derivative/Disabled'
 * '<S585>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Filter/Disabled'
 * '<S586>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Filter ICs/Disabled'
 * '<S587>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/I Gain/External Parameters'
 * '<S588>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Ideal P Gain/Passthrough'
 * '<S589>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Ideal P Gain Fdbk/Disabled'
 * '<S590>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Integrator/Continuous'
 * '<S591>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Integrator ICs/Internal IC'
 * '<S592>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/N Copy/Disabled wSignal Specification'
 * '<S593>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/N Gain/Disabled'
 * '<S594>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/P Copy/Disabled'
 * '<S595>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Parallel P Gain/External Parameters'
 * '<S596>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Reset Signal/External Reset'
 * '<S597>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Saturation/External'
 * '<S598>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Saturation/External/Saturation Dynamic'
 * '<S599>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Saturation Fdbk/Disabled'
 * '<S600>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Sum/Sum_PI'
 * '<S601>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Sum Fdbk/Disabled'
 * '<S602>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Tracking Mode/Disabled'
 * '<S603>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Tracking Mode Sum/Passthrough'
 * '<S604>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Tsamp - Integral/TsSignalSpecification'
 * '<S605>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/Tsamp - Ngain/Passthrough'
 * '<S606>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/postSat Signal/Forward_Path'
 * '<S607>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/preInt Signal/Internal PreInt'
 * '<S608>' : 'windEmulatorStep4_WECSim/sidCtrl/PID Controller/preSat Signal/Forward_Path'
 * '<S609>' : 'windEmulatorStep4_WECSim/sidCtrl/setpointGenerator/fromFile'
 * '<S610>' : 'windEmulatorStep4_WECSim/writeAcs800Pdos/ACS800ConvertREF1'
 * '<S611>' : 'windEmulatorStep4_WECSim/writeAcs800Pdos/ACS800ConvertREF2'
 * '<S612>' : 'windEmulatorStep4_WECSim/writeAcs880Pdos/ACS880ConvertREF2'
 */
#endif                                 /* windEmulatorStep4_WECSim_h_ */
