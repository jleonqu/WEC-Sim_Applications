/*
 * windEmulatorStep4.h
 *
 * Code generation for model "windEmulatorStep4".
 *
 * Model version              : 10.2
 * Simulink Coder version : 25.2 (R2025b) 28-Jul-2025
 * C++ source code generated on : Fri Jun 26 11:42:45 2026
 *
 * Target selection: speedgoat.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef windEmulatorStep4_h_
#define windEmulatorStep4_h_
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
#include "windEmulatorStep4_1e9c788f_1_gateway.h"
#include "nesl_rtw.h"
#include "windEmulatorStep4_types.h"
#include "abbState.h"
#include "expType.h"

extern "C"
{

#include "rt_nonfinite.h"

}

#include <stddef.h>

extern "C"
{

#include "rtGetInf.h"

}

#include <cstring>
#include "windEmulatorStep4_cal.h"
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

#ifndef rtmGetOdeF
#define rtmGetOdeF(rtm)                ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
#define rtmSetOdeF(rtm, val)           ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
#define rtmGetOdeY(rtm)                ((rtm)->odeY)
#endif

#ifndef rtmSetOdeY
#define rtmSetOdeY(rtm, val)           ((rtm)->odeY = (val))
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

/* Block signals (default storage) */
struct B_windEmulatorStep4_T {
  invPowerBus BusAssignment;           /* '<S9>/Bus Assignment' */
  invPowerBus BusAssignment_l;         /* '<S11>/Bus Assignment' */
  hptoSignalBus BusAssignment_n;       /* '<S6>/Bus Assignment' */
  acs880SignalBus BusAssignment_a;     /* '<S10>/Bus Assignment' */
  acs800SignalBus BusAssignment_h;     /* '<S8>/Bus Assignment' */
  shaftSignalBus BusAssignment_g;      /* '<S12>/Bus Assignment' */
  sidInfoBus BusAssignment_j;          /* '<S243>/Bus Assignment' */
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
  real_T rampValue;                    /* '<S243>/Switch' */
  real_T Switch2;                      /* '<S297>/Switch2' */
  real_T torqueSlewRate;               /* '<S297>/torqueSlewRate' */
  real_T MultiportSwitch;              /* '<S243>/Multiport Switch' */
  real_T Product_g;                    /* '<S243>/Product' */
  real_T Switch1;                      /* '<S297>/Switch1' */
  real_T speedSlewRate;                /* '<S297>/speedSlewRate' */
  real_T MultiportSwitch1;             /* '<S243>/Multiport Switch1' */
  real_T Product1;                     /* '<S243>/Product1' */
  real_T Sum;                          /* '<S13>/Sum' */
  real_T PProdOut;                     /* '<S283>/PProd Out' */
  real_T Integrator;                   /* '<S278>/Integrator' */
  real_T Sum_k;                        /* '<S288>/Sum' */
  real_T Switch;                       /* '<S286>/Switch' */
  real_T Switch2_h;                    /* '<S286>/Switch2' */
  real_T RateLimiter;                  /* '<S7>/Rate Limiter' */
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
  real_T Add;                          /* '<S64>/Add' */
  real_T ControlSignal31;              /* '<S119>/Product' */
  real_T RateLimiter1;                 /* '<S61>/Rate Limiter1' */
  real_T Add1;                         /* '<S64>/Add1' */
  real_T ControlSignal31_o;            /* '<S120>/Product' */
  real_T Switch_g;                     /* '<S64>/Switch' */
  real_T Gain_l;                       /* '<S64>/Gain' */
  real_T RateLimiter_b;                /* '<S61>/Rate Limiter' */
  real_T wError;                       /* '<S63>/Sum' */
  real_T ProportionalGain;             /* '<S105>/Proportional Gain' */
  real_T Integrator_b;                 /* '<S100>/Integrator' */
  real_T DerivativeGain;               /* '<S93>/Derivative Gain' */
  real_T Filter;                       /* '<S95>/Filter' */
  real_T SumD;                         /* '<S95>/SumD' */
  real_T FilterCoefficient;            /* '<S103>/Filter Coefficient' */
  real_T Sum_c;                        /* '<S110>/Sum' */
  real_T Switch_i;                     /* '<S108>/Switch' */
  real_T Switch2_k;                    /* '<S108>/Switch2' */
  real_T ContolTorque;                 /* '<S63>/Gain2' */
  real_T Step;                         /* '<S58>/Step' */
  real_T Clock;                        /* '<S58>/Clock' */
  real_T Sum_m;                        /* '<S58>/Sum' */
  real_T Product_a;                    /* '<S58>/Product' */
  real_T Output;                       /* '<S58>/Output' */
  real_T Saturation;                   /* '<S53>/Saturation' */
  real_T kDampingNow;                  /* '<S59>/kDampingNow' */
  real_T RateLimiter_a;                /* '<S186>/Rate Limiter' */
  real_T INPUT_1_1_1[4];               /* '<S229>/INPUT_1_1_1' */
  real_T Internal;                     /* '<S219>/Internal' */
  real_T INPUT_2_1_1[4];               /* '<S229>/INPUT_2_1_1' */
  real_T Internal_j;                   /* '<S233>/Internal' */
  real_T Gain_lp;                      /* '<S198>/Gain' */
  real_T INPUT_4_1_1[4];               /* '<S229>/INPUT_4_1_1' */
  real_T Gain_d;                       /* '<S54>/Gain' */
  real_T INPUT_5_1_1[4];               /* '<S229>/INPUT_5_1_1' */
  real_T Internal_h;                   /* '<S230>/Internal' */
  real_T Gain_g;                       /* '<S197>/Gain' */
  real_T INPUT_3_1_1[4];               /* '<S229>/INPUT_3_1_1' */
  real_T RTP_1;                        /* '<S196>/RTP_1' */
  real_T STATE_1[37];                  /* '<S229>/STATE_1' */
  real_T OUTPUT_1_0[11];               /* '<S229>/OUTPUT_1_0' */
  real_T Product_mo;                   /* '<S59>/Product' */
  real_T kSpringNow;                   /* '<S59>/kSpringNow' */
  real_T Product1_g;                   /* '<S59>/Product1' */
  real_T ForceD;                       /* '<S59>/Add' */
  real_T Product2;                     /* '<S55>/Product2' */
  real_T Gain_n;                       /* '<S55>/Gain' */
  real_T TorqueInputRef;               /* '<S53>/Product' */
  real_T Abs;                          /* '<S53>/Abs' */
  real_T Add_m;                        /* '<S127>/Add' */
  real_T ControlSignal31_d;            /* '<S182>/Product' */
  real_T RateLimiter1_m;               /* '<S124>/Rate Limiter1' */
  real_T Add1_f;                       /* '<S127>/Add1' */
  real_T ControlSignal31_m;            /* '<S183>/Product' */
  real_T Switch_n;                     /* '<S127>/Switch' */
  real_T Gain_f;                       /* '<S127>/Gain' */
  real_T RateLimiter_aw;               /* '<S124>/Rate Limiter' */
  real_T wError_c;                     /* '<S126>/Sum' */
  real_T ProportionalGain_h;           /* '<S168>/Proportional Gain' */
  real_T Integrator_l;                 /* '<S163>/Integrator' */
  real_T DerivativeGain_g;             /* '<S156>/Derivative Gain' */
  real_T Filter_j;                     /* '<S158>/Filter' */
  real_T SumD_c;                       /* '<S158>/SumD' */
  real_T FilterCoefficient_g;          /* '<S166>/Filter Coefficient' */
  real_T Sum_n;                        /* '<S173>/Sum' */
  real_T Switch_c;                     /* '<S171>/Switch' */
  real_T Switch2_m;                    /* '<S171>/Switch2' */
  real_T ContolTorque_f;               /* '<S126>/Gain2' */
  real_T ControlTorqueLoad;            /* '<S53>/Switch' */
  real_T Pressure;                     /* '<S192>/Gain' */
  real_T psibar;                       /* '<S6>/psi -> bar' */
  real_T ShaftSpeedPump;               /* '<S199>/Gain' */
  real_T Sum_f;                        /* '<S125>/Sum' */
  real_T DiscreteTimeIntegrator;       /* '<S125>/Discrete-Time Integrator' */
  real_T Gain1_b;                      /* '<S125>/Gain1' */
  real_T ControlSignal1;               /* '<S53>/Switch' */
  real_T Abs2;                         /* '<S53>/Abs2' */
  real_T PressureRef;                  /* '<S53>/Switch2' */
  real_T Sum_me;                       /* '<S60>/Sum' */
  real_T DiscreteTimeIntegrator_i;     /* '<S60>/Discrete-Time Integrator' */
  real_T Gain1_j;                      /* '<S60>/Gain1' */
  real_T Sum_d;                        /* '<S123>/Sum' */
  real_T DiscreteTimeIntegrator_e;     /* '<S123>/Discrete-Time Integrator' */
  real_T Gain1_h;                      /* '<S123>/Gain1' */
  real_T ControlSignal2;               /* '<S53>/Switch' */
  real_T FlowMotor1;                   /* '<S189>/Gain' */
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
  real_T CastToDouble3_b;              /* '<S240>/Cast To Double3' */
  real_T encoderCountsToRad;           /* '<S240>/encoderCountsToRad' */
  real_T CastToDouble1_m;              /* '<S240>/Cast To Double1' */
  real_T Gain_b;                       /* '<S240>/Gain' */
  real_T Add2;                         /* '<S240>/Add2' */
  real_T TSamp;                        /* '<S241>/TSamp' */
  real_T Uk1;                          /* '<S241>/UD' */
  real_T Diff;                         /* '<S241>/Diff' */
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
  real_T Switch3;                      /* '<S297>/Switch3' */
  real_T caseCounterSignalsNow;        /* '<S297>/caseCounterSignalsNow' */
  real_T DataTypeConversion_f;         /* '<S29>/Data Type Conversion' */
  real_T time_s;                       /* '<S29>/time_s' */
  real_T Product_k;                    /* '<S299>/Product' */
  real_T Product_d;                    /* '<S298>/Product' */
  real_T Product_l;                    /* '<S300>/Product' */
  real_T IntegralGain;                 /* '<S97>/Integral Gain' */
  real_T IntegralGain_h;               /* '<S160>/Integral Gain' */
  real_T SwitchLogic;                  /* '<S53>/Switch1' */
  real_T FlowPump1;                    /* '<S187>/m3toL' */
  real_T Sum_a;                        /* '<S210>/Sum' */
  real_T FlowAccumulator;              /* '<S188>/Gain' */
  real_T IProdOut;                     /* '<S275>/IProd Out' */
  real_T CastToDouble_g;               /* '<S297>/Cast To Double' */
  real_T CastToDouble1_i;              /* '<S297>/Cast To Double1' */
  real_T Divide1;                      /* '<S297>/Divide1' */
  real_T vecPercent;                   /* '<S297>/vecPercent' */
  real_T fromFileTorqueNow_Nm;         /* '<S297>/fromFileTorqueNow_Nm' */
  real_T fromFileSpeedNow_rpm;         /* '<S297>/fromFileSpeedNow_rpm' */
  real_T Gain_k;                       /* '<S53>/Gain' */
  real_T Abs1;                         /* '<S53>/Abs1' */
  real_T Switch_f;                     /* '<S124>/Switch' */
  real_T Switch_cn;                    /* '<S61>/Switch' */
  real_T Gain_j;                       /* '<S123>/Gain' */
  real_T Add_k;                        /* '<S123>/Add' */
  real_T Saturation_f;                 /* '<S123>/Saturation' */
  real_T Gain_ge;                      /* '<S60>/Gain' */
  real_T Add_i;                        /* '<S60>/Add' */
  real_T Saturation_a;                 /* '<S60>/Saturation' */
  real_T Gain2_n;                      /* '<S125>/Gain2' */
  real_T Divide;                       /* '<S125>/Divide' */
  real_T Product_aq;                   /* '<S125>/Product' */
  real_T Gain_a;                       /* '<S125>/Gain' */
  real_T Add_iu;                       /* '<S125>/Add' */
  real_T Add1_l;                       /* '<S125>/Add1' */
  real_T Saturation_j;                 /* '<S125>/Saturation' */
  real_T Switch_ne;                    /* '<S62>/Switch' */
  real_T Switch1_g;                    /* '<S127>/Switch1' */
  real_T ControlSignal3;               /* '<S183>/Switch' */
  real_T Saturation_p;                 /* '<S183>/Saturation' */
  real_T ControlSignal3_h;             /* '<S182>/Switch' */
  real_T Saturation_af;                /* '<S182>/Saturation' */
  real_T Switch1_o;                    /* '<S64>/Switch1' */
  real_T ControlSignal3_e;             /* '<S120>/Switch' */
  real_T Saturation_k;                 /* '<S120>/Saturation' */
  real_T ControlSignal3_f;             /* '<S119>/Switch' */
  real_T Saturation_a3;                /* '<S119>/Saturation' */
  real_T time_c;                       /* '<S4>/FexcRamp' */
  real_T ramp_l;                       /* '<S4>/FexcRamp' */
  uint32_T MultiportSwitch_h;          /* '<S297>/Multiport Switch' */
  uint32_T Mod;                        /* '<S297>/Mod' */
  uint32_T runCounter;                 /* '<S31>/runCounter' */
  uint32_T stepCounter;                /* '<S31>/stepCounter' */
  uint32_T readEncoderCounter;         /* '<S12>/readEncoderCounter' */
  uint32_T absEncoderCounts;           /* '<S39>/absEncoderCounts' */
  uint32_T CastTouint32;               /* '<S243>/Cast To uint32' */
  uint32_T Memory;                     /* '<S29>/Memory' */
  uint32_T Sum_e;                      /* '<S29>/Sum' */
  uint32_T loopCounter;                /* '<S29>/loopCounter' */
  uint32_T vecIndex;                   /* '<S297>/vecIndex' */
  uint32_T fileSamples;                /* '<S297>/fileSamples' */
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
  int32_T lastRawCounts;               /* '<S240>/lastRawCounts' */
  int32_T CastToDouble_p;              /* '<S240>/Cast To Double' */
  int32_T Add_d;                       /* '<S240>/Add' */
  int32_T Abs_a;                       /* '<S240>/Abs' */
  int32_T Switch_g0;                   /* '<S240>/Switch' */
  int32_T lastTurn;                    /* '<S240>/lastTurn' */
  int32_T Add1_p;                      /* '<S240>/Add1' */
  int32_T absEncoderTurns;             /* '<S39>/absEncoderTurns' */
  int32_T Sign;                        /* '<S240>/Sign' */
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
  boolean_T Outofbounds;               /* '<S297>/Out of bounds' */
  boolean_T LowerRelop1;               /* '<S286>/LowerRelop1' */
  boolean_T UpperRelop;                /* '<S286>/UpperRelop' */
  boolean_T RelationalOperator;        /* '<S119>/Relational Operator' */
  boolean_T Memory_f;                  /* '<S121>/Memory' */
  boolean_T Logic[2];                  /* '<S121>/Logic' */
  boolean_T RelationalOperator_k;      /* '<S120>/Relational Operator' */
  boolean_T Memory_in;                 /* '<S122>/Memory' */
  boolean_T Logic_g[2];                /* '<S122>/Logic' */
  boolean_T LowerRelop1_g;             /* '<S108>/LowerRelop1' */
  boolean_T UpperRelop_g;              /* '<S108>/UpperRelop' */
  boolean_T RelationalOperator_e;      /* '<S182>/Relational Operator' */
  boolean_T Memory_o;                  /* '<S184>/Memory' */
  boolean_T Logic_c[2];                /* '<S184>/Logic' */
  boolean_T RelationalOperator_g;      /* '<S183>/Relational Operator' */
  boolean_T Memory_a;                  /* '<S185>/Memory' */
  boolean_T Logic_p[2];                /* '<S185>/Logic' */
  boolean_T LowerRelop1_h;             /* '<S171>/LowerRelop1' */
  boolean_T UpperRelop_m;              /* '<S171>/UpperRelop' */
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
  boolean_T RelationalOperator_l;      /* '<S240>/Relational Operator' */
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
  boolean_T resetHilIntegrator_h;      /* '<S4>/FexcRamp' */
  boolean_T resetSidIntegrator_h;      /* '<S4>/FexcRamp' */
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
struct DW_windEmulatorStep4_T {
  real_T Integrator_DSTATE;            /* '<S100>/Integrator' */
  real_T Filter_DSTATE;                /* '<S95>/Filter' */
  real_T INPUT_1_1_1_Discrete_142920204[2];/* '<S229>/INPUT_1_1_1' */
  real_T INPUT_2_1_1_Discrete_1327804636[2];/* '<S229>/INPUT_2_1_1' */
  real_T INPUT_4_1_1_Discrete_3227796860[2];/* '<S229>/INPUT_4_1_1' */
  real_T INPUT_5_1_1_Discrete_4244925644[2];/* '<S229>/INPUT_5_1_1' */
  real_T INPUT_3_1_1_Discrete_1917098348[2];/* '<S229>/INPUT_3_1_1' */
  real_T STATE_1_Discrete_1041191992[22];/* '<S229>/STATE_1' */
  real_T Integrator_DSTATE_e;          /* '<S163>/Integrator' */
  real_T Filter_DSTATE_b;              /* '<S158>/Filter' */
  real_T DiscreteTimeIntegrator_DSTATE;/* '<S125>/Discrete-Time Integrator' */
  real_T DiscreteTimeIntegrator_DSTATE_n;/* '<S60>/Discrete-Time Integrator' */
  real_T DiscreteTimeIntegrator_DSTATE_l;/* '<S123>/Discrete-Time Integrator' */
  real_T UD_DSTATE;                    /* '<S241>/UD' */
  real_T PrevY;                        /* '<S297>/torqueSlewRate' */
  real_T LastMajorTime;                /* '<S297>/torqueSlewRate' */
  real_T PrevY_a;                      /* '<S297>/speedSlewRate' */
  real_T LastMajorTime_j;              /* '<S297>/speedSlewRate' */
  real_T PrevY_f;                      /* '<S7>/Rate Limiter' */
  real_T PrevY_m;                      /* '<S61>/Rate Limiter1' */
  real_T PrevY_l;                      /* '<S61>/Rate Limiter' */
  real_T LastMajorTime_d;              /* '<S61>/Rate Limiter' */
  real_T PrevY_k;                      /* '<S186>/Rate Limiter' */
  real_T STATE_1_ZcValueStore;         /* '<S229>/STATE_1' */
  real_T OUTPUT_1_0_Discrete;          /* '<S229>/OUTPUT_1_0' */
  real_T OUTPUT_1_0_ZcValueStore;      /* '<S229>/OUTPUT_1_0' */
  real_T PrevY_fq;                     /* '<S124>/Rate Limiter1' */
  real_T PrevY_e;                      /* '<S124>/Rate Limiter' */
  real_T LastMajorTime_k;              /* '<S124>/Rate Limiter' */
  real_T PrevY_b;                      /* '<S2>/acs880RateLim' */
  real_T LastMajorTime_a;              /* '<S2>/acs880RateLim' */
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

  struct {
    void *AQHandles;
    void *SLRTSigHandles;
  } TAQSigLogging_InsertedFor_acs88;   /* synthesized block */

  void* RTP_1_RtpManager;              /* '<S196>/RTP_1' */
  void* STATE_1_Simulator;             /* '<S229>/STATE_1' */
  void* STATE_1_SimData;               /* '<S229>/STATE_1' */
  void* STATE_1_DiagMgr;               /* '<S229>/STATE_1' */
  void* STATE_1_ZcLogger;              /* '<S229>/STATE_1' */
  void* STATE_1_TsInfo;                /* '<S229>/STATE_1' */
  void* OUTPUT_1_0_Simulator;          /* '<S229>/OUTPUT_1_0' */
  void* OUTPUT_1_0_SimData;            /* '<S229>/OUTPUT_1_0' */
  void* OUTPUT_1_0_DiagMgr;            /* '<S229>/OUTPUT_1_0' */
  void* OUTPUT_1_0_ZcLogger;           /* '<S229>/OUTPUT_1_0' */
  void* OUTPUT_1_0_TsInfo;             /* '<S229>/OUTPUT_1_0' */
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

  struct {
    void *LoggedData;
  } Scope_PWORK;                       /* '<S187>/Scope' */

  struct {
    void *LoggedData;
  } Scope_PWORK_e;                     /* '<S192>/Scope' */

  struct {
    void *LoggedData;
  } Scope_PWORK_b;                     /* '<S188>/Scope' */

  int32_T lastRawCounts_PreviousInput; /* '<S240>/lastRawCounts' */
  int32_T lastTurn_PreviousInput;      /* '<S240>/lastTurn' */
  int32_T sfEvent;                     /* '<S4>/FexcRamp' */
  int32_T sfEvent_l;                   /* '<S18>/ABB Fieldbus Control' */
  int32_T sfEvent_d;                   /* '<S16>/ABB Fieldbus Control' */
  uint32_T Memory_PreviousInput;       /* '<S29>/Memory' */
  uint32_T is_c3_windEmulatorStep4;    /* '<S4>/FexcRamp' */
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
  int_T STATE_1_Modes[15];             /* '<S229>/STATE_1' */
  int_T OUTPUT_1_0_Modes;              /* '<S229>/OUTPUT_1_0' */
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
  int8_T Integrator_PrevResetState;    /* '<S100>/Integrator' */
  int8_T Filter_PrevResetState;        /* '<S95>/Filter' */
  int8_T Integrator_PrevResetState_g;  /* '<S163>/Integrator' */
  int8_T Filter_PrevResetState_g;      /* '<S158>/Filter' */
  uint8_T STATE_1_ZcSignalDir;         /* '<S229>/STATE_1' */
  uint8_T STATE_1_ZcStateStore;        /* '<S229>/STATE_1' */
  uint8_T OUTPUT_1_0_ZcSignalDir;      /* '<S229>/OUTPUT_1_0' */
  uint8_T OUTPUT_1_0_ZcStateStore;     /* '<S229>/OUTPUT_1_0' */
  uint8_T is_active_c3_windEmulatorStep4;/* '<S4>/FexcRamp' */
  uint8_T is_active_c7_windEmulatorStep4;/* '<S18>/ABB Fieldbus Control' */
  uint8_T is_active_UpdateStateMachine;/* '<S18>/ABB Fieldbus Control' */
  uint8_T is_active_UpdateControlWord; /* '<S18>/ABB Fieldbus Control' */
  uint8_T temporalCounter_i1_l;        /* '<S18>/ABB Fieldbus Control' */
  uint8_T is_active_c9_windEmulatorStep4;/* '<S16>/ABB Fieldbus Control' */
  uint8_T is_active_UpdateStateMachine_a;/* '<S16>/ABB Fieldbus Control' */
  uint8_T is_active_UpdateControlWord_p;/* '<S16>/ABB Fieldbus Control' */
  uint8_T temporalCounter_i1_g;        /* '<S16>/ABB Fieldbus Control' */
  boolean_T Memory_PreviousInput_i;    /* '<S2>/Memory' */
  boolean_T Memory1_PreviousInput;     /* '<S2>/Memory1' */
  boolean_T Memory2_PreviousInput;     /* '<S2>/Memory2' */
  boolean_T Memory_PreviousInput_k;    /* '<S4>/Memory' */
  boolean_T Memory1_PreviousInput_d;   /* '<S4>/Memory1' */
  boolean_T Memory2_PreviousInput_l;   /* '<S4>/Memory2' */
  boolean_T PrevLimited;               /* '<S297>/torqueSlewRate' */
  boolean_T PrevLimited_d;             /* '<S297>/speedSlewRate' */
  boolean_T Memory_PreviousInput_kk;   /* '<S121>/Memory' */
  boolean_T Memory_PreviousInput_h;    /* '<S122>/Memory' */
  boolean_T PrevLimited_g;             /* '<S61>/Rate Limiter' */
  boolean_T RTP_1_SetParametersNeeded; /* '<S196>/RTP_1' */
  boolean_T STATE_1_FirstOutput;       /* '<S229>/STATE_1' */
  boolean_T OUTPUT_1_0_FirstOutput;    /* '<S229>/OUTPUT_1_0' */
  boolean_T Memory_PreviousInput_g;    /* '<S184>/Memory' */
  boolean_T Memory_PreviousInput_n;    /* '<S185>/Memory' */
  boolean_T PrevLimited_dz;            /* '<S124>/Rate Limiter' */
  boolean_T PrevLimited_o;             /* '<S2>/acs880RateLim' */
  boolean_T Memory_PreviousInput_d;    /* '<S1>/Memory' */
  boolean_T Memory1_PreviousInput_p;   /* '<S1>/Memory1' */
  boolean_T Memory2_PreviousInput_h;   /* '<S1>/Memory2' */
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
struct X_windEmulatorStep4_T {
  real_T Integrator_CSTATE;            /* '<S278>/Integrator' */
  real_T Internal_CSTATE[3];           /* '<S219>/Internal' */
  real_T Internal_CSTATE_j;            /* '<S233>/Internal' */
  real_T Internal_CSTATE_a;            /* '<S230>/Internal' */
};

/* State derivatives (default storage) */
struct XDot_windEmulatorStep4_T {
  real_T Integrator_CSTATE;            /* '<S278>/Integrator' */
  real_T Internal_CSTATE[3];           /* '<S219>/Internal' */
  real_T Internal_CSTATE_j;            /* '<S233>/Internal' */
  real_T Internal_CSTATE_a;            /* '<S230>/Internal' */
};

/* State disabled  */
struct XDis_windEmulatorStep4_T {
  boolean_T Integrator_CSTATE;         /* '<S278>/Integrator' */
  boolean_T Internal_CSTATE[3];        /* '<S219>/Internal' */
  boolean_T Internal_CSTATE_j;         /* '<S233>/Internal' */
  boolean_T Internal_CSTATE_a;         /* '<S230>/Internal' */
};

/* Zero-crossing (trigger) state */
struct PrevZCX_windEmulatorStep4_T {
  ZCSigState Integrator_Reset_ZCE;     /* '<S278>/Integrator' */
};

#ifndef ODE4_INTG
#define ODE4_INTG

/* ODE4 Integration Data */
struct ODE4_IntgData {
  real_T *y;                           /* output */
  real_T *f[4];                        /* derivatives */
};

#endif

/* External inputs (root inport signals with default storage) */
struct ExtU_windEmulatorStep4_T {
  real_T inportSpeed_rpm;              /* '<Root>/inportSpeed_rpm' */
  real_T inportTorque_Nm;              /* '<Root>/inportTorque_Nm' */
  real_T inportCaseCounter;            /* '<Root>/inportCaseCounter' */
};

/* Real-time Model Data Structure */
struct tag_RTM_windEmulatorStep4_T {
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

  X_windEmulatorStep4_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  XDis_windEmulatorStep4_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[6];
  real_T odeF[4][6];
  ODE4_IntgData intgData;

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

  extern struct B_windEmulatorStep4_T windEmulatorStep4_B;

#ifdef __cplusplus

}

#endif

/* Continuous states (default storage) */
extern X_windEmulatorStep4_T windEmulatorStep4_X;

/* Disabled states (default storage) */
extern XDis_windEmulatorStep4_T windEmulatorStep4_XDis;

/* Block states (default storage) */
extern struct DW_windEmulatorStep4_T windEmulatorStep4_DW;

/* Zero-crossing (trigger) state */
extern PrevZCX_windEmulatorStep4_T windEmulatorStep4_PrevZCX;

#ifdef __cplusplus

extern "C"
{

#endif

  /* External inputs (root inport signals with default storage) */
  extern struct ExtU_windEmulatorStep4_T windEmulatorStep4_U;

#ifdef __cplusplus

}

#endif

#ifdef __cplusplus

extern "C"
{

#endif

  /* Model entry point functions */
  extern void windEmulatorStep4_initialize(void);
  extern void windEmulatorStep4_step(void);
  extern void windEmulatorStep4_terminate(void);

#ifdef __cplusplus

}

#endif

/* Real-time Model object */
#ifdef __cplusplus

extern "C"
{

#endif

  extern RT_MODEL_windEmulatorStep4_T *const windEmulatorStep4_M;

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
 * '<Root>' : 'windEmulatorStep4'
 * '<S1>'   : 'windEmulatorStep4/acs800Ctrl'
 * '<S2>'   : 'windEmulatorStep4/acs880Ctrl'
 * '<S3>'   : 'windEmulatorStep4/ctrlSignalSelector'
 * '<S4>'   : 'windEmulatorStep4/expCtrl'
 * '<S5>'   : 'windEmulatorStep4/fileAndUI'
 * '<S6>'   : 'windEmulatorStep4/hptoSim'
 * '<S7>'   : 'windEmulatorStep4/processHptoInputs'
 * '<S8>'   : 'windEmulatorStep4/readAcs800Pdos'
 * '<S9>'   : 'windEmulatorStep4/readAcs800Power'
 * '<S10>'  : 'windEmulatorStep4/readAcs880Pdos'
 * '<S11>'  : 'windEmulatorStep4/readAcs880Power'
 * '<S12>'  : 'windEmulatorStep4/readShaftSignals'
 * '<S13>'  : 'windEmulatorStep4/sidCtrl'
 * '<S14>'  : 'windEmulatorStep4/writeAcs800Pdos'
 * '<S15>'  : 'windEmulatorStep4/writeAcs880Pdos'
 * '<S16>'  : 'windEmulatorStep4/acs800Ctrl/ACS880FieldbusControl'
 * '<S17>'  : 'windEmulatorStep4/acs800Ctrl/ACS880FieldbusControl/ABB Fieldbus Control'
 * '<S18>'  : 'windEmulatorStep4/acs880Ctrl/ACS880FieldbusControl'
 * '<S19>'  : 'windEmulatorStep4/acs880Ctrl/ACS880FieldbusControl/ABB Fieldbus Control'
 * '<S20>'  : 'windEmulatorStep4/expCtrl/FexcRamp'
 * '<S21>'  : 'windEmulatorStep4/fileAndUI/acs800CtrlSignalsToFile'
 * '<S22>'  : 'windEmulatorStep4/fileAndUI/acs800CtrlSignalsUI'
 * '<S23>'  : 'windEmulatorStep4/fileAndUI/acs800SignalsToFile'
 * '<S24>'  : 'windEmulatorStep4/fileAndUI/acs800SignalsUI'
 * '<S25>'  : 'windEmulatorStep4/fileAndUI/acs880CtrlSignalsToFile'
 * '<S26>'  : 'windEmulatorStep4/fileAndUI/acs880CtrlSignalsUI'
 * '<S27>'  : 'windEmulatorStep4/fileAndUI/acs880SignalsToFile'
 * '<S28>'  : 'windEmulatorStep4/fileAndUI/acs880SignalsUI'
 * '<S29>'  : 'windEmulatorStep4/fileAndUI/diagnostics'
 * '<S30>'  : 'windEmulatorStep4/fileAndUI/expCtrlSignalsToFile'
 * '<S31>'  : 'windEmulatorStep4/fileAndUI/expCtrlSignalsUI'
 * '<S32>'  : 'windEmulatorStep4/fileAndUI/hptoCtrlSignals'
 * '<S33>'  : 'windEmulatorStep4/fileAndUI/hptoCtrlSignalsUI'
 * '<S34>'  : 'windEmulatorStep4/fileAndUI/hptoSignalsToFile'
 * '<S35>'  : 'windEmulatorStep4/fileAndUI/hptoSignalsUI'
 * '<S36>'  : 'windEmulatorStep4/fileAndUI/invPowerAcs800ToFile'
 * '<S37>'  : 'windEmulatorStep4/fileAndUI/invPowerAcs880ToFile'
 * '<S38>'  : 'windEmulatorStep4/fileAndUI/shaftSignalsToFile'
 * '<S39>'  : 'windEmulatorStep4/fileAndUI/shaftSignalsUI'
 * '<S40>'  : 'windEmulatorStep4/fileAndUI/sidInfoSignalsToFile'
 * '<S41>'  : 'windEmulatorStep4/fileAndUI/acs800CtrlSignalsUI/ctrlWordDetail'
 * '<S42>'  : 'windEmulatorStep4/fileAndUI/acs800CtrlSignalsUI/ctrlWordDetail/parseCtrlWord'
 * '<S43>'  : 'windEmulatorStep4/fileAndUI/acs800SignalsUI/powerCals'
 * '<S44>'  : 'windEmulatorStep4/fileAndUI/acs800SignalsUI/statusWordDetail'
 * '<S45>'  : 'windEmulatorStep4/fileAndUI/acs800SignalsUI/statusWordDetail/Parse Status Word'
 * '<S46>'  : 'windEmulatorStep4/fileAndUI/acs880CtrlSignalsUI/ctrlWordDetail'
 * '<S47>'  : 'windEmulatorStep4/fileAndUI/acs880CtrlSignalsUI/ctrlWordDetail/parseCtrlWord'
 * '<S48>'  : 'windEmulatorStep4/fileAndUI/acs880SignalsUI/powerCals'
 * '<S49>'  : 'windEmulatorStep4/fileAndUI/acs880SignalsUI/statusWordDetail'
 * '<S50>'  : 'windEmulatorStep4/fileAndUI/acs880SignalsUI/statusWordDetail/Parse Status Word'
 * '<S51>'  : 'windEmulatorStep4/fileAndUI/hptoSignalsUI/powerCalcs'
 * '<S52>'  : 'windEmulatorStep4/hptoSim/hptoModel'
 * '<S53>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm'
 * '<S54>'  : 'windEmulatorStep4/hptoSim/hptoModel/HPTO'
 * '<S55>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/ForceToTorque'
 * '<S56>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin'
 * '<S57>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin'
 * '<S58>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/Ramp'
 * '<S59>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/ForceToTorque/DesiredForce'
 * '<S60>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/PressureControl1'
 * '<S61>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1'
 * '<S62>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/TorqueControl1'
 * '<S63>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem'
 * '<S64>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem1'
 * '<S65>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller'
 * '<S66>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Anti-windup'
 * '<S67>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/D Gain'
 * '<S68>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/External Derivative'
 * '<S69>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Filter'
 * '<S70>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Filter ICs'
 * '<S71>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/I Gain'
 * '<S72>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Ideal P Gain'
 * '<S73>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Ideal P Gain Fdbk'
 * '<S74>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Integrator'
 * '<S75>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Integrator ICs'
 * '<S76>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/N Copy'
 * '<S77>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/N Gain'
 * '<S78>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/P Copy'
 * '<S79>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Parallel P Gain'
 * '<S80>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Reset Signal'
 * '<S81>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Saturation'
 * '<S82>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Saturation Fdbk'
 * '<S83>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Sum'
 * '<S84>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Sum Fdbk'
 * '<S85>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Tracking Mode'
 * '<S86>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Tracking Mode Sum'
 * '<S87>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Tsamp - Integral'
 * '<S88>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Tsamp - Ngain'
 * '<S89>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/postSat Signal'
 * '<S90>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/preInt Signal'
 * '<S91>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/preSat Signal'
 * '<S92>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Anti-windup/Passthrough'
 * '<S93>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/D Gain/Internal Parameters'
 * '<S94>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/External Derivative/Error'
 * '<S95>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Filter/Disc. Forward Euler Filter'
 * '<S96>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Filter ICs/Internal IC - Filter'
 * '<S97>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/I Gain/Internal Parameters'
 * '<S98>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Ideal P Gain/Passthrough'
 * '<S99>'  : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Ideal P Gain Fdbk/Disabled'
 * '<S100>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Integrator/Discrete'
 * '<S101>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Integrator ICs/Internal IC'
 * '<S102>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/N Copy/Disabled'
 * '<S103>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/N Gain/Internal Parameters'
 * '<S104>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/P Copy/Disabled'
 * '<S105>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Parallel P Gain/Internal Parameters'
 * '<S106>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Reset Signal/External Reset'
 * '<S107>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Saturation/External'
 * '<S108>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Saturation/External/Saturation Dynamic'
 * '<S109>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Saturation Fdbk/Disabled'
 * '<S110>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Sum/Sum_PID'
 * '<S111>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Sum Fdbk/Disabled'
 * '<S112>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Tracking Mode/Disabled'
 * '<S113>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Tracking Mode Sum/Passthrough'
 * '<S114>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Tsamp - Integral/TsSignalSpecification'
 * '<S115>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/Tsamp - Ngain/Passthrough'
 * '<S116>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/postSat Signal/Forward_Path'
 * '<S117>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/preInt Signal/Internal PreInt'
 * '<S118>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem/PID Controller/preSat Signal/Forward_Path'
 * '<S119>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem1/Subsystem'
 * '<S120>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem1/Subsystem1'
 * '<S121>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem1/Subsystem/S-R Flip-Flop'
 * '<S122>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureAbovePmin/SpeedControl1/Subsystem1/Subsystem1/S-R Flip-Flop'
 * '<S123>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/PressureControl1'
 * '<S124>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1'
 * '<S125>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/TorqueControl1'
 * '<S126>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem'
 * '<S127>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem1'
 * '<S128>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller'
 * '<S129>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Anti-windup'
 * '<S130>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/D Gain'
 * '<S131>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/External Derivative'
 * '<S132>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Filter'
 * '<S133>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Filter ICs'
 * '<S134>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/I Gain'
 * '<S135>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Ideal P Gain'
 * '<S136>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Ideal P Gain Fdbk'
 * '<S137>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Integrator'
 * '<S138>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Integrator ICs'
 * '<S139>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/N Copy'
 * '<S140>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/N Gain'
 * '<S141>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/P Copy'
 * '<S142>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Parallel P Gain'
 * '<S143>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Reset Signal'
 * '<S144>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Saturation'
 * '<S145>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Saturation Fdbk'
 * '<S146>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Sum'
 * '<S147>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Sum Fdbk'
 * '<S148>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Tracking Mode'
 * '<S149>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Tracking Mode Sum'
 * '<S150>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Tsamp - Integral'
 * '<S151>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Tsamp - Ngain'
 * '<S152>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/postSat Signal'
 * '<S153>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/preInt Signal'
 * '<S154>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/preSat Signal'
 * '<S155>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Anti-windup/Passthrough'
 * '<S156>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/D Gain/Internal Parameters'
 * '<S157>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/External Derivative/Error'
 * '<S158>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Filter/Disc. Forward Euler Filter'
 * '<S159>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Filter ICs/Internal IC - Filter'
 * '<S160>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/I Gain/Internal Parameters'
 * '<S161>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Ideal P Gain/Passthrough'
 * '<S162>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Ideal P Gain Fdbk/Disabled'
 * '<S163>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Integrator/Discrete'
 * '<S164>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Integrator ICs/Internal IC'
 * '<S165>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/N Copy/Disabled'
 * '<S166>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/N Gain/Internal Parameters'
 * '<S167>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/P Copy/Disabled'
 * '<S168>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Parallel P Gain/Internal Parameters'
 * '<S169>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Reset Signal/External Reset'
 * '<S170>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Saturation/External'
 * '<S171>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Saturation/External/Saturation Dynamic'
 * '<S172>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Saturation Fdbk/Disabled'
 * '<S173>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Sum/Sum_PID'
 * '<S174>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Sum Fdbk/Disabled'
 * '<S175>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Tracking Mode/Disabled'
 * '<S176>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Tracking Mode Sum/Passthrough'
 * '<S177>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Tsamp - Integral/TsSignalSpecification'
 * '<S178>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/Tsamp - Ngain/Passthrough'
 * '<S179>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/postSat Signal/Forward_Path'
 * '<S180>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/preInt Signal/Internal PreInt'
 * '<S181>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem/PID Controller/preSat Signal/Forward_Path'
 * '<S182>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem1/Subsystem'
 * '<S183>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem1/Subsystem1'
 * '<S184>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem1/Subsystem/S-R Flip-Flop'
 * '<S185>' : 'windEmulatorStep4/hptoSim/hptoModel/ControlAlgorithm/PressureBelowPmin/SpeedControl1/Subsystem1/Subsystem1/S-R Flip-Flop'
 * '<S186>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/ChargeCircuit'
 * '<S187>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/FlowSensor'
 * '<S188>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/FlowSensor1'
 * '<S189>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/FlowSensor2'
 * '<S190>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/LinearVelInput'
 * '<S191>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/PS-Simulink Converter'
 * '<S192>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/PressureSensor'
 * '<S193>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/Simulink-PS Converter'
 * '<S194>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/Simulink-PS Converter1'
 * '<S195>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/Simulink-PS Converter2'
 * '<S196>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/Solver Configuration'
 * '<S197>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/Subsystem1'
 * '<S198>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/Subsystem2'
 * '<S199>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/Subsystem3'
 * '<S200>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/ChargeCircuit/Simulink-PS Converter1'
 * '<S201>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/ChargeCircuit/Simulink-PS Converter1/EVAL_KEY'
 * '<S202>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/FlowSensor/PS-Simulink Converter'
 * '<S203>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/FlowSensor/PS-Simulink Converter/EVAL_KEY'
 * '<S204>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/FlowSensor1/PS-Simulink Converter'
 * '<S205>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/FlowSensor1/PS-Simulink Converter/EVAL_KEY'
 * '<S206>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/FlowSensor2/PS-Simulink Converter'
 * '<S207>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/FlowSensor2/PS-Simulink Converter/EVAL_KEY'
 * '<S208>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/LinearVelInput/PS-Simulink Converter'
 * '<S209>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/LinearVelInput/Subsystem'
 * '<S210>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/LinearVelInput/WavebotModel'
 * '<S211>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/LinearVelInput/PS-Simulink Converter/EVAL_KEY'
 * '<S212>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/LinearVelInput/Subsystem/PS-Simulink Converter'
 * '<S213>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/LinearVelInput/Subsystem/PS-Simulink Converter1'
 * '<S214>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/LinearVelInput/Subsystem/PS-Simulink Converter2'
 * '<S215>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/LinearVelInput/Subsystem/PS-Simulink Converter/EVAL_KEY'
 * '<S216>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/LinearVelInput/Subsystem/PS-Simulink Converter1/EVAL_KEY'
 * '<S217>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/LinearVelInput/Subsystem/PS-Simulink Converter2/EVAL_KEY'
 * '<S218>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/LinearVelInput/WavebotModel/Simulink-PS Converter'
 * '<S219>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/LinearVelInput/WavebotModel/WaveBotTF'
 * '<S220>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/LinearVelInput/WavebotModel/Simulink-PS Converter/EVAL_KEY'
 * '<S221>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/LinearVelInput/WavebotModel/WaveBotTF/Input Delay'
 * '<S222>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/LinearVelInput/WavebotModel/WaveBotTF/Output Delay'
 * '<S223>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/PS-Simulink Converter/EVAL_KEY'
 * '<S224>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/PressureSensor/PS-Simulink Converter'
 * '<S225>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/PressureSensor/PS-Simulink Converter/EVAL_KEY'
 * '<S226>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/Simulink-PS Converter/EVAL_KEY'
 * '<S227>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/Simulink-PS Converter1/EVAL_KEY'
 * '<S228>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/Simulink-PS Converter2/EVAL_KEY'
 * '<S229>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/Solver Configuration/EVAL_KEY'
 * '<S230>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/Subsystem1/LTI System2'
 * '<S231>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/Subsystem1/LTI System2/Input Delay'
 * '<S232>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/Subsystem1/LTI System2/Output Delay'
 * '<S233>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/Subsystem2/LTI System2'
 * '<S234>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/Subsystem2/LTI System2/Input Delay'
 * '<S235>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/Subsystem2/LTI System2/Output Delay'
 * '<S236>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/Subsystem3/PS-Simulink Converter'
 * '<S237>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/Subsystem3/PS-Simulink Converter1'
 * '<S238>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/Subsystem3/PS-Simulink Converter/EVAL_KEY'
 * '<S239>' : 'windEmulatorStep4/hptoSim/hptoModel/HPTO/Subsystem3/PS-Simulink Converter1/EVAL_KEY'
 * '<S240>' : 'windEmulatorStep4/readShaftSignals/countsToRads'
 * '<S241>' : 'windEmulatorStep4/readShaftSignals/posToVel'
 * '<S242>' : 'windEmulatorStep4/sidCtrl/PID Controller'
 * '<S243>' : 'windEmulatorStep4/sidCtrl/setpointGenerator'
 * '<S244>' : 'windEmulatorStep4/sidCtrl/PID Controller/Anti-windup'
 * '<S245>' : 'windEmulatorStep4/sidCtrl/PID Controller/D Gain'
 * '<S246>' : 'windEmulatorStep4/sidCtrl/PID Controller/External Derivative'
 * '<S247>' : 'windEmulatorStep4/sidCtrl/PID Controller/Filter'
 * '<S248>' : 'windEmulatorStep4/sidCtrl/PID Controller/Filter ICs'
 * '<S249>' : 'windEmulatorStep4/sidCtrl/PID Controller/I Gain'
 * '<S250>' : 'windEmulatorStep4/sidCtrl/PID Controller/Ideal P Gain'
 * '<S251>' : 'windEmulatorStep4/sidCtrl/PID Controller/Ideal P Gain Fdbk'
 * '<S252>' : 'windEmulatorStep4/sidCtrl/PID Controller/Integrator'
 * '<S253>' : 'windEmulatorStep4/sidCtrl/PID Controller/Integrator ICs'
 * '<S254>' : 'windEmulatorStep4/sidCtrl/PID Controller/N Copy'
 * '<S255>' : 'windEmulatorStep4/sidCtrl/PID Controller/N Gain'
 * '<S256>' : 'windEmulatorStep4/sidCtrl/PID Controller/P Copy'
 * '<S257>' : 'windEmulatorStep4/sidCtrl/PID Controller/Parallel P Gain'
 * '<S258>' : 'windEmulatorStep4/sidCtrl/PID Controller/Reset Signal'
 * '<S259>' : 'windEmulatorStep4/sidCtrl/PID Controller/Saturation'
 * '<S260>' : 'windEmulatorStep4/sidCtrl/PID Controller/Saturation Fdbk'
 * '<S261>' : 'windEmulatorStep4/sidCtrl/PID Controller/Sum'
 * '<S262>' : 'windEmulatorStep4/sidCtrl/PID Controller/Sum Fdbk'
 * '<S263>' : 'windEmulatorStep4/sidCtrl/PID Controller/Tracking Mode'
 * '<S264>' : 'windEmulatorStep4/sidCtrl/PID Controller/Tracking Mode Sum'
 * '<S265>' : 'windEmulatorStep4/sidCtrl/PID Controller/Tsamp - Integral'
 * '<S266>' : 'windEmulatorStep4/sidCtrl/PID Controller/Tsamp - Ngain'
 * '<S267>' : 'windEmulatorStep4/sidCtrl/PID Controller/postSat Signal'
 * '<S268>' : 'windEmulatorStep4/sidCtrl/PID Controller/preInt Signal'
 * '<S269>' : 'windEmulatorStep4/sidCtrl/PID Controller/preSat Signal'
 * '<S270>' : 'windEmulatorStep4/sidCtrl/PID Controller/Anti-windup/Passthrough'
 * '<S271>' : 'windEmulatorStep4/sidCtrl/PID Controller/D Gain/Disabled'
 * '<S272>' : 'windEmulatorStep4/sidCtrl/PID Controller/External Derivative/Disabled'
 * '<S273>' : 'windEmulatorStep4/sidCtrl/PID Controller/Filter/Disabled'
 * '<S274>' : 'windEmulatorStep4/sidCtrl/PID Controller/Filter ICs/Disabled'
 * '<S275>' : 'windEmulatorStep4/sidCtrl/PID Controller/I Gain/External Parameters'
 * '<S276>' : 'windEmulatorStep4/sidCtrl/PID Controller/Ideal P Gain/Passthrough'
 * '<S277>' : 'windEmulatorStep4/sidCtrl/PID Controller/Ideal P Gain Fdbk/Disabled'
 * '<S278>' : 'windEmulatorStep4/sidCtrl/PID Controller/Integrator/Continuous'
 * '<S279>' : 'windEmulatorStep4/sidCtrl/PID Controller/Integrator ICs/Internal IC'
 * '<S280>' : 'windEmulatorStep4/sidCtrl/PID Controller/N Copy/Disabled wSignal Specification'
 * '<S281>' : 'windEmulatorStep4/sidCtrl/PID Controller/N Gain/Disabled'
 * '<S282>' : 'windEmulatorStep4/sidCtrl/PID Controller/P Copy/Disabled'
 * '<S283>' : 'windEmulatorStep4/sidCtrl/PID Controller/Parallel P Gain/External Parameters'
 * '<S284>' : 'windEmulatorStep4/sidCtrl/PID Controller/Reset Signal/External Reset'
 * '<S285>' : 'windEmulatorStep4/sidCtrl/PID Controller/Saturation/External'
 * '<S286>' : 'windEmulatorStep4/sidCtrl/PID Controller/Saturation/External/Saturation Dynamic'
 * '<S287>' : 'windEmulatorStep4/sidCtrl/PID Controller/Saturation Fdbk/Disabled'
 * '<S288>' : 'windEmulatorStep4/sidCtrl/PID Controller/Sum/Sum_PI'
 * '<S289>' : 'windEmulatorStep4/sidCtrl/PID Controller/Sum Fdbk/Disabled'
 * '<S290>' : 'windEmulatorStep4/sidCtrl/PID Controller/Tracking Mode/Disabled'
 * '<S291>' : 'windEmulatorStep4/sidCtrl/PID Controller/Tracking Mode Sum/Passthrough'
 * '<S292>' : 'windEmulatorStep4/sidCtrl/PID Controller/Tsamp - Integral/TsSignalSpecification'
 * '<S293>' : 'windEmulatorStep4/sidCtrl/PID Controller/Tsamp - Ngain/Passthrough'
 * '<S294>' : 'windEmulatorStep4/sidCtrl/PID Controller/postSat Signal/Forward_Path'
 * '<S295>' : 'windEmulatorStep4/sidCtrl/PID Controller/preInt Signal/Internal PreInt'
 * '<S296>' : 'windEmulatorStep4/sidCtrl/PID Controller/preSat Signal/Forward_Path'
 * '<S297>' : 'windEmulatorStep4/sidCtrl/setpointGenerator/fromFile'
 * '<S298>' : 'windEmulatorStep4/writeAcs800Pdos/ACS800ConvertREF1'
 * '<S299>' : 'windEmulatorStep4/writeAcs800Pdos/ACS800ConvertREF2'
 * '<S300>' : 'windEmulatorStep4/writeAcs880Pdos/ACS880ConvertREF2'
 */
#endif                                 /* windEmulatorStep4_h_ */
