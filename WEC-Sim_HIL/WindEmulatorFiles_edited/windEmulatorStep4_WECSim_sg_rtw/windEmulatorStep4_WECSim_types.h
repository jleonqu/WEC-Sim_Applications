/*
 * windEmulatorStep4_WECSim_types.h
 *
 * Code generation for model "windEmulatorStep4_WECSim".
 *
 * Model version              : 10.5
 * Simulink Coder version : 25.2 (R2025b) 28-Jul-2025
 * C++ source code generated on : Mon Jul  6 11:39:01 2026
 *
 * Target selection: speedgoat.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef windEmulatorStep4_WECSim_types_h_
#define windEmulatorStep4_WECSim_types_h_
#include "rtwtypes.h"
#include "expType.h"
#include "sidType.h"
#include "abbState.h"
#ifndef DEFINED_TYPEDEF_FOR_bus_body2_hydroForce_hf_fExt_
#define DEFINED_TYPEDEF_FOR_bus_body2_hydroForce_hf_fExt_

struct bus_body2_hydroForce_hf_fExt
{
  real_T re[6];
  real_T im[6];
  real_T md[6];
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_bus_body2_hydroForce_hf_storage_
#define DEFINED_TYPEDEF_FOR_bus_body2_hydroForce_hf_storage_

struct bus_body2_hydroForce_hf_storage
{
  real_T fAddedMass[36];
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_bus_body2_hydroForce_hf_
#define DEFINED_TYPEDEF_FOR_bus_body2_hydroForce_hf_

struct bus_body2_hydroForce_hf
{
  real_T linearHydroRestCoef[36];
  real_T userDefinedFe;
  real_T quadDrag[36];
  real_T linearDamping[36];
  real_T volume;
  real_T centerBuoyancy[3];
  real_T centerGravity[3];
  bus_body2_hydroForce_hf_fExt fExt;
  real_T fAddedMass[36];
  real_T fDamping[36];
  real_T totDOF[36];
  bus_body2_hydroForce_hf_storage storage;
  real_T mass;
  real_T adjustedMass;
  real_T adjustedInertia[3];
  real_T adjustedInertiaProducts[3];
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_bus_body2_hydroForce_
#define DEFINED_TYPEDEF_FOR_bus_body2_hydroForce_

struct bus_body2_hydroForce
{
  bus_body2_hydroForce_hf hf1;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_bus_body1_hydroForce_hf_fExt_
#define DEFINED_TYPEDEF_FOR_bus_body1_hydroForce_hf_fExt_

struct bus_body1_hydroForce_hf_fExt
{
  real_T re[6];
  real_T im[6];
  real_T md[6];
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_bus_body1_hydroForce_hf_storage_
#define DEFINED_TYPEDEF_FOR_bus_body1_hydroForce_hf_storage_

struct bus_body1_hydroForce_hf_storage
{
  real_T fAddedMass[36];
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_bus_body1_hydroForce_hf_
#define DEFINED_TYPEDEF_FOR_bus_body1_hydroForce_hf_

struct bus_body1_hydroForce_hf
{
  real_T linearHydroRestCoef[36];
  real_T userDefinedFe;
  real_T quadDrag[36];
  real_T linearDamping[36];
  real_T volume;
  real_T centerBuoyancy[3];
  real_T centerGravity[3];
  bus_body1_hydroForce_hf_fExt fExt;
  real_T fAddedMass[36];
  real_T fDamping[36];
  real_T totDOF[36];
  bus_body1_hydroForce_hf_storage storage;
  real_T mass;
  real_T adjustedMass;
  real_T adjustedInertia[3];
  real_T adjustedInertiaProducts[3];
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_bus_body1_hydroForce_
#define DEFINED_TYPEDEF_FOR_bus_body1_hydroForce_

struct bus_body1_hydroForce
{
  bus_body1_hydroForce_hf hf1;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_acs800CtrlBus_
#define DEFINED_TYPEDEF_FOR_acs800CtrlBus_

struct acs800CtrlBus
{
  uint16_T ctrlWord;
  int32_T state;
  real_T torqueSetpoint_Nm;
  real_T torqueSetpoint_percent;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_acs880CtrlBus_
#define DEFINED_TYPEDEF_FOR_acs880CtrlBus_

struct acs880CtrlBus
{
  uint16_T ctrlWord;
  int32_T state;
  real_T torqueSetpoint_Nm;
  real_T torqueSetpoint_percent;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_expCtrlBus_
#define DEFINED_TYPEDEF_FOR_expCtrlBus_

struct expCtrlBus
{
  uint16_T expType;
  boolean_T runHil;
  boolean_T runSid;
  boolean_T resetHilIntegrator;
  boolean_T resetSidIntegrator;
  real_T time;
  real_T ramp;
  uint32_T runCounter;
  uint32_T stepCounter;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_hptoSignalBus_
#define DEFINED_TYPEDEF_FOR_hptoSignalBus_

struct hptoSignalBus
{
  real_T genTorqueCmd_Nm;
  real_T pressure_bar;
  real_T hmOutputShafTorque_Nm;
  real_T genShaftSpeed_rpm;
  real_T excShaftSpeed_rpm;
  real_T excShaftTorque_Nm;
  real_T ctrlSignal1;
  real_T ctrlSignal2;
  real_T genPumpFlow_lpm;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_sidCtrlBus_
#define DEFINED_TYPEDEF_FOR_sidCtrlBus_

struct sidCtrlBus
{
  real_T acs800Torque_Nm;
  real_T acs880Torque_Nm;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_acs880SignalBus_
#define DEFINED_TYPEDEF_FOR_acs880SignalBus_

struct acs880SignalBus
{
  real_T motorVoltage_V;
  real_T motorCurrent_A;
  real_T frequency_Hz;
  real_T motorSpeed_rpm;
  real_T motorTorque_Nm;
  real_T shaftPower_W;
  uint16_T statusWord;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_acs800SignalBus_
#define DEFINED_TYPEDEF_FOR_acs800SignalBus_

struct acs800SignalBus
{
  real_T dcBusVoltage_V;
  real_T frequency_Hz;
  real_T temperature;
  real_T motorSpeed_rpm;
  real_T motorTorque_Nm;
  real_T shaftPower_W;
  uint16_T statusWord;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_hptoCtrlBus_
#define DEFINED_TYPEDEF_FOR_hptoCtrlBus_

struct hptoCtrlBus
{
  real_T speedRef_rpm;
  real_T excForce_N;
  real_T genSpeedActual;
  boolean_T speedCtrlReset;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_shaftSignalBus_
#define DEFINED_TYPEDEF_FOR_shaftSignalBus_

struct shaftSignalBus
{
  real_T torqueActual_Nm;
  uint32_T absEncoderCounts;
  int32_T absEncoderTurns;
  real_T absEncoderPosition_rad;
  real_T absEncoderSpeed_rpm;
  uint8_T absEncoderStatus1;
  uint8_T absEncoderStatus2;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_invPowerBus_
#define DEFINED_TYPEDEF_FOR_invPowerBus_

struct invPowerBus
{
  boolean_T L1InaccurateU;
  boolean_T L1InaccurateI;
  real_T L1Voltage;
  real_T L1Current;
  real_T L1PowFactor;
  real_T L1ActivePow;
  real_T L1THDu;
  real_T L1THDi;
  boolean_T L2InaccurateU;
  boolean_T L2InaccurateI;
  real_T L2Voltage;
  real_T L2Current;
  real_T L2PowFactor;
  real_T L2ActivePow;
  real_T L2THDu;
  real_T L2THDi;
  boolean_T L3InaccurateU;
  boolean_T L3InaccurateI;
  real_T L3Voltage;
  real_T L3Current;
  real_T L3PowFactor;
  real_T L3ActivePow;
  real_T L3THDu;
  real_T L3THDi;
  real_T Frequency;
  real_T TotalPowFactor;
  real_T TotalActivePow;
  real_T L1L2Voltage;
  real_T L2L3Voltage;
  real_T L3L1Volage;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_sidInfoBus_
#define DEFINED_TYPEDEF_FOR_sidInfoBus_

struct sidInfoBus
{
  real_T acs800TorqueSetpoint_Nm;
  real_T acs880SpeedSetpoint_rpm;
  uint32_T sidType;
  real_T fromFileCaseCounter;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_KY1U3Kyrwv5e6VnUIBWG5G_
#define DEFINED_TYPEDEF_FOR_struct_KY1U3Kyrwv5e6VnUIBWG5G_

struct struct_KY1U3Kyrwv5e6VnUIBWG5G
{
  real_T PG;
  real_T IG;
};

#endif

#ifndef struct_cell_wrap_windEmulatorStep4_W_T
#define struct_cell_wrap_windEmulatorStep4_W_T

struct cell_wrap_windEmulatorStep4_W_T
{
  uint32_T f1[8];
};

#endif                              /* struct_cell_wrap_windEmulatorStep4_W_T */

#ifndef struct_dsp_simulink_MovingAverage_wi_T
#define struct_dsp_simulink_MovingAverage_wi_T

struct dsp_simulink_MovingAverage_wi_T
{
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
  boolean_T TunablePropsChanged;
  cell_wrap_windEmulatorStep4_W_T inputVarSize;
  int32_T NumChannels;
  int32_T FrameLength;
  real_T pCumSum;
  real_T pCumSumRev[2499];
  real_T pCumRevIndex;
  real_T pModValueRev;
};

#endif                              /* struct_dsp_simulink_MovingAverage_wi_T */

/* Forward declaration for rtModel */
typedef struct tag_RTM_windEmulatorStep4_WECSim_T
  RT_MODEL_windEmulatorStep4_WECSim_T;

#endif                                 /* windEmulatorStep4_WECSim_types_h_ */
