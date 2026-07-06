/*
 * windEmulatorStep4_WECSim_private.h
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

#ifndef windEmulatorStep4_WECSim_private_h_
#define windEmulatorStep4_WECSim_private_h_
#include "rtwtypes.h"
#include "multiword_types.h"
#include "zero_crossing_types.h"
#include "windEmulatorStep4_WECSim.h"
#include "windEmulatorStep4_WECSim_cal.h"
#include "windEmulatorStep4_WECSim_types.h"

/* Private macros used by the generated code to access rtModel */
#ifndef rtmIsMajorTimeStep
#define rtmIsMajorTimeStep(rtm)        (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
#define rtmIsMinorTimeStep(rtm)        (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmSetTFinal
#define rtmSetTFinal(rtm, val)         ((rtm)->Timing.tFinal = (val))
#endif

#ifndef rtmSetTPtr
#define rtmSetTPtr(rtm, val)           ((rtm)->Timing.t = (val))
#endif

#ifndef UCHAR_MAX
#include <limits.h>
#endif

#if ( UCHAR_MAX != (0xFFU) ) || ( SCHAR_MAX != (0x7F) )
#error Code was generated for compiler with different sized uchar/char. \
Consider adjusting Test hardware word size settings on the \
Hardware Implementation pane to match your compiler word sizes as \
defined in limits.h of the compiler. Alternatively, you can \
select the Test hardware is the same as production hardware option and \
select the Enable portable word sizes option on the Code Generation > \
Verification pane for ERT based targets, which will disable the \
preprocessor word size checks.
#endif

#if ( USHRT_MAX != (0xFFFFU) ) || ( SHRT_MAX != (0x7FFF) )
#error Code was generated for compiler with different sized ushort/short. \
Consider adjusting Test hardware word size settings on the \
Hardware Implementation pane to match your compiler word sizes as \
defined in limits.h of the compiler. Alternatively, you can \
select the Test hardware is the same as production hardware option and \
select the Enable portable word sizes option on the Code Generation > \
Verification pane for ERT based targets, which will disable the \
preprocessor word size checks.
#endif

#if ( UINT_MAX != (0xFFFFFFFFU) ) || ( INT_MAX != (0x7FFFFFFF) )
#error Code was generated for compiler with different sized uint/int. \
Consider adjusting Test hardware word size settings on the \
Hardware Implementation pane to match your compiler word sizes as \
defined in limits.h of the compiler. Alternatively, you can \
select the Test hardware is the same as production hardware option and \
select the Enable portable word sizes option on the Code Generation > \
Verification pane for ERT based targets, which will disable the \
preprocessor word size checks.
#endif

/* Skipping ulong/long check: insufficient preprocessor integer range. */

/* Skipping ulong_long/long_long check: insufficient preprocessor integer range. */
extern unsigned int xmlecatArr_0_count;
extern unsigned char xmlecatArr_0[];
extern int_T slrtEcatDCM[8];           // From master shift controller
namespace slrealtime
{
  namespace tracing
  {
    struct IamRoot;
  }
}

extern real_T rt_atan2d_snf(real_T u0, real_T u1);
extern void* slrtRegisterSignalToLoggingService(uintptr_t sigAddr);
real_T rt_TDelayInterpolate(
  real_T tMinusDelay,                 /* tMinusDelay = currentSimTime - delay */
  real_T tStart,
  real_T *uBuf,
  int_T bufSz,
  int_T *lastIdx,
  int_T oldestIdx,
  int_T newIdx,
  real_T initOutput,
  boolean_T discrete,
  boolean_T minorStepAndTAtLastMajorOutput)
  ;
extern "C" void slrealtimeenablelogging(SimStruct *rts);
extern void windEmulator_parseCtrlWord_Init(DW_parseCtrlWord_windEmulator_T
  *localDW);
extern void windEmulatorStep4_parseCtrlWord(uint16_T rtu_ctrlWord,
  B_parseCtrlWord_windEmulatorS_T *localB, DW_parseCtrlWord_windEmulator_T
  *localDW);
extern void windEmulator_MovingAverage_Init(DW_MovingAverage_windEmulator_T
  *localDW);
extern void windEmulatorStep4_MovingAverage(real_T rtu_0,
  B_MovingAverage_windEmulatorS_T *localB, DW_MovingAverage_windEmulator_T
  *localDW);
extern void windEmulat_ParseStatusWord_Init(DW_ParseStatusWord_windEmulat_T
  *localDW);
extern void windEmulatorSte_ParseStatusWord(uint16_T rtu_StatusWord,
  B_ParseStatusWord_windEmulato_T *localB, DW_ParseStatusWord_windEmulat_T
  *localDW);
extern void windEmulat_MovingAverage_e_Init(DW_MovingAverage_windEmulat_f_T
  *localDW);
extern void windEmulatorSte_MovingAverage_p(real_T rtu_0,
  B_MovingAverage_windEmulato_c_T *localB, DW_MovingAverage_windEmulat_f_T
  *localDW);
extern void windEmul_NonlinearWaveElevation(const real_T rtu_Displacement[3],
  const real_T rtu_Displacement_h[3], B_NonlinearWaveElevation_wind_T *localB,
  NonlinearWaveElevation_cal_type *windEmulatorS_PageSwitching_arg);
extern void windEmul_quaternion2EulXYZ_Init(DW_quaternion2EulXYZ_windEmul_T
  *localDW);
extern void windEmulatorS_quaternion2EulXYZ(const real_T rtu_Q[4],
  B_quaternion2EulXYZ_windEmula_T *localB, DW_quaternion2EulXYZ_windEmul_T
  *localDW);
extern void windEmulat_MATLABFunction1_Init(DW_MATLABFunction1_windEmulat_T
  *localDW);
extern void windEmulatorSte_MATLABFunction1(const real_T rtu_disp[2], real_T
  rtu_enable, real_T rtu_wavenumber, real_T rtu_direction,
  B_MATLABFunction1_windEmulato_T *localB, DW_MATLABFunction1_windEmulat_T
  *localDW);
extern void windEmu_YawForceTransforms_Init(DW_YawForceTransforms_windEmu_T
  *localDW);
extern void windEmulator_YawForceTransforms(real_T rtu_yaw, const real_T
  rtu_dispGlobal[3], const real_T rtu_dispGlobal_e[3], const real_T
  rtu_F_RadiationDampingLocal[6], const real_T rtu_F_AddedMassLocal[6], const
  real_T rtu_F_ExcitationLocal[6], const real_T rtu_F_RestoringLocal[6],
  B_YawForceTransforms_windEmul_T *localB, DW_YawForceTransforms_windEmu_T
  *localDW);
extern void win_YawKinematicTransforms_Init(DW_YawKinematicTransforms_win_T
  *localDW);
extern void windEmul_YawKinematicTransforms(real_T rtu_yaw, const real_T
  rtu_dispGlobal[3], const real_T rtu_dispGlobal_k[3], const real_T
  rtu_velGlobal[3], const real_T rtu_velGlobal_i[3], const real_T rtu_accGlobal
  [6], B_YawKinematicTransforms_wind_T *localB, DW_YawKinematicTransforms_win_T *
  localDW);
extern void windEmulator_MovingAverage_Term(DW_MovingAverage_windEmulator_T
  *localDW);
extern void windEmulat_MovingAverage_g_Term(DW_MovingAverage_windEmulat_f_T
  *localDW);

/* private model entry point functions */
extern void windEmulatorStep4_WECSim_derivatives(void);
extern void windEmulatorStep4_WECSim_forcingfunction(void);
extern void windEmulatorStep4_WECSim_massmatrix(void);

#endif                                 /* windEmulatorStep4_WECSim_private_h_ */
