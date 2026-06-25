/*
 * windEmulatorStep4_private.h
 *
 * Code generation for model "windEmulatorStep4".
 *
 * Model version              : 4.1
 * Simulink Coder version : 9.8 (R2022b) 13-May-2022
 * C++ source code generated on : Thu Jun 25 12:29:59 2026
 *
 * Target selection: slrealtime.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_windEmulatorStep4_private_h_
#define RTW_HEADER_windEmulatorStep4_private_h_
#include "rtwtypes.h"
#include "multiword_types.h"
#include "zero_crossing_types.h"
#include "windEmulatorStep4.h"
#include "windEmulatorStep4_types.h"

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

extern void* slrtRegisterSignalToLoggingService(uintptr_t sigAddr);
extern "C" void slrealtimeenablelogging(SimStruct *rts);
extern void windEmulatorStep4_parseCtrlWord(uint16_T rtu_ctrlWord,
  B_parseCtrlWord_windEmulatorS_T *localB);
extern void windEmulator_MovingAverage_Init(DW_MovingAverage_windEmulator_T
  *localDW);
extern void windEmulatorStep4_MovingAverage(real_T rtu_0,
  B_MovingAverage_windEmulatorS_T *localB, DW_MovingAverage_windEmulator_T
  *localDW);
extern void windEmulatorSte_ParseStatusWord(uint16_T rtu_StatusWord,
  B_ParseStatusWord_windEmulato_T *localB);
extern void windEmulator_MovingAverage_Term(DW_MovingAverage_windEmulator_T
  *localDW);

/* private model entry point functions */
extern void windEmulatorStep4_derivatives(void);

#endif                             /* RTW_HEADER_windEmulatorStep4_private_h_ */
