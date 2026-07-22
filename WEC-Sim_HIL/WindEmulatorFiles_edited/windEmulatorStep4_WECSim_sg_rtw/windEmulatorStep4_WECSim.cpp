/*
 * windEmulatorStep4_WECSim.cpp
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

#include "windEmulatorStep4_WECSim.h"
#include "rtwtypes.h"
#include "windEmulatorStep4_WECSim_types.h"
#include "windEmulatorStep4_WECSim_cal.h"
#include "windEmulatorStep4_WECSim_private.h"
#include <cstring>
#include <cmath>

extern "C"
{

#include "rt_nonfinite.h"

}

#include "rte_windEmulatorStep4_WECSim_parameters.h"
#include "abbState.h"
#include "sidType.h"
#include <stddef.h>
#include "expType.h"
#include "zero_crossing_types.h"
#include "rt_defines.h"

/* Named constants for MATLAB Function: '<S41>/parseCtrlWord' */
const int32_T windEmulatorStep4_WE_CALL_EVENT = -1;

/* Named constants for MATLAB Function: '<S44>/Parse Status Word' */
const int32_T windEmulatorStep4__CALL_EVENT_n = -1;

/* Named constants for MATLAB Function: '<S81>/quaternion2EulXYZ' */
const int32_T windEmulatorStep4__CALL_EVENT_g = -1;

/* Named constants for MATLAB Function: '<S126>/MATLAB Function1' */
const int32_T windEmulatorStep4__CALL_EVENT_b = -1;

/* Named constants for MATLAB Function: '<S70>/Yaw Force Transforms' */
const int32_T windEmulatorStep4__CALL_EVENT_c = -1;

/* Named constants for MATLAB Function: '<S133>/Yaw Kinematic Transforms' */
const int32_T windEmulatorStep4__CALL_EVENT_d = -1;

/* Named constants for Chart: '<S16>/ABB Fieldbus Control' */
const uint32_T windEmula_IN_notReadyToSwitchOn = 3U;
const uint32_T windEmulat_IN_operationDisabled = 4U;
const uint32_T windEmulato_IN_operationEnabled = 5U;
const uint32_T windEmulatorStep4_IN_initialize = 2U;
const int32_T windEmulatorStep4__CALL_EVENT_k = -1;
const uint32_T windEmulatorStep4__IN_DelayOFF1 = 1U;
const uint8_T windEmulator_IN_NO_ACTIVE_CHILD = 0U;
const uint32_T windEmulator_IN_readyToSwitchOn = 6U;

/* Named constants for Chart: '<S4>/FexcRamp' */
const uint32_T windEmulatorStep4_WECSi_IN_idle = 1U;
const uint32_T windEmulatorStep4_WECSi_IN_init = 2U;
const uint32_T windEmulatorStep4_WEC_IN_rampup = 4U;
const uint32_T windEmulatorStep4_WEC_IN_runing = 5U;
const uint32_T windEmulatorStep4_W_IN_rampdown = 3U;
const real_T windEmulatorStep4_WECSim_period = 0.004;

/* Block signals (default storage) */
B_windEmulatorStep4_WECSim_T windEmulatorStep4_WECSim_B;

/* Continuous states */
X_windEmulatorStep4_WECSim_T windEmulatorStep4_WECSim_X;

/* Disabled State Vector */
XDis_windEmulatorStep4_WECSim_T windEmulatorStep4_WECSim_XDis;

/* Block states (default storage) */
DW_windEmulatorStep4_WECSim_T windEmulatorStep4_WECSim_DW;

/* Previous zero-crossings (trigger) states */
PrevZCX_windEmulatorStep4_WECSim_T windEmulatorStep4_WECSim_PrevZCX;

/* External inputs (root inport signals with default storage) */
ExtU_windEmulatorStep4_WECSim_T windEmulatorStep4_WECSim_U;

/* Mass Matrices */
MassMatrix_windEmulatorStep4_WECSim_T windEmulatorStep4_WECSim_MassMatrix;

/* Real-time model */
RT_MODEL_windEmulatorStep4_WECSim_T windEmulatorStep4_WECSim_M_ =
  RT_MODEL_windEmulatorStep4_WECSim_T();
RT_MODEL_windEmulatorStep4_WECSim_T *const windEmulatorStep4_WECSim_M =
  &windEmulatorStep4_WECSim_M_;

/* Forward declaration for local functions */
static void windEmulatorSt_SystemCore_setup(dsp_simulink_MovingAverage_wi_T *obj);

/* Forward declaration for local functions */
static void windEmulator_SystemCore_setup_n(dsp_simulink_MovingAverage_wi_T *obj);

/* Forward declaration for local functions */
static void windEmulato_swParseStatusWord_m(void);
static void windEmulat_cwBuildControlWord_o(void);
static void windEmulatorStep_cwInitialize_a(void);
static void windEmulatorS_swParseStatusWord(void);
static void windEmulator_cwBuildControlWord(void);
static void windEmulatorStep4__cwInitialize(void);
void Root_EtherCATInit_callback(void * const ptr_rtm )
{
  int_T status = 1;
  static const char_T *errMsg;
  int_T j;
  static char_T msg[256];
  std::string logfile( "" );           //logFileName );
  std::string DeviceType( "I8254x" );
  mwStateClear( 0 );
  LOG(info,0) << "EtherCAT going to state 8";
  status = slrtEcatInit(0,
                        DeviceType.c_str(),
                        "1",
                        1,
                        (unsigned char *)xmlecatArr_0,
                        xmlecatArr_0_count,
                        1,
                        0,
                        logfile.c_str(),
                        0.004,
                        2,
                        8 );
  if (status != XPC_ECAT_OK) {
    if ((((uint32_T)status >> 16) & 0xffff) == 0xffff ) {
      // Our error conditions, negative numbers.
      switch ( status )
      {
       case -10:        // very rare, sg_getEthercatInterface can't be executed!
        errMsg =
          "Speedgoat library files for EtherCAT port identification are not properly installed on the target";
        break;

       case -11:          // rare, sg_getEthercatInterface didn't create eciface
        errMsg = "Ethernet port mapping failed";// eciface didn't get created
        break;

       case -12:
              // common, User entered a 1 based port number that isn't reserved.
        errMsg = "EtherCAT port 1 is not reserved for EtherCAT";
        break;
      }
    } else {
      if ((uint32_T)status == 0x9811000C )
        errMsg =
          "Network port 1 is not accessible to EtherCAT.\nIt is either non-existant or not configured for EtherCAT.";
      else
        errMsg = xpcPrintEtherCATError(0, 0);
    }

    rtmSetErrorStatus(windEmulatorStep4_WECSim_M, errMsg);
    return;
  }
}

/*
 * Time delay interpolation routine
 *
 * The linear interpolation is performed using the formula:
 *
 * (t2 - tMinusDelay)         (tMinusDelay - t1)
 * u(t)  =  ----------------- * u1  +  ------------------- * u2
 * (t2 - t1)                  (t2 - t1)
 */
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
{
  int_T i;
  real_T yout, t1, t2, u1, u2;
  real_T* tBuf = uBuf + bufSz;

  /*
   * If there is only one data point in the buffer, this data point must be
   * the t= 0 and tMinusDelay > t0, it ask for something unknown. The best
   * guess if initial output as well
   */
  if ((newIdx == 0) && (oldestIdx ==0 ) && (tMinusDelay > tStart))
    return initOutput;

  /*
   * If tMinusDelay is less than zero, should output initial value
   */
  if (tMinusDelay <= tStart)
    return initOutput;

  /* For fixed buffer extrapolation:
   * if tMinusDelay is small than the time at oldestIdx, if discrete, output
   * tailptr value,  else use tailptr and tailptr+1 value to extrapolate
   * It is also for fixed buffer. Note: The same condition can happen for transport delay block where
   * use tStart and and t[tail] other than using t[tail] and t[tail+1].
   * See below
   */
  if ((tMinusDelay <= tBuf[oldestIdx] ) ) {
    if (discrete) {
      return(uBuf[oldestIdx]);
    } else {
      int_T tempIdx= oldestIdx + 1;
      if (oldestIdx == bufSz-1)
        tempIdx = 0;
      t1= tBuf[oldestIdx];
      t2= tBuf[tempIdx];
      u1= uBuf[oldestIdx];
      u2= uBuf[tempIdx];
      if (t2 == t1) {
        if (tMinusDelay >= t2) {
          yout = u2;
        } else {
          yout = u1;
        }
      } else {
        real_T f1 = (t2-tMinusDelay) / (t2-t1);
        real_T f2 = 1.0 - f1;

        /*
         * Use Lagrange's interpolation formula.  Exact outputs at t1, t2.
         */
        yout = f1*u1 + f2*u2;
      }

      return yout;
    }
  }

  /*
   * When block does not have direct feedthrough, we use the table of
   * values to extrapolate off the end of the table for delays that are less
   * than 0 (less then step size).  This is not completely accurate.  The
   * chain of events is as follows for a given time t.  Major output - look
   * in table.  Update - add entry to table.  Now, if we call the output at
   * time t again, there is a new entry in the table. For very small delays,
   * this means that we will have a different answer from the previous call
   * to the output fcn at the same time t.  The following code prevents this
   * from happening.
   */
  if (minorStepAndTAtLastMajorOutput) {
    /* pretend that the new entry has not been added to table */
    if (newIdx != 0) {
      if (*lastIdx == newIdx) {
        (*lastIdx)--;
      }

      newIdx--;
    } else {
      if (*lastIdx == newIdx) {
        *lastIdx = bufSz-1;
      }

      newIdx = bufSz - 1;
    }
  }

  i = *lastIdx;
  if (tBuf[i] < tMinusDelay) {
    /* Look forward starting at last index */
    while (tBuf[i] < tMinusDelay) {
      /* May occur if the delay is less than step-size - extrapolate */
      if (i == newIdx)
        break;
      i = ( i < (bufSz-1) ) ? (i+1) : 0;/* move through buffer */
    }
  } else {
    /*
     * Look backwards starting at last index which can happen when the
     * delay time increases.
     */
    while (tBuf[i] >= tMinusDelay) {
      /*
       * Due to the entry condition at top of function, we
       * should never hit the end.
       */
      i = (i > 0) ? i-1 : (bufSz-1);   /* move through buffer */
    }

    i = ( i < (bufSz-1) ) ? (i+1) : 0;
  }

  *lastIdx = i;
  if (discrete) {
    /*
     * tempEps = 128 * eps;
     * localEps = max(tempEps, tempEps*fabs(tBuf[i]))/2;
     */
    double tempEps = (DBL_EPSILON) * 128.0;
    double localEps = tempEps * std::abs(tBuf[i]);
    if (tempEps > localEps) {
      localEps = tempEps;
    }

    localEps = localEps / 2.0;
    if (tMinusDelay >= (tBuf[i] - localEps)) {
      yout = uBuf[i];
    } else {
      if (i == 0) {
        yout = uBuf[bufSz-1];
      } else {
        yout = uBuf[i-1];
      }
    }
  } else {
    if (i == 0) {
      t1 = tBuf[bufSz-1];
      u1 = uBuf[bufSz-1];
    } else {
      t1 = tBuf[i-1];
      u1 = uBuf[i-1];
    }

    t2 = tBuf[i];
    u2 = uBuf[i];
    if (t2 == t1) {
      if (tMinusDelay >= t2) {
        yout = u2;
      } else {
        yout = u1;
      }
    } else {
      real_T f1 = (t2-tMinusDelay) / (t2-t1);
      real_T f2 = 1.0 - f1;

      /*
       * Use Lagrange's interpolation formula.  Exact outputs at t1, t2.
       */
      yout = f1*u1 + f2*u2;
    }
  }

  return(yout);
}

/* Projection for root system: '<Root>' */
void windEmulatorStep4_WECSim_projection(void)
{
  NeslSimulationData *simulationData;
  NeslSimulator *simulator;
  NeuDiagnosticManager *diagnosticManager;
  NeuDiagnosticTree *diagnosticTree;
  char *msg;
  real_T tmp_1[52];
  real_T tmp_5[12];
  real_T time;
  real_T time_0;
  real_T time_tmp;
  int32_T tmp_3;
  int_T tmp_2[14];
  int_T tmp_6[4];
  boolean_T tmp;
  boolean_T tmp_0;
  boolean_T tmp_4;

  /* Projection for SimscapeExecutionBlock: '<S216>/STATE_1' incorporates:
   *  SimscapeExecutionBlock: '<S332>/STATE_1'
   */
  simulationData = static_cast<NeslSimulationData *>
    (windEmulatorStep4_WECSim_DW.STATE_1_SimData_h);
  time_tmp = windEmulatorStep4_WECSim_M->Timing.t[0];
  time = time_tmp;
  simulationData->mData->mTime.mN = 1;
  simulationData->mData->mTime.mX = &time;
  simulationData->mData->mContStates.mN = 2;
  simulationData->mData->mContStates.mX =
    &windEmulatorStep4_WECSim_X.windEmulatorStep4_WECSimhptoSim[0];
  simulationData->mData->mDiscStates.mN = 0;
  simulationData->mData->mDiscStates.mX =
    &windEmulatorStep4_WECSim_DW.STATE_1_Discrete;
  simulationData->mData->mModeVector.mN = 0;
  simulationData->mData->mModeVector.mX =
    &windEmulatorStep4_WECSim_DW.STATE_1_Modes_j;
  tmp = false;
  simulationData->mData->mFoundZcEvents = tmp;
  simulationData->mData->mHadEvents = false;
  tmp = rtmIsMajorTimeStep(windEmulatorStep4_WECSim_M);
  simulationData->mData->mIsMajorTimeStep = tmp;
  tmp_0 = false;
  simulationData->mData->mIsSolverAssertCheck = tmp_0;
  simulationData->mData->mIsSolverCheckingCIC = false;
  tmp_0 = rtsiIsSolverComputingJacobian(&windEmulatorStep4_WECSim_M->solverInfo);
  simulationData->mData->mIsComputingJacobian = tmp_0;
  simulationData->mData->mIsEvaluatingF0 = false;
  simulationData->mData->mIsSolverRequestingReset = false;
  tmp_0 = rtsiIsModeUpdateTimeStep(&windEmulatorStep4_WECSim_M->solverInfo);
  simulationData->mData->mIsModeUpdateTimeStep = tmp_0;
  tmp_2[0] = 0;
  tmp_1[0] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[0];
  tmp_1[1] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[1];
  tmp_1[2] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[2];
  tmp_1[3] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[3];
  tmp_2[1] = 4;
  tmp_1[4] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[0];
  tmp_1[5] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[1];
  tmp_1[6] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[2];
  tmp_1[7] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[3];
  tmp_2[2] = 8;
  tmp_1[8] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[0];
  tmp_1[9] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[1];
  tmp_1[10] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[2];
  tmp_1[11] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[3];
  tmp_2[3] = 12;
  tmp_1[12] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[0];
  tmp_1[13] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[1];
  tmp_1[14] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[2];
  tmp_1[15] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[3];
  tmp_2[4] = 16;
  tmp_1[16] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[0];
  tmp_1[17] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[1];
  tmp_1[18] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[2];
  tmp_1[19] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[3];
  tmp_2[5] = 20;
  tmp_1[20] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[0];
  tmp_1[21] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[1];
  tmp_1[22] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[2];
  tmp_1[23] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[3];
  tmp_2[6] = 24;
  tmp_1[24] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[0];
  tmp_1[25] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[1];
  tmp_1[26] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[2];
  tmp_1[27] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[3];
  tmp_2[7] = 28;
  tmp_1[28] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[0];
  tmp_1[29] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[1];
  tmp_1[30] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[2];
  tmp_1[31] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[3];
  tmp_2[8] = 32;
  tmp_1[32] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[0];
  tmp_1[33] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[1];
  tmp_1[34] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[2];
  tmp_1[35] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[3];
  tmp_2[9] = 36;
  tmp_1[36] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[0];
  tmp_1[37] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[1];
  tmp_1[38] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[2];
  tmp_1[39] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[3];
  tmp_2[10] = 40;
  tmp_1[40] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[0];
  tmp_1[41] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[1];
  tmp_1[42] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[2];
  tmp_1[43] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[3];
  tmp_2[11] = 44;
  tmp_1[44] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[0];
  tmp_1[45] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[1];
  tmp_1[46] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[2];
  tmp_1[47] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[3];
  tmp_2[12] = 48;
  tmp_1[48] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[0];
  tmp_1[49] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[1];
  tmp_1[50] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[2];
  tmp_1[51] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[3];
  tmp_2[13] = 52;
  simulationData->mData->mInputValues.mN = 52;
  simulationData->mData->mInputValues.mX = &tmp_1[0];
  simulationData->mData->mInputOffsets.mN = 14;
  simulationData->mData->mInputOffsets.mX = &tmp_2[0];
  simulator = static_cast<NeslSimulator *>
    (windEmulatorStep4_WECSim_DW.STATE_1_Simulator_f);
  diagnosticManager = static_cast<NeuDiagnosticManager *>
    (windEmulatorStep4_WECSim_DW.STATE_1_DiagMgr_o);
  diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
  tmp_3 = ne_simulator_method(simulator, NESL_SIM_PROJECTION, simulationData,
    diagnosticManager);
  if (tmp_3 != 0) {
    tmp_4 = error_buffer_is_empty(rtmGetErrorStatus(windEmulatorStep4_WECSim_M));
    if (tmp_4) {
      msg = rtw_diagnostics_msg(diagnosticTree);
      rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
    }
  }

  /* End of Projection for SimscapeExecutionBlock: '<S216>/STATE_1' */

  /* Projection for SimscapeExecutionBlock: '<S332>/STATE_1' */
  simulationData = static_cast<NeslSimulationData *>
    (windEmulatorStep4_WECSim_DW.STATE_1_SimData_a);
  time_0 = time_tmp;
  simulationData->mData->mTime.mN = 1;
  simulationData->mData->mTime.mX = &time_0;
  simulationData->mData->mContStates.mN = 35;
  simulationData->mData->mContStates.mX =
    &windEmulatorStep4_WECSim_X.windEmulatorStep4_WECSimhptoS_h[0];
  simulationData->mData->mDiscStates.mN = 6;
  simulationData->mData->mDiscStates.mX =
    &windEmulatorStep4_WECSim_DW.STATE_1_Discrete_208214823[0];
  simulationData->mData->mModeVector.mN = 21;
  simulationData->mData->mModeVector.mX =
    &windEmulatorStep4_WECSim_DW.STATE_1_Modes_i[0];
  tmp_4 = false;
  simulationData->mData->mFoundZcEvents = tmp_4;
  simulationData->mData->mHadEvents = false;
  simulationData->mData->mIsMajorTimeStep = tmp;
  tmp = false;
  simulationData->mData->mIsSolverAssertCheck = tmp;
  simulationData->mData->mIsSolverCheckingCIC = false;
  tmp = rtsiIsSolverComputingJacobian(&windEmulatorStep4_WECSim_M->solverInfo);
  simulationData->mData->mIsComputingJacobian = tmp;
  simulationData->mData->mIsEvaluatingF0 = false;
  simulationData->mData->mIsSolverRequestingReset = false;
  simulationData->mData->mIsModeUpdateTimeStep = tmp_0;
  tmp_6[0] = 0;
  tmp_5[0] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[0];
  tmp_5[1] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[1];
  tmp_5[2] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[2];
  tmp_5[3] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[3];
  tmp_6[1] = 4;
  tmp_5[4] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[0];
  tmp_5[5] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[1];
  tmp_5[6] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[2];
  tmp_5[7] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[3];
  tmp_6[2] = 8;
  tmp_5[8] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[0];
  tmp_5[9] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[1];
  tmp_5[10] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[2];
  tmp_5[11] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[3];
  tmp_6[3] = 12;
  simulationData->mData->mInputValues.mN = 12;
  simulationData->mData->mInputValues.mX = &tmp_5[0];
  simulationData->mData->mInputOffsets.mN = 4;
  simulationData->mData->mInputOffsets.mX = &tmp_6[0];
  simulator = static_cast<NeslSimulator *>
    (windEmulatorStep4_WECSim_DW.STATE_1_Simulator_i);
  diagnosticManager = static_cast<NeuDiagnosticManager *>
    (windEmulatorStep4_WECSim_DW.STATE_1_DiagMgr_g);
  diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
  tmp_3 = ne_simulator_method(simulator, NESL_SIM_PROJECTION, simulationData,
    diagnosticManager);
  if (tmp_3 != 0) {
    tmp = error_buffer_is_empty(rtmGetErrorStatus(windEmulatorStep4_WECSim_M));
    if (tmp) {
      msg = rtw_diagnostics_msg(diagnosticTree);
      rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
    }
  }
}

/* ForcingFunction for root system: '<Root>' */
void windEmulatorStep4_WECSim_forcingfunction(void)
{
  NeslSimulationData *simulationData;
  NeslSimulator *simulator;
  NeuDiagnosticManager *diagnosticManager;
  NeuDiagnosticTree *diagnosticTree;
  XDot_windEmulatorStep4_WECSim_T *_rtXdot;
  char *msg;
  real_T tmp_1[52];
  real_T tmp_4[12];
  real_T time;
  real_T time_0;
  real_T time_tmp;
  int_T tmp_2[14];
  int_T tmp_5[4];
  int_T is;
  uint32_T ri;
  boolean_T tmp;
  boolean_T tmp_0;
  boolean_T tmp_3;
  _rtXdot = ((XDot_windEmulatorStep4_WECSim_T *)
             windEmulatorStep4_WECSim_M->derivs);

  /* ForcingFunction for Integrator: '<S590>/Integrator' */
  if (!windEmulatorStep4_WECSim_B.BusAssignment_b.resetSidIntegrator) {
    _rtXdot->Integrator_CSTATE = windEmulatorStep4_WECSim_B.IProdOut;
  } else {
    /* level reset is active */
    _rtXdot->Integrator_CSTATE = 0.0;
  }

  /* End of ForcingFunction for Integrator: '<S590>/Integrator' */

  /* ForcingFunction for StateSpace: '<S531>/Internal' */
  _rtXdot->Internal_CSTATE[0] = 0.0;
  _rtXdot->Internal_CSTATE[1] = 0.0;
  _rtXdot->Internal_CSTATE[2] = 0.0;
  for (ri = windEmulatorStep4_WECSim_cal->Internal_A_jc[0U]; ri <
       windEmulatorStep4_WECSim_cal->Internal_A_jc[1U]; ri++) {
    _rtXdot->Internal_CSTATE[windEmulatorStep4_WECSim_cal->Internal_A_ir[ri]] +=
      windEmulatorStep4_WECSim_cal->Internal_A_pr[ri] *
      windEmulatorStep4_WECSim_X.Internal_CSTATE[0U];
  }

  for (ri = windEmulatorStep4_WECSim_cal->Internal_A_jc[1U]; ri <
       windEmulatorStep4_WECSim_cal->Internal_A_jc[2U]; ri++) {
    _rtXdot->Internal_CSTATE[windEmulatorStep4_WECSim_cal->Internal_A_ir[ri]] +=
      windEmulatorStep4_WECSim_cal->Internal_A_pr[ri] *
      windEmulatorStep4_WECSim_X.Internal_CSTATE[1U];
  }

  for (ri = windEmulatorStep4_WECSim_cal->Internal_A_jc[2U]; ri <
       windEmulatorStep4_WECSim_cal->Internal_A_jc[3U]; ri++) {
    _rtXdot->Internal_CSTATE[windEmulatorStep4_WECSim_cal->Internal_A_ir[ri]] +=
      windEmulatorStep4_WECSim_cal->Internal_A_pr[ri] *
      windEmulatorStep4_WECSim_X.Internal_CSTATE[2U];
  }

  for (ri = windEmulatorStep4_WECSim_cal->Internal_B_jc[0U]; ri <
       windEmulatorStep4_WECSim_cal->Internal_B_jc[1U]; ri++) {
    _rtXdot->Internal_CSTATE[windEmulatorStep4_WECSim_cal->Internal_B_ir] +=
      windEmulatorStep4_WECSim_cal->Internal_B_pr *
      windEmulatorStep4_WECSim_B.Sum_a;
  }

  /* End of ForcingFunction for StateSpace: '<S531>/Internal' */

  /* ForcingFunction for StateSpace: '<S545>/Internal' */
  _rtXdot->Internal_CSTATE_j = 0.0;
  for (ri = windEmulatorStep4_WECSim_cal->Internal_A_jc_i[0U]; ri <
       windEmulatorStep4_WECSim_cal->Internal_A_jc_i[1U]; ri++) {
    _rtXdot->Internal_CSTATE_j += windEmulatorStep4_WECSim_cal->Internal_A_pr_j *
      windEmulatorStep4_WECSim_X.Internal_CSTATE_j;
  }

  for (ri = windEmulatorStep4_WECSim_cal->Internal_B_jc_i[0U]; ri <
       windEmulatorStep4_WECSim_cal->Internal_B_jc_i[1U]; ri++) {
    _rtXdot->Internal_CSTATE_j += windEmulatorStep4_WECSim_cal->Internal_B_pr_g *
      windEmulatorStep4_WECSim_B.ControlSignal2;
  }

  /* End of ForcingFunction for StateSpace: '<S545>/Internal' */

  /* ForcingFunction for StateSpace: '<S542>/Internal' */
  _rtXdot->Internal_CSTATE_a = 0.0;
  for (ri = windEmulatorStep4_WECSim_cal->Internal_A_jc_e[0U]; ri <
       windEmulatorStep4_WECSim_cal->Internal_A_jc_e[1U]; ri++) {
    _rtXdot->Internal_CSTATE_a += windEmulatorStep4_WECSim_cal->Internal_A_pr_e *
      windEmulatorStep4_WECSim_X.Internal_CSTATE_a;
  }

  for (ri = windEmulatorStep4_WECSim_cal->Internal_B_jc_h[0U]; ri <
       windEmulatorStep4_WECSim_cal->Internal_B_jc_h[1U]; ri++) {
    _rtXdot->Internal_CSTATE_a += windEmulatorStep4_WECSim_cal->Internal_B_pr_k *
      windEmulatorStep4_WECSim_B.ControlSignal1;
  }

  /* End of ForcingFunction for StateSpace: '<S542>/Internal' */

  /* ForcingFunction for SimscapeExecutionBlock: '<S216>/STATE_1' incorporates:
   *  SimscapeExecutionBlock: '<S332>/STATE_1'
   */
  simulationData = static_cast<NeslSimulationData *>
    (windEmulatorStep4_WECSim_DW.STATE_1_SimData_h);
  time_tmp = windEmulatorStep4_WECSim_M->Timing.t[0];
  time = time_tmp;
  simulationData->mData->mTime.mN = 1;
  simulationData->mData->mTime.mX = &time;
  simulationData->mData->mContStates.mN = 2;
  simulationData->mData->mContStates.mX =
    &windEmulatorStep4_WECSim_X.windEmulatorStep4_WECSimhptoSim[0];
  simulationData->mData->mDiscStates.mN = 0;
  simulationData->mData->mDiscStates.mX =
    &windEmulatorStep4_WECSim_DW.STATE_1_Discrete;
  simulationData->mData->mModeVector.mN = 0;
  simulationData->mData->mModeVector.mX =
    &windEmulatorStep4_WECSim_DW.STATE_1_Modes_j;
  tmp = false;
  simulationData->mData->mFoundZcEvents = tmp;
  simulationData->mData->mHadEvents = false;
  tmp = rtmIsMajorTimeStep(windEmulatorStep4_WECSim_M);
  simulationData->mData->mIsMajorTimeStep = tmp;
  tmp_0 = false;
  simulationData->mData->mIsSolverAssertCheck = tmp_0;
  simulationData->mData->mIsSolverCheckingCIC = false;
  tmp_0 = rtsiIsSolverComputingJacobian(&windEmulatorStep4_WECSim_M->solverInfo);
  simulationData->mData->mIsComputingJacobian = tmp_0;
  simulationData->mData->mIsEvaluatingF0 = false;
  simulationData->mData->mIsSolverRequestingReset = false;
  tmp_0 = rtsiIsModeUpdateTimeStep(&windEmulatorStep4_WECSim_M->solverInfo);
  simulationData->mData->mIsModeUpdateTimeStep = tmp_0;
  tmp_2[0] = 0;
  tmp_1[0] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[0];
  tmp_1[1] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[1];
  tmp_1[2] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[2];
  tmp_1[3] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[3];
  tmp_2[1] = 4;
  tmp_1[4] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[0];
  tmp_1[5] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[1];
  tmp_1[6] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[2];
  tmp_1[7] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[3];
  tmp_2[2] = 8;
  tmp_1[8] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[0];
  tmp_1[9] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[1];
  tmp_1[10] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[2];
  tmp_1[11] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[3];
  tmp_2[3] = 12;
  tmp_1[12] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[0];
  tmp_1[13] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[1];
  tmp_1[14] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[2];
  tmp_1[15] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[3];
  tmp_2[4] = 16;
  tmp_1[16] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[0];
  tmp_1[17] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[1];
  tmp_1[18] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[2];
  tmp_1[19] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[3];
  tmp_2[5] = 20;
  tmp_1[20] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[0];
  tmp_1[21] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[1];
  tmp_1[22] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[2];
  tmp_1[23] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[3];
  tmp_2[6] = 24;
  tmp_1[24] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[0];
  tmp_1[25] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[1];
  tmp_1[26] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[2];
  tmp_1[27] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[3];
  tmp_2[7] = 28;
  tmp_1[28] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[0];
  tmp_1[29] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[1];
  tmp_1[30] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[2];
  tmp_1[31] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[3];
  tmp_2[8] = 32;
  tmp_1[32] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[0];
  tmp_1[33] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[1];
  tmp_1[34] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[2];
  tmp_1[35] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[3];
  tmp_2[9] = 36;
  tmp_1[36] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[0];
  tmp_1[37] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[1];
  tmp_1[38] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[2];
  tmp_1[39] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[3];
  tmp_2[10] = 40;
  tmp_1[40] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[0];
  tmp_1[41] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[1];
  tmp_1[42] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[2];
  tmp_1[43] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[3];
  tmp_2[11] = 44;
  tmp_1[44] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[0];
  tmp_1[45] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[1];
  tmp_1[46] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[2];
  tmp_1[47] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[3];
  tmp_2[12] = 48;
  tmp_1[48] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[0];
  tmp_1[49] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[1];
  tmp_1[50] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[2];
  tmp_1[51] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[3];
  tmp_2[13] = 52;
  simulationData->mData->mInputValues.mN = 52;
  simulationData->mData->mInputValues.mX = &tmp_1[0];
  simulationData->mData->mInputOffsets.mN = 14;
  simulationData->mData->mInputOffsets.mX = &tmp_2[0];
  simulationData->mData->mDx.mN = 2;
  simulationData->mData->mDx.mX = &_rtXdot->windEmulatorStep4_WECSimhptoSim[0];
  simulator = static_cast<NeslSimulator *>
    (windEmulatorStep4_WECSim_DW.STATE_1_Simulator_f);
  diagnosticManager = static_cast<NeuDiagnosticManager *>
    (windEmulatorStep4_WECSim_DW.STATE_1_DiagMgr_o);
  diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
  is = ne_simulator_method(simulator, NESL_SIM_DERIVATIVES, simulationData,
    diagnosticManager);
  if (is != 0) {
    tmp_3 = error_buffer_is_empty(rtmGetErrorStatus(windEmulatorStep4_WECSim_M));
    if (tmp_3) {
      msg = rtw_diagnostics_msg(diagnosticTree);
      rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
    }
  }

  /* End of ForcingFunction for SimscapeExecutionBlock: '<S216>/STATE_1' */

  /* ForcingFunction for SimscapeInputBlock: '<S332>/INPUT_3_1_1' */
  _rtXdot->windEmulatorStep4_WECSimhptoS_n =
    (windEmulatorStep4_WECSim_B.velocity[4] -
     windEmulatorStep4_WECSim_X.windEmulatorStep4_WECSimhptoS_n) * 1000.0;

  /* ForcingFunction for SimscapeExecutionBlock: '<S332>/STATE_1' */
  simulationData = static_cast<NeslSimulationData *>
    (windEmulatorStep4_WECSim_DW.STATE_1_SimData_a);
  time_0 = time_tmp;
  simulationData->mData->mTime.mN = 1;
  simulationData->mData->mTime.mX = &time_0;
  simulationData->mData->mContStates.mN = 35;
  simulationData->mData->mContStates.mX =
    &windEmulatorStep4_WECSim_X.windEmulatorStep4_WECSimhptoS_h[0];
  simulationData->mData->mDiscStates.mN = 6;
  simulationData->mData->mDiscStates.mX =
    &windEmulatorStep4_WECSim_DW.STATE_1_Discrete_208214823[0];
  simulationData->mData->mModeVector.mN = 21;
  simulationData->mData->mModeVector.mX =
    &windEmulatorStep4_WECSim_DW.STATE_1_Modes_i[0];
  tmp_3 = false;
  simulationData->mData->mFoundZcEvents = tmp_3;
  simulationData->mData->mHadEvents = false;
  simulationData->mData->mIsMajorTimeStep = tmp;
  tmp = false;
  simulationData->mData->mIsSolverAssertCheck = tmp;
  simulationData->mData->mIsSolverCheckingCIC = false;
  tmp = rtsiIsSolverComputingJacobian(&windEmulatorStep4_WECSim_M->solverInfo);
  simulationData->mData->mIsComputingJacobian = tmp;
  simulationData->mData->mIsEvaluatingF0 = false;
  simulationData->mData->mIsSolverRequestingReset = false;
  simulationData->mData->mIsModeUpdateTimeStep = tmp_0;
  tmp_5[0] = 0;
  tmp_4[0] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[0];
  tmp_4[1] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[1];
  tmp_4[2] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[2];
  tmp_4[3] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[3];
  tmp_5[1] = 4;
  tmp_4[4] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[0];
  tmp_4[5] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[1];
  tmp_4[6] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[2];
  tmp_4[7] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[3];
  tmp_5[2] = 8;
  tmp_4[8] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[0];
  tmp_4[9] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[1];
  tmp_4[10] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[2];
  tmp_4[11] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[3];
  tmp_5[3] = 12;
  simulationData->mData->mInputValues.mN = 12;
  simulationData->mData->mInputValues.mX = &tmp_4[0];
  simulationData->mData->mInputOffsets.mN = 4;
  simulationData->mData->mInputOffsets.mX = &tmp_5[0];
  simulationData->mData->mDx.mN = 35;
  simulationData->mData->mDx.mX = &_rtXdot->windEmulatorStep4_WECSimhptoS_h[0];
  simulator = static_cast<NeslSimulator *>
    (windEmulatorStep4_WECSim_DW.STATE_1_Simulator_i);
  diagnosticManager = static_cast<NeuDiagnosticManager *>
    (windEmulatorStep4_WECSim_DW.STATE_1_DiagMgr_g);
  diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
  is = ne_simulator_method(simulator, NESL_SIM_FORCINGFUNCTION, simulationData,
    diagnosticManager);
  if (is != 0) {
    tmp = error_buffer_is_empty(rtmGetErrorStatus(windEmulatorStep4_WECSim_M));
    if (tmp) {
      msg = rtw_diagnostics_msg(diagnosticTree);
      rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
    }
  }
}

/* MassMatrix for root system: '<Root>' */
void windEmulatorStep4_WECSim_massmatrix(void)
{
  NeslSimulationData *simulationData;
  NeslSimulator *simulator;
  NeuDiagnosticManager *diagnosticManager;
  NeuDiagnosticTree *diagnosticTree;
  char *msg;
  real_T tmp_0[12];
  real_T time;
  real_T *tmp_2;
  real_T *tmp_3;
  int32_T tmp_4;
  int_T tmp_1[4];
  boolean_T tmp;

  /* MassMatrix for SimscapeExecutionBlock: '<S332>/STATE_1' */
  simulationData = static_cast<NeslSimulationData *>
    (windEmulatorStep4_WECSim_DW.STATE_1_SimData_a);
  time = windEmulatorStep4_WECSim_M->Timing.t[0];
  simulationData->mData->mTime.mN = 1;
  simulationData->mData->mTime.mX = &time;
  simulationData->mData->mContStates.mN = 35;
  simulationData->mData->mContStates.mX =
    &windEmulatorStep4_WECSim_X.windEmulatorStep4_WECSimhptoS_h[0];
  simulationData->mData->mDiscStates.mN = 6;
  simulationData->mData->mDiscStates.mX =
    &windEmulatorStep4_WECSim_DW.STATE_1_Discrete_208214823[0];
  simulationData->mData->mModeVector.mN = 21;
  simulationData->mData->mModeVector.mX =
    &windEmulatorStep4_WECSim_DW.STATE_1_Modes_i[0];
  tmp = false;
  simulationData->mData->mFoundZcEvents = tmp;
  simulationData->mData->mHadEvents = false;
  simulationData->mData->mIsMajorTimeStep = rtmIsMajorTimeStep
    (windEmulatorStep4_WECSim_M);
  tmp = false;
  simulationData->mData->mIsSolverAssertCheck = tmp;
  simulationData->mData->mIsSolverCheckingCIC = false;
  tmp = rtsiIsSolverComputingJacobian(&windEmulatorStep4_WECSim_M->solverInfo);
  simulationData->mData->mIsComputingJacobian = tmp;
  simulationData->mData->mIsEvaluatingF0 = false;
  simulationData->mData->mIsSolverRequestingReset = false;
  simulationData->mData->mIsModeUpdateTimeStep = rtsiIsModeUpdateTimeStep
    (&windEmulatorStep4_WECSim_M->solverInfo);
  tmp_1[0] = 0;
  tmp_0[0] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[0];
  tmp_0[1] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[1];
  tmp_0[2] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[2];
  tmp_0[3] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[3];
  tmp_1[1] = 4;
  tmp_0[4] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[0];
  tmp_0[5] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[1];
  tmp_0[6] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[2];
  tmp_0[7] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[3];
  tmp_1[2] = 8;
  tmp_0[8] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[0];
  tmp_0[9] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[1];
  tmp_0[10] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[2];
  tmp_0[11] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[3];
  tmp_1[3] = 12;
  simulationData->mData->mInputValues.mN = 12;
  simulationData->mData->mInputValues.mX = &tmp_0[0];
  simulationData->mData->mInputOffsets.mN = 4;
  simulationData->mData->mInputOffsets.mX = &tmp_1[0];
  tmp_2 = windEmulatorStep4_WECSim_M->massMatrixPr;
  tmp_3 = double_pointer_shift(tmp_2,
    windEmulatorStep4_WECSim_DW.STATE_1_MASS_MATRIX_PR);
  simulationData->mData->mMassMatrixPr.mN = 17;
  simulationData->mData->mMassMatrixPr.mX = tmp_3;
  simulator = static_cast<NeslSimulator *>
    (windEmulatorStep4_WECSim_DW.STATE_1_Simulator_i);
  diagnosticManager = static_cast<NeuDiagnosticManager *>
    (windEmulatorStep4_WECSim_DW.STATE_1_DiagMgr_g);
  diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
  tmp_4 = ne_simulator_method(simulator, NESL_SIM_MASSMATRIX, simulationData,
    diagnosticManager);
  if (tmp_4 != 0) {
    tmp = error_buffer_is_empty(rtmGetErrorStatus(windEmulatorStep4_WECSim_M));
    if (tmp) {
      msg = rtw_diagnostics_msg(diagnosticTree);
      rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
    }
  }

  /* End of MassMatrix for SimscapeExecutionBlock: '<S332>/STATE_1' */
}

void local_evaluateMassMatrix(RTWSolverInfo *si, real_T *Mdest )
{
  /* Refresh global mass matrix */
  windEmulatorStep4_WECSim_massmatrix();

  /* Copy the mass matrix from system to the destination, if needed. */
  if (Mdest != rtsiGetSolverMassMatrixPr(si)) {
    real_T *Msrc = rtsiGetSolverMassMatrixPr(si);
    int_T nzmax = rtsiGetSolverMassMatrixNzMax(si);
    (void) std::memcpy(Mdest, Msrc,
                       static_cast<uint_T>(nzmax)*sizeof(real_T));
  }
}

void local_evaluateFminusMv(RTWSolverInfo *si, const real_T *v, real_T *fminusMv
  )
{
  /* Refresh forcing function */
  rtsiSetdX(si,fminusMv);
  windEmulatorStep4_WECSim_forcingfunction();

  /* Refresh global mass matrix */
  windEmulatorStep4_WECSim_massmatrix();

  /* Form f - M*v */
  {
    real_T *elptr = rtsiGetSolverMassMatrixPr(si);
    int_T *iptr = rtsiGetSolverMassMatrixIr(si);
    int_T *jc = rtsiGetSolverMassMatrixJc(si);
    int_T nx = 44;
    int_T col,row;
    for (col = 0; col < nx; col++) {
      for (row = jc[col]; row < jc[col+1]; row++) {
        fminusMv[*iptr++] -= (*v) * (*elptr++);
      }

      v++;
    }
  }
}

/* Simplified version of numjac.cpp, for use with RTW. */
void local_numjac( RTWSolverInfo *si, real_T *y, const real_T *v, const real_T
                  *Fty, real_T *fac, real_T *dFdy )
{
  /* constants */
  real_T THRESH = 1e-6;
  real_T EPS = 2.2e-16;                /* utGetEps(); */
  real_T BL = std::pow(EPS, 0.75);
  real_T BU = std::pow(EPS, 0.25);
  real_T FACMIN = std::pow(EPS, 0.78);
  real_T FACMAX = 0.1;
  int_T nx = 44;
  real_T *x = rtsiGetContStates(si);
  boolean_T *xdis = rtsiGetContStateDisabledPtr(si);
  real_T del;
  real_T difmax;
  real_T FdelRowmax;
  real_T temp;
  real_T Fdiff;
  real_T maybe;
  real_T xscale;
  real_T fscale;
  real_T *p;
  int_T rowmax;
  int_T i,j;
  if (x != y)
    (void) std::memcpy(x, y,
                       static_cast<uint_T>(nx)*sizeof(real_T));
  rtsiSetSolverComputingJacobian(si,true);
  for (p = dFdy, j = 0; j < nx; j++, p += nx) {
    /* Zero column j of dFdy if state j is currently disabled. */
    if (xdis[j]) {
      (void) std::memset(p, 0,
                         (uint_T)nx*sizeof(p[0]));
      continue;
    }

    /* Select an increment del for a difference approximation to
       column j of dFdy.  The vector fac accounts for experience
       gained in previous calls to numjac. */
    xscale = std::abs(x[j]);
    if (xscale < THRESH)
      xscale = THRESH;
    temp = (x[j] + fac[j]*xscale);
    del = temp - y[j];
    while (del == 0.0) {
      if (fac[j] < FACMAX) {
        fac[j] *= 100.0;
        if (fac[j] > FACMAX)
          fac[j] = FACMAX;
        temp = (x[j] + fac[j]*xscale);
        del = temp - x[j];
      } else {
        del = THRESH;                  /* thresh is nonzero */
        break;
      }
    }

    /* Keep del pointing into region. */
    if (Fty[j] >= 0.0)
      del = std::abs(del);
    else
      del = -std::abs(del);

    /* Form a difference approximation to column j of dFdy. */
    temp = x[j];
    x[j] += del;
    windEmulatorStep4_WECSim_step();
    local_evaluateFminusMv(si,v,p );
    x[j] = temp;
    difmax = 0.0;
    rowmax = 0;
    FdelRowmax = p[0];
    temp = 1.0 / del;
    for (i = 0; i < nx; i++) {
      Fdiff = p[i] - Fty[i];
      maybe = std::abs(Fdiff);
      if (maybe > difmax) {
        difmax = maybe;
        rowmax = i;
        FdelRowmax = p[i];
      }

      p[i] = temp * Fdiff;
    }

    /* Adjust fac for next call to numjac. */
    if (((FdelRowmax != 0.0) && (Fty[rowmax] != 0.0)) || (difmax == 0.0)) {
      fscale = std::abs(FdelRowmax);
      if (fscale < std::abs(Fty[rowmax]))
        fscale = std::abs(Fty[rowmax]);
      if (difmax <= BL*fscale) {
        /* The difference is small, so increase the increment. */
        fac[j] *= 10.0;
        if (fac[j] > FACMAX)
          fac[j] = FACMAX;
      } else if (difmax > BU*fscale) {
        /* The difference is large, so reduce the increment. */
        fac[j] *= 0.1;
        if (fac[j] < FACMIN)
          fac[j] = FACMIN;
      }
    }
  }

  rtsiSetSolverComputingJacobian(si,false);
}                                      /* end local_numjac */

/*
 * This function updates continuous states using the ODE14X fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  /* Solver Matrices */
  static int_T rt_ODE14x_N[4] = { 12, 8, 6, 4 };

  time_T t0 = rtsiGetT(si);
  time_T t1 = t0;
  time_T h = rtsiGetStepSize(si);
  real_T *x1 = rtsiGetContStates(si);
  int_T order = rtsiGetSolverExtrapolationOrder(si);
  int_T numIter = rtsiGetSolverNumberNewtonIterations(si);
  ODE14X_IntgData *id = static_cast<ODE14X_IntgData *>(rtsiGetSolverData(si));
  real_T *x0 = id->x0;
  real_T *f0 = id->f0;
  real_T *x1start = id->x1start;
  real_T *f1 = id->f1;
  real_T *Delta = id->Delta;
  real_T *E = id->E;
  real_T *fac = id->fac;
  real_T *dfdx = id->DFDX;
  real_T *W = id->W;
  int_T *pivots = id->pivots;
  real_T *xtmp = id->xtmp;
  real_T *ztmp = id->ztmp;
  boolean_T *xdis = rtsiGetContStateDisabledPtr(si);
  int_T *Mpattern_ir = rtsiGetSolverMassMatrixIr(si);
  int_T *Mpattern_jc = rtsiGetSolverMassMatrixJc(si);
  real_T *M = id->M;
  real_T *M1 = id->M1;
  real_T *xdot = id->xdot;
  real_T *Edot = id->Edot;
  real_T *fminusMxdot = id->fminusMxdot;
  int_T col,row,rowidx;
  int_T *N = &(rt_ODE14x_N[0]);
  int_T i,j,k,iter;
  int_T nx = 44;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) std::memcpy(x0, x1,
                     static_cast<uint_T>(nx)*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  local_evaluateMassMatrix(si,M );
  rtsiSetdX(si, xdot);
  windEmulatorStep4_WECSim_derivatives();

  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  windEmulatorStep4_WECSim_forcingfunction();

  /* Form fminusMxdot = f(x) - M(x)*xdot, d(fminusMxdot)/dx = df/dx - d(Mv)/dx */
  (void) std::memcpy(fminusMxdot, f0,
                     static_cast<uint_T>(nx)*sizeof(real_T));
  for (col = 0; col < nx; col++) {
    for (rowidx = Mpattern_jc[col]; rowidx < Mpattern_jc[col+1]; rowidx++) {
      real_T m_row_col = M[rowidx];
      row = Mpattern_ir[rowidx];
      fminusMxdot[row] -= m_row_col*xdot[col];
    }
  }

  local_numjac(si,x0,xdot,fminusMxdot,fac,dfdx );
  for (j = 0; j < order; j++) {
    real_T *p;
    real_T hN = h/N[j];

    /* Get the iteration matrix and solution at t0 */

    /* [L,U] = lu(M - hN*J) */
    (void) std::memcpy(W, dfdx,
                       static_cast<uint_T>(nx)*static_cast<uint_T>(nx)*sizeof
                       (real_T));
    for (p = W, i = 0; i < nx*nx; i++, p++) {
      *p *= (-hN);
    }

    for (col = 0, p = W; col < nx; col++, p += nx) {
      if (xdis[col]) {
        (void) std::memset(p, 0,
                           static_cast<uint_T>(nx)*sizeof(p[0]));
        p[col] = 1.0;
      } else {
        for (rowidx = Mpattern_jc[col]; rowidx < Mpattern_jc[col+1]; rowidx++) {
          real_T m_row_col = M[rowidx];
          row = Mpattern_ir[rowidx];
          p[row] += m_row_col;
        }
      }
    }

    rt_lu_real(W, nx,
               pivots);

    /* First Newton's iteration at t0. */
    /* rhs = hN*f0 */
    for (i = 0; i < nx; i++) {
      Delta[i] = hN*f0[i];
    }

    /* Delta = (U \ (L \ rhs)) */
    rt_ForwardSubstitutionRR_Dbl(W, Delta,
      f1, nx,
      1, pivots,
      1);
    rt_BackwardSubstitutionRR_Dbl(W+nx*nx-1, f1+nx-1,
      Delta, nx,
      1, 0);

    /* ytmp = y0 + Delta
       ztmp = (ytmp-y0)/h
     */
    (void) std::memcpy(x1, x0,
                       static_cast<uint_T>(nx)*sizeof(real_T));
    for (i = 0; i < nx; i++) {
      x1[i] += Delta[i];
      ztmp[i] = Delta[i]/hN;
    }

    /* Additional Newton's iterations, if desired.
       for iter = 2:NewtIter
       rhs = hN*feval(odefun,tn,ytmp,extraArgs{:}) - M*(ytmp - yn);
       if statedepM   % only for state dep. Mdel ~= 0
       Mdel = M - feval(massfun,tn,ytmp);
       rhs = rhs + Mdel*ztmp*h;
       end
       Delta = ( U \ ( L \ rhs ) );
       ytmp = ytmp + Delta;
       ztmp = (ytmp - yn)/h
       end
     */
    rtsiSetT(si, t0);
    rtsiSetdX(si, f1);
    for (iter = 1; iter < numIter; iter++) {
      windEmulatorStep4_WECSim_step();
      windEmulatorStep4_WECSim_forcingfunction();
      for (i = 0; i < nx; i++) {
        Delta[i] = hN*f1[i];
        xtmp[i] = x1[i] - x0[i];
      }

      /* rhs = hN*f(tn,ytmp) - M*(ytmp-yn) */
      for (col = 0; col < nx; col++) {
        for (rowidx = Mpattern_jc[col]; rowidx < Mpattern_jc[col+1]; rowidx++) {
          real_T m_row_col = M[rowidx];
          row = Mpattern_ir[rowidx];
          Delta[row] -= m_row_col*xtmp[col];
        }
      }

      /* rhs = rhs - (Mtmp - M)*ztmp*h */
      local_evaluateMassMatrix(si,M1 );
      for (i = 0; i < rtsiGetSolverMassMatrixNzMax(si); i++) {
        M1[i] -= M[i];
      }

      for (col = 0; col < nx; col++) {
        for (rowidx = Mpattern_jc[col]; rowidx < Mpattern_jc[col+1]; rowidx++) {
          real_T m_row_col = M1[rowidx];
          row = Mpattern_ir[rowidx];
          Delta[row] -= hN*m_row_col*ztmp[col];
        }
      }

      rt_ForwardSubstitutionRR_Dbl(W, Delta,
        f1, nx,
        1, pivots,
        1);
      rt_BackwardSubstitutionRR_Dbl(W+nx*nx-1, f1+nx-1,
        Delta, nx,
        1, 0);

      /* ytmp = ytmp + delta
         ztmp = (ytmp - yn)/h
       */
      for (i = 0; i < nx; i++) {
        x1[i] += Delta[i];
        ztmp[i] = (x1[i] - x0[i])/hN;
      }
    }

    /* Steps from t0+hN to t1 -- subintegration of N(j) steps for extrapolation
       ttmp = t0;
       for i = 2:N(j)
       ttmp = ttmp + hN
       ytmp0 = ytmp;
       for iter = 1:NewtIter
       rhs = (ytmp0 - ytmp) + hN*feval(odefun,ttmp,ytmp,extraArgs{:});
       Delta = ( U \ ( L \ rhs ) );
       ytmp = ytmp + Delta;
       end
       end
     */
    for (k = 1; k < N[j]; k++) {
      t1 = t0 + k*hN;
      (void) std::memcpy(x1start, x1,
                         static_cast<uint_T>(nx)*sizeof(real_T));
      rtsiSetT(si, t1);
      rtsiSetdX(si, f1);
      for (iter = 0; iter < numIter; iter++) {
        windEmulatorStep4_WECSim_step();
        windEmulatorStep4_WECSim_forcingfunction();
        if (iter == 0) {
          for (i = 0; i < nx; i++) {
            Delta[i] = hN*f1[i];
          }
        } else {
          for (i = 0; i < nx; i++) {
            Delta[i] = hN*f1[i];
            xtmp[i] = (x1[i]-x1start[i]);
          }

          /* rhs = hN*f(tn,ytmp) - M*(ytmp-yn) */
          for (col = 0; col < nx; col++) {
            for (rowidx = Mpattern_jc[col]; rowidx < Mpattern_jc[col+1]; rowidx
                 ++) {
              real_T m_row_col = M[rowidx];
              row = Mpattern_ir[rowidx];
              Delta[row] -= m_row_col*xtmp[col];
            }
          }
        }

        /* For state-dep.,  Mdel = M(ttmp,ytmp) - M */
        windEmulatorStep4_WECSim_step();
        local_evaluateMassMatrix(si,M1 );
        for (i = 0; i < rtsiGetSolverMassMatrixNzMax(si); i++) {
          M1[i] -= M[i];
        }

        /* rhs = rhs - Mdel*ztmp*h */
        for (col = 0; col < nx; col++) {
          for (rowidx = Mpattern_jc[col]; rowidx < Mpattern_jc[col+1]; rowidx++)
          {
            real_T m_row_col = M1[rowidx];
            row = Mpattern_ir[rowidx];
            Delta[row] -= hN*m_row_col*ztmp[col];
          }
        }

        rt_ForwardSubstitutionRR_Dbl(W, Delta,
          f1, nx,
          1, pivots,
          1);
        rt_BackwardSubstitutionRR_Dbl(W+nx*nx-1, f1+nx-1,
          Delta, nx,
          1, 0);

        /* ytmp = ytmp + Delta
           ztmp = (ytmp - ytmp0)/h
         */
        for (i = 0; i < nx; i++) {
          x1[i] += Delta[i];
          ztmp[i] = (x1[i] - x1start[i])/hN;
        }
      }
    }

    /* Extrapolate to order j
       E(:,j) = ytmp
       for k = j:-1:2
       coef = N(k-1)/(N(j) - N(k-1))
       E(:,k-1) = E(:,k) + coef*( E(:,k) - E(:,k-1) )
       end
     */
    (void) std::memcpy(&(E[nx*j]), x1,
                       static_cast<uint_T>(nx)*sizeof(real_T));
    for (k = j; k > 0; k--) {
      real_T coef = static_cast<real_T>(N[k-1]) / (N[j]-N[k-1]);
      for (i = 0; i < nx; i++) {
        x1[i] = E[nx*k+i] + coef*(E[nx*k+i] - E[nx*(k-1)+i]);
      }

      (void) std::memcpy(&(E[nx*(k-1)]), x1,
                         static_cast<uint_T>(nx)*sizeof(real_T));
    }

    /* Extrapolate the derivative */
    for (i = 0; i < nx; i++) {
      xdot[i] = (x1[i] - x1start[i])/hN;
    }

    (void) std::memcpy(&(Edot[nx*j]), xdot,
                       static_cast<uint_T>(nx)*sizeof(real_T));
    for (k = j; k > 0; k--) {
      real_T coef = static_cast<real_T>(N[k-1]) / (N[j]-N[k-1]);
      for (i = 0; i < nx; i++) {
        xdot[i] = Edot[nx*k+i] + coef*(Edot[nx*k+i] - Edot[nx*(k-1)+i]);
      }

      (void) std::memcpy(&(Edot[nx*(k-1)]), xdot,
                         static_cast<uint_T>(nx)*sizeof(real_T));
    }
  }

  /* x1 = E(:,1); */
  (void) std::memcpy(x1, E,
                     static_cast<uint_T>(nx)*sizeof(real_T));

  /* Extrapolated xdot */
  (void) std::memcpy(xdot, Edot,
                     static_cast<uint_T>(nx)*sizeof(real_T));

  /* t1 = t0 + h; */
  rtsiSetT(si,rtsiGetSolverStopTime(si));
  windEmulatorStep4_WECSim_step();
  windEmulatorStep4_WECSim_projection();
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/*
 * System initialize for atomic system:
 *    '<S41>/parseCtrlWord'
 *    '<S46>/parseCtrlWord'
 */
void windEmulator_parseCtrlWord_Init(DW_parseCtrlWord_windEmulator_T *localDW)
{
  localDW->sfEvent = windEmulatorStep4_WE_CALL_EVENT;
}

/*
 * Output and update for atomic system:
 *    '<S41>/parseCtrlWord'
 *    '<S46>/parseCtrlWord'
 */
void windEmulatorStep4_parseCtrlWord(uint16_T rtu_ctrlWord,
  B_parseCtrlWord_windEmulatorS_T *localB, DW_parseCtrlWord_windEmulator_T
  *localDW)
{
  localDW->sfEvent = windEmulatorStep4_WE_CALL_EVENT;
  localB->off1Ctrl = ((rtu_ctrlWord & 1U) != 0U);
  localB->off2Ctrl = ((rtu_ctrlWord & 2U) != 0U);
  localB->off3Ctrl = ((rtu_ctrlWord & 4U) != 0U);
  localB->enableOperation = ((rtu_ctrlWord & 8U) != 0U);
  localB->rampOutZero = ((rtu_ctrlWord & 16U) != 0U);
  localB->rampHold = ((rtu_ctrlWord & 32U) != 0U);
  localB->rampInZero = ((rtu_ctrlWord & 64U) != 0U);
  localB->reset = ((rtu_ctrlWord & 128U) != 0U);
  localB->inching1 = ((rtu_ctrlWord & 256U) != 0U);
  localB->inching2 = ((rtu_ctrlWord & 512U) != 0U);
  localB->remoteCmd = ((rtu_ctrlWord & 1024U) != 0U);
  localB->extCtrlLoc = ((rtu_ctrlWord & 2048U) != 0U);
}

static void windEmulatorSt_SystemCore_setup(dsp_simulink_MovingAverage_wi_T *obj)
{
  obj->isInitialized = 1;

  /* Start for MATLABSystem: '<S43>/Moving Average' */
  obj->NumChannels = 1;
  obj->FrameLength = 1;
  obj->pCumSum = 0.0;
  std::memset(&obj->pCumSumRev[0], 0, 2499U * sizeof(real_T));
  obj->pCumRevIndex = 1.0;
  obj->pModValueRev = 0.0;
  obj->isSetupComplete = true;
  obj->TunablePropsChanged = false;
}

/* System initialize for atomic system: */
void windEmulator_MovingAverage_Init(DW_MovingAverage_windEmulator_T *localDW)
{
  /* Start for MATLABSystem: '<S43>/Moving Average' */
  localDW->obj.isInitialized = 0;
  localDW->obj.NumChannels = -1;
  localDW->obj.FrameLength = -1;
  localDW->obj.matlabCodegenIsDeleted = false;
  localDW->objisempty = true;
  windEmulatorSt_SystemCore_setup(&localDW->obj);

  /* InitializeConditions for MATLABSystem: '<S43>/Moving Average' */
  localDW->obj.pCumSum = 0.0;
  std::memset(&localDW->obj.pCumSumRev[0], 0, 2499U * sizeof(real_T));
  localDW->obj.pCumRevIndex = 1.0;
  localDW->obj.pModValueRev = 0.0;
}

/* Output and update for atomic system: */
void windEmulatorStep4_MovingAverage(real_T rtu_0,
  B_MovingAverage_windEmulatorS_T *localB, DW_MovingAverage_windEmulator_T
  *localDW)
{
  real_T csum;
  real_T cumRevIndex;
  real_T modValueRev;
  real_T tmp;
  real_T z;

  /* MATLABSystem: '<S43>/Moving Average' */
  if (localDW->obj.TunablePropsChanged) {
    localDW->obj.TunablePropsChanged = false;
  }

  cumRevIndex = localDW->obj.pCumRevIndex;
  csum = localDW->obj.pCumSum;
  std::memcpy(&localB->csumrev[0], &localDW->obj.pCumSumRev[0], 2499U * sizeof
              (real_T));
  modValueRev = localDW->obj.pModValueRev;
  z = 0.0;
  tmp = 0.0;
  csum += rtu_0;
  if (modValueRev == 0.0) {
    z = localB->csumrev[static_cast<int32_T>(cumRevIndex) - 1] + csum;
  }

  localB->csumrev[static_cast<int32_T>(cumRevIndex) - 1] = rtu_0;
  if (cumRevIndex != 2499.0) {
    cumRevIndex++;
  } else {
    cumRevIndex = 1.0;
    csum = 0.0;
    for (int32_T i = 2497; i >= 0; i--) {
      localB->csumrev[i] += localB->csumrev[i + 1];
    }
  }

  if (modValueRev == 0.0) {
    tmp = z / 2500.0;
  }

  if (modValueRev > 0.0) {
    modValueRev--;
  } else {
    modValueRev = 0.0;
  }

  localDW->obj.pCumSum = csum;
  std::memcpy(&localDW->obj.pCumSumRev[0], &localB->csumrev[0], 2499U * sizeof
              (real_T));
  localDW->obj.pCumRevIndex = cumRevIndex;
  localDW->obj.pModValueRev = modValueRev;

  /* MATLABSystem: '<S43>/Moving Average' */
  localB->MovingAverage = tmp;
}

/* Termination for atomic system: */
void windEmulator_MovingAverage_Term(DW_MovingAverage_windEmulator_T *localDW)
{
  /* Terminate for MATLABSystem: '<S43>/Moving Average' */
  if (!localDW->obj.matlabCodegenIsDeleted) {
    localDW->obj.matlabCodegenIsDeleted = true;
    if ((localDW->obj.isInitialized == 1) && localDW->obj.isSetupComplete) {
      localDW->obj.NumChannels = -1;
      localDW->obj.FrameLength = -1;
    }
  }

  /* End of Terminate for MATLABSystem: '<S43>/Moving Average' */
}

/*
 * System initialize for atomic system:
 *    '<S44>/Parse Status Word'
 *    '<S49>/Parse Status Word'
 */
void windEmulat_ParseStatusWord_Init(DW_ParseStatusWord_windEmulat_T *localDW)
{
  localDW->sfEvent = windEmulatorStep4__CALL_EVENT_n;
}

/*
 * Output and update for atomic system:
 *    '<S44>/Parse Status Word'
 *    '<S49>/Parse Status Word'
 */
void windEmulatorSte_ParseStatusWord(uint16_T rtu_StatusWord,
  B_ParseStatusWord_windEmulato_T *localB, DW_ParseStatusWord_windEmulat_T
  *localDW)
{
  localDW->sfEvent = windEmulatorStep4__CALL_EVENT_n;
  localB->rdy_on = ((rtu_StatusWord & 1U) != 0U);
  localB->rdy_run = ((rtu_StatusWord & 2U) != 0U);
  localB->rdy_ref = ((rtu_StatusWord & 4U) != 0U);
  localB->tripped = ((rtu_StatusWord & 8U) != 0U);
  localB->off2 = ((rtu_StatusWord & 16U) != 0U);
  localB->off3 = ((rtu_StatusWord & 32U) != 0U);
  localB->swc_on_inhib = ((rtu_StatusWord & 64U) != 0U);
  localB->alarm = ((rtu_StatusWord & 128U) != 0U);
  localB->at_setpoint = ((rtu_StatusWord & 256U) != 0U);
  localB->remote = ((rtu_StatusWord & 512U) != 0U);
  localB->above_limit = ((rtu_StatusWord & 1024U) != 0U);
  localB->ext_ctrl_loc = ((rtu_StatusWord & 2048U) != 0U);
  localB->ext_run_enable = ((rtu_StatusWord & 4096U) != 0U);
  localB->msw_b13 = ((rtu_StatusWord & 8192U) != 0U);
  localB->msw_b14 = ((rtu_StatusWord & 16384U) != 0U);
  localB->comm_err = ((rtu_StatusWord & 32768U) != 0U);
}

static void windEmulator_SystemCore_setup_n(dsp_simulink_MovingAverage_wi_T *obj)
{
  obj->isInitialized = 1;

  /* Start for MATLABSystem: '<S51>/Moving Average' */
  obj->NumChannels = 1;
  obj->FrameLength = 1;
  obj->pCumSum = 0.0;
  std::memset(&obj->pCumSumRev[0], 0, 2499U * sizeof(real_T));
  obj->pCumRevIndex = 1.0;
  obj->pModValueRev = 0.0;
  obj->isSetupComplete = true;
  obj->TunablePropsChanged = false;
}

/* System initialize for atomic system: */
void windEmulat_MovingAverage_e_Init(DW_MovingAverage_windEmulat_f_T *localDW)
{
  /* Start for MATLABSystem: '<S51>/Moving Average' */
  localDW->obj.isInitialized = 0;
  localDW->obj.NumChannels = -1;
  localDW->obj.FrameLength = -1;
  localDW->obj.matlabCodegenIsDeleted = false;
  localDW->objisempty = true;
  windEmulator_SystemCore_setup_n(&localDW->obj);

  /* InitializeConditions for MATLABSystem: '<S51>/Moving Average' */
  localDW->obj.pCumSum = 0.0;
  std::memset(&localDW->obj.pCumSumRev[0], 0, 2499U * sizeof(real_T));
  localDW->obj.pCumRevIndex = 1.0;
  localDW->obj.pModValueRev = 0.0;
}

/* Output and update for atomic system: */
void windEmulatorSte_MovingAverage_p(real_T rtu_0,
  B_MovingAverage_windEmulato_c_T *localB, DW_MovingAverage_windEmulat_f_T
  *localDW)
{
  real_T csum;
  real_T cumRevIndex;
  real_T modValueRev;
  real_T tmp;
  real_T z;

  /* MATLABSystem: '<S51>/Moving Average' */
  if (localDW->obj.TunablePropsChanged) {
    localDW->obj.TunablePropsChanged = false;
  }

  cumRevIndex = localDW->obj.pCumRevIndex;
  csum = localDW->obj.pCumSum;
  std::memcpy(&localB->csumrev[0], &localDW->obj.pCumSumRev[0], 2499U * sizeof
              (real_T));
  modValueRev = localDW->obj.pModValueRev;
  z = 0.0;
  tmp = 0.0;
  csum += rtu_0;
  if (modValueRev == 0.0) {
    z = localB->csumrev[static_cast<int32_T>(cumRevIndex) - 1] + csum;
  }

  localB->csumrev[static_cast<int32_T>(cumRevIndex) - 1] = rtu_0;
  if (cumRevIndex != 2499.0) {
    cumRevIndex++;
  } else {
    cumRevIndex = 1.0;
    csum = 0.0;
    for (int32_T i = 2497; i >= 0; i--) {
      localB->csumrev[i] += localB->csumrev[i + 1];
    }
  }

  if (modValueRev == 0.0) {
    tmp = z / 2500.0;
  }

  if (modValueRev > 0.0) {
    modValueRev--;
  } else {
    modValueRev = 0.0;
  }

  localDW->obj.pCumSum = csum;
  std::memcpy(&localDW->obj.pCumSumRev[0], &localB->csumrev[0], 2499U * sizeof
              (real_T));
  localDW->obj.pCumRevIndex = cumRevIndex;
  localDW->obj.pModValueRev = modValueRev;

  /* MATLABSystem: '<S51>/Moving Average' */
  localB->MovingAverage = tmp;
}

/* Termination for atomic system: */
void windEmulat_MovingAverage_g_Term(DW_MovingAverage_windEmulat_f_T *localDW)
{
  /* Terminate for MATLABSystem: '<S51>/Moving Average' */
  if (!localDW->obj.matlabCodegenIsDeleted) {
    localDW->obj.matlabCodegenIsDeleted = true;
    if ((localDW->obj.isInitialized == 1) && localDW->obj.isSetupComplete) {
      localDW->obj.NumChannels = -1;
      localDW->obj.FrameLength = -1;
    }
  }

  /* End of Terminate for MATLABSystem: '<S51>/Moving Average' */
}

/*
 * Output and update for atomic system:
 *    '<S60>/Nonlinear Wave Elevation'
 *    '<S139>/Nonlinear Wave Elevation'
 */
void windEmul_NonlinearWaveElevation(const real_T rtu_Displacement[3], const
  real_T rtu_Displacement_h[3], B_NonlinearWaveElevation_wind_T *localB,
  NonlinearWaveElevation_cal_type *windEmulatorS_PageSwitching_arg)
{
  /* Sum: '<S65>/Add' incorporates:
   *  Constant: '<S65>/Center of Gravity'
   *  Constant: '<S65>/Constant'
   */
  localB->x_cg[0] = rtu_Displacement[0] -
    windEmulatorS_PageSwitching_arg->CenterofGravity_Value[0];
  localB->x_cg[3] = rtu_Displacement_h[0] -
    windEmulatorS_PageSwitching_arg->Constant_Value[0];
  localB->x_cg[1] = rtu_Displacement[1] -
    windEmulatorS_PageSwitching_arg->CenterofGravity_Value[1];
  localB->x_cg[4] = rtu_Displacement_h[1] -
    windEmulatorS_PageSwitching_arg->Constant_Value[1];
  localB->x_cg[2] = rtu_Displacement[2] -
    windEmulatorS_PageSwitching_arg->CenterofGravity_Value[2];
  localB->x_cg[5] = rtu_Displacement_h[2] -
    windEmulatorS_PageSwitching_arg->Constant_Value[2];

  /* Constant: '<S80>/zero' */
  localB->zero = windEmulatorS_PageSwitching_arg->zero_Value;
}

real_T rt_atan2d_snf(real_T u0, real_T u1)
{
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    int32_T tmp;
    int32_T tmp_0;
    if (u1 > 0.0) {
      tmp = 1;
    } else {
      tmp = -1;
    }

    if (u0 > 0.0) {
      tmp_0 = 1;
    } else {
      tmp_0 = -1;
    }

    y = std::atan2(static_cast<real_T>(tmp_0), static_cast<real_T>(tmp));
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = std::atan2(u0, u1);
  }

  return y;
}

/*
 * System initialize for atomic system:
 *    '<S81>/quaternion2EulXYZ'
 *    '<S160>/quaternion2EulXYZ'
 */
void windEmul_quaternion2EulXYZ_Init(DW_quaternion2EulXYZ_windEmul_T *localDW)
{
  localDW->sfEvent = windEmulatorStep4__CALL_EVENT_g;
}

/*
 * Output and update for atomic system:
 *    '<S81>/quaternion2EulXYZ'
 *    '<S160>/quaternion2EulXYZ'
 */
void windEmulatorS_quaternion2EulXYZ(const real_T rtu_Q[4],
  B_quaternion2EulXYZ_windEmula_T *localB, DW_quaternion2EulXYZ_windEmul_T
  *localDW)
{
  real_T E_tmp;
  real_T absxk;
  real_T scale;
  real_T t;
  real_T y;
  int32_T exitg1;
  int32_T k;
  boolean_T b;
  localDW->sfEvent = windEmulatorStep4__CALL_EVENT_g;
  scale = 3.3121686421112381E-170;
  absxk = std::abs(rtu_Q[0]);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    y = t * t;
  }

  absxk = std::abs(rtu_Q[1]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  absxk = std::abs(rtu_Q[2]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  absxk = std::abs(rtu_Q[3]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  y = scale * std::sqrt(y);
  b = rtIsNaN(y);
  if (b) {
    k = 0;
    do {
      exitg1 = 0;
      if (k < 4) {
        if (rtIsNaN(rtu_Q[k])) {
          exitg1 = 1;
        } else {
          k++;
        }
      } else {
        y = (rtInf);
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  scale = rtu_Q[0] / y;
  absxk = rtu_Q[1] / y;
  t = rtu_Q[2] / y;
  y = rtu_Q[3] / y;
  E_tmp = t * t;
  localB->E[0] = rt_atan2d_snf((scale * absxk + t * y) * 2.0, 1.0 - (absxk *
    absxk + E_tmp) * 2.0);
  localB->E[1] = std::asin((scale * t - absxk * y) * 2.0);
  localB->E[2] = rt_atan2d_snf((scale * y + absxk * t) * 2.0, 1.0 - (y * y +
    E_tmp) * 2.0);
}

/*
 * System initialize for atomic system:
 *    '<S126>/MATLAB Function1'
 *    '<S205>/MATLAB Function1'
 */
void windEmulat_MATLABFunction1_Init(DW_MATLABFunction1_windEmulat_T *localDW)
{
  localDW->sfEvent = windEmulatorStep4__CALL_EVENT_b;
}

/*
 * Output and update for atomic system:
 *    '<S126>/MATLAB Function1'
 *    '<S205>/MATLAB Function1'
 */
void windEmulatorSte_MATLABFunction1(const real_T rtu_disp[2], real_T rtu_enable,
  real_T rtu_wavenumber, real_T rtu_direction, B_MATLABFunction1_windEmulato_T
  *localB, DW_MATLABFunction1_windEmulat_T *localDW)
{
  localDW->sfEvent = windEmulatorStep4__CALL_EVENT_b;
  localB->dispPhase = 0.0;
  if (rtu_enable == 1.0) {
    real_T dispPhase_tmp;
    dispPhase_tmp = rtu_direction * 3.1415926535897931 / 180.0;
    localB->dispPhase = (std::cos(dispPhase_tmp) * rtu_disp[0] + std::sin
                         (dispPhase_tmp) * rtu_disp[1]) * -rtu_wavenumber;
  }
}

/*
 * System initialize for atomic system:
 *    '<S70>/Yaw Force Transforms'
 *    '<S149>/Yaw Force Transforms'
 */
void windEmu_YawForceTransforms_Init(DW_YawForceTransforms_windEmu_T *localDW)
{
  localDW->sfEvent = windEmulatorStep4__CALL_EVENT_c;
}

/*
 * Output and update for atomic system:
 *    '<S70>/Yaw Force Transforms'
 *    '<S149>/Yaw Force Transforms'
 */
void windEmulator_YawForceTransforms(real_T rtu_yaw, const real_T
  rtu_dispGlobal[3], const real_T rtu_dispGlobal_e[3], const real_T
  rtu_F_RadiationDampingLocal[6], const real_T rtu_F_AddedMassLocal[6], const
  real_T rtu_F_ExcitationLocal[6], const real_T rtu_F_RestoringLocal[6],
  B_YawForceTransforms_windEmul_T *localB, DW_YawForceTransforms_windEmu_T
  *localDW)
{
  real_T rotMatYaw[9];
  real_T rotMatYaw_0[3];

  /* SignalConversion generated from: '<S132>/ SFunction ' */
  localB->TmpSignalConversionAtSFunctionI[0] = rtu_dispGlobal[0];
  localB->TmpSignalConversionAtSFunctionI[3] = rtu_dispGlobal_e[0];
  localB->TmpSignalConversionAtSFunctionI[1] = rtu_dispGlobal[1];
  localB->TmpSignalConversionAtSFunctionI[4] = rtu_dispGlobal_e[1];
  localB->TmpSignalConversionAtSFunctionI[2] = rtu_dispGlobal[2];
  localB->TmpSignalConversionAtSFunctionI[5] = rtu_dispGlobal_e[2];
  localDW->sfEvent = windEmulatorStep4__CALL_EVENT_c;
  if (rtu_yaw == 1.0) {
    real_T rotMatYaw_1;
    real_T rotMatYaw_2;
    real_T rotMatYaw_3;
    real_T rotMatYaw_4;
    real_T rotMatYaw_tmp;
    real_T rotMatYaw_tmp_0;
    real_T rotMatYaw_tmp_1;
    real_T rotMatYaw_tmp_tmp;
    real_T rotMatYaw_tmp_tmp_0;
    real_T tmp;
    rotMatYaw_tmp = std::sin(localB->TmpSignalConversionAtSFunctionI[5]);
    rotMatYaw_tmp_0 = std::cos(localB->TmpSignalConversionAtSFunctionI[5]);
    rotMatYaw[0] = rotMatYaw_tmp_0;
    rotMatYaw[3] = -rotMatYaw_tmp;
    rotMatYaw[6] = 0.0;
    rotMatYaw[1] = rotMatYaw_tmp;
    rotMatYaw[4] = rotMatYaw_tmp_0;
    rotMatYaw[7] = 0.0;
    rotMatYaw[2] = 0.0;
    rotMatYaw[5] = 0.0;
    rotMatYaw[8] = 1.0;
    rotMatYaw_tmp = 0.0;
    rotMatYaw_tmp_0 = 0.0;
    rotMatYaw_4 = 0.0;
    for (int32_T i = 0; i < 3; i++) {
      tmp = rtu_F_RadiationDampingLocal[i];
      rotMatYaw_tmp += rotMatYaw[3 * i] * tmp;
      rotMatYaw_tmp_0 += rotMatYaw[3 * i + 1] * tmp;
      rotMatYaw_4 += rotMatYaw[3 * i + 2] * tmp;
    }

    rotMatYaw_0[2] = rotMatYaw_4;
    rotMatYaw_0[1] = rotMatYaw_tmp_0;
    rotMatYaw_0[0] = rotMatYaw_tmp;
    rotMatYaw_tmp = 0.0;
    rotMatYaw_tmp_0 = 0.0;
    rotMatYaw_4 = 0.0;
    rotMatYaw_1 = 0.0;
    rotMatYaw_2 = 0.0;
    rotMatYaw_3 = 0.0;
    for (int32_T i = 0; i < 3; i++) {
      tmp = rtu_F_RadiationDampingLocal[i + 3];
      rotMatYaw_tmp_tmp = rotMatYaw[3 * i];
      rotMatYaw_tmp += rotMatYaw_tmp_tmp * tmp;
      rotMatYaw_tmp_tmp_0 = rotMatYaw[3 * i + 1];
      rotMatYaw_tmp_0 += rotMatYaw_tmp_tmp_0 * tmp;
      rotMatYaw_tmp_1 = rotMatYaw[3 * i + 2];
      rotMatYaw_4 += rotMatYaw_tmp_1 * tmp;
      localB->F_RadiationDamping[i] = rotMatYaw_0[i];
      tmp = rtu_F_AddedMassLocal[i];
      rotMatYaw_1 += rotMatYaw_tmp_tmp * tmp;
      rotMatYaw_2 += rotMatYaw_tmp_tmp_0 * tmp;
      rotMatYaw_3 += rotMatYaw_tmp_1 * tmp;
    }

    localB->F_RadiationDamping[3] = rotMatYaw_tmp;
    localB->F_RadiationDamping[4] = rotMatYaw_tmp_0;
    localB->F_RadiationDamping[5] = rotMatYaw_4;
    rotMatYaw_0[2] = rotMatYaw_3;
    rotMatYaw_0[1] = rotMatYaw_2;
    rotMatYaw_0[0] = rotMatYaw_1;
    rotMatYaw_tmp = 0.0;
    rotMatYaw_tmp_0 = 0.0;
    rotMatYaw_4 = 0.0;
    rotMatYaw_1 = 0.0;
    rotMatYaw_2 = 0.0;
    rotMatYaw_3 = 0.0;
    for (int32_T i = 0; i < 3; i++) {
      tmp = rtu_F_AddedMassLocal[i + 3];
      rotMatYaw_tmp_tmp = rotMatYaw[3 * i];
      rotMatYaw_tmp += rotMatYaw_tmp_tmp * tmp;
      rotMatYaw_tmp_tmp_0 = rotMatYaw[3 * i + 1];
      rotMatYaw_tmp_0 += rotMatYaw_tmp_tmp_0 * tmp;
      rotMatYaw_tmp_1 = rotMatYaw[3 * i + 2];
      rotMatYaw_4 += rotMatYaw_tmp_1 * tmp;
      localB->F_AddedMass[i] = rotMatYaw_0[i];
      tmp = rtu_F_ExcitationLocal[i];
      rotMatYaw_1 += rotMatYaw_tmp_tmp * tmp;
      rotMatYaw_2 += rotMatYaw_tmp_tmp_0 * tmp;
      rotMatYaw_3 += rotMatYaw_tmp_1 * tmp;
    }

    localB->F_AddedMass[3] = rotMatYaw_tmp;
    localB->F_AddedMass[4] = rotMatYaw_tmp_0;
    localB->F_AddedMass[5] = rotMatYaw_4;
    rotMatYaw_0[2] = rotMatYaw_3;
    rotMatYaw_0[1] = rotMatYaw_2;
    rotMatYaw_0[0] = rotMatYaw_1;
    rotMatYaw_tmp = 0.0;
    rotMatYaw_tmp_0 = 0.0;
    rotMatYaw_4 = 0.0;
    rotMatYaw_1 = 0.0;
    rotMatYaw_2 = 0.0;
    rotMatYaw_3 = 0.0;
    for (int32_T i = 0; i < 3; i++) {
      tmp = rtu_F_ExcitationLocal[i + 3];
      rotMatYaw_tmp_tmp = rotMatYaw[3 * i];
      rotMatYaw_tmp += rotMatYaw_tmp_tmp * tmp;
      rotMatYaw_tmp_tmp_0 = rotMatYaw[3 * i + 1];
      rotMatYaw_tmp_0 += rotMatYaw_tmp_tmp_0 * tmp;
      rotMatYaw_tmp_1 = rotMatYaw[3 * i + 2];
      rotMatYaw_4 += rotMatYaw_tmp_1 * tmp;
      localB->F_Excitation[i] = rotMatYaw_0[i];
      tmp = rtu_F_RestoringLocal[i];
      rotMatYaw_1 += rotMatYaw_tmp_tmp * tmp;
      rotMatYaw_2 += rotMatYaw_tmp_tmp_0 * tmp;
      rotMatYaw_3 += rotMatYaw_tmp_1 * tmp;
    }

    localB->F_Excitation[3] = rotMatYaw_tmp;
    localB->F_Excitation[4] = rotMatYaw_tmp_0;
    localB->F_Excitation[5] = rotMatYaw_4;
    rotMatYaw_0[2] = rotMatYaw_3;
    rotMatYaw_0[1] = rotMatYaw_2;
    rotMatYaw_0[0] = rotMatYaw_1;
    rotMatYaw_tmp = 0.0;
    rotMatYaw_tmp_0 = 0.0;
    rotMatYaw_4 = 0.0;
    for (int32_T i = 0; i < 3; i++) {
      tmp = rtu_F_RestoringLocal[i + 3];
      rotMatYaw_tmp += rotMatYaw[3 * i] * tmp;
      rotMatYaw_tmp_0 += rotMatYaw[3 * i + 1] * tmp;
      rotMatYaw_4 += rotMatYaw[3 * i + 2] * tmp;
      localB->F_Restoring[i] = rotMatYaw_0[i];
    }

    localB->F_Restoring[3] = rotMatYaw_tmp;
    localB->F_Restoring[4] = rotMatYaw_tmp_0;
    localB->F_Restoring[5] = rotMatYaw_4;
  } else {
    for (int32_T i = 0; i < 6; i++) {
      localB->F_RadiationDamping[i] = rtu_F_RadiationDampingLocal[i];
      localB->F_AddedMass[i] = rtu_F_AddedMassLocal[i];
      localB->F_Excitation[i] = rtu_F_ExcitationLocal[i];
      localB->F_Restoring[i] = rtu_F_RestoringLocal[i];
    }
  }
}

/*
 * System initialize for atomic system:
 *    '<S133>/Yaw Kinematic Transforms'
 *    '<S212>/Yaw Kinematic Transforms'
 */
void win_YawKinematicTransforms_Init(DW_YawKinematicTransforms_win_T *localDW)
{
  localDW->sfEvent = windEmulatorStep4__CALL_EVENT_d;
}

/*
 * Output and update for atomic system:
 *    '<S133>/Yaw Kinematic Transforms'
 *    '<S212>/Yaw Kinematic Transforms'
 */
void windEmul_YawKinematicTransforms(real_T rtu_yaw, const real_T
  rtu_dispGlobal[3], const real_T rtu_dispGlobal_k[3], const real_T
  rtu_velGlobal[3], const real_T rtu_velGlobal_i[3], const real_T rtu_accGlobal
  [6], B_YawKinematicTransforms_wind_T *localB, DW_YawKinematicTransforms_win_T *
  localDW)
{
  real_T rotMat[9];
  real_T rotMatYawTranspose[9];
  real_T rotMat_tmp_3[9];
  real_T rotMatYawTranspose_0[3];
  real_T rotMatYawTranspose_tmp;
  real_T rotMatYawTranspose_tmp_0;
  real_T rotMat_tmp;
  real_T rotMat_tmp_0;
  real_T rotMat_tmp_1;
  real_T rotMat_tmp_2;
  real_T rotMat_tmp_tmp;
  int32_T i;
  int32_T i_0;
  int32_T rotMat_tmp_4;
  int32_T rotMat_tmp_5;
  static const int8_T b[3] = { 0, 0, 1 };

  real_T rotMat_tmp_tmp_0;
  real_T rotMat_tmp_tmp_1;
  real_T rotMat_tmp_tmp_tmp;

  /* SignalConversion generated from: '<S134>/ SFunction ' */
  localB->TmpSignalConversionAtSFunctionI[0] = rtu_dispGlobal[0];
  localB->TmpSignalConversionAtSFunctionI[3] = rtu_dispGlobal_k[0];

  /* SignalConversion generated from: '<S134>/ SFunction ' */
  localB->TmpSignalConversionAtSFunctio_e[0] = rtu_velGlobal[0];
  localB->TmpSignalConversionAtSFunctio_e[3] = rtu_velGlobal_i[0];

  /* SignalConversion generated from: '<S134>/ SFunction ' */
  localB->TmpSignalConversionAtSFunctionI[1] = rtu_dispGlobal[1];
  localB->TmpSignalConversionAtSFunctionI[4] = rtu_dispGlobal_k[1];

  /* SignalConversion generated from: '<S134>/ SFunction ' */
  localB->TmpSignalConversionAtSFunctio_e[1] = rtu_velGlobal[1];
  localB->TmpSignalConversionAtSFunctio_e[4] = rtu_velGlobal_i[1];

  /* SignalConversion generated from: '<S134>/ SFunction ' */
  localB->TmpSignalConversionAtSFunctionI[2] = rtu_dispGlobal[2];
  localB->TmpSignalConversionAtSFunctionI[5] = rtu_dispGlobal_k[2];

  /* SignalConversion generated from: '<S134>/ SFunction ' */
  localB->TmpSignalConversionAtSFunctio_e[2] = rtu_velGlobal[2];
  localB->TmpSignalConversionAtSFunctio_e[5] = rtu_velGlobal_i[2];
  localDW->sfEvent = windEmulatorStep4__CALL_EVENT_d;
  if (rtu_yaw == 1.0) {
    rotMatYawTranspose_tmp = std::sin(localB->TmpSignalConversionAtSFunctionI[5]);
    rotMatYawTranspose_tmp_0 = std::cos(localB->TmpSignalConversionAtSFunctionI
      [5]);
    rotMatYawTranspose[0] = rotMatYawTranspose_tmp_0;
    rotMatYawTranspose[1] = -rotMatYawTranspose_tmp;
    rotMatYawTranspose[2] = 0.0;
    rotMatYawTranspose[3] = rotMatYawTranspose_tmp;
    rotMatYawTranspose[4] = rotMatYawTranspose_tmp_0;
    rotMatYawTranspose[5] = 0.0;
    rotMat_tmp = std::cos(localB->TmpSignalConversionAtSFunctionI[4]);
    rotMat_tmp_0 = std::sin(localB->TmpSignalConversionAtSFunctionI[4]);
    rotMat_tmp_1 = std::cos(localB->TmpSignalConversionAtSFunctionI[3]);
    rotMat_tmp_2 = std::sin(localB->TmpSignalConversionAtSFunctionI[3]);
    rotMat_tmp_3[0] = rotMat_tmp * rotMatYawTranspose_tmp_0;
    rotMat_tmp_3[3] = -rotMat_tmp * rotMatYawTranspose_tmp;
    rotMat_tmp_3[6] = rotMat_tmp_0;
    rotMat_tmp_tmp = rotMat_tmp_2 * rotMat_tmp_0;
    rotMat_tmp_3[1] = rotMat_tmp_tmp * rotMatYawTranspose_tmp_0 + rotMat_tmp_1 *
      rotMatYawTranspose_tmp;
    rotMat_tmp_3[4] = rotMat_tmp_1 * rotMatYawTranspose_tmp_0 - rotMat_tmp_tmp *
      rotMatYawTranspose_tmp;
    rotMat_tmp_3[7] = -rotMat_tmp_2 * rotMat_tmp;
    rotMat_tmp_tmp = rotMat_tmp_1 * rotMat_tmp_0;
    rotMat_tmp_3[2] = rotMat_tmp_2 * rotMatYawTranspose_tmp - rotMat_tmp_tmp *
      rotMatYawTranspose_tmp_0;
    rotMat_tmp_3[5] = rotMat_tmp_tmp * rotMatYawTranspose_tmp + rotMat_tmp_2 *
      rotMatYawTranspose_tmp_0;
    rotMat_tmp_3[8] = rotMat_tmp_1 * rotMat_tmp;
    for (i = 0; i < 3; i++) {
      rotMatYawTranspose[i + 6] = b[i];
      rotMat[3 * i] = 0.0;
      rotMat[3 * i + 1] = 0.0;
      rotMat[3 * i + 2] = 0.0;
    }

    rotMatYawTranspose_tmp_0 = 0.0;
    rotMat_tmp = 0.0;
    rotMat_tmp_0 = 0.0;
    rotMatYawTranspose_0[0] = 0.0;
    rotMatYawTranspose_0[1] = 0.0;
    rotMatYawTranspose_0[2] = 0.0;
    for (i = 0; i < 3; i++) {
      rotMat_tmp_1 = rotMat[3 * i];
      rotMat_tmp_5 = 3 * i + 1;
      rotMat_tmp_2 = rotMat[rotMat_tmp_5];
      rotMat_tmp_4 = 3 * i + 2;
      rotMat_tmp_tmp = rotMat[rotMat_tmp_4];
      for (i_0 = 0; i_0 < 3; i_0++) {
        rotMatYawTranspose_tmp = rotMat_tmp_3[3 * i + i_0];
        rotMat_tmp_1 += rotMatYawTranspose[3 * i_0] * rotMatYawTranspose_tmp;
        rotMat_tmp_2 += rotMatYawTranspose[3 * i_0 + 1] * rotMatYawTranspose_tmp;
        rotMat_tmp_tmp += rotMatYawTranspose[3 * i_0 + 2] *
          rotMatYawTranspose_tmp;
      }

      rotMat[rotMat_tmp_4] = rotMat_tmp_tmp;
      rotMat[rotMat_tmp_5] = rotMat_tmp_2;
      rotMat[3 * i] = rotMat_tmp_1;
      rotMatYawTranspose_tmp = localB->TmpSignalConversionAtSFunctionI[i];
      rotMat_tmp_1 = rotMatYawTranspose[3 * i];
      rotMatYawTranspose_tmp_0 += rotMat_tmp_1 * rotMatYawTranspose_tmp;
      rotMat_tmp_2 = rotMatYawTranspose[rotMat_tmp_5];
      rotMat_tmp += rotMat_tmp_2 * rotMatYawTranspose_tmp;
      rotMat_tmp_tmp = rotMatYawTranspose[rotMat_tmp_4];
      rotMat_tmp_0 += rotMat_tmp_tmp * rotMatYawTranspose_tmp;
      rotMatYawTranspose_tmp = localB->TmpSignalConversionAtSFunctio_e[i];
      rotMatYawTranspose_0[0] += rotMat_tmp_1 * rotMatYawTranspose_tmp;
      rotMatYawTranspose_0[1] += rotMat_tmp_2 * rotMatYawTranspose_tmp;
      rotMatYawTranspose_0[2] += rotMat_tmp_tmp * rotMatYawTranspose_tmp;
    }

    rotMatYawTranspose_tmp = rt_atan2d_snf(-rotMat[7], rotMat[8]);
    rotMat_tmp_1 = std::asin(rotMat[6]);
    rotMat_tmp_2 = rt_atan2d_snf(-rotMat[3], rotMat[0]);
    localB->dispLoc[3] = rotMatYawTranspose_tmp;
    localB->dispLoc[4] = rotMat_tmp_1;
    localB->dispLoc[5] = rotMat_tmp_2;
    localB->dispLoc[0] = rotMatYawTranspose_tmp_0;
    localB->dispLoc[1] = rotMat_tmp;
    localB->dispLoc[2] = rotMat_tmp_0;
    rotMat_tmp_1 = 0.0;
    rotMat_tmp_2 = 0.0;
    rotMat_tmp_tmp = 0.0;
    rotMatYawTranspose_tmp_0 = 0.0;
    rotMat_tmp = 0.0;
    rotMat_tmp_0 = 0.0;
    for (i = 0; i < 3; i++) {
      rotMatYawTranspose_tmp = localB->TmpSignalConversionAtSFunctio_e[i + 3];
      rotMat_tmp_tmp_0 = rotMatYawTranspose[3 * i];
      rotMat_tmp_1 += rotMat_tmp_tmp_0 * rotMatYawTranspose_tmp;
      rotMat_tmp_tmp_1 = rotMatYawTranspose[3 * i + 1];
      rotMat_tmp_2 += rotMat_tmp_tmp_1 * rotMatYawTranspose_tmp;
      rotMat_tmp_tmp_tmp = rotMatYawTranspose[3 * i + 2];
      rotMat_tmp_tmp += rotMat_tmp_tmp_tmp * rotMatYawTranspose_tmp;
      localB->velLoc[i] = rotMatYawTranspose_0[i];
      rotMatYawTranspose_tmp = rtu_accGlobal[i];
      rotMatYawTranspose_tmp_0 += rotMat_tmp_tmp_0 * rotMatYawTranspose_tmp;
      rotMat_tmp += rotMat_tmp_tmp_1 * rotMatYawTranspose_tmp;
      rotMat_tmp_0 += rotMat_tmp_tmp_tmp * rotMatYawTranspose_tmp;
    }

    localB->velLoc[3] = rotMat_tmp_1;
    localB->velLoc[4] = rotMat_tmp_2;
    localB->velLoc[5] = rotMat_tmp_tmp;
    rotMatYawTranspose_0[2] = rotMat_tmp_0;
    rotMatYawTranspose_0[1] = rotMat_tmp;
    rotMatYawTranspose_0[0] = rotMatYawTranspose_tmp_0;
    rotMat_tmp_1 = 0.0;
    rotMat_tmp_2 = 0.0;
    rotMat_tmp_tmp = 0.0;
    for (i = 0; i < 3; i++) {
      rotMatYawTranspose_tmp = rtu_accGlobal[i + 3];
      rotMat_tmp_1 += rotMatYawTranspose[3 * i] * rotMatYawTranspose_tmp;
      rotMat_tmp_2 += rotMatYawTranspose[3 * i + 1] * rotMatYawTranspose_tmp;
      rotMat_tmp_tmp += rotMatYawTranspose[3 * i + 2] * rotMatYawTranspose_tmp;
      localB->accLoc[i] = rotMatYawTranspose_0[i];
    }

    localB->accLoc[3] = rotMat_tmp_1;
    localB->accLoc[4] = rotMat_tmp_2;
    localB->accLoc[5] = rotMat_tmp_tmp;
  } else {
    for (i = 0; i < 6; i++) {
      localB->dispLoc[i] = localB->TmpSignalConversionAtSFunctionI[i];
      localB->velLoc[i] = localB->TmpSignalConversionAtSFunctio_e[i];
      localB->accLoc[i] = rtu_accGlobal[i];
    }
  }
}

/* Function for Chart: '<S18>/ABB Fieldbus Control' */
static void windEmulato_swParseStatusWord_m(void)
{
  uint16_T sw;
  sw = windEmulatorStep4_WECSim_B.BusAssignment_a.statusWord;
  windEmulatorStep4_WECSim_DW.swRDY_ON = 0.0;
  windEmulatorStep4_WECSim_DW.swRDY_RUN = 0.0;
  windEmulatorStep4_WECSim_DW.swRDY_REF = 0.0;
  windEmulatorStep4_WECSim_DW.swTRIPPED = 0.0;
  windEmulatorStep4_WECSim_DW.swOFF_2_STA = 0.0;
  windEmulatorStep4_WECSim_DW.swOFF_3_STA = 0.0;
  windEmulatorStep4_WECSim_DW.swSWC_ON_INHIB = 0.0;
  windEmulatorStep4_WECSim_DW.swWARNING = 0.0;
  windEmulatorStep4_WECSim_DW.swAT_SETPOINT = 0.0;
  windEmulatorStep4_WECSim_DW.swREMOTE = 0.0;
  windEmulatorStep4_WECSim_DW.swABOVE_LIMIT = 0.0;
  windEmulatorStep4_WECSim_DW.swEXT_CTRL_LOC = 0.0;
  windEmulatorStep4_WECSim_DW.swEXT_RUN_ENABLE = 0.0;
  windEmulatorStep4_WECSim_DW.swMSW_B13 = 0.0;
  windEmulatorStep4_WECSim_DW.swMSW_B14 = 0.0;
  windEmulatorStep4_WECSim_DW.swCOMM_ERR = 0.0;
  if (windEmulatorStep4_WECSim_B.BusAssignment_a.statusWord >= 32768) {
    uint32_T q0;
    windEmulatorStep4_WECSim_DW.swCOMM_ERR = 1.0;
    q0 = windEmulatorStep4_WECSim_B.BusAssignment_a.statusWord;
    q0 -= 32768U;
    sw = static_cast<uint16_T>(q0);
  }

  if (sw >= 16384) {
    windEmulatorStep4_WECSim_DW.swMSW_B14 = 1.0;
    sw = static_cast<uint16_T>(sw - 16384);
  }

  if (sw >= 8192) {
    windEmulatorStep4_WECSim_DW.swMSW_B13 = 1.0;
    sw = static_cast<uint16_T>(sw - 8192);
  }

  if (sw >= 4096) {
    windEmulatorStep4_WECSim_DW.swEXT_RUN_ENABLE = 1.0;
    sw = static_cast<uint16_T>(sw - 4096);
  }

  if (sw >= 2048) {
    windEmulatorStep4_WECSim_DW.swEXT_CTRL_LOC = 1.0;
    sw = static_cast<uint16_T>(sw - 2048);
  }

  if (sw >= 1024) {
    windEmulatorStep4_WECSim_DW.swABOVE_LIMIT = 1.0;
    sw = static_cast<uint16_T>(sw - 1024);
  }

  if (sw >= 512) {
    windEmulatorStep4_WECSim_DW.swREMOTE = 1.0;
    sw = static_cast<uint16_T>(sw - 512);
  }

  if (sw >= 256) {
    windEmulatorStep4_WECSim_DW.swAT_SETPOINT = 1.0;
    sw = static_cast<uint16_T>(sw - 256);
  }

  if (sw >= 128) {
    windEmulatorStep4_WECSim_DW.swWARNING = 1.0;
    sw = static_cast<uint16_T>(sw - 128);
  }

  if (sw >= 64) {
    windEmulatorStep4_WECSim_DW.swSWC_ON_INHIB = 1.0;
    sw = static_cast<uint16_T>(sw - 64);
  }

  if (sw >= 32) {
    windEmulatorStep4_WECSim_DW.swOFF_3_STA = 1.0;
    sw = static_cast<uint16_T>(sw - 32);
  }

  if (sw >= 16) {
    windEmulatorStep4_WECSim_DW.swOFF_2_STA = 1.0;
    sw = static_cast<uint16_T>(sw - 16);
  }

  if (sw >= 8) {
    windEmulatorStep4_WECSim_DW.swTRIPPED = 1.0;
    sw = static_cast<uint16_T>(sw - 8);
  }

  if (sw >= 4) {
    windEmulatorStep4_WECSim_DW.swRDY_REF = 1.0;
    sw = static_cast<uint16_T>(sw - 4);
  }

  if (sw >= 2) {
    windEmulatorStep4_WECSim_DW.swRDY_RUN = 1.0;
    sw = static_cast<uint16_T>(sw - 2);
  }

  if (sw >= 1) {
    windEmulatorStep4_WECSim_DW.swRDY_ON = 1.0;
  }
}

/* Function for Chart: '<S18>/ABB Fieldbus Control' */
static void windEmulat_cwBuildControlWord_o(void)
{
  windEmulatorStep4_WECSim_B.ControlWord = 0U;
  if (windEmulatorStep4_WECSim_DW.cwOFF1_CONTROL != 0.0) {
    windEmulatorStep4_WECSim_B.ControlWord = 1U;
  }

  if (windEmulatorStep4_WECSim_DW.cwOFF2_CONTROL != 0.0) {
    windEmulatorStep4_WECSim_B.ControlWord = static_cast<uint16_T>
      (windEmulatorStep4_WECSim_B.ControlWord | 2);
  }

  if (windEmulatorStep4_WECSim_DW.cwOFF3_CONTROL != 0.0) {
    windEmulatorStep4_WECSim_B.ControlWord = static_cast<uint16_T>
      (windEmulatorStep4_WECSim_B.ControlWord | 4);
  }

  if (windEmulatorStep4_WECSim_DW.cwENABLE_OPERATION != 0.0) {
    windEmulatorStep4_WECSim_B.ControlWord = static_cast<uint16_T>
      (windEmulatorStep4_WECSim_B.ControlWord | 8);
  }

  if (windEmulatorStep4_WECSim_DW.cwRAMP_OUT_ZERO != 0.0) {
    windEmulatorStep4_WECSim_B.ControlWord = static_cast<uint16_T>
      (windEmulatorStep4_WECSim_B.ControlWord | 16);
  }

  if (windEmulatorStep4_WECSim_DW.cwRAMP_HOLD != 0.0) {
    windEmulatorStep4_WECSim_B.ControlWord = static_cast<uint16_T>
      (windEmulatorStep4_WECSim_B.ControlWord | 32);
  }

  if (windEmulatorStep4_WECSim_DW.cwRAMP_IN_ZERO != 0.0) {
    windEmulatorStep4_WECSim_B.ControlWord = static_cast<uint16_T>
      (windEmulatorStep4_WECSim_B.ControlWord | 64);
  }

  if (windEmulatorStep4_WECSim_DW.cwRESET != 0.0) {
    windEmulatorStep4_WECSim_B.ControlWord = static_cast<uint16_T>
      (windEmulatorStep4_WECSim_B.ControlWord | 128);
  }

  if (windEmulatorStep4_WECSim_DW.cwREMOTE_CMD != 0.0) {
    windEmulatorStep4_WECSim_B.ControlWord = static_cast<uint16_T>
      (windEmulatorStep4_WECSim_B.ControlWord | 1024);
  }

  if (windEmulatorStep4_WECSim_B.ACS880CtrlMode) {
    windEmulatorStep4_WECSim_B.ControlWord = static_cast<uint16_T>
      (windEmulatorStep4_WECSim_B.ControlWord | 2048);
  }
}

/* Function for Chart: '<S18>/ABB Fieldbus Control' */
static void windEmulatorStep_cwInitialize_a(void)
{
  windEmulatorStep4_WECSim_DW.cwOFF1_CONTROL = 0.0;
  windEmulatorStep4_WECSim_DW.cwOFF2_CONTROL = 1.0;
  windEmulatorStep4_WECSim_DW.cwOFF3_CONTROL = 1.0;
  windEmulatorStep4_WECSim_DW.cwENABLE_OPERATION = 0.0;
  windEmulatorStep4_WECSim_DW.cwRAMP_OUT_ZERO = 1.0;
  windEmulatorStep4_WECSim_DW.cwRAMP_HOLD = 1.0;
  windEmulatorStep4_WECSim_DW.cwRAMP_IN_ZERO = 1.0;
  windEmulatorStep4_WECSim_DW.cwRESET = 0.0;
  windEmulatorStep4_WECSim_DW.cwREMOTE_CMD = 1.0;
}

/* Function for Chart: '<S16>/ABB Fieldbus Control' */
static void windEmulatorS_swParseStatusWord(void)
{
  uint16_T sw;
  sw = windEmulatorStep4_WECSim_B.BusAssignment_h.statusWord;
  windEmulatorStep4_WECSim_DW.swRDY_ON_f = 0.0;
  windEmulatorStep4_WECSim_DW.swRDY_RUN_f = 0.0;
  windEmulatorStep4_WECSim_DW.swRDY_REF_a = 0.0;
  windEmulatorStep4_WECSim_DW.swTRIPPED_e = 0.0;
  windEmulatorStep4_WECSim_DW.swOFF_2_STA_c = 0.0;
  windEmulatorStep4_WECSim_DW.swOFF_3_STA_i = 0.0;
  windEmulatorStep4_WECSim_DW.swSWC_ON_INHIB_l = 0.0;
  windEmulatorStep4_WECSim_DW.swWARNING_k = 0.0;
  windEmulatorStep4_WECSim_DW.swAT_SETPOINT_b = 0.0;
  windEmulatorStep4_WECSim_DW.swREMOTE_j = 0.0;
  windEmulatorStep4_WECSim_DW.swABOVE_LIMIT_a = 0.0;
  windEmulatorStep4_WECSim_DW.swEXT_CTRL_LOC_c = 0.0;
  windEmulatorStep4_WECSim_DW.swEXT_RUN_ENABLE_h = 0.0;
  windEmulatorStep4_WECSim_DW.swMSW_B13_l = 0.0;
  windEmulatorStep4_WECSim_DW.swMSW_B14_i = 0.0;
  windEmulatorStep4_WECSim_DW.swCOMM_ERR_p = 0.0;
  if (windEmulatorStep4_WECSim_B.BusAssignment_h.statusWord >= 32768) {
    uint32_T q0;
    windEmulatorStep4_WECSim_DW.swCOMM_ERR_p = 1.0;
    q0 = windEmulatorStep4_WECSim_B.BusAssignment_h.statusWord;
    q0 -= 32768U;
    sw = static_cast<uint16_T>(q0);
  }

  if (sw >= 16384) {
    windEmulatorStep4_WECSim_DW.swMSW_B14_i = 1.0;
    sw = static_cast<uint16_T>(sw - 16384);
  }

  if (sw >= 8192) {
    windEmulatorStep4_WECSim_DW.swMSW_B13_l = 1.0;
    sw = static_cast<uint16_T>(sw - 8192);
  }

  if (sw >= 4096) {
    windEmulatorStep4_WECSim_DW.swEXT_RUN_ENABLE_h = 1.0;
    sw = static_cast<uint16_T>(sw - 4096);
  }

  if (sw >= 2048) {
    windEmulatorStep4_WECSim_DW.swEXT_CTRL_LOC_c = 1.0;
    sw = static_cast<uint16_T>(sw - 2048);
  }

  if (sw >= 1024) {
    windEmulatorStep4_WECSim_DW.swABOVE_LIMIT_a = 1.0;
    sw = static_cast<uint16_T>(sw - 1024);
  }

  if (sw >= 512) {
    windEmulatorStep4_WECSim_DW.swREMOTE_j = 1.0;
    sw = static_cast<uint16_T>(sw - 512);
  }

  if (sw >= 256) {
    windEmulatorStep4_WECSim_DW.swAT_SETPOINT_b = 1.0;
    sw = static_cast<uint16_T>(sw - 256);
  }

  if (sw >= 128) {
    windEmulatorStep4_WECSim_DW.swWARNING_k = 1.0;
    sw = static_cast<uint16_T>(sw - 128);
  }

  if (sw >= 64) {
    windEmulatorStep4_WECSim_DW.swSWC_ON_INHIB_l = 1.0;
    sw = static_cast<uint16_T>(sw - 64);
  }

  if (sw >= 32) {
    windEmulatorStep4_WECSim_DW.swOFF_3_STA_i = 1.0;
    sw = static_cast<uint16_T>(sw - 32);
  }

  if (sw >= 16) {
    windEmulatorStep4_WECSim_DW.swOFF_2_STA_c = 1.0;
    sw = static_cast<uint16_T>(sw - 16);
  }

  if (sw >= 8) {
    windEmulatorStep4_WECSim_DW.swTRIPPED_e = 1.0;
    sw = static_cast<uint16_T>(sw - 8);
  }

  if (sw >= 4) {
    windEmulatorStep4_WECSim_DW.swRDY_REF_a = 1.0;
    sw = static_cast<uint16_T>(sw - 4);
  }

  if (sw >= 2) {
    windEmulatorStep4_WECSim_DW.swRDY_RUN_f = 1.0;
    sw = static_cast<uint16_T>(sw - 2);
  }

  if (sw >= 1) {
    windEmulatorStep4_WECSim_DW.swRDY_ON_f = 1.0;
  }
}

/* Function for Chart: '<S16>/ABB Fieldbus Control' */
static void windEmulator_cwBuildControlWord(void)
{
  windEmulatorStep4_WECSim_B.ControlWord_l = 0U;
  if (windEmulatorStep4_WECSim_DW.cwOFF1_CONTROL_i != 0.0) {
    windEmulatorStep4_WECSim_B.ControlWord_l = 1U;
  }

  if (windEmulatorStep4_WECSim_DW.cwOFF2_CONTROL_b != 0.0) {
    windEmulatorStep4_WECSim_B.ControlWord_l = static_cast<uint16_T>
      (windEmulatorStep4_WECSim_B.ControlWord_l | 2);
  }

  if (windEmulatorStep4_WECSim_DW.cwOFF3_CONTROL_n != 0.0) {
    windEmulatorStep4_WECSim_B.ControlWord_l = static_cast<uint16_T>
      (windEmulatorStep4_WECSim_B.ControlWord_l | 4);
  }

  if (windEmulatorStep4_WECSim_DW.cwENABLE_OPERATION_l != 0.0) {
    windEmulatorStep4_WECSim_B.ControlWord_l = static_cast<uint16_T>
      (windEmulatorStep4_WECSim_B.ControlWord_l | 8);
  }

  if (windEmulatorStep4_WECSim_DW.cwRAMP_OUT_ZERO_m != 0.0) {
    windEmulatorStep4_WECSim_B.ControlWord_l = static_cast<uint16_T>
      (windEmulatorStep4_WECSim_B.ControlWord_l | 16);
  }

  if (windEmulatorStep4_WECSim_DW.cwRAMP_HOLD_j != 0.0) {
    windEmulatorStep4_WECSim_B.ControlWord_l = static_cast<uint16_T>
      (windEmulatorStep4_WECSim_B.ControlWord_l | 32);
  }

  if (windEmulatorStep4_WECSim_DW.cwRAMP_IN_ZERO_m != 0.0) {
    windEmulatorStep4_WECSim_B.ControlWord_l = static_cast<uint16_T>
      (windEmulatorStep4_WECSim_B.ControlWord_l | 64);
  }

  if (windEmulatorStep4_WECSim_DW.cwRESET_e != 0.0) {
    windEmulatorStep4_WECSim_B.ControlWord_l = static_cast<uint16_T>
      (windEmulatorStep4_WECSim_B.ControlWord_l | 128);
  }

  if (windEmulatorStep4_WECSim_DW.cwREMOTE_CMD_c != 0.0) {
    windEmulatorStep4_WECSim_B.ControlWord_l = static_cast<uint16_T>
      (windEmulatorStep4_WECSim_B.ControlWord_l | 1024);
  }

  if (windEmulatorStep4_WECSim_B.ACS800CtrlMode) {
    windEmulatorStep4_WECSim_B.ControlWord_l = static_cast<uint16_T>
      (windEmulatorStep4_WECSim_B.ControlWord_l | 2048);
  }
}

/* Function for Chart: '<S16>/ABB Fieldbus Control' */
static void windEmulatorStep4__cwInitialize(void)
{
  windEmulatorStep4_WECSim_DW.cwOFF1_CONTROL_i = 0.0;
  windEmulatorStep4_WECSim_DW.cwOFF2_CONTROL_b = 1.0;
  windEmulatorStep4_WECSim_DW.cwOFF3_CONTROL_n = 1.0;
  windEmulatorStep4_WECSim_DW.cwENABLE_OPERATION_l = 0.0;
  windEmulatorStep4_WECSim_DW.cwRAMP_OUT_ZERO_m = 1.0;
  windEmulatorStep4_WECSim_DW.cwRAMP_HOLD_j = 1.0;
  windEmulatorStep4_WECSim_DW.cwRAMP_IN_ZERO_m = 1.0;
  windEmulatorStep4_WECSim_DW.cwRESET_e = 0.0;
  windEmulatorStep4_WECSim_DW.cwREMOTE_CMD_c = 1.0;
}

/* Model step function */
void windEmulatorStep4_WECSim_step(void)
{
  if (rtmIsMajorTimeStep(windEmulatorStep4_WECSim_M)) {
    /* set solver stop time */
    if (!(windEmulatorStep4_WECSim_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&windEmulatorStep4_WECSim_M->solverInfo,
                            ((windEmulatorStep4_WECSim_M->Timing.clockTickH0 + 1)
        * windEmulatorStep4_WECSim_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&windEmulatorStep4_WECSim_M->solverInfo,
                            ((windEmulatorStep4_WECSim_M->Timing.clockTick0 + 1)
        * windEmulatorStep4_WECSim_M->Timing.stepSize0 +
        windEmulatorStep4_WECSim_M->Timing.clockTickH0 *
        windEmulatorStep4_WECSim_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(windEmulatorStep4_WECSim_M)) {
    windEmulatorStep4_WECSim_M->Timing.t[0] = rtsiGetT
      (&windEmulatorStep4_WECSim_M->solverInfo);
  }

  {
    NeParameterBundle expl_temp;
    NeslRtpManager *rtpManager;
    NeslSimulationData *simulationData;
    NeslSimulator *simulator;
    NeuDiagnosticManager *diag;
    NeuDiagnosticTree *diagTree;
    char *msg;
    struct_KY1U3Kyrwv5e6VnUIBWG5G *tmp_m;
    struct_KY1U3Kyrwv5e6VnUIBWG5G *tmp_s;
    struct_KY1U3Kyrwv5e6VnUIBWG5G *tmp_t;
    real_T tmp_a[74];
    real_T tmp_2[57];
    real_T tmp_6[54];
    real_T tmp_c[54];
    real_T tmp_4[52];
    real_T tmp_0[20];
    real_T tmp_8[12];
    real_T tmp_e[6];
    real_T tmp_f[6];
    real_T tmp[3];
    const real_T *tmp_i;
    real_T Clock_tmp;
    real_T deltaT_tmp;
    real_T rateLimiterRate;
    real_T riseValLimit;
    real_T time;
    real_T time_0;
    real_T time_1;
    real_T time_2;
    real_T time_3;
    real_T time_4;
    real_T time_5;
    real_T time_6;
    real_T time_7;
    real_T time_8;
    real_T time_9;
    real_T time_a;
    real_T time_b;
    real_T time_c;
    real_T tmp_10;
    real_T tmp_j;
    real_T tmp_k;
    real_T tmp_l;
    real_T tmp_n;
    real_T tmp_o;
    real_T tmp_p;
    real_T tmp_q;
    real_T tmp_r;
    real_T tmp_v;
    real_T tmp_w;
    real_T tmp_x;
    real_T tmp_y;
    real_T tmp_z;
    real_T u1;
    real_T *parameterBundle_mRealParameters;
    int32_T i;
    int32_T i_0;
    int32_T isHit;
    int32_T isHit_0;
    int_T tmp_7[15];
    int_T tmp_d[15];
    int_T tmp_5[14];
    int_T tmp_3[7];
    int_T tmp_1[6];
    int_T tmp_b[5];
    int_T tmp_9[4];
    uint32_T q0;
    uint32_T qY;
    boolean_T f;
    boolean_T tmp_11;
    boolean_T tmp_g;
    boolean_T tmp_h;
    boolean_T tmp_u;

    /* Constant: '<S1>/ACS800CtrlMode' */
    tmp_11 = *get_ctrlModeTorque();

    /* Gain: '<S12>/rad//s->rpm' */
    tmp_10 = *get_radps2rpm();

    /* Gain: '<S230>/Gain' */
    tmp_z = *get_m3persecond2lpm();

    /* Gain: '<S245>/Gain' */
    tmp_y = *get_pa2psi();

    /* Gain: '<S439>/Gain' */
    tmp_x = *get_TorqueLoadMax();

    /* Switch: '<S598>/Switch' incorporates:
     *  Constant: '<S13>/Constant2'
     */
    Clock_tmp = *get_acs880SpeedPILimLo();

    /* Switch: '<S598>/Switch2' incorporates:
     *  Constant: '<S13>/Constant1'
     */
    tmp_w = *get_acs880SpeedPILimUp();

    /* Switch generated from: '<S365>/Switch' */
    tmp_v = *get_minTorqueRef_Nm();

    /* Chart: '<S4>/FexcRamp' */
    u1 = *get_rampTime();

    /* Switch: '<S373>/Switch' incorporates:
     *  Constant: '<S373>/DeadBandController'
     *  Switch generated from: '<S365>/Switch'
     */
    tmp_u = *get_deadBandController();

    /* Gain: '<S437>/Gain' incorporates:
     *  Switch generated from: '<S365>/Switch'
     */
    tmp_t = get_TorqueInputControl();

    /* Gain: '<S372>/Gain' incorporates:
     *  Switch generated from: '<S365>/Switch'
     */
    tmp_s = get_PressureControl();

    /* RateLimiter: '<S609>/torqueSlewRate' */
    tmp_r = *get_fromFileTorqueSlewRate();

    /* RateLimiter: '<S609>/speedSlewRate' */
    tmp_q = *get_fromFileSpeedSlewRate();

    /* RateLimiter: '<S498>/Rate Limiter' incorporates:
     *  Constant: '<S498>/Constant1'
     */
    tmp_p = *get_rpm2radps();

    /* Gain: '<S365>/Gain' incorporates:
     *  Switch: '<S365>/Switch2'
     */
    tmp_o = *get_Dm_max();

    /* Gain: '<S1>/Nm -> %' */
    tmp_n = *get_acs880RatedTorque();

    /* Product: '<S494>/Product' incorporates:
     *  Constant: '<S494>/Constant1'
     */
    tmp_m = get_SpeedControl();

    /* Product: '<S495>/Product' incorporates:
     *  Constant: '<S495>/Constant1'
     */
    tmp_l = *get_belowMinPGain();

    /* RateLimiter: '<S436>/Rate Limiter' */
    tmp_k = *get_deadbandTorqueSlewRate();

    /* Switch: '<S483>/Switch' incorporates:
     *  Constant: '<S438>/Constant1'
     */
    tmp_j = *get_genMaxTorque();

    {
      /* user code (Output function Header) */
      {
        /*------------ S-Function Block: <Root>/EtherCAT Init Process Received Frames ------------*/
        unsigned int data[6]= { 0 };

        int32_T msdata[4] = { 0 };

        xpcEtherCATReadProcessData(0,NULL);
        mwErrorGet((int_T)0,
                   &data[0], &data[1], &data[2], (int *)&data[3], &data[4], (int
                    *)&data[5]);
        memcpy(&windEmulatorStep4_WECSim_B.EtherCATInit[0], data,6*sizeof
               (int32_T));
        mwErrorClear( (int_T)0 );

        // Clear all momentary triggered values
      }

      /* SimscapeInputBlock: '<S541>/INPUT_1_1_1' incorporates:
       *  SimscapeInputBlock: '<S541>/INPUT_5_1_1'
       */
      tmp_g = rtmIsMajorTimeStep(windEmulatorStep4_WECSim_M);
      if (tmp_g) {
        /* S-Function (slecatpdorx): '<S10>/EtherCAT PDO Receive1' */
        {
          /*------------ S-Function Block: <S10>/EtherCAT PDO Receive1 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.ACS880MotorVoltage;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 440;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 4 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (16 == 8) && (2 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (16 == 8) && (2 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                 sigIdx*2, 2);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 16;
          }
        }

        /* DataTypeConversion: '<S10>/Cast To Double' */
        windEmulatorStep4_WECSim_B.CastToDouble =
          windEmulatorStep4_WECSim_B.ACS880MotorVoltage;

        /* Gain: '<S10>/Gain' */
        windEmulatorStep4_WECSim_B.Gain = *get_acs880MotorVoltsScaling() *
          windEmulatorStep4_WECSim_B.CastToDouble;

        /* S-Function (slecatpdorx): '<S10>/EtherCAT PDO Receive2' */
        {
          /*------------ S-Function Block: <S10>/EtherCAT PDO Receive2 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.ACS880MotorCurrent;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 424;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 4 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (16 == 8) && (2 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (16 == 8) && (2 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                 sigIdx*2, 2);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 16;
          }
        }

        /* DataTypeConversion: '<S10>/Cast To Double1' */
        windEmulatorStep4_WECSim_B.CastToDouble1 =
          windEmulatorStep4_WECSim_B.ACS880MotorCurrent;

        /* Gain: '<S10>/Gain2' */
        windEmulatorStep4_WECSim_B.Gain2 = *get_acs880MotorCurrentScaling() *
          windEmulatorStep4_WECSim_B.CastToDouble1;

        /* S-Function (slecatpdorx): '<S10>/EtherCAT PDO Receive3' */
        {
          /*------------ S-Function Block: <S10>/EtherCAT PDO Receive3 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.ACS880OutputFreq;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 408;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 4 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (16 == 8) && (2 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (16 == 8) && (2 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                 sigIdx*2, 2);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 16;
          }
        }

        /* DataTypeConversion: '<S10>/Cast To Double2' */
        windEmulatorStep4_WECSim_B.CastToDouble2 =
          windEmulatorStep4_WECSim_B.ACS880OutputFreq;

        /* Gain: '<S10>/Gain3' */
        windEmulatorStep4_WECSim_B.Gain3 = *get_acs880FreqScaling() *
          windEmulatorStep4_WECSim_B.CastToDouble2;

        /* S-Function (slecatpdorx): '<S10>/EtherCAT PDO Receive4' */
        {
          /*------------ S-Function Block: <S10>/EtherCAT PDO Receive4 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.ACS880MotorSpeed;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 392;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 4 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (16 == 8) && (2 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (16 == 8) && (2 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                 sigIdx*2, 2);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 16;
          }
        }

        /* DataTypeConversion: '<S10>/Cast To Double3' */
        windEmulatorStep4_WECSim_B.CastToDouble3 =
          windEmulatorStep4_WECSim_B.ACS880MotorSpeed;

        /* Gain: '<S10>/Gain4' */
        windEmulatorStep4_WECSim_B.Gain4 = *get_acs880SpeedScaling() *
          windEmulatorStep4_WECSim_B.CastToDouble3;

        /* S-Function (slecatpdorx): '<S10>/EtherCAT PDO Receive5' */
        {
          /*------------ S-Function Block: <S10>/EtherCAT PDO Receive5 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.ACS880MotorTorque;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 456;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 4 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (16 == 8) && (2 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (16 == 8) && (2 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                 sigIdx*2, 2);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 16;
          }
        }

        /* DataTypeConversion: '<S10>/Data Type Conversion1' */
        windEmulatorStep4_WECSim_B.DataTypeConversion1 =
          windEmulatorStep4_WECSim_B.ACS880MotorTorque;

        /* Gain: '<S10>/Gain1' */
        windEmulatorStep4_WECSim_B.Gain1 = *get_acs880TorqueScaling() *
          windEmulatorStep4_WECSim_B.DataTypeConversion1;

        /* S-Function (slecatpdorx): '<S10>/EtherCAT PDO Receive6' */
        {
          /*------------ S-Function Block: <S10>/EtherCAT PDO Receive6 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.ACS880MotorShaftPower;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 472;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 4 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (16 == 8) && (2 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (16 == 8) && (2 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                 sigIdx*2, 2);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 16;
          }
        }

        /* DataTypeConversion: '<S10>/Data Type Conversion2' */
        windEmulatorStep4_WECSim_B.DataTypeConversion2 =
          windEmulatorStep4_WECSim_B.ACS880MotorShaftPower;

        /* Gain: '<S10>/Gain5' */
        windEmulatorStep4_WECSim_B.Gain5 = *get_acs880PowerScaling() *
          windEmulatorStep4_WECSim_B.DataTypeConversion2;

        /* S-Function (slecatpdorx): '<S10>/EtherCAT PDO Receive' */
        {
          /*------------ S-Function Block: <S10>/EtherCAT PDO Receive PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.EtherCATPDOReceive;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 344;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 5 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (16 == 8) && (2 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (16 == 8) && (2 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                 sigIdx*2, 2);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 16;
          }
        }

        /* BusAssignment: '<S10>/Bus Assignment' */
        windEmulatorStep4_WECSim_B.BusAssignment_a.motorVoltage_V =
          windEmulatorStep4_WECSim_B.Gain;
        windEmulatorStep4_WECSim_B.BusAssignment_a.motorCurrent_A =
          windEmulatorStep4_WECSim_B.Gain2;
        windEmulatorStep4_WECSim_B.BusAssignment_a.frequency_Hz =
          windEmulatorStep4_WECSim_B.Gain3;
        windEmulatorStep4_WECSim_B.BusAssignment_a.motorSpeed_rpm =
          windEmulatorStep4_WECSim_B.Gain4;
        windEmulatorStep4_WECSim_B.BusAssignment_a.motorTorque_Nm =
          windEmulatorStep4_WECSim_B.Gain1;
        windEmulatorStep4_WECSim_B.BusAssignment_a.shaftPower_W =
          windEmulatorStep4_WECSim_B.Gain5;
        windEmulatorStep4_WECSim_B.BusAssignment_a.statusWord =
          windEmulatorStep4_WECSim_B.EtherCATPDOReceive;

        /* ToAsyncQueueBlock generated from: '<S27>/acs880Signals' */
        slrtLogSignal
          (windEmulatorStep4_WECSim_DW.TAQSigLogging_InsertedFor_acs88.SLRTSigHandles,
           (((windEmulatorStep4_WECSim_M->Timing.clockTick1+
              windEmulatorStep4_WECSim_M->Timing.clockTickH1* 4294967296.0)) *
            0.004));

        /* Gain: '<S48>/rpm -> rad//s' */
        windEmulatorStep4_WECSim_B.rpmrads = tmp_p *
          windEmulatorStep4_WECSim_B.BusAssignment_a.motorSpeed_rpm;

        /* Product: '<S48>/Product' */
        windEmulatorStep4_WECSim_B.Product = windEmulatorStep4_WECSim_B.rpmrads *
          windEmulatorStep4_WECSim_B.BusAssignment_a.motorTorque_Nm;
        windEmulatorStep4_MovingAverage(windEmulatorStep4_WECSim_B.Product,
          &windEmulatorStep4_WECSim_B.MovingAverage_p,
          &windEmulatorStep4_WECSim_DW.MovingAverage_p);

        /* Gain: '<S48>/shaftPowerAverage_W' */
        windEmulatorStep4_WECSim_B.shaftPowerAverage_W =
          windEmulatorStep4_WECSim_cal->shaftPowerAverage_W_Gain *
          windEmulatorStep4_WECSim_B.MovingAverage_p.MovingAverage;

        /* Gain: '<S48>/shaftPower_W' */
        windEmulatorStep4_WECSim_B.shaftPower_W =
          windEmulatorStep4_WECSim_cal->shaftPower_W_Gain *
          windEmulatorStep4_WECSim_B.Product;

        /* MATLAB Function: '<S49>/Parse Status Word' */
        windEmulatorSte_ParseStatusWord
          (windEmulatorStep4_WECSim_B.BusAssignment_a.statusWord,
           &windEmulatorStep4_WECSim_B.sf_ParseStatusWord_h,
           &windEmulatorStep4_WECSim_DW.sf_ParseStatusWord_h);

        /* Logic: '<S49>/aboveLimit' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_WECSim_B.aboveLimit =
          (windEmulatorStep4_WECSim_cal->Constant_Value_ov &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord_h.above_limit != 0));

        /* Logic: '<S49>/alarm' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_WECSim_B.alarm =
          (windEmulatorStep4_WECSim_cal->Constant_Value_ov &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord_h.alarm != 0));

        /* Logic: '<S49>/atSetpoint' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_WECSim_B.atSetpoint =
          (windEmulatorStep4_WECSim_cal->Constant_Value_ov &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord_h.at_setpoint != 0));

        /* Logic: '<S49>/commErr' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_WECSim_B.commErr =
          (windEmulatorStep4_WECSim_cal->Constant_Value_ov &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord_h.comm_err != 0));

        /* Logic: '<S49>/extCtrlLoc' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_WECSim_B.extCtrlLoc =
          (windEmulatorStep4_WECSim_cal->Constant_Value_ov &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord_h.ext_ctrl_loc != 0));

        /* Logic: '<S49>/extRunEnable' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_WECSim_B.extRunEnable =
          (windEmulatorStep4_WECSim_cal->Constant_Value_ov &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord_h.ext_run_enable != 0));

        /* Logic: '<S49>/mswB13' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_WECSim_B.mswB13 =
          (windEmulatorStep4_WECSim_cal->Constant_Value_ov &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord_h.msw_b13 != 0));

        /* Logic: '<S49>/mswB14' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_WECSim_B.mswB14 =
          (windEmulatorStep4_WECSim_cal->Constant_Value_ov &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord_h.msw_b14 != 0));

        /* Logic: '<S49>/off2' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_WECSim_B.off2 =
          (windEmulatorStep4_WECSim_cal->Constant_Value_ov &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord_h.off2 != 0));

        /* Logic: '<S49>/off3' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_WECSim_B.off3 =
          (windEmulatorStep4_WECSim_cal->Constant_Value_ov &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord_h.off3 != 0));

        /* Logic: '<S49>/rdyOn' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_WECSim_B.rdyOn =
          (windEmulatorStep4_WECSim_cal->Constant_Value_ov &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord_h.rdy_on != 0));

        /* Logic: '<S49>/rdyRef' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_WECSim_B.rdyRef =
          (windEmulatorStep4_WECSim_cal->Constant_Value_ov &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord_h.rdy_ref != 0));

        /* Logic: '<S49>/rdyRun' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_WECSim_B.rdyRun =
          (windEmulatorStep4_WECSim_cal->Constant_Value_ov &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord_h.rdy_run != 0));

        /* Logic: '<S49>/remote' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_WECSim_B.remote =
          (windEmulatorStep4_WECSim_cal->Constant_Value_ov &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord_h.remote != 0));

        /* Logic: '<S49>/switchOnInhibit' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_WECSim_B.switchOnInhibit =
          (windEmulatorStep4_WECSim_cal->Constant_Value_ov &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord_h.swc_on_inhib != 0));

        /* Logic: '<S49>/tripped' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_WECSim_B.tripped =
          (windEmulatorStep4_WECSim_cal->Constant_Value_ov &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord_h.tripped != 0));

        /* Gain: '<S28>/frequency_Hz' */
        windEmulatorStep4_WECSim_B.frequency_Hz =
          windEmulatorStep4_WECSim_cal->frequency_Hz_Gain *
          windEmulatorStep4_WECSim_B.BusAssignment_a.frequency_Hz;

        /* Gain: '<S28>/motorCurrent_A' */
        windEmulatorStep4_WECSim_B.motorCurrent_A =
          windEmulatorStep4_WECSim_cal->motorCurrent_A_Gain *
          windEmulatorStep4_WECSim_B.BusAssignment_a.motorCurrent_A;

        /* Gain: '<S28>/motorSpeed_rpm' */
        windEmulatorStep4_WECSim_B.motorSpeed_rpm =
          windEmulatorStep4_WECSim_cal->motorSpeed_rpm_Gain *
          windEmulatorStep4_WECSim_B.BusAssignment_a.motorSpeed_rpm;

        /* Gain: '<S28>/motorTorque_Nm' */
        windEmulatorStep4_WECSim_B.motorTorque_Nm =
          windEmulatorStep4_WECSim_cal->motorTorque_Nm_Gain *
          windEmulatorStep4_WECSim_B.BusAssignment_a.motorTorque_Nm;

        /* Gain: '<S28>/motorVoltage_V' */
        windEmulatorStep4_WECSim_B.motorVoltage_V =
          windEmulatorStep4_WECSim_cal->motorVoltage_V_Gain *
          windEmulatorStep4_WECSim_B.BusAssignment_a.motorVoltage_V;

        /* Gain: '<S28>/shaftPower_W' */
        windEmulatorStep4_WECSim_B.shaftPower_W_a =
          windEmulatorStep4_WECSim_cal->shaftPower_W_Gain_e *
          windEmulatorStep4_WECSim_B.BusAssignment_a.shaftPower_W;

        /* Bias: '<S28>/statusWord' */
        windEmulatorStep4_WECSim_B.statusWord = static_cast<uint16_T>
          (windEmulatorStep4_WECSim_B.BusAssignment_a.statusWord +
           windEmulatorStep4_WECSim_cal->statusWord_Bias);

        /* Memory: '<S2>/Memory' */
        windEmulatorStep4_WECSim_B.Memory_b =
          windEmulatorStep4_WECSim_DW.Memory_PreviousInput_i;

        /* RelationalOperator: '<S2>/NotEqual' incorporates:
         *  Constant: '<S2>/powerUpButton'
         */
        windEmulatorStep4_WECSim_B.NotEqual =
          (windEmulatorStep4_WECSim_cal->powerUpButton_Value !=
           windEmulatorStep4_WECSim_B.Memory_b);

        /* Memory: '<S2>/Memory1' */
        windEmulatorStep4_WECSim_B.Memory1 =
          windEmulatorStep4_WECSim_DW.Memory1_PreviousInput;

        /* RelationalOperator: '<S2>/NotEqual1' incorporates:
         *  Constant: '<S2>/powerDownButton'
         */
        windEmulatorStep4_WECSim_B.NotEqual1 =
          (windEmulatorStep4_WECSim_cal->powerDownButton_Value !=
           windEmulatorStep4_WECSim_B.Memory1);

        /* Memory: '<S2>/Memory2' */
        windEmulatorStep4_WECSim_B.Memory2 =
          windEmulatorStep4_WECSim_DW.Memory2_PreviousInput;

        /* RelationalOperator: '<S2>/NotEqual2' incorporates:
         *  Constant: '<S2>/resetFaultButton'
         */
        windEmulatorStep4_WECSim_B.NotEqual2 =
          (windEmulatorStep4_WECSim_cal->resetFaultButton_Value !=
           windEmulatorStep4_WECSim_B.Memory2);

        /* DataTypeConversion: '<S2>/Cast To Double' */
        windEmulatorStep4_WECSim_B.CastToDouble_o =
          windEmulatorStep4_WECSim_B.NotEqual2;

        /* Memory: '<S4>/Memory' */
        windEmulatorStep4_WECSim_B.Memory_i =
          windEmulatorStep4_WECSim_DW.Memory_PreviousInput_k;

        /* RelationalOperator: '<S4>/NotEqual' incorporates:
         *  Constant: '<S4>/eStopButton'
         */
        windEmulatorStep4_WECSim_B.NotEqual_i =
          (windEmulatorStep4_WECSim_cal->eStopButton_Value !=
           windEmulatorStep4_WECSim_B.Memory_i);

        /* Constant: '<S2>/ACS880CtrlMode' */
        windEmulatorStep4_WECSim_B.ACS880CtrlMode = tmp_11;

        /* Chart: '<S18>/ABB Fieldbus Control' */
        if (windEmulatorStep4_WECSim_DW.temporalCounter_i1_l < 31) {
          windEmulatorStep4_WECSim_DW.temporalCounter_i1_l = static_cast<uint8_T>
            (windEmulatorStep4_WECSim_DW.temporalCounter_i1_l + 1);
        }

        windEmulatorStep4_WECSim_DW.sfEvent_l = windEmulatorStep4__CALL_EVENT_k;
        if (windEmulatorStep4_WECSim_DW.is_active_c7_windEmulatorStep4_ == 0) {
          windEmulatorStep4_WECSim_DW.is_active_c7_windEmulatorStep4_ = 1U;
          windEmulatorStep4_WECSim_DW.is_active_UpdateStateMachine = 1U;
          windEmulatorStep4_WECSim_DW.is_UpdateStateMachine =
            windEmulatorStep4_IN_initialize;
          windEmulatorStep_cwInitialize_a();
          windEmulatorStep4_WECSim_B.state_e = abbStateEnum_init;
          windEmulatorStep4_WECSim_DW.is_active_UpdateControlWord = 1U;
        } else {
          windEmulato_swParseStatusWord_m();
          windEmulatorStep4_WECSim_DW.cwRESET =
            windEmulatorStep4_WECSim_B.CastToDouble_o;
          switch (windEmulatorStep4_WECSim_DW.is_UpdateStateMachine) {
           case windEmulatorStep4__IN_DelayOFF1:
            if (windEmulatorStep4_WECSim_DW.temporalCounter_i1_l >= 25) {
              windEmulatorStep4_WECSim_DW.is_UpdateStateMachine =
                windEmula_IN_notReadyToSwitchOn;
              windEmulatorStep4_WECSim_DW.cwOFF1_CONTROL = 0.0;
            } else {
              windEmulatorStep4_WECSim_B.state_e = abbStateEnum_delayOff1;
            }
            break;

           case windEmulatorStep4_IN_initialize:
            windEmulatorStep4_WECSim_DW.is_UpdateStateMachine =
              windEmula_IN_notReadyToSwitchOn;
            windEmulatorStep4_WECSim_DW.cwOFF1_CONTROL = 0.0;
            break;

           case windEmula_IN_notReadyToSwitchOn:
            f = ((windEmulatorStep4_WECSim_DW.swRDY_ON == 1.0) &&
                 (windEmulatorStep4_WECSim_DW.swWARNING == 0.0));
            if (f) {
              windEmulatorStep4_WECSim_DW.is_UpdateStateMachine =
                windEmulator_IN_readyToSwitchOn;
              windEmulatorStep4_WECSim_DW.cwOFF1_CONTROL = 0.0;
            } else {
              windEmulatorStep4_WECSim_B.state_e =
                abbStateEnum_notReadyToSwitchOn;
            }
            break;

           case windEmulat_IN_operationDisabled:
            f = (windEmulatorStep4_WECSim_B.NotEqual &&
                 (windEmulatorStep4_WECSim_DW.swRDY_RUN == 1.0));
            if (f) {
              windEmulatorStep4_WECSim_DW.is_UpdateStateMachine =
                windEmulato_IN_operationEnabled;
              windEmulatorStep4_WECSim_DW.cwOFF1_CONTROL = 1.0;
              windEmulatorStep4_WECSim_DW.cwENABLE_OPERATION = 1.0;
            } else {
              f = (windEmulatorStep4_WECSim_B.NotEqual1 ||
                   windEmulatorStep4_WECSim_B.NotEqual_i ||
                   (windEmulatorStep4_WECSim_DW.swTRIPPED == 1.0) ||
                   (windEmulatorStep4_WECSim_DW.swSWC_ON_INHIB == 1.0) ||
                   (windEmulatorStep4_WECSim_DW.swWARNING == 1.0));
              if (f) {
                windEmulatorStep4_WECSim_DW.temporalCounter_i1_l = 0U;
                windEmulatorStep4_WECSim_DW.is_UpdateStateMachine =
                  windEmulatorStep4__IN_DelayOFF1;
              } else {
                windEmulatorStep4_WECSim_B.state_e =
                  abbStateEnum_operationDisabled;
              }
            }
            break;

           case windEmulato_IN_operationEnabled:
            if (windEmulatorStep4_WECSim_B.NotEqual1) {
              windEmulatorStep4_WECSim_DW.cwENABLE_OPERATION = 0.0;
              windEmulatorStep4_WECSim_DW.is_UpdateStateMachine =
                windEmulat_IN_operationDisabled;
              windEmulatorStep4_WECSim_DW.cwOFF1_CONTROL = 1.0;
            } else {
              f = (windEmulatorStep4_WECSim_B.NotEqual_i ||
                   (windEmulatorStep4_WECSim_DW.swTRIPPED == 1.0) ||
                   (windEmulatorStep4_WECSim_DW.swSWC_ON_INHIB == 1.0) ||
                   (windEmulatorStep4_WECSim_DW.swWARNING == 1.0));
              if (f) {
                windEmulatorStep4_WECSim_DW.cwENABLE_OPERATION = 0.0;
                windEmulatorStep4_WECSim_DW.temporalCounter_i1_l = 0U;
                windEmulatorStep4_WECSim_DW.is_UpdateStateMachine =
                  windEmulatorStep4__IN_DelayOFF1;
              } else {
                windEmulatorStep4_WECSim_B.state_e =
                  abbStateEnum_operationEnabled;
              }
            }
            break;

           default:
            /* case IN_readyToSwitchOn: */
            f = (windEmulatorStep4_WECSim_B.NotEqual &&
                 (windEmulatorStep4_WECSim_DW.swREMOTE == 1.0));
            if (f) {
              windEmulatorStep4_WECSim_DW.is_UpdateStateMachine =
                windEmulat_IN_operationDisabled;
              windEmulatorStep4_WECSim_DW.cwOFF1_CONTROL = 1.0;
            } else {
              windEmulatorStep4_WECSim_B.state_e = abbStateEnum_readyToSwitchOn;
            }
            break;
          }

          windEmulat_cwBuildControlWord_o();
        }

        /* End of Chart: '<S18>/ABB Fieldbus Control' */

        /* DataTypeConversion: '<S2>/Cast To Double1' */
        windEmulatorStep4_WECSim_B.CastToDouble1_a =
          windEmulatorStep4_WECSim_B.state_e;

        /* Constant: '<S4>/expType' */
        windEmulatorStep4_WECSim_B.expType_a =
          windEmulatorStep4_WECSim_cal->expType_Value;

        /* DataTypeConversion: '<S4>/ToUint16' incorporates:
         *  Constant: '<S4>/expType'
         */
        windEmulatorStep4_WECSim_B.ToUint16 =
          windEmulatorStep4_WECSim_B.expType_a;

        /* RelationalOperator: '<S4>/Equal' incorporates:
         *  Constant: '<S4>/expModeHil'
         *  Constant: '<S4>/expType'
         */
        windEmulatorStep4_WECSim_B.Equal = (windEmulatorStep4_WECSim_B.expType_a
          == windEmulatorStep4_WECSim_cal->expModeHil_Value);

        /* RelationalOperator: '<S4>/Equal1' incorporates:
         *  Constant: '<S4>/expModeSid'
         *  Constant: '<S4>/expType'
         */
        windEmulatorStep4_WECSim_B.Equal1 =
          (windEmulatorStep4_WECSim_B.expType_a ==
           windEmulatorStep4_WECSim_cal->expModeSid_Value);

        /* Memory: '<S4>/Memory1' */
        windEmulatorStep4_WECSim_B.Memory1_n =
          windEmulatorStep4_WECSim_DW.Memory1_PreviousInput_d;

        /* RelationalOperator: '<S4>/NotEqual1' incorporates:
         *  Constant: '<S4>/startButton'
         */
        windEmulatorStep4_WECSim_B.NotEqual1_a =
          (windEmulatorStep4_WECSim_cal->startButton_Value !=
           windEmulatorStep4_WECSim_B.Memory1_n);

        /* Memory: '<S4>/Memory2' */
        windEmulatorStep4_WECSim_B.Memory2_c =
          windEmulatorStep4_WECSim_DW.Memory2_PreviousInput_l;

        /* RelationalOperator: '<S4>/NotEqual2' incorporates:
         *  Constant: '<S4>/stopButton'
         */
        windEmulatorStep4_WECSim_B.NotEqual2_c =
          (windEmulatorStep4_WECSim_cal->stopButton_Value !=
           windEmulatorStep4_WECSim_B.Memory2_c);

        /* Constant: '<S4>/expRunTime' */
        windEmulatorStep4_WECSim_B.expRunTime =
          windEmulatorStep4_WECSim_cal->expRunTime_Value;

        /* Chart: '<S4>/FexcRamp' incorporates:
         *  Constant: '<S4>/expType'
         */
        if (windEmulatorStep4_WECSim_DW.temporalCounter_i1 < MAX_uint32_T) {
          windEmulatorStep4_WECSim_DW.temporalCounter_i1++;
        }

        windEmulatorStep4_WECSim_DW.sfEvent = windEmulatorStep4__CALL_EVENT_k;
        if (windEmulatorStep4_WECSim_DW.is_active_c3_windEmulatorStep4_ == 0) {
          windEmulatorStep4_WECSim_DW.is_active_c3_windEmulatorStep4_ = 1U;
          windEmulatorStep4_WECSim_DW.is_c3_windEmulatorStep4_WECSim =
            windEmulatorStep4_WECSi_IN_init;
          windEmulatorStep4_WECSim_B.runCounter_b = 0U;
          windEmulatorStep4_WECSim_B.stepCounter_n = 1U;
          windEmulatorStep4_WECSim_B.ramp_l = 0.0;
          windEmulatorStep4_WECSim_DW.rampLast = 0.0;
          windEmulatorStep4_WECSim_B.time_c = 0.0;
          windEmulatorStep4_WECSim_DW.rampUpTime = 0.0;
          windEmulatorStep4_WECSim_DW.runTime = 0.0;
          windEmulatorStep4_WECSim_B.resetHilIntegrator_h = true;
          windEmulatorStep4_WECSim_B.resetSidIntegrator_h = true;
        } else {
          switch (windEmulatorStep4_WECSim_DW.is_c3_windEmulatorStep4_WECSim) {
           case windEmulatorStep4_WECSi_IN_idle:
            if (windEmulatorStep4_WECSim_B.NotEqual1_a) {
              windEmulatorStep4_WECSim_DW.temporalCounter_i1 = 0U;
              windEmulatorStep4_WECSim_DW.is_c3_windEmulatorStep4_WECSim =
                windEmulatorStep4_WEC_IN_rampup;
              q0 = windEmulatorStep4_WECSim_B.runCounter_b;
              qY = q0 + 1U;
              if (qY < q0) {
                qY = MAX_uint32_T;
              }

              windEmulatorStep4_WECSim_B.runCounter_b = qY;
              windEmulatorStep4_WECSim_B.resetHilIntegrator_h =
                (windEmulatorStep4_WECSim_B.expType_a != expTypeEnum_hil);
              windEmulatorStep4_WECSim_B.resetSidIntegrator_h =
                (windEmulatorStep4_WECSim_B.expType_a != expTypeEnum_sid);
            } else {
              windEmulatorStep4_WECSim_B.ramp_l = 0.0;
              windEmulatorStep4_WECSim_DW.rampLast = 0.0;
              windEmulatorStep4_WECSim_B.stepCounter_n = 1U;
              windEmulatorStep4_WECSim_B.time_c = 0.0;
              windEmulatorStep4_WECSim_DW.rampUpTime = 0.0;
              windEmulatorStep4_WECSim_DW.runTime = 0.0;
              windEmulatorStep4_WECSim_B.resetHilIntegrator_h = true;
              windEmulatorStep4_WECSim_B.resetSidIntegrator_h = true;
            }
            break;

           case windEmulatorStep4_WECSi_IN_init:
            windEmulatorStep4_WECSim_DW.is_c3_windEmulatorStep4_WECSim =
              windEmulatorStep4_WECSi_IN_idle;
            break;

           case windEmulatorStep4_W_IN_rampdown:
            if (windEmulatorStep4_WECSim_B.ramp_l <= 0.0) {
              windEmulatorStep4_WECSim_DW.is_c3_windEmulatorStep4_WECSim =
                windEmulatorStep4_WECSi_IN_idle;
            } else {
              u1 = static_cast<real_T>
                (windEmulatorStep4_WECSim_DW.temporalCounter_i1) * 0.004 / u1;
              if ((u1 <= 0.0) || rtIsNaN(u1)) {
                u1 = 0.0;
              }

              windEmulatorStep4_WECSim_B.ramp_l =
                windEmulatorStep4_WECSim_DW.rampLast - u1;
              windEmulatorStep4_WECSim_B.time_c =
                (windEmulatorStep4_WECSim_DW.rampUpTime +
                 windEmulatorStep4_WECSim_DW.runTime) + static_cast<real_T>
                (windEmulatorStep4_WECSim_DW.temporalCounter_i1) * 0.004;
              q0 = windEmulatorStep4_WECSim_B.stepCounter_n;
              qY = q0 + 1U;
              if (qY < q0) {
                qY = MAX_uint32_T;
              }

              windEmulatorStep4_WECSim_B.stepCounter_n = qY;
              windEmulatorStep4_WECSim_B.resetHilIntegrator_h =
                (windEmulatorStep4_WECSim_B.expType_a != expTypeEnum_hil);
              windEmulatorStep4_WECSim_B.resetSidIntegrator_h =
                (windEmulatorStep4_WECSim_B.expType_a != expTypeEnum_sid);
            }
            break;

           case windEmulatorStep4_WEC_IN_rampup:
            if (windEmulatorStep4_WECSim_B.NotEqual2_c) {
              windEmulatorStep4_WECSim_DW.temporalCounter_i1 = 0U;
              windEmulatorStep4_WECSim_DW.is_c3_windEmulatorStep4_WECSim =
                windEmulatorStep4_W_IN_rampdown;
            } else if (windEmulatorStep4_WECSim_B.ramp_l >= 1.0) {
              windEmulatorStep4_WECSim_DW.temporalCounter_i1 = 0U;
              windEmulatorStep4_WECSim_DW.is_c3_windEmulatorStep4_WECSim =
                windEmulatorStep4_WEC_IN_runing;
            } else {
              u1 = static_cast<real_T>
                (windEmulatorStep4_WECSim_DW.temporalCounter_i1) * 0.004 / u1;
              if ((u1 >= 1.0) || rtIsNaN(u1)) {
                windEmulatorStep4_WECSim_B.ramp_l = 1.0;
              } else {
                windEmulatorStep4_WECSim_B.ramp_l = u1;
              }

              q0 = windEmulatorStep4_WECSim_B.stepCounter_n;
              qY = q0 + 1U;
              if (qY < q0) {
                qY = MAX_uint32_T;
              }

              windEmulatorStep4_WECSim_B.stepCounter_n = qY;
              windEmulatorStep4_WECSim_DW.rampLast =
                windEmulatorStep4_WECSim_B.ramp_l;
              windEmulatorStep4_WECSim_DW.rampUpTime = static_cast<real_T>
                (windEmulatorStep4_WECSim_DW.temporalCounter_i1) * 0.004;
              windEmulatorStep4_WECSim_B.time_c =
                windEmulatorStep4_WECSim_DW.rampUpTime;
            }
            break;

           default:
            /* case IN_runing: */
            f = (windEmulatorStep4_WECSim_B.NotEqual2_c ||
                 (windEmulatorStep4_WECSim_B.time_c >=
                  windEmulatorStep4_WECSim_B.expRunTime));
            if (f) {
              windEmulatorStep4_WECSim_DW.temporalCounter_i1 = 0U;
              windEmulatorStep4_WECSim_DW.is_c3_windEmulatorStep4_WECSim =
                windEmulatorStep4_W_IN_rampdown;
            } else {
              windEmulatorStep4_WECSim_B.ramp_l = 1.0;
              windEmulatorStep4_WECSim_DW.rampLast =
                windEmulatorStep4_WECSim_B.ramp_l;
              windEmulatorStep4_WECSim_DW.runTime = static_cast<real_T>
                (windEmulatorStep4_WECSim_DW.temporalCounter_i1) * 0.004;
              windEmulatorStep4_WECSim_B.time_c =
                windEmulatorStep4_WECSim_DW.rampUpTime +
                windEmulatorStep4_WECSim_DW.runTime;
              q0 = windEmulatorStep4_WECSim_B.stepCounter_n;
              qY = q0 + 1U;
              if (qY < q0) {
                qY = MAX_uint32_T;
              }

              windEmulatorStep4_WECSim_B.stepCounter_n = qY;
              windEmulatorStep4_WECSim_B.resetHilIntegrator_h =
                (windEmulatorStep4_WECSim_B.expType_a != expTypeEnum_hil);
              windEmulatorStep4_WECSim_B.resetSidIntegrator_h =
                (windEmulatorStep4_WECSim_B.expType_a != expTypeEnum_sid);
            }
            break;
          }
        }

        /* BusAssignment: '<S4>/Bus Assignment' */
        windEmulatorStep4_WECSim_B.BusAssignment_b.expType =
          windEmulatorStep4_WECSim_B.ToUint16;
        windEmulatorStep4_WECSim_B.BusAssignment_b.runHil =
          windEmulatorStep4_WECSim_B.Equal;
        windEmulatorStep4_WECSim_B.BusAssignment_b.runSid =
          windEmulatorStep4_WECSim_B.Equal1;
        windEmulatorStep4_WECSim_B.BusAssignment_b.time =
          windEmulatorStep4_WECSim_B.time_c;
        windEmulatorStep4_WECSim_B.BusAssignment_b.ramp =
          windEmulatorStep4_WECSim_B.ramp_l;
        windEmulatorStep4_WECSim_B.BusAssignment_b.runCounter =
          windEmulatorStep4_WECSim_B.runCounter_b;
        windEmulatorStep4_WECSim_B.BusAssignment_b.stepCounter =
          windEmulatorStep4_WECSim_B.stepCounter_n;
        windEmulatorStep4_WECSim_B.BusAssignment_b.resetHilIntegrator =
          windEmulatorStep4_WECSim_B.resetHilIntegrator_h;
        windEmulatorStep4_WECSim_B.BusAssignment_b.resetSidIntegrator =
          windEmulatorStep4_WECSim_B.resetSidIntegrator_h;

        /* DataTypeConversion: '<S3>/toExpTypeEnum' */
        windEmulatorStep4_WECSim_B.toExpTypeEnum =
          windEmulatorStep4_WECSim_B.BusAssignment_b.expType;

        /* Switch: '<S555>/Switch' */
        if (windEmulatorStep4_WECSim_B.BusAssignment_b.runSid) {
          /* Switch: '<S555>/Switch' */
          windEmulatorStep4_WECSim_B.rampValue =
            windEmulatorStep4_WECSim_B.BusAssignment_b.ramp;
        } else {
          /* Switch: '<S555>/Switch' incorporates:
           *  Constant: '<S555>/Constant1'
           */
          windEmulatorStep4_WECSim_B.rampValue =
            windEmulatorStep4_WECSim_cal->Constant1_Value_e;
        }

        /* End of Switch: '<S555>/Switch' */

        /* MultiPortSwitch: '<S609>/Multiport Switch' incorporates:
         *  Constant: '<S555>/sidType'
         */
        switch (windEmulatorStep4_WECSim_cal->sidType_Value) {
         case sidTypeEnum_off:
          /* MultiPortSwitch: '<S609>/Multiport Switch' incorporates:
           *  Constant: '<S609>/Constant'
           */
          windEmulatorStep4_WECSim_B.MultiportSwitch_h =
            windEmulatorStep4_WECSim_cal->Constant_Value_im;
          break;

         case sidTypeEnum_manual:
          /* MultiPortSwitch: '<S609>/Multiport Switch' incorporates:
           *  Constant: '<S609>/Constant'
           */
          windEmulatorStep4_WECSim_B.MultiportSwitch_h =
            windEmulatorStep4_WECSim_cal->Constant_Value_im;
          break;

         case sidTypeEnum_fromFile:
          /* MultiPortSwitch: '<S609>/Multiport Switch' */
          windEmulatorStep4_WECSim_B.MultiportSwitch_h =
            windEmulatorStep4_WECSim_B.BusAssignment_b.stepCounter;
          break;

         default:
          /* MultiPortSwitch: '<S609>/Multiport Switch' incorporates:
           *  Constant: '<S609>/Constant'
           */
          windEmulatorStep4_WECSim_B.MultiportSwitch_h =
            windEmulatorStep4_WECSim_cal->Constant_Value_im;
          break;
        }

        /* End of MultiPortSwitch: '<S609>/Multiport Switch' */

        /* Math: '<S609>/Mod' incorporates:
         *  Constant: '<S609>/Length of input'
         */
        q0 = windEmulatorStep4_WECSim_B.MultiportSwitch_h;
        qY = windEmulatorStep4_WECSim_cal->Lengthofinput_Value;
        if (qY == 0U) {
          /* Math: '<S609>/Mod' */
          windEmulatorStep4_WECSim_B.Mod = q0;
        } else {
          /* Math: '<S609>/Mod' */
          windEmulatorStep4_WECSim_B.Mod = q0 % qY;
        }

        /* End of Math: '<S609>/Mod' */

        /* RelationalOperator: '<S609>/Out of bounds' incorporates:
         *  Constant: '<S609>/Length of input'
         */
        windEmulatorStep4_WECSim_B.Outofbounds = (windEmulatorStep4_WECSim_B.Mod
          > windEmulatorStep4_WECSim_cal->Lengthofinput_Value);
      }

      /* Switch: '<S609>/Switch2' */
      if (windEmulatorStep4_WECSim_B.Outofbounds) {
        /* Switch: '<S609>/Switch2' incorporates:
         *  Constant: '<S609>/Set bound'
         */
        windEmulatorStep4_WECSim_B.Switch2 =
          windEmulatorStep4_WECSim_cal->Setbound_Value;
      } else {
        /* Switch: '<S609>/Switch2' incorporates:
         *  Inport: '<Root>/inportTorque_Nm'
         */
        windEmulatorStep4_WECSim_B.Switch2 =
          windEmulatorStep4_WECSim_U.inportTorque_Nm;
      }

      /* End of Switch: '<S609>/Switch2' */

      /* RateLimiter: '<S609>/torqueSlewRate' */
      if (windEmulatorStep4_WECSim_DW.LastMajorTime == (rtInf)) {
        /* RateLimiter: '<S609>/torqueSlewRate' */
        windEmulatorStep4_WECSim_B.torqueSlewRate =
          windEmulatorStep4_WECSim_B.Switch2;
      } else {
        deltaT_tmp = windEmulatorStep4_WECSim_M->Timing.t[0];
        u1 = deltaT_tmp - windEmulatorStep4_WECSim_DW.LastMajorTime;
        if (windEmulatorStep4_WECSim_DW.LastMajorTime == deltaT_tmp) {
          if (windEmulatorStep4_WECSim_DW.PrevLimited) {
            /* RateLimiter: '<S609>/torqueSlewRate' */
            windEmulatorStep4_WECSim_B.torqueSlewRate =
              windEmulatorStep4_WECSim_DW.PrevY;
          } else {
            /* RateLimiter: '<S609>/torqueSlewRate' */
            windEmulatorStep4_WECSim_B.torqueSlewRate =
              windEmulatorStep4_WECSim_B.Switch2;
          }
        } else {
          riseValLimit = u1 * tmp_r;
          rateLimiterRate = windEmulatorStep4_WECSim_B.Switch2 -
            windEmulatorStep4_WECSim_DW.PrevY;
          if (rateLimiterRate > riseValLimit) {
            /* RateLimiter: '<S609>/torqueSlewRate' */
            windEmulatorStep4_WECSim_B.torqueSlewRate =
              windEmulatorStep4_WECSim_DW.PrevY + riseValLimit;
            f = true;
          } else {
            riseValLimit = -tmp_r;
            u1 *= riseValLimit;
            if (rateLimiterRate < u1) {
              /* RateLimiter: '<S609>/torqueSlewRate' */
              windEmulatorStep4_WECSim_B.torqueSlewRate =
                windEmulatorStep4_WECSim_DW.PrevY + u1;
              f = true;
            } else {
              /* RateLimiter: '<S609>/torqueSlewRate' */
              windEmulatorStep4_WECSim_B.torqueSlewRate =
                windEmulatorStep4_WECSim_B.Switch2;
              f = false;
            }
          }

          if (rtsiIsModeUpdateTimeStep(&windEmulatorStep4_WECSim_M->solverInfo))
          {
            windEmulatorStep4_WECSim_DW.PrevLimited = f;
          }
        }
      }

      /* MultiPortSwitch: '<S555>/Multiport Switch' incorporates:
       *  Constant: '<S555>/sidType'
       */
      switch (windEmulatorStep4_WECSim_cal->sidType_Value) {
       case sidTypeEnum_off:
        /* MultiPortSwitch: '<S555>/Multiport Switch' incorporates:
         *  Constant: '<S555>/Constant3'
         */
        windEmulatorStep4_WECSim_B.MultiportSwitch =
          windEmulatorStep4_WECSim_cal->Constant3_Value;
        break;

       case sidTypeEnum_manual:
        /* MultiPortSwitch: '<S555>/Multiport Switch' incorporates:
         *  Constant: '<S555>/manualTorqueSetpoint_Nm'
         */
        windEmulatorStep4_WECSim_B.MultiportSwitch =
          windEmulatorStep4_WECSim_cal->manualTorqueSetpoint_Nm_Value;
        break;

       case sidTypeEnum_fromFile:
        /* Gain: '<S609>/fromFileTorqueNow_Nm' */
        windEmulatorStep4_WECSim_B.fromFileTorqueNow_Nm =
          windEmulatorStep4_WECSim_cal->fromFileTorqueNow_Nm_Gain *
          windEmulatorStep4_WECSim_B.torqueSlewRate;

        /* MultiPortSwitch: '<S555>/Multiport Switch' */
        windEmulatorStep4_WECSim_B.MultiportSwitch =
          windEmulatorStep4_WECSim_B.fromFileTorqueNow_Nm;
        break;

       default:
        /* MultiPortSwitch: '<S555>/Multiport Switch' incorporates:
         *  Constant: '<S555>/Constant3'
         */
        windEmulatorStep4_WECSim_B.MultiportSwitch =
          windEmulatorStep4_WECSim_cal->Constant3_Value;
        break;
      }

      /* End of MultiPortSwitch: '<S555>/Multiport Switch' */

      /* Product: '<S555>/Product' */
      windEmulatorStep4_WECSim_B.Product_g =
        windEmulatorStep4_WECSim_B.rampValue *
        windEmulatorStep4_WECSim_B.MultiportSwitch;

      /* Switch: '<S609>/Switch1' */
      if (windEmulatorStep4_WECSim_B.Outofbounds) {
        /* Switch: '<S609>/Switch1' incorporates:
         *  Constant: '<S609>/Set bound'
         */
        windEmulatorStep4_WECSim_B.Switch1 =
          windEmulatorStep4_WECSim_cal->Setbound_Value;
      } else {
        /* Switch: '<S609>/Switch1' incorporates:
         *  Inport: '<Root>/inportSpeed_rpm'
         */
        windEmulatorStep4_WECSim_B.Switch1 =
          windEmulatorStep4_WECSim_U.inportSpeed_rpm;
      }

      /* End of Switch: '<S609>/Switch1' */

      /* RateLimiter: '<S609>/speedSlewRate' */
      if (windEmulatorStep4_WECSim_DW.LastMajorTime_j == (rtInf)) {
        /* RateLimiter: '<S609>/speedSlewRate' */
        windEmulatorStep4_WECSim_B.speedSlewRate =
          windEmulatorStep4_WECSim_B.Switch1;
      } else {
        deltaT_tmp = windEmulatorStep4_WECSim_M->Timing.t[0];
        u1 = deltaT_tmp - windEmulatorStep4_WECSim_DW.LastMajorTime_j;
        if (windEmulatorStep4_WECSim_DW.LastMajorTime_j == deltaT_tmp) {
          if (windEmulatorStep4_WECSim_DW.PrevLimited_d) {
            /* RateLimiter: '<S609>/speedSlewRate' */
            windEmulatorStep4_WECSim_B.speedSlewRate =
              windEmulatorStep4_WECSim_DW.PrevY_a;
          } else {
            /* RateLimiter: '<S609>/speedSlewRate' */
            windEmulatorStep4_WECSim_B.speedSlewRate =
              windEmulatorStep4_WECSim_B.Switch1;
          }
        } else {
          riseValLimit = u1 * tmp_q;
          rateLimiterRate = windEmulatorStep4_WECSim_B.Switch1 -
            windEmulatorStep4_WECSim_DW.PrevY_a;
          if (rateLimiterRate > riseValLimit) {
            /* RateLimiter: '<S609>/speedSlewRate' */
            windEmulatorStep4_WECSim_B.speedSlewRate =
              windEmulatorStep4_WECSim_DW.PrevY_a + riseValLimit;
            f = true;
          } else {
            riseValLimit = -tmp_q;
            u1 *= riseValLimit;
            if (rateLimiterRate < u1) {
              /* RateLimiter: '<S609>/speedSlewRate' */
              windEmulatorStep4_WECSim_B.speedSlewRate =
                windEmulatorStep4_WECSim_DW.PrevY_a + u1;
              f = true;
            } else {
              /* RateLimiter: '<S609>/speedSlewRate' */
              windEmulatorStep4_WECSim_B.speedSlewRate =
                windEmulatorStep4_WECSim_B.Switch1;
              f = false;
            }
          }

          if (rtsiIsModeUpdateTimeStep(&windEmulatorStep4_WECSim_M->solverInfo))
          {
            windEmulatorStep4_WECSim_DW.PrevLimited_d = f;
          }
        }
      }

      /* MultiPortSwitch: '<S555>/Multiport Switch1' incorporates:
       *  Constant: '<S555>/sidType'
       */
      switch (windEmulatorStep4_WECSim_cal->sidType_Value) {
       case sidTypeEnum_off:
        /* MultiPortSwitch: '<S555>/Multiport Switch1' incorporates:
         *  Constant: '<S555>/Constant3'
         */
        windEmulatorStep4_WECSim_B.MultiportSwitch1 =
          windEmulatorStep4_WECSim_cal->Constant3_Value;
        break;

       case sidTypeEnum_manual:
        /* MultiPortSwitch: '<S555>/Multiport Switch1' incorporates:
         *  Constant: '<S555>/manualSpeedSetpoint_rpm'
         */
        windEmulatorStep4_WECSim_B.MultiportSwitch1 =
          windEmulatorStep4_WECSim_cal->manualSpeedSetpoint_rpm_Value;
        break;

       case sidTypeEnum_fromFile:
        /* Gain: '<S609>/fromFileSpeedNow_rpm' */
        windEmulatorStep4_WECSim_B.fromFileSpeedNow_rpm =
          windEmulatorStep4_WECSim_cal->fromFileSpeedNow_rpm_Gain *
          windEmulatorStep4_WECSim_B.speedSlewRate;

        /* MultiPortSwitch: '<S555>/Multiport Switch1' */
        windEmulatorStep4_WECSim_B.MultiportSwitch1 =
          windEmulatorStep4_WECSim_B.fromFileSpeedNow_rpm;
        break;

       default:
        /* MultiPortSwitch: '<S555>/Multiport Switch1' incorporates:
         *  Constant: '<S555>/Constant3'
         */
        windEmulatorStep4_WECSim_B.MultiportSwitch1 =
          windEmulatorStep4_WECSim_cal->Constant3_Value;
        break;
      }

      /* End of MultiPortSwitch: '<S555>/Multiport Switch1' */

      /* Product: '<S555>/Product1' */
      windEmulatorStep4_WECSim_B.Product1 = windEmulatorStep4_WECSim_B.rampValue
        * windEmulatorStep4_WECSim_B.MultiportSwitch1;

      /* Sum: '<S13>/Sum' */
      windEmulatorStep4_WECSim_B.Sum = windEmulatorStep4_WECSim_B.Product1 -
        windEmulatorStep4_WECSim_B.BusAssignment_a.motorSpeed_rpm;

      /* Product: '<S595>/PProd Out' incorporates:
       *  Constant: '<S13>/acs880SpeedPGain'
       */
      windEmulatorStep4_WECSim_B.PProdOut = windEmulatorStep4_WECSim_B.Sum *
        windEmulatorStep4_WECSim_cal->acs880SpeedPGain_Value;

      /* Integrator: '<S590>/Integrator' incorporates:
       *  RateLimiter: '<S2>/acs880RateLim'
       *  RateLimiter: '<S373>/Rate Limiter'
       *  RateLimiter: '<S436>/Rate Limiter'
       *  SimscapeExecutionBlock: '<S216>/OUTPUT_1_0'
       *  SimscapeExecutionBlock: '<S216>/OUTPUT_1_1'
       *  SimscapeExecutionBlock: '<S216>/STATE_1'
       *  SimscapeExecutionBlock: '<S332>/OUTPUT_1_0'
       *  SimscapeExecutionBlock: '<S332>/STATE_1'
       */
      tmp_h = rtsiIsModeUpdateTimeStep(&windEmulatorStep4_WECSim_M->solverInfo);
      if (tmp_h) {
        f = (((windEmulatorStep4_WECSim_PrevZCX.Integrator_Reset_ZCE ==
               POS_ZCSIG) !=
              windEmulatorStep4_WECSim_B.BusAssignment_b.resetSidIntegrator) &&
             (windEmulatorStep4_WECSim_PrevZCX.Integrator_Reset_ZCE !=
              UNINITIALIZED_ZCSIG));
        windEmulatorStep4_WECSim_PrevZCX.Integrator_Reset_ZCE =
          windEmulatorStep4_WECSim_B.BusAssignment_b.resetSidIntegrator;

        /* evaluate zero-crossings and the level of the reset signal */
        if (f || windEmulatorStep4_WECSim_B.BusAssignment_b.resetSidIntegrator)
        {
          windEmulatorStep4_WECSim_X.Integrator_CSTATE =
            windEmulatorStep4_WECSim_cal->PIDController_InitialConditio_o;
        }
      }

      /* Integrator: '<S590>/Integrator' */
      windEmulatorStep4_WECSim_B.Integrator =
        windEmulatorStep4_WECSim_X.Integrator_CSTATE;

      /* Sum: '<S600>/Sum' */
      windEmulatorStep4_WECSim_B.Sum_k = windEmulatorStep4_WECSim_B.PProdOut +
        windEmulatorStep4_WECSim_B.Integrator;

      /* RelationalOperator: '<S598>/LowerRelop1' incorporates:
       *  Constant: '<S13>/Constant1'
       */
      windEmulatorStep4_WECSim_B.LowerRelop1 = (windEmulatorStep4_WECSim_B.Sum_k
        > tmp_w);

      /* RelationalOperator: '<S598>/UpperRelop' incorporates:
       *  Constant: '<S13>/Constant2'
       */
      windEmulatorStep4_WECSim_B.UpperRelop = (windEmulatorStep4_WECSim_B.Sum_k <
        Clock_tmp);

      /* Switch: '<S598>/Switch' */
      if (windEmulatorStep4_WECSim_B.UpperRelop) {
        /* Switch: '<S598>/Switch' incorporates:
         *  Constant: '<S13>/Constant2'
         */
        windEmulatorStep4_WECSim_B.Switch = Clock_tmp;
      } else {
        /* Switch: '<S598>/Switch' */
        windEmulatorStep4_WECSim_B.Switch = windEmulatorStep4_WECSim_B.Sum_k;
      }

      /* Switch: '<S598>/Switch2' */
      if (windEmulatorStep4_WECSim_B.LowerRelop1) {
        /* Switch: '<S598>/Switch2' incorporates:
         *  Constant: '<S13>/Constant1'
         */
        windEmulatorStep4_WECSim_B.Switch2_h = tmp_w;
      } else {
        /* Switch: '<S598>/Switch2' */
        windEmulatorStep4_WECSim_B.Switch2_h = windEmulatorStep4_WECSim_B.Switch;
      }

      /* BusAssignment: '<S13>/Bus Assignment' */
      windEmulatorStep4_WECSim_B.BusAssignment_f.acs800Torque_Nm =
        windEmulatorStep4_WECSim_B.Product_g;
      windEmulatorStep4_WECSim_B.BusAssignment_f.acs880Torque_Nm =
        windEmulatorStep4_WECSim_B.Switch2_h;
      if (tmp_g) {
        /* RateLimiter: '<S498>/Rate Limiter' incorporates:
         *  Constant: '<S498>/Constant1'
         */
        riseValLimit = 2000.0 * tmp_p;
        rateLimiterRate = riseValLimit - windEmulatorStep4_WECSim_DW.PrevY_k;
        if (rateLimiterRate >
            windEmulatorStep4_WECSim_cal->RateLimiter_RisingLim *
            windEmulatorStep4_WECSim_period) {
          /* RateLimiter: '<S498>/Rate Limiter' */
          windEmulatorStep4_WECSim_B.RateLimiter =
            windEmulatorStep4_WECSim_cal->RateLimiter_RisingLim *
            windEmulatorStep4_WECSim_period +
            windEmulatorStep4_WECSim_DW.PrevY_k;
        } else if (rateLimiterRate <
                   windEmulatorStep4_WECSim_cal->RateLimiter_FallingLim *
                   windEmulatorStep4_WECSim_period) {
          /* RateLimiter: '<S498>/Rate Limiter' */
          windEmulatorStep4_WECSim_B.RateLimiter =
            windEmulatorStep4_WECSim_cal->RateLimiter_FallingLim *
            windEmulatorStep4_WECSim_period +
            windEmulatorStep4_WECSim_DW.PrevY_k;
        } else {
          /* RateLimiter: '<S498>/Rate Limiter' */
          windEmulatorStep4_WECSim_B.RateLimiter = 2000.0 * tmp_p;
        }

        windEmulatorStep4_WECSim_DW.PrevY_k =
          windEmulatorStep4_WECSim_B.RateLimiter;

        /* SimscapeInputBlock: '<S541>/INPUT_1_1_1' */
        windEmulatorStep4_WECSim_B.INPUT_1_1_1[0] =
          windEmulatorStep4_WECSim_B.RateLimiter;
        windEmulatorStep4_WECSim_B.INPUT_1_1_1[1] = 0.0;
        windEmulatorStep4_WECSim_B.INPUT_1_1_1[2] = 0.0;
        windEmulatorStep4_WECSim_DW.INPUT_1_1_1_Discrete_142920204[0] =
          !(windEmulatorStep4_WECSim_B.INPUT_1_1_1[0] ==
            windEmulatorStep4_WECSim_DW.INPUT_1_1_1_Discrete_142920204[1]);
        windEmulatorStep4_WECSim_DW.INPUT_1_1_1_Discrete_142920204[1] =
          windEmulatorStep4_WECSim_B.INPUT_1_1_1[0];
        windEmulatorStep4_WECSim_B.INPUT_1_1_1[0] =
          windEmulatorStep4_WECSim_DW.INPUT_1_1_1_Discrete_142920204[1];
        windEmulatorStep4_WECSim_B.INPUT_1_1_1[3] =
          windEmulatorStep4_WECSim_DW.INPUT_1_1_1_Discrete_142920204[0];
      }

      /* StateSpace: '<S531>/Internal' */
      windEmulatorStep4_WECSim_B.Internal = 0.0;

      /* StateSpace: '<S531>/Internal' */
      for (q0 = windEmulatorStep4_WECSim_cal->Internal_C_jc[0U]; q0 <
           windEmulatorStep4_WECSim_cal->Internal_C_jc[1U]; q0++) {
        /* StateSpace: '<S531>/Internal' */
        windEmulatorStep4_WECSim_B.Internal +=
          windEmulatorStep4_WECSim_cal->Internal_C_pr *
          windEmulatorStep4_WECSim_X.Internal_CSTATE[0U];
      }

      for (q0 = windEmulatorStep4_WECSim_cal->Internal_C_jc[1U]; q0 <
           windEmulatorStep4_WECSim_cal->Internal_C_jc[2U]; q0++) {
        /* StateSpace: '<S531>/Internal' */
        windEmulatorStep4_WECSim_B.Internal +=
          windEmulatorStep4_WECSim_cal->Internal_C_pr *
          windEmulatorStep4_WECSim_X.Internal_CSTATE[1U];
      }

      for (q0 = windEmulatorStep4_WECSim_cal->Internal_C_jc[2U]; q0 <
           windEmulatorStep4_WECSim_cal->Internal_C_jc[3U]; q0++) {
        /* StateSpace: '<S531>/Internal' */
        windEmulatorStep4_WECSim_B.Internal +=
          windEmulatorStep4_WECSim_cal->Internal_C_pr *
          windEmulatorStep4_WECSim_X.Internal_CSTATE[2U];
      }

      /* SimscapeInputBlock: '<S541>/INPUT_2_1_1' */
      windEmulatorStep4_WECSim_B.INPUT_2_1_1[0] =
        windEmulatorStep4_WECSim_B.Internal;
      windEmulatorStep4_WECSim_B.INPUT_2_1_1[1] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_2_1_1[2] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_2_1_1[3] = 0.0;

      /* StateSpace: '<S545>/Internal' */
      windEmulatorStep4_WECSim_B.Internal_j = 0.0;

      /* StateSpace: '<S545>/Internal' */
      for (q0 = windEmulatorStep4_WECSim_cal->Internal_C_jc_k[0U]; q0 <
           windEmulatorStep4_WECSim_cal->Internal_C_jc_k[1U]; q0++) {
        /* StateSpace: '<S545>/Internal' */
        windEmulatorStep4_WECSim_B.Internal_j +=
          windEmulatorStep4_WECSim_cal->Internal_C_pr_a *
          windEmulatorStep4_WECSim_X.Internal_CSTATE_j;
      }

      /* Gain: '<S510>/Gain' */
      windEmulatorStep4_WECSim_B.Gain_l =
        windEmulatorStep4_WECSim_cal->Gain_Gain *
        windEmulatorStep4_WECSim_B.Internal_j;

      /* SimscapeInputBlock: '<S541>/INPUT_4_1_1' */
      windEmulatorStep4_WECSim_B.INPUT_4_1_1[0] =
        windEmulatorStep4_WECSim_B.Gain_l;
      windEmulatorStep4_WECSim_B.INPUT_4_1_1[1] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_4_1_1[2] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_4_1_1[3] = 0.0;
      if (tmp_g) {
        /* RateLimiter: '<S7>/Rate Limiter' incorporates:
         *  Constant: '<S7>/speedReference'
         */
        rateLimiterRate = windEmulatorStep4_WECSim_cal->speedReference_Value -
          windEmulatorStep4_WECSim_DW.PrevY_f;
        if (rateLimiterRate >
            windEmulatorStep4_WECSim_cal->RateLimiter_RisingLim_c *
            windEmulatorStep4_WECSim_period) {
          /* RateLimiter: '<S7>/Rate Limiter' */
          windEmulatorStep4_WECSim_B.RateLimiter_j =
            windEmulatorStep4_WECSim_cal->RateLimiter_RisingLim_c *
            windEmulatorStep4_WECSim_period +
            windEmulatorStep4_WECSim_DW.PrevY_f;
        } else if (rateLimiterRate <
                   windEmulatorStep4_WECSim_cal->RateLimiter_FallingLim_e *
                   windEmulatorStep4_WECSim_period) {
          /* RateLimiter: '<S7>/Rate Limiter' */
          windEmulatorStep4_WECSim_B.RateLimiter_j =
            windEmulatorStep4_WECSim_cal->RateLimiter_FallingLim_e *
            windEmulatorStep4_WECSim_period +
            windEmulatorStep4_WECSim_DW.PrevY_f;
        } else {
          /* RateLimiter: '<S7>/Rate Limiter' */
          windEmulatorStep4_WECSim_B.RateLimiter_j =
            windEmulatorStep4_WECSim_cal->speedReference_Value;
        }

        windEmulatorStep4_WECSim_DW.PrevY_f =
          windEmulatorStep4_WECSim_B.RateLimiter_j;

        /* End of RateLimiter: '<S7>/Rate Limiter' */

        /* Gain: '<S7>/f->w' incorporates:
         *  Constant: '<S7>/excForceFreq_Hz'
         */
        windEmulatorStep4_WECSim_B.fw = windEmulatorStep4_WECSim_cal->fw_Gain *
          windEmulatorStep4_WECSim_cal->excForceFreq_Hz_Value;

        /* Product: '<S7>/Product4' */
        windEmulatorStep4_WECSim_B.Product4 = windEmulatorStep4_WECSim_B.fw *
          windEmulatorStep4_WECSim_B.BusAssignment_b.time;

        /* Trigonometry: '<S7>/Sin2' */
        windEmulatorStep4_WECSim_B.Sin2 = std::sin
          (windEmulatorStep4_WECSim_B.Product4);

        /* Gain: '<S7>/excForceAmpNow_N' incorporates:
         *  Constant: '<S7>/excForceAmp_N'
         */
        windEmulatorStep4_WECSim_B.excForceAmpNow_N =
          windEmulatorStep4_WECSim_cal->excForceAmpNow_N_Gain *
          windEmulatorStep4_WECSim_cal->excForceAmp_N_Value;

        /* Product: '<S7>/ExcitationForce_N' */
        windEmulatorStep4_WECSim_B.waveBotExcitationForce_N =
          windEmulatorStep4_WECSim_B.Sin2 *
          windEmulatorStep4_WECSim_B.excForceAmpNow_N;

        /* Switch: '<S7>/Switch' */
        if (windEmulatorStep4_WECSim_B.BusAssignment_b.runHil) {
          /* Switch: '<S7>/Switch' */
          windEmulatorStep4_WECSim_B.rampValue_a =
            windEmulatorStep4_WECSim_B.BusAssignment_b.ramp;
        } else {
          /* Switch: '<S7>/Switch' incorporates:
           *  Constant: '<S7>/Constant1'
           */
          windEmulatorStep4_WECSim_B.rampValue_a =
            windEmulatorStep4_WECSim_cal->Constant1_Value;
        }

        /* End of Switch: '<S7>/Switch' */

        /* Product: '<S7>/Product' */
        windEmulatorStep4_WECSim_B.Product_m =
          windEmulatorStep4_WECSim_B.waveBotExcitationForce_N *
          windEmulatorStep4_WECSim_B.rampValue_a;

        /* S-Function (slecatpdorx): '<S8>/EtherCAT PDO Receive8' */
        {
          /*------------ S-Function Block: <S8>/EtherCAT PDO Receive8 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.ACS800DcBusVoltage;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 616;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 4 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (16 == 8) && (2 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (16 == 8) && (2 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                 sigIdx*2, 2);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 16;
          }
        }

        /* DataTypeConversion: '<S8>/Data Type Conversion4' */
        windEmulatorStep4_WECSim_B.DataTypeConversion4 =
          windEmulatorStep4_WECSim_B.ACS800DcBusVoltage;

        /* Gain: '<S8>/Gain3' */
        windEmulatorStep4_WECSim_B.Gain3_g = *get_acs800DcBusVoltsScaling() *
          windEmulatorStep4_WECSim_B.DataTypeConversion4;

        /* S-Function (slecatpdorx): '<S8>/EtherCAT PDO Receive9' */
        {
          /*------------ S-Function Block: <S8>/EtherCAT PDO Receive9 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.ACS800Frequency;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 600;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 4 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (16 == 8) && (2 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (16 == 8) && (2 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                 sigIdx*2, 2);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 16;
          }
        }

        /* DataTypeConversion: '<S8>/Data Type Conversion3' */
        windEmulatorStep4_WECSim_B.DataTypeConversion3 =
          windEmulatorStep4_WECSim_B.ACS800Frequency;

        /* Gain: '<S8>/Gain4' */
        windEmulatorStep4_WECSim_B.Gain4_p = *get_acs800FreqScaling() *
          windEmulatorStep4_WECSim_B.DataTypeConversion3;

        /* S-Function (slecatpdorx): '<S8>/EtherCAT PDO Receive10' */
        {
          /*------------ S-Function Block: <S8>/EtherCAT PDO Receive10 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.ACS800Temperature;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 632;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 4 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (16 == 8) && (2 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (16 == 8) && (2 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                 sigIdx*2, 2);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 16;
          }
        }

        /* DataTypeConversion: '<S8>/Data Type Conversion5' */
        windEmulatorStep4_WECSim_B.DataTypeConversion5 =
          windEmulatorStep4_WECSim_B.ACS800Temperature;

        /* Gain: '<S8>/Gain5' */
        windEmulatorStep4_WECSim_B.Gain5_g = *get_acs800TempScaling() *
          windEmulatorStep4_WECSim_B.DataTypeConversion5;

        /* S-Function (slecatpdorx): '<S8>/EtherCAT PDO Receive11' */
        {
          /*------------ S-Function Block: <S8>/EtherCAT PDO Receive11 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.ACS800ActualFeedback;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 536;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 4 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (16 == 8) && (2 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (16 == 8) && (2 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                 sigIdx*2, 2);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 16;
          }
        }

        /* DataTypeConversion: '<S8>/Data Type Conversion' */
        windEmulatorStep4_WECSim_B.DataTypeConversion =
          windEmulatorStep4_WECSim_B.ACS800ActualFeedback;

        /* Gain: '<S8>/Gain' */
        windEmulatorStep4_WECSim_B.Gain_c = *get_acs800SpeedScaling() *
          windEmulatorStep4_WECSim_B.DataTypeConversion;

        /* S-Function (slecatpdorx): '<S8>/EtherCAT PDO Receive12' */
        {
          /*------------ S-Function Block: <S8>/EtherCAT PDO Receive12 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.ACS800Torque;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 568;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 4 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (16 == 8) && (2 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (16 == 8) && (2 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                 sigIdx*2, 2);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 16;
          }
        }

        /* DataTypeConversion: '<S8>/Data Type Conversion1' */
        windEmulatorStep4_WECSim_B.DataTypeConversion1_e =
          windEmulatorStep4_WECSim_B.ACS800Torque;

        /* Gain: '<S8>/Gain1' */
        windEmulatorStep4_WECSim_B.Gain1_f = *get_acs800TorqueScaling() *
          windEmulatorStep4_WECSim_B.DataTypeConversion1_e;

        /* S-Function (slecatpdorx): '<S8>/EtherCAT PDO Receive13' */
        {
          /*------------ S-Function Block: <S8>/EtherCAT PDO Receive13 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.ACS800Power;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 584;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 4 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (16 == 8) && (2 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (16 == 8) && (2 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                 sigIdx*2, 2);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 16;
          }
        }

        /* DataTypeConversion: '<S8>/Data Type Conversion2' */
        windEmulatorStep4_WECSim_B.DataTypeConversion2_i =
          windEmulatorStep4_WECSim_B.ACS800Power;

        /* Gain: '<S8>/Gain2' */
        windEmulatorStep4_WECSim_B.Gain2_d = *get_acs800PowerScaling() *
          windEmulatorStep4_WECSim_B.DataTypeConversion2_i;

        /* S-Function (slecatpdorx): '<S8>/EtherCAT PDO Receive7' */
        {
          /*------------ S-Function Block: <S8>/EtherCAT PDO Receive7 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.ACS800StatusWord;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 520;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 5 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (16 == 8) && (2 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (16 == 8) && (2 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                 sigIdx*2, 2);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (16 == 16) && (2 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (16 == 32) && (2 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 16, sigOutputPtr+
                                   sigIdx*2, 2);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 16;
          }
        }

        /* BusAssignment: '<S8>/Bus Assignment' */
        windEmulatorStep4_WECSim_B.BusAssignment_h.dcBusVoltage_V =
          windEmulatorStep4_WECSim_B.Gain3_g;
        windEmulatorStep4_WECSim_B.BusAssignment_h.frequency_Hz =
          windEmulatorStep4_WECSim_B.Gain4_p;
        windEmulatorStep4_WECSim_B.BusAssignment_h.temperature =
          windEmulatorStep4_WECSim_B.Gain5_g;
        windEmulatorStep4_WECSim_B.BusAssignment_h.motorSpeed_rpm =
          windEmulatorStep4_WECSim_B.Gain_c;
        windEmulatorStep4_WECSim_B.BusAssignment_h.motorTorque_Nm =
          windEmulatorStep4_WECSim_B.Gain1_f;
        windEmulatorStep4_WECSim_B.BusAssignment_h.shaftPower_W =
          windEmulatorStep4_WECSim_B.Gain2_d;
        windEmulatorStep4_WECSim_B.BusAssignment_h.statusWord =
          windEmulatorStep4_WECSim_B.ACS800StatusWord;

        /* BusAssignment: '<S7>/Bus Assignment' */
        windEmulatorStep4_WECSim_B.BusAssignment_c.speedRef_rpm =
          windEmulatorStep4_WECSim_B.RateLimiter_j;
        windEmulatorStep4_WECSim_B.BusAssignment_c.excForce_N =
          windEmulatorStep4_WECSim_B.Product_m;
        windEmulatorStep4_WECSim_B.BusAssignment_c.genSpeedActual =
          windEmulatorStep4_WECSim_B.BusAssignment_h.motorSpeed_rpm;
        windEmulatorStep4_WECSim_B.BusAssignment_c.speedCtrlReset =
          windEmulatorStep4_WECSim_B.BusAssignment_b.resetHilIntegrator;

        /* Gain: '<S366>/Gain' */
        windEmulatorStep4_WECSim_B.Gain_d = tmp_p *
          windEmulatorStep4_WECSim_B.BusAssignment_c.genSpeedActual;

        /* SimscapeInputBlock: '<S541>/INPUT_5_1_1' */
        windEmulatorStep4_WECSim_B.INPUT_5_1_1[0] =
          windEmulatorStep4_WECSim_B.Gain_d;
        windEmulatorStep4_WECSim_B.INPUT_5_1_1[1] = 0.0;
        windEmulatorStep4_WECSim_B.INPUT_5_1_1[2] = 0.0;
        windEmulatorStep4_WECSim_DW.INPUT_5_1_1_Discrete_4244925644[0] =
          !(windEmulatorStep4_WECSim_B.INPUT_5_1_1[0] ==
            windEmulatorStep4_WECSim_DW.INPUT_5_1_1_Discrete_4244925644[1]);
        windEmulatorStep4_WECSim_DW.INPUT_5_1_1_Discrete_4244925644[1] =
          windEmulatorStep4_WECSim_B.INPUT_5_1_1[0];
        windEmulatorStep4_WECSim_B.INPUT_5_1_1[0] =
          windEmulatorStep4_WECSim_DW.INPUT_5_1_1_Discrete_4244925644[1];
        windEmulatorStep4_WECSim_B.INPUT_5_1_1[3] =
          windEmulatorStep4_WECSim_DW.INPUT_5_1_1_Discrete_4244925644[0];
      }

      /* StateSpace: '<S542>/Internal' */
      windEmulatorStep4_WECSim_B.Internal_h = 0.0;

      /* StateSpace: '<S542>/Internal' */
      for (q0 = windEmulatorStep4_WECSim_cal->Internal_C_jc_h[0U]; q0 <
           windEmulatorStep4_WECSim_cal->Internal_C_jc_h[1U]; q0++) {
        /* StateSpace: '<S542>/Internal' */
        windEmulatorStep4_WECSim_B.Internal_h +=
          windEmulatorStep4_WECSim_cal->Internal_C_pr_m *
          windEmulatorStep4_WECSim_X.Internal_CSTATE_a;
      }

      /* Gain: '<S509>/Gain' */
      windEmulatorStep4_WECSim_B.Gain_g =
        windEmulatorStep4_WECSim_cal->Gain_Gain_j *
        windEmulatorStep4_WECSim_B.Internal_h;

      /* SimscapeInputBlock: '<S541>/INPUT_3_1_1' */
      windEmulatorStep4_WECSim_B.INPUT_3_1_1[0] =
        windEmulatorStep4_WECSim_B.Gain_g;
      windEmulatorStep4_WECSim_B.INPUT_3_1_1[1] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_3_1_1[2] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_3_1_1[3] = 0.0;
      if (tmp_g) {
        /* SimscapeRtp: '<S508>/RTP_1' incorporates:
         *  Constant: '<S366>/Subsystem_around_RTP_D290B913_fluid_volume'
         *  Constant: '<S498>/Subsystem_around_RTP_D2E1D090_liquid_pressure'
         *  Constant: '<S498>/Subsystem_around_RTP_D2E1D090_liquid_volume'
         */
        if (windEmulatorStep4_WECSim_DW.RTP_1_SetParametersNeeded) {
          tmp[0] = windEmulatorStep4_WECSim_cal->RTP_D290B913_fluid_volume_Value;
          tmp[1] = windEmulatorStep4_WECSim_cal->RTP_D2E1D090_liquid_pressure_Va;
          tmp[2] = windEmulatorStep4_WECSim_cal->RTP_D2E1D090_liquid_volume_Valu;
          parameterBundle_mRealParameters = &tmp[0];
          rtpManager = static_cast<NeslRtpManager *>
            (windEmulatorStep4_WECSim_DW.RTP_1_RtpManager);
          diag = rtw_create_diagnostics();
          diagTree = neu_diagnostic_manager_get_initial_tree(diag);
          expl_temp.mRealParameters.mN = 3;
          expl_temp.mRealParameters.mX = parameterBundle_mRealParameters;
          expl_temp.mLogicalParameters.mN = 0;
          expl_temp.mLogicalParameters.mX = NULL;
          expl_temp.mIntegerParameters.mN = 0;
          expl_temp.mIntegerParameters.mX = NULL;
          expl_temp.mIndexParameters.mN = 0;
          expl_temp.mIndexParameters.mX = NULL;
          f = nesl_rtp_manager_set_rtps(rtpManager,
            windEmulatorStep4_WECSim_M->Timing.t[0], expl_temp, diag);
          if (!f) {
            f = error_buffer_is_empty(rtmGetErrorStatus
              (windEmulatorStep4_WECSim_M));
            if (f) {
              msg = rtw_diagnostics_msg(diagTree);
              rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
            }
          }
        }

        windEmulatorStep4_WECSim_DW.RTP_1_SetParametersNeeded = false;

        /* End of SimscapeRtp: '<S508>/RTP_1' */

        /* SimscapeExecutionBlock: '<S541>/STATE_1' incorporates:
         *  SimscapeExecutionBlock: '<S541>/OUTPUT_1_0'
         */
        simulationData = static_cast<NeslSimulationData *>
          (windEmulatorStep4_WECSim_DW.STATE_1_SimData);
        u1 = windEmulatorStep4_WECSim_M->Timing.t[0];
        time = u1;
        simulationData->mData->mTime.mN = 1;
        simulationData->mData->mTime.mX = &time;
        simulationData->mData->mContStates.mN = 0;
        simulationData->mData->mContStates.mX = NULL;
        simulationData->mData->mDiscStates.mN = 22;
        simulationData->mData->mDiscStates.mX =
          &windEmulatorStep4_WECSim_DW.STATE_1_Discrete_1041191992[0];
        simulationData->mData->mModeVector.mN = 15;
        simulationData->mData->mModeVector.mX =
          &windEmulatorStep4_WECSim_DW.STATE_1_Modes[0];
        f = false;
        simulationData->mData->mFoundZcEvents = f;
        simulationData->mData->mHadEvents = false;
        simulationData->mData->mIsMajorTimeStep = true;
        f = false;
        simulationData->mData->mIsSolverAssertCheck = f;
        simulationData->mData->mIsSolverCheckingCIC = false;
        simulationData->mData->mIsComputingJacobian = false;
        simulationData->mData->mIsEvaluatingF0 = false;
        simulationData->mData->mIsSolverRequestingReset = false;
        simulationData->mData->mIsModeUpdateTimeStep = true;
        tmp_1[0] = 0;
        tmp_0[0] = windEmulatorStep4_WECSim_B.INPUT_1_1_1[0];
        tmp_0[1] = windEmulatorStep4_WECSim_B.INPUT_1_1_1[1];
        tmp_0[2] = windEmulatorStep4_WECSim_B.INPUT_1_1_1[2];
        tmp_0[3] = windEmulatorStep4_WECSim_B.INPUT_1_1_1[3];
        tmp_1[1] = 4;
        tmp_0[4] = windEmulatorStep4_WECSim_B.INPUT_2_1_1[0];
        tmp_0[5] = windEmulatorStep4_WECSim_B.INPUT_2_1_1[1];
        tmp_0[6] = windEmulatorStep4_WECSim_B.INPUT_2_1_1[2];
        tmp_0[7] = windEmulatorStep4_WECSim_B.INPUT_2_1_1[3];
        tmp_1[2] = 8;
        tmp_0[8] = windEmulatorStep4_WECSim_B.INPUT_4_1_1[0];
        tmp_0[9] = windEmulatorStep4_WECSim_B.INPUT_4_1_1[1];
        tmp_0[10] = windEmulatorStep4_WECSim_B.INPUT_4_1_1[2];
        tmp_0[11] = windEmulatorStep4_WECSim_B.INPUT_4_1_1[3];
        tmp_1[3] = 12;
        tmp_0[12] = windEmulatorStep4_WECSim_B.INPUT_5_1_1[0];
        tmp_0[13] = windEmulatorStep4_WECSim_B.INPUT_5_1_1[1];
        tmp_0[14] = windEmulatorStep4_WECSim_B.INPUT_5_1_1[2];
        tmp_0[15] = windEmulatorStep4_WECSim_B.INPUT_5_1_1[3];
        tmp_1[4] = 16;
        tmp_0[16] = windEmulatorStep4_WECSim_B.INPUT_3_1_1[0];
        tmp_0[17] = windEmulatorStep4_WECSim_B.INPUT_3_1_1[1];
        tmp_0[18] = windEmulatorStep4_WECSim_B.INPUT_3_1_1[2];
        tmp_0[19] = windEmulatorStep4_WECSim_B.INPUT_3_1_1[3];
        tmp_1[5] = 20;
        simulationData->mData->mInputValues.mN = 20;
        simulationData->mData->mInputValues.mX = &tmp_0[0];
        simulationData->mData->mInputOffsets.mN = 6;
        simulationData->mData->mInputOffsets.mX = &tmp_1[0];
        simulationData->mData->mOutputs.mN = 37;
        simulationData->mData->mOutputs.mX =
          &windEmulatorStep4_WECSim_B.STATE_1[0];
        simulationData->mData->mTolerances.mN = 0;
        simulationData->mData->mTolerances.mX = NULL;
        simulationData->mData->mCstateHasChanged = false;
        simulationData->mData->mDstateHasChanged = false;
        deltaT_tmp = windEmulatorStep4_WECSim_M->Timing.t[1];
        time_0 = deltaT_tmp;
        simulationData->mData->mTime.mN = 1;
        simulationData->mData->mTime.mX = &time_0;
        isHit = 0;
        simulationData->mData->mSampleHits.mN = 1;
        simulationData->mData->mSampleHits.mX = &isHit;
        simulationData->mData->mIsFundamentalSampleHit = true;
        simulationData->mData->mHadEvents = false;
        simulator = static_cast<NeslSimulator *>
          (windEmulatorStep4_WECSim_DW.STATE_1_Simulator);
        diag = static_cast<NeuDiagnosticManager *>
          (windEmulatorStep4_WECSim_DW.STATE_1_DiagMgr);
        diagTree = neu_diagnostic_manager_get_initial_tree(diag);
        i = ne_simulator_method(simulator, NESL_SIM_OUTPUTS, simulationData,
          diag);
        if (i != 0) {
          f = error_buffer_is_empty(rtmGetErrorStatus(windEmulatorStep4_WECSim_M));
          if (f) {
            msg = rtw_diagnostics_msg(diagTree);
            rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
          }
        }

        /* End of SimscapeExecutionBlock: '<S541>/STATE_1' */

        /* SimscapeExecutionBlock: '<S541>/OUTPUT_1_0' */
        simulationData = static_cast<NeslSimulationData *>
          (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_SimData);
        time_1 = u1;
        simulationData->mData->mTime.mN = 1;
        simulationData->mData->mTime.mX = &time_1;
        simulationData->mData->mContStates.mN = 0;
        simulationData->mData->mContStates.mX = NULL;
        simulationData->mData->mDiscStates.mN = 0;
        simulationData->mData->mDiscStates.mX =
          &windEmulatorStep4_WECSim_DW.OUTPUT_1_0_Discrete;
        simulationData->mData->mModeVector.mN = 0;
        simulationData->mData->mModeVector.mX =
          &windEmulatorStep4_WECSim_DW.OUTPUT_1_0_Modes;
        f = false;
        simulationData->mData->mFoundZcEvents = f;
        simulationData->mData->mHadEvents = false;
        simulationData->mData->mIsMajorTimeStep = true;
        f = false;
        simulationData->mData->mIsSolverAssertCheck = f;
        simulationData->mData->mIsSolverCheckingCIC = false;
        simulationData->mData->mIsComputingJacobian = false;
        simulationData->mData->mIsEvaluatingF0 = false;
        simulationData->mData->mIsSolverRequestingReset = false;
        simulationData->mData->mIsModeUpdateTimeStep = true;
        tmp_3[0] = 0;
        tmp_2[0] = windEmulatorStep4_WECSim_B.INPUT_1_1_1[0];
        tmp_2[1] = windEmulatorStep4_WECSim_B.INPUT_1_1_1[1];
        tmp_2[2] = windEmulatorStep4_WECSim_B.INPUT_1_1_1[2];
        tmp_2[3] = windEmulatorStep4_WECSim_B.INPUT_1_1_1[3];
        tmp_3[1] = 4;
        tmp_2[4] = windEmulatorStep4_WECSim_B.INPUT_2_1_1[0];
        tmp_2[5] = windEmulatorStep4_WECSim_B.INPUT_2_1_1[1];
        tmp_2[6] = windEmulatorStep4_WECSim_B.INPUT_2_1_1[2];
        tmp_2[7] = windEmulatorStep4_WECSim_B.INPUT_2_1_1[3];
        tmp_3[2] = 8;
        tmp_2[8] = windEmulatorStep4_WECSim_B.INPUT_4_1_1[0];
        tmp_2[9] = windEmulatorStep4_WECSim_B.INPUT_4_1_1[1];
        tmp_2[10] = windEmulatorStep4_WECSim_B.INPUT_4_1_1[2];
        tmp_2[11] = windEmulatorStep4_WECSim_B.INPUT_4_1_1[3];
        tmp_3[3] = 12;
        tmp_2[12] = windEmulatorStep4_WECSim_B.INPUT_5_1_1[0];
        tmp_2[13] = windEmulatorStep4_WECSim_B.INPUT_5_1_1[1];
        tmp_2[14] = windEmulatorStep4_WECSim_B.INPUT_5_1_1[2];
        tmp_2[15] = windEmulatorStep4_WECSim_B.INPUT_5_1_1[3];
        tmp_3[4] = 16;
        tmp_2[16] = windEmulatorStep4_WECSim_B.INPUT_3_1_1[0];
        tmp_2[17] = windEmulatorStep4_WECSim_B.INPUT_3_1_1[1];
        tmp_2[18] = windEmulatorStep4_WECSim_B.INPUT_3_1_1[2];
        tmp_2[19] = windEmulatorStep4_WECSim_B.INPUT_3_1_1[3];
        tmp_3[5] = 20;
        std::memcpy(&tmp_2[20], &windEmulatorStep4_WECSim_B.STATE_1[0], 37U *
                    sizeof(real_T));
        tmp_3[6] = 57;
        simulationData->mData->mInputValues.mN = 57;
        simulationData->mData->mInputValues.mX = &tmp_2[0];
        simulationData->mData->mInputOffsets.mN = 7;
        simulationData->mData->mInputOffsets.mX = &tmp_3[0];
        simulationData->mData->mOutputs.mN = 11;
        simulationData->mData->mOutputs.mX =
          &windEmulatorStep4_WECSim_B.OUTPUT_1_0[0];
        simulationData->mData->mTolerances.mN = 0;
        simulationData->mData->mTolerances.mX = NULL;
        simulationData->mData->mCstateHasChanged = false;
        simulationData->mData->mDstateHasChanged = false;
        time_2 = deltaT_tmp;
        simulationData->mData->mTime.mN = 1;
        simulationData->mData->mTime.mX = &time_2;
        isHit_0 = 0;
        simulationData->mData->mSampleHits.mN = 1;
        simulationData->mData->mSampleHits.mX = &isHit_0;
        simulationData->mData->mIsFundamentalSampleHit = true;
        simulationData->mData->mHadEvents = false;
        simulator = static_cast<NeslSimulator *>
          (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_Simulator);
        diag = static_cast<NeuDiagnosticManager *>
          (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_DiagMgr);
        diagTree = neu_diagnostic_manager_get_initial_tree(diag);
        i = ne_simulator_method(simulator, NESL_SIM_OUTPUTS, simulationData,
          diag);
        if (i != 0) {
          f = error_buffer_is_empty(rtmGetErrorStatus(windEmulatorStep4_WECSim_M));
          if (f) {
            msg = rtw_diagnostics_msg(diagTree);
            rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
          }
        }

        /* Gain: '<S504>/Gain' */
        windEmulatorStep4_WECSim_B.Pressure =
          windEmulatorStep4_WECSim_cal->Gain_Gain_c *
          windEmulatorStep4_WECSim_B.OUTPUT_1_0[8];

        /* Gain: '<S6>/psi -> bar' */
        windEmulatorStep4_WECSim_B.psibar = *get_psi2bar() *
          windEmulatorStep4_WECSim_B.Pressure;

        /* Gain: '<S511>/Gain' */
        windEmulatorStep4_WECSim_B.ShaftSpeedPump = tmp_10 *
          windEmulatorStep4_WECSim_B.OUTPUT_1_0[9];
      }

      /* Step: '<S370>/Step' incorporates:
       *  RateLimiter: '<S2>/acs880RateLim'
       *  RateLimiter: '<S373>/Rate Limiter'
       *  RateLimiter: '<S436>/Rate Limiter'
       *  SimscapeExecutionBlock: '<S216>/OUTPUT_1_0'
       *  SimscapeExecutionBlock: '<S216>/OUTPUT_1_1'
       *  SimscapeExecutionBlock: '<S216>/STATE_1'
       *  SimscapeExecutionBlock: '<S332>/OUTPUT_1_0'
       *  SimscapeExecutionBlock: '<S332>/STATE_1'
       *  Step: '<S58>/Step'
       */
      deltaT_tmp = windEmulatorStep4_WECSim_M->Timing.t[0];
      if (deltaT_tmp < windEmulatorStep4_WECSim_cal->Ramp_start) {
        /* Step: '<S370>/Step' */
        windEmulatorStep4_WECSim_B.Step = windEmulatorStep4_WECSim_cal->Step_Y0;
      } else {
        /* Step: '<S370>/Step' */
        windEmulatorStep4_WECSim_B.Step =
          windEmulatorStep4_WECSim_cal->Ramp_slope;
      }

      /* End of Step: '<S370>/Step' */

      /* Clock: '<S370>/Clock' incorporates:
       *  Clock: '<S57>/Clock'
       *  SimscapeExecutionBlock: '<S216>/STATE_1'
       */
      Clock_tmp = windEmulatorStep4_WECSim_M->Timing.t[0];

      /* Clock: '<S370>/Clock' */
      windEmulatorStep4_WECSim_B.Clock = Clock_tmp;

      /* Sum: '<S370>/Sum' incorporates:
       *  Constant: '<S370>/Constant'
       */
      windEmulatorStep4_WECSim_B.Sum_m = windEmulatorStep4_WECSim_B.Clock -
        windEmulatorStep4_WECSim_cal->Ramp_start;

      /* Product: '<S370>/Product' */
      windEmulatorStep4_WECSim_B.Product_a = windEmulatorStep4_WECSim_B.Step *
        windEmulatorStep4_WECSim_B.Sum_m;

      /* Sum: '<S370>/Output' incorporates:
       *  Constant: '<S370>/Constant1'
       */
      windEmulatorStep4_WECSim_B.Output = windEmulatorStep4_WECSim_B.Product_a +
        windEmulatorStep4_WECSim_cal->Ramp_InitialOutput;

      /* Saturate: '<S365>/Saturation' */
      riseValLimit = windEmulatorStep4_WECSim_B.Output;
      u1 = windEmulatorStep4_WECSim_cal->Saturation_LowerSat_em;
      rateLimiterRate = windEmulatorStep4_WECSim_cal->Saturation_UpperSat_da;
      if (riseValLimit > rateLimiterRate) {
        /* Saturate: '<S365>/Saturation' */
        windEmulatorStep4_WECSim_B.Saturation = rateLimiterRate;
      } else if (riseValLimit < u1) {
        /* Saturate: '<S365>/Saturation' */
        windEmulatorStep4_WECSim_B.Saturation = u1;
      } else {
        /* Saturate: '<S365>/Saturation' */
        windEmulatorStep4_WECSim_B.Saturation = riseValLimit;
      }

      /* End of Saturate: '<S365>/Saturation' */
      if (tmp_g) {
        /* Gain: '<S371>/kDampingNow' incorporates:
         *  Constant: '<S371>/kDamping'
         */
        windEmulatorStep4_WECSim_B.kDampingNow =
          windEmulatorStep4_WECSim_cal->kDampingNow_Gain *
          windEmulatorStep4_WECSim_cal->kDamping_Value;

        /* Product: '<S371>/Product' */
        windEmulatorStep4_WECSim_B.Product_mo =
          windEmulatorStep4_WECSim_B.kDampingNow *
          windEmulatorStep4_WECSim_B.OUTPUT_1_0[7];

        /* Gain: '<S371>/kSpringNow' incorporates:
         *  Constant: '<S371>/kSpring'
         */
        windEmulatorStep4_WECSim_B.kSpringNow =
          windEmulatorStep4_WECSim_cal->kSpringNow_Gain *
          windEmulatorStep4_WECSim_cal->kSpring_Value;

        /* Product: '<S371>/Product1' */
        windEmulatorStep4_WECSim_B.Product1_g =
          windEmulatorStep4_WECSim_B.OUTPUT_1_0[6] *
          windEmulatorStep4_WECSim_B.kSpringNow;

        /* Sum: '<S371>/Add' */
        windEmulatorStep4_WECSim_B.ForceD =
          windEmulatorStep4_WECSim_B.Product_mo +
          windEmulatorStep4_WECSim_B.Product1_g;

        /* Product: '<S367>/Product2' incorporates:
         *  Constant: '<S367>/Constant2'
         */
        riseValLimit = *get_A() * 2.0 * 3.1415926535897931 / *get_T() /
          *get_targetOmega();

        /* Product: '<S367>/Product2' */
        windEmulatorStep4_WECSim_B.Product2 = windEmulatorStep4_WECSim_B.ForceD *
          riseValLimit;

        /* Gain: '<S367>/Gain' */
        windEmulatorStep4_WECSim_B.Gain_n =
          windEmulatorStep4_WECSim_cal->Gain_Gain_cc *
          windEmulatorStep4_WECSim_B.Product2;
      }

      /* Product: '<S365>/Product' */
      windEmulatorStep4_WECSim_B.TorqueInputRef =
        windEmulatorStep4_WECSim_B.Saturation *
        windEmulatorStep4_WECSim_B.Gain_n;

      /* Abs: '<S365>/Abs' */
      windEmulatorStep4_WECSim_B.Abs = std::abs
        (windEmulatorStep4_WECSim_B.TorqueInputRef);

      /* Sum: '<S437>/Sum' */
      windEmulatorStep4_WECSim_B.Sum_f = windEmulatorStep4_WECSim_B.OUTPUT_1_0
        [10] - windEmulatorStep4_WECSim_B.TorqueInputRef;
      if (tmp_g) {
        /* DiscreteIntegrator: '<S437>/Discrete-Time Integrator' */
        windEmulatorStep4_WECSim_B.DiscreteTimeIntegrator =
          windEmulatorStep4_WECSim_DW.DiscreteTimeIntegrator_DSTATE;

        /* Gain: '<S437>/Gain1' */
        windEmulatorStep4_WECSim_B.Gain1_b = tmp_t->IG *
          windEmulatorStep4_WECSim_B.DiscreteTimeIntegrator;
      }

      /* Switch generated from: '<S365>/Switch' */
      if (windEmulatorStep4_WECSim_B.Abs >= tmp_v) {
        /* Switch: '<S374>/Switch' */
        if (windEmulatorStep4_WECSim_B.TorqueInputRef >=
            windEmulatorStep4_WECSim_cal->Switch_Threshold) {
          /* Switch: '<S374>/Switch' incorporates:
           *  Constant: '<S374>/Constant'
           */
          windEmulatorStep4_WECSim_B.Switch_ne =
            windEmulatorStep4_WECSim_cal->Constant_Value_ly;
        } else {
          /* Switch: '<S374>/Switch' incorporates:
           *  Constant: '<S374>/Constant1'
           */
          windEmulatorStep4_WECSim_B.Switch_ne =
            windEmulatorStep4_WECSim_cal->Constant1_Value_b;
        }

        /* End of Switch: '<S374>/Switch' */

        /* Switch generated from: '<S365>/Switch' */
        windEmulatorStep4_WECSim_B.ControlSignal1 =
          windEmulatorStep4_WECSim_B.Switch_ne;
      } else {
        /* Gain: '<S437>/Gain2' */
        windEmulatorStep4_WECSim_B.Gain2_n =
          windEmulatorStep4_WECSim_cal->Gain2_Gain *
          windEmulatorStep4_WECSim_B.TorqueInputRef;

        /* Product: '<S437>/Divide' */
        windEmulatorStep4_WECSim_B.Divide = windEmulatorStep4_WECSim_B.Gain2_n /
          windEmulatorStep4_WECSim_B.Pressure;

        /* Product: '<S437>/Product' incorporates:
         *  Constant: '<S437>/UnitsConversion'
         */
        riseValLimit = 6.283185307179586E+6 / (tmp_o * 6894.75);

        /* Product: '<S437>/Product' */
        windEmulatorStep4_WECSim_B.Product_aq =
          windEmulatorStep4_WECSim_B.Divide * riseValLimit;

        /* Gain: '<S437>/Gain' */
        windEmulatorStep4_WECSim_B.Gain_a = tmp_t->PG *
          windEmulatorStep4_WECSim_B.Sum_f;

        /* Sum: '<S437>/Add' */
        windEmulatorStep4_WECSim_B.Add_iu = windEmulatorStep4_WECSim_B.Gain_a +
          windEmulatorStep4_WECSim_B.Gain1_b;

        /* Sum: '<S437>/Add1' */
        windEmulatorStep4_WECSim_B.Add1_lh = windEmulatorStep4_WECSim_B.Add_iu +
          windEmulatorStep4_WECSim_B.Product_aq;

        /* Saturate: '<S437>/Saturation' */
        riseValLimit = windEmulatorStep4_WECSim_B.Add1_lh;
        u1 = windEmulatorStep4_WECSim_cal->Saturation_LowerSat_b;
        rateLimiterRate = windEmulatorStep4_WECSim_cal->Saturation_UpperSat_mj;
        if (riseValLimit > rateLimiterRate) {
          /* Saturate: '<S437>/Saturation' */
          windEmulatorStep4_WECSim_B.Saturation_j = rateLimiterRate;
        } else if (riseValLimit < u1) {
          /* Saturate: '<S437>/Saturation' */
          windEmulatorStep4_WECSim_B.Saturation_j = u1;
        } else {
          /* Saturate: '<S437>/Saturation' */
          windEmulatorStep4_WECSim_B.Saturation_j = riseValLimit;
        }

        /* End of Saturate: '<S437>/Saturation' */

        /* Switch generated from: '<S365>/Switch' */
        windEmulatorStep4_WECSim_B.ControlSignal1 =
          windEmulatorStep4_WECSim_B.Saturation_j;
      }

      /* Abs: '<S365>/Abs2' */
      windEmulatorStep4_WECSim_B.Abs2 = std::abs
        (windEmulatorStep4_WECSim_B.TorqueInputRef);

      /* Switch: '<S365>/Switch2' */
      if (windEmulatorStep4_WECSim_B.Abs2 >= tmp_v) {
        /* Gain: '<S365>/Gain' */
        riseValLimit = 1.0 / (tmp_o * 1.0E-6 * 0.15915494309189535) / 6894.75;

        /* Gain: '<S365>/Gain' */
        windEmulatorStep4_WECSim_B.Gain_k = riseValLimit *
          windEmulatorStep4_WECSim_B.TorqueInputRef;

        /* Abs: '<S365>/Abs1' */
        windEmulatorStep4_WECSim_B.Abs1_a = std::abs
          (windEmulatorStep4_WECSim_B.Gain_k);

        /* Switch: '<S365>/Switch2' */
        windEmulatorStep4_WECSim_B.PressureRef =
          windEmulatorStep4_WECSim_B.Abs1_a;
      } else {
        /* Switch: '<S365>/Switch2' incorporates:
         *  Constant: '<S53>/Constant1'
         */
        windEmulatorStep4_WECSim_B.PressureRef = *get_minPressureRef_psi();
      }

      /* Sum: '<S372>/Sum' */
      windEmulatorStep4_WECSim_B.Sum_me = windEmulatorStep4_WECSim_B.Pressure -
        windEmulatorStep4_WECSim_B.PressureRef;
      if (tmp_g) {
        /* DiscreteIntegrator: '<S372>/Discrete-Time Integrator' */
        windEmulatorStep4_WECSim_B.DiscreteTimeIntegrator_i =
          windEmulatorStep4_WECSim_DW.DiscreteTimeIntegrator_DSTATE_n;

        /* Gain: '<S372>/Gain1' */
        windEmulatorStep4_WECSim_B.Gain1_j = tmp_s->IG *
          windEmulatorStep4_WECSim_B.DiscreteTimeIntegrator_i;

        /* DiscreteIntegrator: '<S435>/Discrete-Time Integrator' */
        windEmulatorStep4_WECSim_B.DiscreteTimeIntegrator_e =
          windEmulatorStep4_WECSim_DW.DiscreteTimeIntegrator_DSTATE_l;

        /* Gain: '<S435>/Gain1' */
        windEmulatorStep4_WECSim_B.Gain1_h = tmp_s->IG *
          windEmulatorStep4_WECSim_B.DiscreteTimeIntegrator_e;
      }

      /* Sum: '<S435>/Sum' */
      windEmulatorStep4_WECSim_B.Sum_d = windEmulatorStep4_WECSim_B.Pressure -
        windEmulatorStep4_WECSim_B.PressureRef;

      /* Switch generated from: '<S365>/Switch' */
      if (windEmulatorStep4_WECSim_B.Abs >= tmp_v) {
        /* Gain: '<S372>/Gain' */
        windEmulatorStep4_WECSim_B.Gain_ge = tmp_s->PG *
          windEmulatorStep4_WECSim_B.Sum_me;

        /* Sum: '<S372>/Add' */
        windEmulatorStep4_WECSim_B.Add_i = windEmulatorStep4_WECSim_B.Gain_ge +
          windEmulatorStep4_WECSim_B.Gain1_j;

        /* Saturate: '<S372>/Saturation' */
        riseValLimit = windEmulatorStep4_WECSim_B.Add_i;
        u1 = windEmulatorStep4_WECSim_cal->Saturation_LowerSat_g;
        rateLimiterRate = windEmulatorStep4_WECSim_cal->Saturation_UpperSat_d;
        if (riseValLimit > rateLimiterRate) {
          /* Saturate: '<S372>/Saturation' */
          windEmulatorStep4_WECSim_B.Saturation_a = rateLimiterRate;
        } else if (riseValLimit < u1) {
          /* Saturate: '<S372>/Saturation' */
          windEmulatorStep4_WECSim_B.Saturation_a = u1;
        } else {
          /* Saturate: '<S372>/Saturation' */
          windEmulatorStep4_WECSim_B.Saturation_a = riseValLimit;
        }

        /* End of Saturate: '<S372>/Saturation' */

        /* Switch generated from: '<S365>/Switch' */
        windEmulatorStep4_WECSim_B.ControlSignal2 =
          windEmulatorStep4_WECSim_B.Saturation_a;
      } else {
        /* Gain: '<S435>/Gain' */
        windEmulatorStep4_WECSim_B.Gain_j = tmp_s->PG *
          windEmulatorStep4_WECSim_B.Sum_d;

        /* Sum: '<S435>/Add' */
        windEmulatorStep4_WECSim_B.Add_k = windEmulatorStep4_WECSim_B.Gain_j +
          windEmulatorStep4_WECSim_B.Gain1_h;

        /* Saturate: '<S435>/Saturation' */
        riseValLimit = windEmulatorStep4_WECSim_B.Add_k;
        u1 = windEmulatorStep4_WECSim_cal->Saturation_LowerSat_e;
        rateLimiterRate = windEmulatorStep4_WECSim_cal->Saturation_UpperSat_p;
        if (riseValLimit > rateLimiterRate) {
          /* Saturate: '<S435>/Saturation' */
          windEmulatorStep4_WECSim_B.Saturation_f = rateLimiterRate;
        } else if (riseValLimit < u1) {
          /* Saturate: '<S435>/Saturation' */
          windEmulatorStep4_WECSim_B.Saturation_f = u1;
        } else {
          /* Saturate: '<S435>/Saturation' */
          windEmulatorStep4_WECSim_B.Saturation_f = riseValLimit;
        }

        /* End of Saturate: '<S435>/Saturation' */

        /* Switch generated from: '<S365>/Switch' */
        windEmulatorStep4_WECSim_B.ControlSignal2 =
          windEmulatorStep4_WECSim_B.Saturation_f;
      }

      if (tmp_g) {
        /* Gain: '<S501>/Gain' */
        windEmulatorStep4_WECSim_B.FlowMotor1 =
          windEmulatorStep4_WECSim_cal->Gain_Gain_f *
          windEmulatorStep4_WECSim_B.OUTPUT_1_0[2];
      }

      /* BusAssignment: '<S6>/Bus Assignment' */
      windEmulatorStep4_WECSim_B.BusAssignment_n.genTorqueCmd_Nm = 0.0;
      windEmulatorStep4_WECSim_B.BusAssignment_n.pressure_bar =
        windEmulatorStep4_WECSim_B.psibar;
      windEmulatorStep4_WECSim_B.BusAssignment_n.hmOutputShafTorque_Nm = 0.0;
      windEmulatorStep4_WECSim_B.BusAssignment_n.genShaftSpeed_rpm =
        windEmulatorStep4_WECSim_B.BusAssignment_c.genSpeedActual;
      windEmulatorStep4_WECSim_B.BusAssignment_n.excShaftSpeed_rpm =
        windEmulatorStep4_WECSim_B.ShaftSpeedPump;
      windEmulatorStep4_WECSim_B.BusAssignment_n.excShaftTorque_Nm =
        windEmulatorStep4_WECSim_B.OUTPUT_1_0[10];
      windEmulatorStep4_WECSim_B.BusAssignment_n.ctrlSignal1 =
        windEmulatorStep4_WECSim_B.ControlSignal1;
      windEmulatorStep4_WECSim_B.BusAssignment_n.ctrlSignal2 =
        windEmulatorStep4_WECSim_B.ControlSignal2;
      windEmulatorStep4_WECSim_B.BusAssignment_n.genPumpFlow_lpm =
        windEmulatorStep4_WECSim_B.FlowMotor1;

      /* MultiPortSwitch: '<S3>/Multiport Switch' incorporates:
       *  DataTypeConversion: '<S3>/toExpTypeEnum'
       */
      switch (windEmulatorStep4_WECSim_B.toExpTypeEnum) {
       case expTypeEnum_off:
        /* MultiPortSwitch: '<S3>/Multiport Switch' incorporates:
         *  Constant: '<S3>/zeroTorque'
         */
        windEmulatorStep4_WECSim_B.MultiportSwitch_p =
          windEmulatorStep4_WECSim_cal->zeroTorque_Value;
        break;

       case expTypeEnum_sid:
        /* MultiPortSwitch: '<S3>/Multiport Switch' */
        windEmulatorStep4_WECSim_B.MultiportSwitch_p =
          windEmulatorStep4_WECSim_B.BusAssignment_f.acs880Torque_Nm;
        break;

       case expTypeEnum_hil:
        /* MultiPortSwitch: '<S3>/Multiport Switch' */
        windEmulatorStep4_WECSim_B.MultiportSwitch_p =
          windEmulatorStep4_WECSim_B.BusAssignment_n.hmOutputShafTorque_Nm;
        break;

       default:
        /* MultiPortSwitch: '<S3>/Multiport Switch' incorporates:
         *  Constant: '<S3>/zeroTorque'
         */
        windEmulatorStep4_WECSim_B.MultiportSwitch_p =
          windEmulatorStep4_WECSim_cal->zeroTorque_Value;
        break;
      }

      /* End of MultiPortSwitch: '<S3>/Multiport Switch' */

      /* RateLimiter: '<S2>/acs880RateLim' */
      if (windEmulatorStep4_WECSim_DW.LastMajorTime_a == (rtInf)) {
        /* RateLimiter: '<S2>/acs880RateLim' */
        windEmulatorStep4_WECSim_B.acs880RateLim =
          windEmulatorStep4_WECSim_B.MultiportSwitch_p;
      } else {
        u1 = deltaT_tmp - windEmulatorStep4_WECSim_DW.LastMajorTime_a;
        if (windEmulatorStep4_WECSim_DW.LastMajorTime_a == deltaT_tmp) {
          if (windEmulatorStep4_WECSim_DW.PrevLimited_o) {
            /* RateLimiter: '<S2>/acs880RateLim' */
            windEmulatorStep4_WECSim_B.acs880RateLim =
              windEmulatorStep4_WECSim_DW.PrevY_b;
          } else {
            /* RateLimiter: '<S2>/acs880RateLim' */
            windEmulatorStep4_WECSim_B.acs880RateLim =
              windEmulatorStep4_WECSim_B.MultiportSwitch_p;
          }
        } else {
          riseValLimit = u1 * *get_acs880RateLimRising();
          rateLimiterRate = windEmulatorStep4_WECSim_B.MultiportSwitch_p -
            windEmulatorStep4_WECSim_DW.PrevY_b;
          if (rateLimiterRate > riseValLimit) {
            /* RateLimiter: '<S2>/acs880RateLim' */
            windEmulatorStep4_WECSim_B.acs880RateLim =
              windEmulatorStep4_WECSim_DW.PrevY_b + riseValLimit;
            f = true;
          } else {
            u1 *= *get_acs880RateLimFalling();
            if (rateLimiterRate < u1) {
              /* RateLimiter: '<S2>/acs880RateLim' */
              windEmulatorStep4_WECSim_B.acs880RateLim =
                windEmulatorStep4_WECSim_DW.PrevY_b + u1;
              f = true;
            } else {
              /* RateLimiter: '<S2>/acs880RateLim' */
              windEmulatorStep4_WECSim_B.acs880RateLim =
                windEmulatorStep4_WECSim_B.MultiportSwitch_p;
              f = false;
            }
          }

          if (tmp_h) {
            windEmulatorStep4_WECSim_DW.PrevLimited_o = f;
          }
        }
      }

      /* Saturate: '<S2>/Saturation' */
      riseValLimit = windEmulatorStep4_WECSim_B.acs880RateLim;
      u1 = *get_acs880SetpointLimLower();
      rateLimiterRate = *get_acs880SetpointLimUpper();
      if (riseValLimit > rateLimiterRate) {
        /* Saturate: '<S2>/Saturation' */
        windEmulatorStep4_WECSim_B.ACS880Setpoint = rateLimiterRate;
      } else if (riseValLimit < u1) {
        /* Saturate: '<S2>/Saturation' */
        windEmulatorStep4_WECSim_B.ACS880Setpoint = u1;
      } else {
        /* Saturate: '<S2>/Saturation' */
        windEmulatorStep4_WECSim_B.ACS880Setpoint = riseValLimit;
      }

      /* End of Saturate: '<S2>/Saturation' */

      /* Gain: '<S2>/Nm -> %' */
      riseValLimit = 1.0 / tmp_n * 100.0;

      /* Gain: '<S2>/Nm -> %' */
      windEmulatorStep4_WECSim_B.Nm = riseValLimit *
        windEmulatorStep4_WECSim_B.ACS880Setpoint;

      /* BusAssignment: '<S2>/Bus Assignment' */
      windEmulatorStep4_WECSim_B.BusAssignment_k.ctrlWord =
        windEmulatorStep4_WECSim_B.ControlWord;
      windEmulatorStep4_WECSim_B.BusAssignment_k.state =
        windEmulatorStep4_WECSim_B.CastToDouble1_a;
      windEmulatorStep4_WECSim_B.BusAssignment_k.torqueSetpoint_Nm =
        windEmulatorStep4_WECSim_B.ACS880Setpoint;
      windEmulatorStep4_WECSim_B.BusAssignment_k.torqueSetpoint_percent =
        windEmulatorStep4_WECSim_B.Nm;
      if (tmp_g) {
        /* ToAsyncQueueBlock generated from: '<S25>/acs880CtrlSignals' */
        slrtLogSignal
          (windEmulatorStep4_WECSim_DW.TAQSigLogging_InsertedFor_acs_a.SLRTSigHandles,
           (((windEmulatorStep4_WECSim_M->Timing.clockTick1+
              windEmulatorStep4_WECSim_M->Timing.clockTickH1* 4294967296.0)) *
            0.004));
      }

      /* MATLAB Function: '<S46>/parseCtrlWord' */
      windEmulatorStep4_parseCtrlWord
        (windEmulatorStep4_WECSim_B.BusAssignment_k.ctrlWord,
         &windEmulatorStep4_WECSim_B.sf_parseCtrlWord_h,
         &windEmulatorStep4_WECSim_DW.sf_parseCtrlWord_h);

      /* Logic: '<S46>/enableOperation' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_WECSim_B.enableOperation =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord_h.enableOperation != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_e);

      /* Logic: '<S46>/extCtrlLoc' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_WECSim_B.extCtrlLoc_f =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord_h.extCtrlLoc != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_e);

      /* Logic: '<S46>/inching1' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_WECSim_B.inching1 =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord_h.inching1 != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_e);

      /* Logic: '<S46>/inching2' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_WECSim_B.inching2 =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord_h.inching2 != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_e);

      /* Logic: '<S46>/off1Ctrl' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_WECSim_B.off1Ctrl =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord_h.off1Ctrl != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_e);

      /* Logic: '<S46>/off2Ctrl' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_WECSim_B.off2Ctrl =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord_h.off2Ctrl != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_e);

      /* Logic: '<S46>/off3Ctrl' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_WECSim_B.off3Ctrl =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord_h.off3Ctrl != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_e);

      /* Logic: '<S46>/rampHold' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_WECSim_B.rampHold =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord_h.rampHold != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_e);

      /* Logic: '<S46>/rampInZero' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_WECSim_B.rampInZero =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord_h.rampInZero != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_e);

      /* Logic: '<S46>/rampOutZero' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_WECSim_B.rampOutZero =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord_h.rampOutZero != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_e);

      /* Logic: '<S46>/remoteCmd' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_WECSim_B.remoteCmd =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord_h.remoteCmd != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_e);

      /* Logic: '<S46>/reset' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_WECSim_B.reset =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord_h.reset != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_e);

      /* Bias: '<S26>/ctrlWord' */
      windEmulatorStep4_WECSim_B.ctrlWord = static_cast<uint16_T>
        (windEmulatorStep4_WECSim_B.BusAssignment_k.ctrlWord +
         windEmulatorStep4_WECSim_cal->ctrlWord_Bias);

      /* Bias: '<S26>/state' */
      windEmulatorStep4_WECSim_B.state =
        windEmulatorStep4_WECSim_B.BusAssignment_k.state +
        windEmulatorStep4_WECSim_cal->state_Bias;

      /* Gain: '<S26>/torqueSetpoint_Nm' */
      windEmulatorStep4_WECSim_B.torqueSetpoint_Nm =
        windEmulatorStep4_WECSim_cal->torqueSetpoint_Nm_Gain *
        windEmulatorStep4_WECSim_B.BusAssignment_k.torqueSetpoint_Nm;

      /* Gain: '<S26>/torqueSetpoint_percent' */
      windEmulatorStep4_WECSim_B.torqueSetpoint_percent =
        windEmulatorStep4_WECSim_cal->torqueSetpoint_percent_Gain *
        windEmulatorStep4_WECSim_B.BusAssignment_k.torqueSetpoint_percent;
      if (tmp_g) {
        /* ToAsyncQueueBlock generated from: '<S23>/acs800Signals' */
        slrtLogSignal
          (windEmulatorStep4_WECSim_DW.TAQSigLogging_InsertedFor_acs80.SLRTSigHandles,
           (((windEmulatorStep4_WECSim_M->Timing.clockTick1+
              windEmulatorStep4_WECSim_M->Timing.clockTickH1* 4294967296.0)) *
            0.004));

        /* Gain: '<S43>/rpm -> rad//s' */
        windEmulatorStep4_WECSim_B.rpmrads_o = tmp_p *
          windEmulatorStep4_WECSim_B.BusAssignment_h.motorSpeed_rpm;

        /* Product: '<S43>/Product' */
        windEmulatorStep4_WECSim_B.Product_o =
          windEmulatorStep4_WECSim_B.rpmrads_o *
          windEmulatorStep4_WECSim_B.BusAssignment_h.motorTorque_Nm;
        windEmulatorStep4_MovingAverage(windEmulatorStep4_WECSim_B.Product_o,
          &windEmulatorStep4_WECSim_B.MovingAverage,
          &windEmulatorStep4_WECSim_DW.MovingAverage);

        /* Gain: '<S43>/shaftPowerAverage_W' */
        windEmulatorStep4_WECSim_B.shaftPowerAverage_W_a =
          windEmulatorStep4_WECSim_cal->shaftPowerAverage_W_Gain_l *
          windEmulatorStep4_WECSim_B.MovingAverage.MovingAverage;

        /* Gain: '<S43>/shaftPower_W' */
        windEmulatorStep4_WECSim_B.shaftPower_W_e =
          windEmulatorStep4_WECSim_cal->shaftPower_W_Gain_b *
          windEmulatorStep4_WECSim_B.Product_o;

        /* MATLAB Function: '<S44>/Parse Status Word' */
        windEmulatorSte_ParseStatusWord
          (windEmulatorStep4_WECSim_B.BusAssignment_h.statusWord,
           &windEmulatorStep4_WECSim_B.sf_ParseStatusWord,
           &windEmulatorStep4_WECSim_DW.sf_ParseStatusWord);

        /* Logic: '<S44>/aboveLimit' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_WECSim_B.aboveLimit_c =
          (windEmulatorStep4_WECSim_cal->Constant_Value_jd &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord.above_limit != 0));

        /* Logic: '<S44>/alarm' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_WECSim_B.alarm_c =
          (windEmulatorStep4_WECSim_cal->Constant_Value_jd &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord.alarm != 0));

        /* Logic: '<S44>/atSetpoint' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_WECSim_B.atSetpoint_a =
          (windEmulatorStep4_WECSim_cal->Constant_Value_jd &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord.at_setpoint != 0));

        /* Logic: '<S44>/commErr' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_WECSim_B.commErr_i =
          (windEmulatorStep4_WECSim_cal->Constant_Value_jd &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord.comm_err != 0));

        /* Logic: '<S44>/extCtrlLoc' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_WECSim_B.extCtrlLoc_e =
          (windEmulatorStep4_WECSim_cal->Constant_Value_jd &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord.ext_ctrl_loc != 0));

        /* Logic: '<S44>/extRunEnable' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_WECSim_B.extRunEnable_k =
          (windEmulatorStep4_WECSim_cal->Constant_Value_jd &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord.ext_run_enable != 0));

        /* Logic: '<S44>/mswB13' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_WECSim_B.mswB13_a =
          (windEmulatorStep4_WECSim_cal->Constant_Value_jd &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord.msw_b13 != 0));

        /* Logic: '<S44>/mswB14' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_WECSim_B.mswB14_m =
          (windEmulatorStep4_WECSim_cal->Constant_Value_jd &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord.msw_b14 != 0));

        /* Logic: '<S44>/off2' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_WECSim_B.off2_n =
          (windEmulatorStep4_WECSim_cal->Constant_Value_jd &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord.off2 != 0));

        /* Logic: '<S44>/off3' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_WECSim_B.off3_p =
          (windEmulatorStep4_WECSim_cal->Constant_Value_jd &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord.off3 != 0));

        /* Logic: '<S44>/rdyOn' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_WECSim_B.rdyOn_b =
          (windEmulatorStep4_WECSim_cal->Constant_Value_jd &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord.rdy_on != 0));

        /* Logic: '<S44>/rdyRef' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_WECSim_B.rdyRef_f =
          (windEmulatorStep4_WECSim_cal->Constant_Value_jd &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord.rdy_ref != 0));

        /* Logic: '<S44>/rdyRun' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_WECSim_B.rdyRun_n =
          (windEmulatorStep4_WECSim_cal->Constant_Value_jd &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord.rdy_run != 0));

        /* Logic: '<S44>/remote' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_WECSim_B.remote_d =
          (windEmulatorStep4_WECSim_cal->Constant_Value_jd &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord.remote != 0));

        /* Logic: '<S44>/switchOnInhibit' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_WECSim_B.switchOnInhibit_p =
          (windEmulatorStep4_WECSim_cal->Constant_Value_jd &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord.swc_on_inhib != 0));

        /* Logic: '<S44>/tripped' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_WECSim_B.tripped_p =
          (windEmulatorStep4_WECSim_cal->Constant_Value_jd &&
           (windEmulatorStep4_WECSim_B.sf_ParseStatusWord.tripped != 0));

        /* Gain: '<S24>/dcBusVoltage_V' */
        windEmulatorStep4_WECSim_B.dcBusVoltage_V =
          windEmulatorStep4_WECSim_cal->dcBusVoltage_V_Gain *
          windEmulatorStep4_WECSim_B.BusAssignment_h.dcBusVoltage_V;

        /* Gain: '<S24>/frequency_Hz' */
        windEmulatorStep4_WECSim_B.frequency_Hz_j =
          windEmulatorStep4_WECSim_cal->frequency_Hz_Gain_i *
          windEmulatorStep4_WECSim_B.BusAssignment_h.frequency_Hz;

        /* Gain: '<S24>/motorSpeed_rpm' */
        windEmulatorStep4_WECSim_B.motorSpeed_rpm_j =
          windEmulatorStep4_WECSim_cal->motorSpeed_rpm_Gain_j *
          windEmulatorStep4_WECSim_B.BusAssignment_h.motorSpeed_rpm;

        /* Gain: '<S24>/motorTorque_Nm' */
        windEmulatorStep4_WECSim_B.motorTorque_Nm_p =
          windEmulatorStep4_WECSim_cal->motorTorque_Nm_Gain_l *
          windEmulatorStep4_WECSim_B.BusAssignment_h.motorTorque_Nm;

        /* Gain: '<S24>/shaftPower_W' */
        windEmulatorStep4_WECSim_B.shaftPower_W_p =
          windEmulatorStep4_WECSim_cal->shaftPower_W_Gain_l *
          windEmulatorStep4_WECSim_B.BusAssignment_h.shaftPower_W;

        /* Bias: '<S24>/statusWord' */
        windEmulatorStep4_WECSim_B.statusWord_m = static_cast<uint16_T>
          (windEmulatorStep4_WECSim_B.BusAssignment_h.statusWord +
           windEmulatorStep4_WECSim_cal->statusWord_Bias_l);

        /* Gain: '<S24>/temperature' */
        windEmulatorStep4_WECSim_B.temperature =
          windEmulatorStep4_WECSim_cal->temperature_Gain *
          windEmulatorStep4_WECSim_B.BusAssignment_h.temperature;

        /* Memory: '<S1>/Memory' */
        windEmulatorStep4_WECSim_B.Memory_j =
          windEmulatorStep4_WECSim_DW.Memory_PreviousInput_d;

        /* RelationalOperator: '<S1>/NotEqual' incorporates:
         *  Constant: '<S1>/powerUpButton'
         */
        windEmulatorStep4_WECSim_B.NotEqual_f =
          (windEmulatorStep4_WECSim_cal->powerUpButton_Value_b !=
           windEmulatorStep4_WECSim_B.Memory_j);

        /* Memory: '<S1>/Memory1' */
        windEmulatorStep4_WECSim_B.Memory1_f =
          windEmulatorStep4_WECSim_DW.Memory1_PreviousInput_p;

        /* RelationalOperator: '<S1>/NotEqual1' incorporates:
         *  Constant: '<S1>/powerDownButton'
         */
        windEmulatorStep4_WECSim_B.NotEqual1_m =
          (windEmulatorStep4_WECSim_cal->powerDownButton_Value_k !=
           windEmulatorStep4_WECSim_B.Memory1_f);

        /* Memory: '<S1>/Memory2' */
        windEmulatorStep4_WECSim_B.Memory2_h =
          windEmulatorStep4_WECSim_DW.Memory2_PreviousInput_h;

        /* RelationalOperator: '<S1>/NotEqual2' incorporates:
         *  Constant: '<S1>/resetFaultButton'
         */
        windEmulatorStep4_WECSim_B.NotEqual2_h =
          (windEmulatorStep4_WECSim_cal->resetFaultButton_Value_n !=
           windEmulatorStep4_WECSim_B.Memory2_h);

        /* DataTypeConversion: '<S1>/Cast To Double' */
        windEmulatorStep4_WECSim_B.CastToDouble_j =
          windEmulatorStep4_WECSim_B.NotEqual2_h;

        /* Constant: '<S1>/ACS800CtrlMode' */
        windEmulatorStep4_WECSim_B.ACS800CtrlMode = tmp_11;

        /* Chart: '<S16>/ABB Fieldbus Control' */
        if (windEmulatorStep4_WECSim_DW.temporalCounter_i1_g < 31) {
          windEmulatorStep4_WECSim_DW.temporalCounter_i1_g = static_cast<uint8_T>
            (windEmulatorStep4_WECSim_DW.temporalCounter_i1_g + 1);
        }

        windEmulatorStep4_WECSim_DW.sfEvent_d = windEmulatorStep4__CALL_EVENT_k;
        if (windEmulatorStep4_WECSim_DW.is_active_c9_windEmulatorStep4_ == 0) {
          windEmulatorStep4_WECSim_DW.is_active_c9_windEmulatorStep4_ = 1U;
          windEmulatorStep4_WECSim_DW.is_active_UpdateStateMachine_a = 1U;
          windEmulatorStep4_WECSim_DW.is_UpdateStateMachine_g =
            windEmulatorStep4_IN_initialize;
          windEmulatorStep4__cwInitialize();
          windEmulatorStep4_WECSim_B.state_ed = abbStateEnum_init;
          windEmulatorStep4_WECSim_DW.is_active_UpdateControlWord_p = 1U;
        } else {
          windEmulatorS_swParseStatusWord();
          windEmulatorStep4_WECSim_DW.cwRESET_e =
            windEmulatorStep4_WECSim_B.CastToDouble_j;
          switch (windEmulatorStep4_WECSim_DW.is_UpdateStateMachine_g) {
           case windEmulatorStep4__IN_DelayOFF1:
            if (windEmulatorStep4_WECSim_DW.temporalCounter_i1_g >= 25) {
              windEmulatorStep4_WECSim_DW.is_UpdateStateMachine_g =
                windEmula_IN_notReadyToSwitchOn;
              windEmulatorStep4_WECSim_DW.cwOFF1_CONTROL_i = 0.0;
            } else {
              windEmulatorStep4_WECSim_B.state_ed = abbStateEnum_delayOff1;
            }
            break;

           case windEmulatorStep4_IN_initialize:
            windEmulatorStep4_WECSim_DW.is_UpdateStateMachine_g =
              windEmula_IN_notReadyToSwitchOn;
            windEmulatorStep4_WECSim_DW.cwOFF1_CONTROL_i = 0.0;
            break;

           case windEmula_IN_notReadyToSwitchOn:
            f = ((windEmulatorStep4_WECSim_DW.swRDY_ON_f == 1.0) &&
                 (windEmulatorStep4_WECSim_DW.swWARNING_k == 0.0));
            if (f) {
              windEmulatorStep4_WECSim_DW.is_UpdateStateMachine_g =
                windEmulator_IN_readyToSwitchOn;
              windEmulatorStep4_WECSim_DW.cwOFF1_CONTROL_i = 0.0;
            } else {
              windEmulatorStep4_WECSim_B.state_ed =
                abbStateEnum_notReadyToSwitchOn;
            }
            break;

           case windEmulat_IN_operationDisabled:
            f = (windEmulatorStep4_WECSim_B.NotEqual_f &&
                 (windEmulatorStep4_WECSim_DW.swRDY_RUN_f == 1.0));
            if (f) {
              windEmulatorStep4_WECSim_DW.is_UpdateStateMachine_g =
                windEmulato_IN_operationEnabled;
              windEmulatorStep4_WECSim_DW.cwOFF1_CONTROL_i = 1.0;
              windEmulatorStep4_WECSim_DW.cwENABLE_OPERATION_l = 1.0;
            } else {
              f = (windEmulatorStep4_WECSim_B.NotEqual1_m ||
                   windEmulatorStep4_WECSim_B.NotEqual_i ||
                   (windEmulatorStep4_WECSim_DW.swTRIPPED_e == 1.0) ||
                   (windEmulatorStep4_WECSim_DW.swSWC_ON_INHIB_l == 1.0) ||
                   (windEmulatorStep4_WECSim_DW.swWARNING_k == 1.0));
              if (f) {
                windEmulatorStep4_WECSim_DW.temporalCounter_i1_g = 0U;
                windEmulatorStep4_WECSim_DW.is_UpdateStateMachine_g =
                  windEmulatorStep4__IN_DelayOFF1;
              } else {
                windEmulatorStep4_WECSim_B.state_ed =
                  abbStateEnum_operationDisabled;
              }
            }
            break;

           case windEmulato_IN_operationEnabled:
            if (windEmulatorStep4_WECSim_B.NotEqual1_m) {
              windEmulatorStep4_WECSim_DW.cwENABLE_OPERATION_l = 0.0;
              windEmulatorStep4_WECSim_DW.is_UpdateStateMachine_g =
                windEmulat_IN_operationDisabled;
              windEmulatorStep4_WECSim_DW.cwOFF1_CONTROL_i = 1.0;
            } else {
              f = (windEmulatorStep4_WECSim_B.NotEqual_i ||
                   (windEmulatorStep4_WECSim_DW.swTRIPPED_e == 1.0) ||
                   (windEmulatorStep4_WECSim_DW.swSWC_ON_INHIB_l == 1.0));
              if (f) {
                windEmulatorStep4_WECSim_DW.cwENABLE_OPERATION_l = 0.0;
                windEmulatorStep4_WECSim_DW.temporalCounter_i1_g = 0U;
                windEmulatorStep4_WECSim_DW.is_UpdateStateMachine_g =
                  windEmulatorStep4__IN_DelayOFF1;
              } else {
                windEmulatorStep4_WECSim_B.state_ed =
                  abbStateEnum_operationEnabled;
              }
            }
            break;

           default:
            /* case IN_readyToSwitchOn: */
            f = (windEmulatorStep4_WECSim_B.NotEqual_f &&
                 (windEmulatorStep4_WECSim_DW.swREMOTE_j == 1.0));
            if (f) {
              windEmulatorStep4_WECSim_DW.is_UpdateStateMachine_g =
                windEmulat_IN_operationDisabled;
              windEmulatorStep4_WECSim_DW.cwOFF1_CONTROL_i = 1.0;
            } else {
              windEmulatorStep4_WECSim_B.state_ed = abbStateEnum_readyToSwitchOn;
            }
            break;
          }

          windEmulator_cwBuildControlWord();
        }

        /* End of Chart: '<S16>/ABB Fieldbus Control' */

        /* DataTypeConversion: '<S1>/Cast To Double1' */
        windEmulatorStep4_WECSim_B.CastToDouble1_ia =
          windEmulatorStep4_WECSim_B.state_ed;
      }

      /* MultiPortSwitch: '<S3>/Multiport Switch1' incorporates:
       *  DataTypeConversion: '<S3>/toExpTypeEnum'
       */
      switch (windEmulatorStep4_WECSim_B.toExpTypeEnum) {
       case expTypeEnum_off:
        /* MultiPortSwitch: '<S3>/Multiport Switch1' incorporates:
         *  Constant: '<S3>/zeroTorque'
         */
        windEmulatorStep4_WECSim_B.MultiportSwitch1_m =
          windEmulatorStep4_WECSim_cal->zeroTorque_Value;
        break;

       case expTypeEnum_sid:
        /* MultiPortSwitch: '<S3>/Multiport Switch1' */
        windEmulatorStep4_WECSim_B.MultiportSwitch1_m =
          windEmulatorStep4_WECSim_B.BusAssignment_f.acs800Torque_Nm;
        break;

       case expTypeEnum_hil:
        /* MultiPortSwitch: '<S3>/Multiport Switch1' */
        windEmulatorStep4_WECSim_B.MultiportSwitch1_m =
          windEmulatorStep4_WECSim_B.BusAssignment_n.genTorqueCmd_Nm;
        break;

       default:
        /* MultiPortSwitch: '<S3>/Multiport Switch1' incorporates:
         *  Constant: '<S3>/zeroTorque'
         */
        windEmulatorStep4_WECSim_B.MultiportSwitch1_m =
          windEmulatorStep4_WECSim_cal->zeroTorque_Value;
        break;
      }

      /* End of MultiPortSwitch: '<S3>/Multiport Switch1' */

      /* Gain: '<S1>/Nm -> %' */
      riseValLimit = 1.0 / tmp_n * 100.0;

      /* Gain: '<S1>/Nm -> %' */
      windEmulatorStep4_WECSim_B.Nm_j = riseValLimit *
        windEmulatorStep4_WECSim_B.MultiportSwitch1_m;

      /* BusAssignment: '<S1>/Bus Assignment' */
      windEmulatorStep4_WECSim_B.BusAssignment_kc.ctrlWord =
        windEmulatorStep4_WECSim_B.ControlWord_l;
      windEmulatorStep4_WECSim_B.BusAssignment_kc.state =
        windEmulatorStep4_WECSim_B.CastToDouble1_ia;
      windEmulatorStep4_WECSim_B.BusAssignment_kc.torqueSetpoint_Nm =
        windEmulatorStep4_WECSim_B.MultiportSwitch1_m;
      windEmulatorStep4_WECSim_B.BusAssignment_kc.torqueSetpoint_percent =
        windEmulatorStep4_WECSim_B.Nm_j;
      if (tmp_g) {
        /* ToAsyncQueueBlock generated from: '<S21>/acs800CtrlSignals' */
        slrtLogSignal
          (windEmulatorStep4_WECSim_DW.TAQSigLogging_InsertedFor_acs_l.SLRTSigHandles,
           (((windEmulatorStep4_WECSim_M->Timing.clockTick1+
              windEmulatorStep4_WECSim_M->Timing.clockTickH1* 4294967296.0)) *
            0.004));
      }

      /* MATLAB Function: '<S41>/parseCtrlWord' */
      windEmulatorStep4_parseCtrlWord
        (windEmulatorStep4_WECSim_B.BusAssignment_kc.ctrlWord,
         &windEmulatorStep4_WECSim_B.sf_parseCtrlWord,
         &windEmulatorStep4_WECSim_DW.sf_parseCtrlWord);

      /* Logic: '<S41>/enableOperation' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_WECSim_B.enableOperation_c =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord.enableOperation != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_ld);

      /* Logic: '<S41>/extCtrlLoc' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_WECSim_B.extCtrlLoc_ec =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord.extCtrlLoc != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_ld);

      /* Logic: '<S41>/inching1' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_WECSim_B.inching1_o =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord.inching1 != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_ld);

      /* Logic: '<S41>/inching2' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_WECSim_B.inching2_j =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord.inching2 != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_ld);

      /* Logic: '<S41>/off1Ctrl' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_WECSim_B.off1Ctrl_a =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord.off1Ctrl != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_ld);

      /* Logic: '<S41>/off2Ctrl' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_WECSim_B.off2Ctrl_h =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord.off2Ctrl != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_ld);

      /* Logic: '<S41>/off3Ctrl' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_WECSim_B.off3Ctrl_a =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord.off3Ctrl != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_ld);

      /* Logic: '<S41>/rampHold' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_WECSim_B.rampHold_p =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord.rampHold != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_ld);

      /* Logic: '<S41>/rampInZero' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_WECSim_B.rampInZero_e =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord.rampInZero != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_ld);

      /* Logic: '<S41>/rampOutZero' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_WECSim_B.rampOutZero_p =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord.rampOutZero != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_ld);

      /* Logic: '<S41>/remoteCmd' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_WECSim_B.remoteCmd_p =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord.remoteCmd != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_ld);

      /* Logic: '<S41>/reset' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_WECSim_B.reset_b =
        ((windEmulatorStep4_WECSim_B.sf_parseCtrlWord.reset != 0) &&
         windEmulatorStep4_WECSim_cal->Constant_Value_ld);

      /* Bias: '<S22>/ctrlWord' */
      windEmulatorStep4_WECSim_B.ctrlWord_g = static_cast<uint16_T>
        (windEmulatorStep4_WECSim_B.BusAssignment_kc.ctrlWord +
         windEmulatorStep4_WECSim_cal->ctrlWord_Bias_a);

      /* Bias: '<S22>/state' */
      windEmulatorStep4_WECSim_B.state_i =
        windEmulatorStep4_WECSim_B.BusAssignment_kc.state +
        windEmulatorStep4_WECSim_cal->state_Bias_l;

      /* Gain: '<S22>/torqueSetpoint_Nm' */
      windEmulatorStep4_WECSim_B.torqueSetpoint_Nm_j =
        windEmulatorStep4_WECSim_cal->torqueSetpoint_Nm_Gain_m *
        windEmulatorStep4_WECSim_B.BusAssignment_kc.torqueSetpoint_Nm;

      /* Gain: '<S22>/torqueSetpoint_percent' */
      windEmulatorStep4_WECSim_B.torqueSetpoint_percent_c =
        windEmulatorStep4_WECSim_cal->torqueSetpoint_percent_Gain_k *
        windEmulatorStep4_WECSim_B.BusAssignment_kc.torqueSetpoint_percent;
      if (tmp_g) {
        /* ToAsyncQueueBlock generated from: '<S34>/hptoSignals' */
        slrtLogSignal
          (windEmulatorStep4_WECSim_DW.TAQSigLogging_InsertedFor_hptoS.SLRTSigHandles,
           (((windEmulatorStep4_WECSim_M->Timing.clockTick1+
              windEmulatorStep4_WECSim_M->Timing.clockTickH1* 4294967296.0)) *
            0.004));
      }

      /* Gain: '<S51>/bar->Pa' */
      windEmulatorStep4_WECSim_B.barPa = *get_bar2pa() *
        windEmulatorStep4_WECSim_B.BusAssignment_n.pressure_bar;

      /* Gain: '<S51>/l//m -> m3//s' */
      windEmulatorStep4_WECSim_B.lmm3s = *get_lpm2m3ps() *
        windEmulatorStep4_WECSim_B.BusAssignment_n.genPumpFlow_lpm;

      /* Product: '<S51>/Product' */
      windEmulatorStep4_WECSim_B.Product_gm = windEmulatorStep4_WECSim_B.barPa *
        windEmulatorStep4_WECSim_B.lmm3s;
      windEmulatorSte_MovingAverage_p(windEmulatorStep4_WECSim_B.Product_gm,
        &windEmulatorStep4_WECSim_B.MovingAverage_pn,
        &windEmulatorStep4_WECSim_DW.MovingAverage_pn);

      /* Gain: '<S51>/rpm -> rad//s' */
      windEmulatorStep4_WECSim_B.rpmrads_oh = tmp_p *
        windEmulatorStep4_WECSim_B.BusAssignment_n.excShaftSpeed_rpm;

      /* Product: '<S51>/Product1' */
      windEmulatorStep4_WECSim_B.Product1_i =
        windEmulatorStep4_WECSim_B.rpmrads_oh *
        windEmulatorStep4_WECSim_B.BusAssignment_n.excShaftTorque_Nm;
      windEmulatorSte_MovingAverage_p(windEmulatorStep4_WECSim_B.Product1_i,
        &windEmulatorStep4_WECSim_B.MovingAverage1,
        &windEmulatorStep4_WECSim_DW.MovingAverage1);

      /* Gain: '<S51>/excShaftPowerAverage_W' */
      windEmulatorStep4_WECSim_B.excShaftPowerAverage_W =
        windEmulatorStep4_WECSim_cal->excShaftPowerAverage_W_Gain *
        windEmulatorStep4_WECSim_B.MovingAverage1.MovingAverage;

      /* Gain: '<S51>/excShaftPower_W' */
      windEmulatorStep4_WECSim_B.excShaftPower_W =
        windEmulatorStep4_WECSim_cal->excShaftPower_W_Gain *
        windEmulatorStep4_WECSim_B.Product1_i;

      /* Gain: '<S51>/hydrPowerAverage_W' */
      windEmulatorStep4_WECSim_B.hydrPowerAverage_W =
        windEmulatorStep4_WECSim_cal->hydrPowerAverage_W_Gain *
        windEmulatorStep4_WECSim_B.MovingAverage_pn.MovingAverage;

      /* Gain: '<S51>/hydrPower_W' */
      windEmulatorStep4_WECSim_B.hydrPower_W =
        windEmulatorStep4_WECSim_cal->hydrPower_W_Gain *
        windEmulatorStep4_WECSim_B.Product_gm;

      /* Gain: '<S35>/ctrlSignal1' */
      windEmulatorStep4_WECSim_B.ctrlSignal1 =
        windEmulatorStep4_WECSim_cal->ctrlSignal1_Gain *
        windEmulatorStep4_WECSim_B.BusAssignment_n.ctrlSignal1;

      /* Gain: '<S35>/ctrlSignal2' */
      windEmulatorStep4_WECSim_B.ctrlSignal2 =
        windEmulatorStep4_WECSim_cal->ctrlSignal2_Gain *
        windEmulatorStep4_WECSim_B.BusAssignment_n.ctrlSignal2;

      /* Gain: '<S35>/excShaftSpeed_rpm' */
      windEmulatorStep4_WECSim_B.excShaftSpeed_rpm =
        windEmulatorStep4_WECSim_cal->excShaftSpeed_rpm_Gain *
        windEmulatorStep4_WECSim_B.BusAssignment_n.excShaftSpeed_rpm;

      /* Gain: '<S35>/excShaftTorque_Nm' */
      windEmulatorStep4_WECSim_B.excShaftTorque_Nm =
        windEmulatorStep4_WECSim_cal->excShaftTorque_Nm_Gain *
        windEmulatorStep4_WECSim_B.BusAssignment_n.excShaftTorque_Nm;

      /* Gain: '<S35>/genPumpFlow_lpm' */
      windEmulatorStep4_WECSim_B.genPumpFlow_lpm =
        windEmulatorStep4_WECSim_cal->genPumpFlow_lpm_Gain *
        windEmulatorStep4_WECSim_B.BusAssignment_n.genPumpFlow_lpm;

      /* Gain: '<S35>/genShaftSpeed_rpm' */
      windEmulatorStep4_WECSim_B.genShaftSpeed_rpm =
        windEmulatorStep4_WECSim_cal->genShaftSpeed_rpm_Gain *
        windEmulatorStep4_WECSim_B.BusAssignment_n.genShaftSpeed_rpm;

      /* Gain: '<S35>/genTorqueCmd_Nm' */
      windEmulatorStep4_WECSim_B.genTorqueCmd_Nm =
        windEmulatorStep4_WECSim_cal->genTorqueCmd_Nm_Gain *
        windEmulatorStep4_WECSim_B.BusAssignment_n.genTorqueCmd_Nm;

      /* Gain: '<S35>/hmOutputShaftTorque_Nm' */
      windEmulatorStep4_WECSim_B.hmOutputShaftTorque_Nm =
        windEmulatorStep4_WECSim_cal->hmOutputShaftTorque_Nm_Gain *
        windEmulatorStep4_WECSim_B.BusAssignment_n.hmOutputShafTorque_Nm;

      /* Gain: '<S35>/pressure_bar' */
      windEmulatorStep4_WECSim_B.pressure_bar =
        windEmulatorStep4_WECSim_cal->pressure_bar_Gain *
        windEmulatorStep4_WECSim_B.BusAssignment_n.pressure_bar;
      if (tmp_g) {
        /* ToAsyncQueueBlock generated from: '<S32>/hptoCtrl' */
        slrtLogSignal
          (windEmulatorStep4_WECSim_DW.TAQSigLogging_InsertedFor_hptoC.SLRTSigHandles,
           (((windEmulatorStep4_WECSim_M->Timing.clockTick1+
              windEmulatorStep4_WECSim_M->Timing.clockTickH1* 4294967296.0)) *
            0.004));

        /* Gain: '<S33>/excForce_N' */
        windEmulatorStep4_WECSim_B.excForce_N =
          windEmulatorStep4_WECSim_cal->excForce_N_Gain *
          windEmulatorStep4_WECSim_B.BusAssignment_c.excForce_N;

        /* Gain: '<S33>/genSpeedActual' */
        windEmulatorStep4_WECSim_B.genSpeedActual =
          windEmulatorStep4_WECSim_cal->genSpeedActual_Gain *
          windEmulatorStep4_WECSim_B.BusAssignment_c.genSpeedActual;

        /* Gain: '<S33>/speedCtrlReset' */
        windEmulatorStep4_WECSim_B.speedCtrlReset = static_cast<uint8_T>
          (windEmulatorStep4_WECSim_B.BusAssignment_c.speedCtrlReset ?
           static_cast<int32_T>
           (windEmulatorStep4_WECSim_cal->speedCtrlReset_Gain) : 0);

        /* Gain: '<S33>/speedRef_rpm' */
        windEmulatorStep4_WECSim_B.speedRef_rpm =
          windEmulatorStep4_WECSim_cal->speedRef_rpm_Gain *
          windEmulatorStep4_WECSim_B.BusAssignment_c.speedRef_rpm;

        /* ToAsyncQueueBlock generated from: '<S30>/expCtrlSignals' */
        slrtLogSignal
          (windEmulatorStep4_WECSim_DW.TAQSigLogging_InsertedFor_expCt.SLRTSigHandles,
           (((windEmulatorStep4_WECSim_M->Timing.clockTick1+
              windEmulatorStep4_WECSim_M->Timing.clockTickH1* 4294967296.0)) *
            0.004));

        /* Bias: '<S31>/expType' */
        windEmulatorStep4_WECSim_B.expType = static_cast<uint16_T>
          (windEmulatorStep4_WECSim_B.BusAssignment_b.expType +
           windEmulatorStep4_WECSim_cal->expType_Bias);

        /* Gain: '<S31>/ramp' */
        windEmulatorStep4_WECSim_B.ramp =
          windEmulatorStep4_WECSim_cal->ramp_Gain *
          windEmulatorStep4_WECSim_B.BusAssignment_b.ramp;

        /* Logic: '<S31>/resetHilIntegrator' incorporates:
         *  Constant: '<S31>/Constant'
         */
        windEmulatorStep4_WECSim_B.resetHilIntegrator =
          (windEmulatorStep4_WECSim_cal->Constant_Value_k &&
           windEmulatorStep4_WECSim_B.BusAssignment_b.resetHilIntegrator);

        /* Logic: '<S31>/resetSidIntegrator' incorporates:
         *  Constant: '<S31>/Constant'
         */
        windEmulatorStep4_WECSim_B.resetSidIntegrator =
          (windEmulatorStep4_WECSim_cal->Constant_Value_k &&
           windEmulatorStep4_WECSim_B.BusAssignment_b.resetSidIntegrator);

        /* Bias: '<S31>/runCounter' */
        windEmulatorStep4_WECSim_B.runCounter =
          windEmulatorStep4_WECSim_B.BusAssignment_b.runCounter +
          windEmulatorStep4_WECSim_cal->runCounter_Bias;

        /* Logic: '<S31>/runHil' incorporates:
         *  Constant: '<S31>/Constant'
         */
        windEmulatorStep4_WECSim_B.runHil =
          (windEmulatorStep4_WECSim_cal->Constant_Value_k &&
           windEmulatorStep4_WECSim_B.BusAssignment_b.runHil);

        /* Logic: '<S31>/runSid' incorporates:
         *  Constant: '<S31>/Constant'
         */
        windEmulatorStep4_WECSim_B.runSid =
          (windEmulatorStep4_WECSim_cal->Constant_Value_k &&
           windEmulatorStep4_WECSim_B.BusAssignment_b.runSid);

        /* Bias: '<S31>/stepCounter' */
        windEmulatorStep4_WECSim_B.stepCounter =
          windEmulatorStep4_WECSim_B.BusAssignment_b.stepCounter +
          windEmulatorStep4_WECSim_cal->stepCounter_Bias;

        /* Gain: '<S31>/time' */
        windEmulatorStep4_WECSim_B.time =
          windEmulatorStep4_WECSim_cal->time_Gain *
          windEmulatorStep4_WECSim_B.BusAssignment_b.time;

        /* S-Function (slecatpdorx): '<S12>/readTorqueInput' */
        {
          /*------------ S-Function Block: <S12>/readTorqueInput PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.readTorqueInput;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 856;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 6 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* DataTypeConversion: '<S12>/Cast To Double' */
        windEmulatorStep4_WECSim_B.CastToDouble_a =
          windEmulatorStep4_WECSim_B.readTorqueInput;

        /* Gain: '<S12>/Gain' */
        windEmulatorStep4_WECSim_B.Gain_nm = *get_futekTorqueScale() *
          windEmulatorStep4_WECSim_B.CastToDouble_a;

        /* S-Function (slecatpdorx): '<S12>/readEncoderCounter' */
        {
          /*------------ S-Function Block: <S12>/readEncoderCounter PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.readEncoderCounter;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 952;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 7 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Memory: '<S552>/lastRawCounts' */
        windEmulatorStep4_WECSim_B.lastRawCounts =
          windEmulatorStep4_WECSim_DW.lastRawCounts_PreviousInput;

        /* DataTypeConversion: '<S552>/Cast To Double' */
        windEmulatorStep4_WECSim_B.CastToDouble_p = static_cast<int32_T>
          (windEmulatorStep4_WECSim_B.readEncoderCounter);

        /* Sum: '<S552>/Add' */
        windEmulatorStep4_WECSim_B.Add_d =
          windEmulatorStep4_WECSim_B.lastRawCounts -
          windEmulatorStep4_WECSim_B.CastToDouble_p;

        /* Abs: '<S552>/Abs' */
        i = windEmulatorStep4_WECSim_B.Add_d;
        if (i < 0) {
          /* Abs: '<S552>/Abs' */
          windEmulatorStep4_WECSim_B.Abs_a = -i;
        } else {
          /* Abs: '<S552>/Abs' */
          windEmulatorStep4_WECSim_B.Abs_a = i;
        }

        /* End of Abs: '<S552>/Abs' */

        /* RelationalOperator: '<S552>/Relational Operator' incorporates:
         *  Constant: '<S552>/Constant'
         */
        windEmulatorStep4_WECSim_B.RelationalOperator =
          (windEmulatorStep4_WECSim_B.Abs_a >=
           windEmulatorStep4_WECSim_cal->Constant_Value_pf);

        /* Switch: '<S552>/Switch' */
        if (windEmulatorStep4_WECSim_B.RelationalOperator) {
          /* Signum: '<S552>/Sign' */
          i = windEmulatorStep4_WECSim_B.Add_d;
          if (i < 0) {
            /* Signum: '<S552>/Sign' */
            windEmulatorStep4_WECSim_B.Sign = -1;
          } else {
            /* Signum: '<S552>/Sign' */
            windEmulatorStep4_WECSim_B.Sign = (i > 0);
          }

          /* End of Signum: '<S552>/Sign' */

          /* Switch: '<S552>/Switch' */
          windEmulatorStep4_WECSim_B.Switch_g0 = windEmulatorStep4_WECSim_B.Sign;
        } else {
          /* Switch: '<S552>/Switch' incorporates:
           *  Constant: '<S552>/Constant1'
           */
          windEmulatorStep4_WECSim_B.Switch_g0 =
            windEmulatorStep4_WECSim_cal->Constant1_Value_jt;
        }

        /* End of Switch: '<S552>/Switch' */

        /* Memory: '<S552>/lastTurn' */
        windEmulatorStep4_WECSim_B.lastTurn =
          windEmulatorStep4_WECSim_DW.lastTurn_PreviousInput;

        /* Sum: '<S552>/Add1' */
        windEmulatorStep4_WECSim_B.Add1_p2 =
          windEmulatorStep4_WECSim_B.Switch_g0 +
          windEmulatorStep4_WECSim_B.lastTurn;

        /* DataTypeConversion: '<S552>/Cast To Double3' */
        windEmulatorStep4_WECSim_B.CastToDouble3_b =
          windEmulatorStep4_WECSim_B.CastToDouble_p;

        /* Gain: '<S552>/encoderCountsToRad' */
        windEmulatorStep4_WECSim_B.encoderCountsToRad =
          *get_absEncoderCountsToRad() *
          windEmulatorStep4_WECSim_B.CastToDouble3_b;

        /* DataTypeConversion: '<S552>/Cast To Double1' */
        windEmulatorStep4_WECSim_B.CastToDouble1_m =
          windEmulatorStep4_WECSim_B.Add1_p2;

        /* Gain: '<S552>/Gain' */
        windEmulatorStep4_WECSim_B.Gain_b =
          windEmulatorStep4_WECSim_cal->Gain_Gain_cf *
          windEmulatorStep4_WECSim_B.CastToDouble1_m;

        /* Sum: '<S552>/Add2' */
        windEmulatorStep4_WECSim_B.Add2 =
          windEmulatorStep4_WECSim_B.encoderCountsToRad +
          windEmulatorStep4_WECSim_B.Gain_b;

        /* SampleTimeMath: '<S553>/TSamp'
         *
         * About '<S553>/TSamp':
         *  y = u * K where K = 1 / ( w * Ts )
         *   */
        windEmulatorStep4_WECSim_B.TSamp = windEmulatorStep4_WECSim_B.Add2 *
          windEmulatorStep4_WECSim_cal->TSamp_WtEt;

        /* UnitDelay: '<S553>/UD' */
        windEmulatorStep4_WECSim_B.Uk1 = windEmulatorStep4_WECSim_DW.UD_DSTATE;

        /* Sum: '<S553>/Diff' */
        windEmulatorStep4_WECSim_B.Diff = windEmulatorStep4_WECSim_B.TSamp -
          windEmulatorStep4_WECSim_B.Uk1;

        /* Gain: '<S12>/rad//s->rpm' */
        windEmulatorStep4_WECSim_B.radsrpm = tmp_10 *
          windEmulatorStep4_WECSim_B.Diff;

        /* S-Function (slecatpdorx): '<S12>/EtherCAT PDO Receive7' */
        {
          /*------------ S-Function Block: <S12>/EtherCAT PDO Receive7 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            windEmulatorStep4_WECSim_B.encoderStatus;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 936;
          for (sigIdx=0; sigIdx < 2; sigIdx++) {
            switch ( 3 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (8 == 8) && (1 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 8, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (8 == 8) && (1 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 8, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 8, sigOutputPtr+sigIdx*
                                 1, 1);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (8 == 16) && (1 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 8, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (8 == 16) && (1 == sizeof(uint16_T))) {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 8, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (8 == 32) && (1 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 8, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (8 == 32) && (1 == sizeof(uint32_T))) {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 8, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 8;
          }
        }

        /* DataTypeConversion: '<S12>/Cast To Double1' */
        windEmulatorStep4_WECSim_B.CastToDouble1_f =
          windEmulatorStep4_WECSim_B.encoderStatus[0];

        /* DataTypeConversion: '<S12>/Cast To Double2' */
        windEmulatorStep4_WECSim_B.CastToDouble2_d =
          windEmulatorStep4_WECSim_B.encoderStatus[1];

        /* BusAssignment: '<S12>/Bus Assignment' */
        windEmulatorStep4_WECSim_B.BusAssignment_g.torqueActual_Nm =
          windEmulatorStep4_WECSim_B.Gain_nm;
        windEmulatorStep4_WECSim_B.BusAssignment_g.absEncoderCounts =
          windEmulatorStep4_WECSim_B.readEncoderCounter;
        windEmulatorStep4_WECSim_B.BusAssignment_g.absEncoderTurns =
          windEmulatorStep4_WECSim_B.Add1_p2;
        windEmulatorStep4_WECSim_B.BusAssignment_g.absEncoderPosition_rad =
          windEmulatorStep4_WECSim_B.Add2;
        windEmulatorStep4_WECSim_B.BusAssignment_g.absEncoderSpeed_rpm =
          windEmulatorStep4_WECSim_B.radsrpm;
        windEmulatorStep4_WECSim_B.BusAssignment_g.absEncoderStatus1 =
          windEmulatorStep4_WECSim_B.CastToDouble1_f;
        windEmulatorStep4_WECSim_B.BusAssignment_g.absEncoderStatus2 =
          windEmulatorStep4_WECSim_B.CastToDouble2_d;

        /* ToAsyncQueueBlock generated from: '<S38>/shaftSignals' */
        slrtLogSignal
          (windEmulatorStep4_WECSim_DW.TAQSigLogging_InsertedFor_shaft.SLRTSigHandles,
           (((windEmulatorStep4_WECSim_M->Timing.clockTick1+
              windEmulatorStep4_WECSim_M->Timing.clockTickH1* 4294967296.0)) *
            0.004));

        /* Bias: '<S39>/absEncoderCounts' */
        windEmulatorStep4_WECSim_B.absEncoderCounts =
          windEmulatorStep4_WECSim_B.BusAssignment_g.absEncoderCounts +
          windEmulatorStep4_WECSim_cal->absEncoderCounts_Bias;

        /* Gain: '<S39>/absEncoderPosition_rad' */
        windEmulatorStep4_WECSim_B.absEncoderPosition_rad =
          windEmulatorStep4_WECSim_cal->absEncoderPosition_rad_Gain *
          windEmulatorStep4_WECSim_B.BusAssignment_g.absEncoderPosition_rad;

        /* Gain: '<S39>/absEncoderSpeed_rpm' */
        windEmulatorStep4_WECSim_B.absEncoderSpeed_rpm =
          windEmulatorStep4_WECSim_cal->absEncoderSpeed_rpm_Gain *
          windEmulatorStep4_WECSim_B.BusAssignment_g.absEncoderSpeed_rpm;

        /* Bias: '<S39>/absEncoderStatus1' */
        windEmulatorStep4_WECSim_B.absEncoderStatus1 = static_cast<uint8_T>
          (windEmulatorStep4_WECSim_B.BusAssignment_g.absEncoderStatus1 +
           windEmulatorStep4_WECSim_cal->absEncoderStatus1_Bias);

        /* Bias: '<S39>/absEncoderStatus2' */
        windEmulatorStep4_WECSim_B.absEncoderStatus2 = static_cast<uint8_T>
          (windEmulatorStep4_WECSim_B.BusAssignment_g.absEncoderStatus2 +
           windEmulatorStep4_WECSim_cal->absEncoderStatus2_Bias);

        /* Bias: '<S39>/absEncoderTurns' */
        windEmulatorStep4_WECSim_B.absEncoderTurns =
          windEmulatorStep4_WECSim_B.BusAssignment_g.absEncoderTurns +
          windEmulatorStep4_WECSim_cal->absEncoderTurns_Bias;

        /* Gain: '<S39>/torqueActual_Nm' */
        windEmulatorStep4_WECSim_B.torqueActual_Nm =
          windEmulatorStep4_WECSim_cal->torqueActual_Nm_Gain *
          windEmulatorStep4_WECSim_B.BusAssignment_g.torqueActual_Nm;

        /* S-Function (slecatpdorx): '<S9>/L1InaccurateURead' */
        {
          /*------------ S-Function Block: <S9>/L1InaccurateURead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L1InaccurateURead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 2699;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 8 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (1 == 8) && (1 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (1 == 8) && (1 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+sigIdx*
                                 1, 1);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(uint16_T))) {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(uint32_T))) {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 1;
          }
        }

        /* Logic: '<S9>/L1InaccurateU' incorporates:
         *  Constant: '<S9>/Constant1'
         */
        windEmulatorStep4_WECSim_B.L1InaccurateU =
          (windEmulatorStep4_WECSim_cal->Constant1_Value_k &&
           windEmulatorStep4_WECSim_B.L1InaccurateURead);

        /* S-Function (slecatpdorx): '<S9>/L1InaccurateIRead' */
        {
          /*------------ S-Function Block: <S9>/L1InaccurateIRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L1InaccurateIRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 2700;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 8 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (1 == 8) && (1 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (1 == 8) && (1 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+sigIdx*
                                 1, 1);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(uint16_T))) {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(uint32_T))) {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 1;
          }
        }

        /* Logic: '<S9>/L1InaccurateI' incorporates:
         *  Constant: '<S9>/Constant1'
         */
        windEmulatorStep4_WECSim_B.L1InaccurateI =
          (windEmulatorStep4_WECSim_cal->Constant1_Value_k &&
           windEmulatorStep4_WECSim_B.L1InaccurateIRead);

        /* S-Function (slecatpdorx): '<S9>/L1VoltageRead' */
        {
          /*------------ S-Function Block: <S9>/L1VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L1VoltageRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 2712;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/L1Voltage' */
        windEmulatorStep4_WECSim_B.L1Voltage = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L1Voltage_Gain) *
          windEmulatorStep4_WECSim_B.L1VoltageRead;

        /* S-Function (slecatpdorx): '<S9>/L1CurrentRead' */
        {
          /*------------ S-Function Block: <S9>/L1CurrentRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L1CurrentRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 2744;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/L1Current' */
        windEmulatorStep4_WECSim_B.L1Current = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L1Current_Gain) *
          windEmulatorStep4_WECSim_B.L1CurrentRead;

        /* S-Function (slecatpdorx): '<S9>/L1PowFactorRead' */
        {
          /*------------ S-Function Block: <S9>/L1PowFactorRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L1PowFactorRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 2872;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/L1PowFactor' */
        windEmulatorStep4_WECSim_B.L1PowFactor = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L1PowFactor_Gain) *
          windEmulatorStep4_WECSim_B.L1PowFactorRead;

        /* S-Function (slecatpdorx): '<S9>/L1ActivePowRead' */
        {
          /*------------ S-Function Block: <S9>/L1ActivePowRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L1ActivePowRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 2776;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/L1ActivePow' */
        windEmulatorStep4_WECSim_B.L1ActivePow = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L1ActivePow_Gain) *
          windEmulatorStep4_WECSim_B.L1ActivePowRead;

        /* S-Function (slecatpdorx): '<S9>/L1THDuRead' */
        {
          /*------------ S-Function Block: <S9>/L1THDuRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L1THDuRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 2984;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/L1THDu' */
        windEmulatorStep4_WECSim_B.L1THDu = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L1THDu_Gain) *
          windEmulatorStep4_WECSim_B.L1THDuRead;

        /* S-Function (slecatpdorx): '<S9>/L1THDiRead' */
        {
          /*------------ S-Function Block: <S9>/L1THDiRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L1THDiRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 3048;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/L1THDi' */
        windEmulatorStep4_WECSim_B.L1THDi = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L1THDi_Gain) *
          windEmulatorStep4_WECSim_B.L1THDiRead;

        /* S-Function (slecatpdorx): '<S9>/L2InaccurateURead' */
        {
          /*------------ S-Function Block: <S9>/L2InaccurateURead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L2InaccurateURead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 3115;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 8 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (1 == 8) && (1 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (1 == 8) && (1 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+sigIdx*
                                 1, 1);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(uint16_T))) {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(uint32_T))) {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 1;
          }
        }

        /* Logic: '<S9>/L2InaccurateU' incorporates:
         *  Constant: '<S9>/Constant2'
         */
        windEmulatorStep4_WECSim_B.L2InaccurateU =
          (windEmulatorStep4_WECSim_cal->Constant2_Value_e &&
           windEmulatorStep4_WECSim_B.L2InaccurateURead);

        /* S-Function (slecatpdorx): '<S9>/L2InaccurateIRead' */
        {
          /*------------ S-Function Block: <S9>/L2InaccurateIRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L2InaccurateIRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 3116;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 8 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (1 == 8) && (1 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (1 == 8) && (1 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+sigIdx*
                                 1, 1);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(uint16_T))) {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(uint32_T))) {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 1;
          }
        }

        /* Logic: '<S9>/L2InaccurateI' incorporates:
         *  Constant: '<S9>/Constant2'
         */
        windEmulatorStep4_WECSim_B.L2InaccurateI =
          (windEmulatorStep4_WECSim_cal->Constant2_Value_e &&
           windEmulatorStep4_WECSim_B.L2InaccurateIRead);

        /* S-Function (slecatpdorx): '<S9>/L2VoltageRead' */
        {
          /*------------ S-Function Block: <S9>/L2VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L2VoltageRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 3128;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/L2Voltage' */
        windEmulatorStep4_WECSim_B.L2Voltage = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L2Voltage_Gain) *
          windEmulatorStep4_WECSim_B.L2VoltageRead;

        /* S-Function (slecatpdorx): '<S9>/L2CurrentRead' */
        {
          /*------------ S-Function Block: <S9>/L2CurrentRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L2CurrentRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 3160;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/L2Current' */
        windEmulatorStep4_WECSim_B.L2Current = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L2Current_Gain) *
          windEmulatorStep4_WECSim_B.L2CurrentRead;

        /* S-Function (slecatpdorx): '<S9>/L2PowFactorRead' */
        {
          /*------------ S-Function Block: <S9>/L2PowFactorRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L2PowFactorRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 3288;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/L2PowFactor' */
        windEmulatorStep4_WECSim_B.L2PowFactor = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L2PowFactor_Gain) *
          windEmulatorStep4_WECSim_B.L2PowFactorRead;

        /* S-Function (slecatpdorx): '<S9>/L2ActivePowRead' */
        {
          /*------------ S-Function Block: <S9>/L2ActivePowRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L2ActivePowRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 3192;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/L2ActivePow' */
        windEmulatorStep4_WECSim_B.L2ActivePow = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L2ActivePow_Gain) *
          windEmulatorStep4_WECSim_B.L2ActivePowRead;

        /* S-Function (slecatpdorx): '<S9>/L2THDuRead' */
        {
          /*------------ S-Function Block: <S9>/L2THDuRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L2THDuRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 3400;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/L2THDu' */
        windEmulatorStep4_WECSim_B.L2THDu = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L2THDu_Gain) *
          windEmulatorStep4_WECSim_B.L2THDuRead;

        /* S-Function (slecatpdorx): '<S9>/L2THDiRead' */
        {
          /*------------ S-Function Block: <S9>/L2THDiRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L2THDiRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 3464;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/L2THDi' */
        windEmulatorStep4_WECSim_B.L2THDi = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L2THDi_Gain) *
          windEmulatorStep4_WECSim_B.L2THDiRead;

        /* S-Function (slecatpdorx): '<S9>/L3InaccurateURead' */
        {
          /*------------ S-Function Block: <S9>/L3InaccurateURead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L3InaccurateURead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 3531;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 8 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (1 == 8) && (1 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (1 == 8) && (1 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+sigIdx*
                                 1, 1);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(uint16_T))) {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(uint32_T))) {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 1;
          }
        }

        /* Logic: '<S9>/L3InaccurateU' incorporates:
         *  Constant: '<S9>/Constant3'
         */
        windEmulatorStep4_WECSim_B.L3InaccurateU =
          (windEmulatorStep4_WECSim_cal->Constant3_Value_dm &&
           windEmulatorStep4_WECSim_B.L3InaccurateURead);

        /* S-Function (slecatpdorx): '<S9>/L3InaccurateIRead' */
        {
          /*------------ S-Function Block: <S9>/L3InaccurateIRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L3InaccurateIRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 3532;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 8 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (1 == 8) && (1 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (1 == 8) && (1 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+sigIdx*
                                 1, 1);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(uint16_T))) {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(uint32_T))) {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 1;
          }
        }

        /* Logic: '<S9>/L3InaccurateI' incorporates:
         *  Constant: '<S9>/Constant3'
         */
        windEmulatorStep4_WECSim_B.L3InaccurateI =
          (windEmulatorStep4_WECSim_cal->Constant3_Value_dm &&
           windEmulatorStep4_WECSim_B.L3InaccurateIRead);

        /* S-Function (slecatpdorx): '<S9>/L3VoltageRead' */
        {
          /*------------ S-Function Block: <S9>/L3VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L3VoltageRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 3544;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/L3Voltage' */
        windEmulatorStep4_WECSim_B.L3Voltage = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L3Voltage_Gain) *
          windEmulatorStep4_WECSim_B.L3VoltageRead;

        /* S-Function (slecatpdorx): '<S9>/L3CurrentRead' */
        {
          /*------------ S-Function Block: <S9>/L3CurrentRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L3CurrentRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 3576;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/L3Current' */
        windEmulatorStep4_WECSim_B.L3Current = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L3Current_Gain) *
          windEmulatorStep4_WECSim_B.L3CurrentRead;

        /* S-Function (slecatpdorx): '<S9>/L3PowFactorRead' */
        {
          /*------------ S-Function Block: <S9>/L3PowFactorRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L3PowFactorRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 3704;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/L3PowFactor' */
        windEmulatorStep4_WECSim_B.L3PowFactor = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L3PowFactor_Gain) *
          windEmulatorStep4_WECSim_B.L3PowFactorRead;

        /* S-Function (slecatpdorx): '<S9>/L3ActivePowRead' */
        {
          /*------------ S-Function Block: <S9>/L3ActivePowRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L3ActivePowRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 3608;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/L3ActivePow' */
        windEmulatorStep4_WECSim_B.L3ActivePow = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L3ActivePow_Gain) *
          windEmulatorStep4_WECSim_B.L3ActivePowRead;

        /* S-Function (slecatpdorx): '<S9>/L3THDuRead' */
        {
          /*------------ S-Function Block: <S9>/L3THDuRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L3THDuRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 3816;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/L3THDu' */
        windEmulatorStep4_WECSim_B.L3THDu = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L3THDu_Gain) *
          windEmulatorStep4_WECSim_B.L3THDuRead;

        /* S-Function (slecatpdorx): '<S9>/L3THDiRead' */
        {
          /*------------ S-Function Block: <S9>/L3THDiRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L3THDiRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 3880;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/L3THDi' */
        windEmulatorStep4_WECSim_B.L3THDi = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L3THDi_Gain) *
          windEmulatorStep4_WECSim_B.L3THDiRead;

        /* S-Function (slecatpdorx): '<S9>/FrequencyRead' */
        {
          /*------------ S-Function Block: <S9>/FrequencyRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.FrequencyRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 3992;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/totalFrequency' */
        windEmulatorStep4_WECSim_B.totalFrequency = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->totalFrequency_Gain) *
          windEmulatorStep4_WECSim_B.FrequencyRead;

        /* S-Function (slecatpdorx): '<S9>/totalPowFactorRead' */
        {
          /*------------ S-Function Block: <S9>/totalPowFactorRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.totalPowFactorRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 4024;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/totalPowFactor' */
        windEmulatorStep4_WECSim_B.totalPowFactor = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->totalPowFactor_Gain) *
          windEmulatorStep4_WECSim_B.totalPowFactorRead;

        /* S-Function (slecatpdorx): '<S9>/totalActivePowRead' */
        {
          /*------------ S-Function Block: <S9>/totalActivePowRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.totalActivePowRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 4088;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/totalActivePow' */
        windEmulatorStep4_WECSim_B.totalActivePow = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->totalActivePow_Gain) *
          windEmulatorStep4_WECSim_B.totalActivePowRead;

        /* S-Function (slecatpdorx): '<S9>/L1L2VoltageRead' */
        {
          /*------------ S-Function Block: <S9>/L1L2VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L1L2VoltageRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 4312;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/L1L2Voltage' */
        windEmulatorStep4_WECSim_B.L1L2Voltage = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L1L2Voltage_Gain) *
          windEmulatorStep4_WECSim_B.L1L2VoltageRead;

        /* S-Function (slecatpdorx): '<S9>/L2L3VoltageRead' */
        {
          /*------------ S-Function Block: <S9>/L2L3VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L2L3VoltageRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 4344;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/L2L3Voltage' */
        windEmulatorStep4_WECSim_B.L2L3Voltage = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L2L3Voltage_Gain) *
          windEmulatorStep4_WECSim_B.L2L3VoltageRead;

        /* S-Function (slecatpdorx): '<S9>/L3L1VoltageRead' */
        {
          /*------------ S-Function Block: <S9>/L3L1VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L3L1VoltageRead;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 4376;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S9>/L3L1Voltage' */
        windEmulatorStep4_WECSim_B.L3L1Voltage = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L3L1Voltage_Gain) *
          windEmulatorStep4_WECSim_B.L3L1VoltageRead;

        /* BusAssignment: '<S9>/Bus Assignment' */
        windEmulatorStep4_WECSim_B.BusAssignment.L1InaccurateU =
          windEmulatorStep4_WECSim_B.L1InaccurateU;
        windEmulatorStep4_WECSim_B.BusAssignment.L1InaccurateI =
          windEmulatorStep4_WECSim_B.L1InaccurateI;
        windEmulatorStep4_WECSim_B.BusAssignment.L1Voltage =
          windEmulatorStep4_WECSim_B.L1Voltage;
        windEmulatorStep4_WECSim_B.BusAssignment.L1Current =
          windEmulatorStep4_WECSim_B.L1Current;
        windEmulatorStep4_WECSim_B.BusAssignment.L1PowFactor =
          windEmulatorStep4_WECSim_B.L1PowFactor;
        windEmulatorStep4_WECSim_B.BusAssignment.L1ActivePow =
          windEmulatorStep4_WECSim_B.L1ActivePow;
        windEmulatorStep4_WECSim_B.BusAssignment.L1THDu =
          windEmulatorStep4_WECSim_B.L1THDu;
        windEmulatorStep4_WECSim_B.BusAssignment.L1THDi =
          windEmulatorStep4_WECSim_B.L1THDi;
        windEmulatorStep4_WECSim_B.BusAssignment.L2InaccurateU =
          windEmulatorStep4_WECSim_B.L2InaccurateU;
        windEmulatorStep4_WECSim_B.BusAssignment.L2InaccurateI =
          windEmulatorStep4_WECSim_B.L2InaccurateI;
        windEmulatorStep4_WECSim_B.BusAssignment.L2Voltage =
          windEmulatorStep4_WECSim_B.L2Voltage;
        windEmulatorStep4_WECSim_B.BusAssignment.L2Current =
          windEmulatorStep4_WECSim_B.L2Current;
        windEmulatorStep4_WECSim_B.BusAssignment.L2PowFactor =
          windEmulatorStep4_WECSim_B.L2PowFactor;
        windEmulatorStep4_WECSim_B.BusAssignment.L2ActivePow =
          windEmulatorStep4_WECSim_B.L2ActivePow;
        windEmulatorStep4_WECSim_B.BusAssignment.L2THDu =
          windEmulatorStep4_WECSim_B.L2THDu;
        windEmulatorStep4_WECSim_B.BusAssignment.L2THDi =
          windEmulatorStep4_WECSim_B.L2THDi;
        windEmulatorStep4_WECSim_B.BusAssignment.L3InaccurateU =
          windEmulatorStep4_WECSim_B.L3InaccurateU;
        windEmulatorStep4_WECSim_B.BusAssignment.L3InaccurateI =
          windEmulatorStep4_WECSim_B.L3InaccurateI;
        windEmulatorStep4_WECSim_B.BusAssignment.L3Voltage =
          windEmulatorStep4_WECSim_B.L3Voltage;
        windEmulatorStep4_WECSim_B.BusAssignment.L3Current =
          windEmulatorStep4_WECSim_B.L3Current;
        windEmulatorStep4_WECSim_B.BusAssignment.L3PowFactor =
          windEmulatorStep4_WECSim_B.L3PowFactor;
        windEmulatorStep4_WECSim_B.BusAssignment.L3ActivePow =
          windEmulatorStep4_WECSim_B.L3ActivePow;
        windEmulatorStep4_WECSim_B.BusAssignment.L3THDu =
          windEmulatorStep4_WECSim_B.L3THDu;
        windEmulatorStep4_WECSim_B.BusAssignment.L3THDi =
          windEmulatorStep4_WECSim_B.L3THDi;
        windEmulatorStep4_WECSim_B.BusAssignment.Frequency =
          windEmulatorStep4_WECSim_B.totalFrequency;
        windEmulatorStep4_WECSim_B.BusAssignment.TotalPowFactor =
          windEmulatorStep4_WECSim_B.totalPowFactor;
        windEmulatorStep4_WECSim_B.BusAssignment.TotalActivePow =
          windEmulatorStep4_WECSim_B.totalActivePow;
        windEmulatorStep4_WECSim_B.BusAssignment.L1L2Voltage =
          windEmulatorStep4_WECSim_B.L1L2Voltage;
        windEmulatorStep4_WECSim_B.BusAssignment.L2L3Voltage =
          windEmulatorStep4_WECSim_B.L2L3Voltage;
        windEmulatorStep4_WECSim_B.BusAssignment.L3L1Volage =
          windEmulatorStep4_WECSim_B.L3L1Voltage;

        /* ToAsyncQueueBlock generated from: '<S36>/invPowerAcs800' */
        slrtLogSignal
          (windEmulatorStep4_WECSim_DW.TAQSigLogging_InsertedFor_invPo.SLRTSigHandles,
           (((windEmulatorStep4_WECSim_M->Timing.clockTick1+
              windEmulatorStep4_WECSim_M->Timing.clockTickH1* 4294967296.0)) *
            0.004));

        /* S-Function (slecatpdorx): '<S11>/L1InaccurateURead' */
        {
          /*------------ S-Function Block: <S11>/L1InaccurateURead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L1InaccurateURead_l;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 987;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 8 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (1 == 8) && (1 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (1 == 8) && (1 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+sigIdx*
                                 1, 1);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(uint16_T))) {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(uint32_T))) {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 1;
          }
        }

        /* Logic: '<S11>/L1InaccurateU' incorporates:
         *  Constant: '<S11>/Constant1'
         */
        windEmulatorStep4_WECSim_B.L1InaccurateU_f =
          (windEmulatorStep4_WECSim_cal->Constant1_Value_cc &&
           windEmulatorStep4_WECSim_B.L1InaccurateURead_l);

        /* S-Function (slecatpdorx): '<S11>/L1InaccurateIRead' */
        {
          /*------------ S-Function Block: <S11>/L1InaccurateIRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L1InaccurateIRead_o;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 988;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 8 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (1 == 8) && (1 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (1 == 8) && (1 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+sigIdx*
                                 1, 1);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(uint16_T))) {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(uint32_T))) {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 1;
          }
        }

        /* Logic: '<S11>/L1InaccurateI' incorporates:
         *  Constant: '<S11>/Constant1'
         */
        windEmulatorStep4_WECSim_B.L1InaccurateI_f =
          (windEmulatorStep4_WECSim_cal->Constant1_Value_cc &&
           windEmulatorStep4_WECSim_B.L1InaccurateIRead_o);

        /* S-Function (slecatpdorx): '<S11>/L1VoltageRead' */
        {
          /*------------ S-Function Block: <S11>/L1VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L1VoltageRead_i;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 1000;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/L1Voltage' */
        windEmulatorStep4_WECSim_B.L1Voltage_f = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L1Voltage_Gain_h) *
          windEmulatorStep4_WECSim_B.L1VoltageRead_i;

        /* S-Function (slecatpdorx): '<S11>/L1CurrentRead' */
        {
          /*------------ S-Function Block: <S11>/L1CurrentRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L1CurrentRead_c;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 1032;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/L1Current' */
        windEmulatorStep4_WECSim_B.L1Current_a = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L1Current_Gain_f) *
          windEmulatorStep4_WECSim_B.L1CurrentRead_c;

        /* S-Function (slecatpdorx): '<S11>/L1PowFactorRead' */
        {
          /*------------ S-Function Block: <S11>/L1PowFactorRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L1PowFactorRead_k;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 1160;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/L1PowFactor' */
        windEmulatorStep4_WECSim_B.L1PowFactor_l = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L1PowFactor_Gain_a) *
          windEmulatorStep4_WECSim_B.L1PowFactorRead_k;

        /* S-Function (slecatpdorx): '<S11>/L1ActivePowRead' */
        {
          /*------------ S-Function Block: <S11>/L1ActivePowRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L1ActivePowRead_h;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 1064;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/L1ActivePow' */
        windEmulatorStep4_WECSim_B.L1ActivePow_m = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L1ActivePow_Gain_o) *
          windEmulatorStep4_WECSim_B.L1ActivePowRead_h;

        /* S-Function (slecatpdorx): '<S11>/L1THDuRead' */
        {
          /*------------ S-Function Block: <S11>/L1THDuRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L1THDuRead_a;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 1272;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/L1THDu' */
        windEmulatorStep4_WECSim_B.L1THDu_l = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L1THDu_Gain_a) *
          windEmulatorStep4_WECSim_B.L1THDuRead_a;

        /* S-Function (slecatpdorx): '<S11>/L1THDiRead' */
        {
          /*------------ S-Function Block: <S11>/L1THDiRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L1THDiRead_a;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 1336;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/L1THDi' */
        windEmulatorStep4_WECSim_B.L1THDi_e = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L1THDi_Gain_f) *
          windEmulatorStep4_WECSim_B.L1THDiRead_a;

        /* S-Function (slecatpdorx): '<S11>/L2InaccurateURead' */
        {
          /*------------ S-Function Block: <S11>/L2InaccurateURead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L2InaccurateURead_h;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 1403;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 8 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (1 == 8) && (1 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (1 == 8) && (1 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+sigIdx*
                                 1, 1);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(uint16_T))) {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(uint32_T))) {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 1;
          }
        }

        /* Logic: '<S11>/L2InaccurateU' incorporates:
         *  Constant: '<S11>/Constant2'
         */
        windEmulatorStep4_WECSim_B.L2InaccurateU_p =
          (windEmulatorStep4_WECSim_cal->Constant2_Value_eh &&
           windEmulatorStep4_WECSim_B.L2InaccurateURead_h);

        /* S-Function (slecatpdorx): '<S11>/L2InaccurateIRead' */
        {
          /*------------ S-Function Block: <S11>/L2InaccurateIRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L2InaccurateIRead_f;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 1404;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 8 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (1 == 8) && (1 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (1 == 8) && (1 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+sigIdx*
                                 1, 1);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(uint16_T))) {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(uint32_T))) {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 1;
          }
        }

        /* Logic: '<S11>/L2InaccurateI' incorporates:
         *  Constant: '<S11>/Constant2'
         */
        windEmulatorStep4_WECSim_B.L2InaccurateI_p =
          (windEmulatorStep4_WECSim_cal->Constant2_Value_eh &&
           windEmulatorStep4_WECSim_B.L2InaccurateIRead_f);

        /* S-Function (slecatpdorx): '<S11>/L2VoltageRead' */
        {
          /*------------ S-Function Block: <S11>/L2VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L2VoltageRead_p;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 1416;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/L2Voltage' */
        windEmulatorStep4_WECSim_B.L2Voltage_j = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L2Voltage_Gain_a) *
          windEmulatorStep4_WECSim_B.L2VoltageRead_p;

        /* S-Function (slecatpdorx): '<S11>/L2CurrentRead' */
        {
          /*------------ S-Function Block: <S11>/L2CurrentRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L2CurrentRead_a;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 1448;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/L2Current' */
        windEmulatorStep4_WECSim_B.L2Current_h = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L2Current_Gain_j) *
          windEmulatorStep4_WECSim_B.L2CurrentRead_a;

        /* S-Function (slecatpdorx): '<S11>/L2PowFactorRead' */
        {
          /*------------ S-Function Block: <S11>/L2PowFactorRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L2PowFactorRead_o;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 1576;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/L2PowFactor' */
        windEmulatorStep4_WECSim_B.L2PowFactor_p = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L2PowFactor_Gain_p) *
          windEmulatorStep4_WECSim_B.L2PowFactorRead_o;

        /* S-Function (slecatpdorx): '<S11>/L2ActivePowRead' */
        {
          /*------------ S-Function Block: <S11>/L2ActivePowRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L2ActivePowRead_e;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 1480;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/L2ActivePow' */
        windEmulatorStep4_WECSim_B.L2ActivePow_b = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L2ActivePow_Gain_b) *
          windEmulatorStep4_WECSim_B.L2ActivePowRead_e;

        /* S-Function (slecatpdorx): '<S11>/L2THDuRead' */
        {
          /*------------ S-Function Block: <S11>/L2THDuRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L2THDuRead_l;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 1688;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/L2THDu' */
        windEmulatorStep4_WECSim_B.L2THDu_i = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L2THDu_Gain_i) *
          windEmulatorStep4_WECSim_B.L2THDuRead_l;

        /* S-Function (slecatpdorx): '<S11>/L2THDiRead' */
        {
          /*------------ S-Function Block: <S11>/L2THDiRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L2THDiRead_j;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 1752;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/L2THDi' */
        windEmulatorStep4_WECSim_B.L2THDi_m = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L2THDi_Gain_j) *
          windEmulatorStep4_WECSim_B.L2THDiRead_j;

        /* S-Function (slecatpdorx): '<S11>/L3InaccurateURead' */
        {
          /*------------ S-Function Block: <S11>/L3InaccurateURead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L3InaccurateURead_o;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 1819;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 8 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (1 == 8) && (1 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (1 == 8) && (1 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+sigIdx*
                                 1, 1);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(uint16_T))) {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(uint32_T))) {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 1;
          }
        }

        /* Logic: '<S11>/L3InaccurateU' incorporates:
         *  Constant: '<S11>/Constant3'
         */
        windEmulatorStep4_WECSim_B.L3InaccurateU_n =
          (windEmulatorStep4_WECSim_cal->Constant3_Value_g &&
           windEmulatorStep4_WECSim_B.L3InaccurateURead_o);

        /* S-Function (slecatpdorx): '<S11>/L3InaccurateIRead' */
        {
          /*------------ S-Function Block: <S11>/L3InaccurateIRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L3InaccurateIRead_p;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 1820;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 8 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (1 == 8) && (1 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (1 == 8) && (1 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+sigIdx*
                                 1, 1);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (1 == 16) && (1 == sizeof(uint16_T))) {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (1 == 32) && (1 == sizeof(uint32_T))) {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 1, sigOutputPtr+
                                   sigIdx*1, 1);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 1;
          }
        }

        /* Logic: '<S11>/L3InaccurateI' incorporates:
         *  Constant: '<S11>/Constant3'
         */
        windEmulatorStep4_WECSim_B.L3InaccurateI_m =
          (windEmulatorStep4_WECSim_cal->Constant3_Value_g &&
           windEmulatorStep4_WECSim_B.L3InaccurateIRead_p);

        /* S-Function (slecatpdorx): '<S11>/L3VoltageRead' */
        {
          /*------------ S-Function Block: <S11>/L3VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L3VoltageRead_h;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 1832;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/L3Voltage' */
        windEmulatorStep4_WECSim_B.L3Voltage_n = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L3Voltage_Gain_h) *
          windEmulatorStep4_WECSim_B.L3VoltageRead_h;

        /* S-Function (slecatpdorx): '<S11>/L3CurrentRead' */
        {
          /*------------ S-Function Block: <S11>/L3CurrentRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L3CurrentRead_h;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 1864;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/L3Current' */
        windEmulatorStep4_WECSim_B.L3Current_h = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L3Current_Gain_i) *
          windEmulatorStep4_WECSim_B.L3CurrentRead_h;

        /* S-Function (slecatpdorx): '<S11>/L3PowFactorRead' */
        {
          /*------------ S-Function Block: <S11>/L3PowFactorRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L3PowFactorRead_o;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 1992;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/L3PowFactor' */
        windEmulatorStep4_WECSim_B.L3PowFactor_o = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L3PowFactor_Gain_k) *
          windEmulatorStep4_WECSim_B.L3PowFactorRead_o;

        /* S-Function (slecatpdorx): '<S11>/L3ActivePowRead' */
        {
          /*------------ S-Function Block: <S11>/L3ActivePowRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L3ActivePowRead_m;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 1896;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/L3ActivePow' */
        windEmulatorStep4_WECSim_B.L3ActivePow_l = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L3ActivePow_Gain_f) *
          windEmulatorStep4_WECSim_B.L3ActivePowRead_m;

        /* S-Function (slecatpdorx): '<S11>/L3THDuRead' */
        {
          /*------------ S-Function Block: <S11>/L3THDuRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L3THDuRead_b;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 2104;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/L3THDu' */
        windEmulatorStep4_WECSim_B.L3THDu_i = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L3THDu_Gain_b) *
          windEmulatorStep4_WECSim_B.L3THDuRead_b;

        /* S-Function (slecatpdorx): '<S11>/L3THDiRead' */
        {
          /*------------ S-Function Block: <S11>/L3THDiRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L3THDiRead_m;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 2168;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/L3THDi' */
        windEmulatorStep4_WECSim_B.L3THDi_j = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L3THDi_Gain_p) *
          windEmulatorStep4_WECSim_B.L3THDiRead_m;

        /* S-Function (slecatpdorx): '<S11>/FrequencyRead' */
        {
          /*------------ S-Function Block: <S11>/FrequencyRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.FrequencyRead_g;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 2280;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/totalFrequency' */
        windEmulatorStep4_WECSim_B.totalFrequency_g = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->totalFrequency_Gain_b) *
          windEmulatorStep4_WECSim_B.FrequencyRead_g;

        /* S-Function (slecatpdorx): '<S11>/totalPowFactorRead' */
        {
          /*------------ S-Function Block: <S11>/totalPowFactorRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.totalPowFactorRead_n;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 2312;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/totalPowFactor' */
        windEmulatorStep4_WECSim_B.totalPowFactor_e = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->totalPowFactor_Gain_f) *
          windEmulatorStep4_WECSim_B.totalPowFactorRead_n;

        /* S-Function (slecatpdorx): '<S11>/totalActivePowRead' */
        {
          /*------------ S-Function Block: <S11>/totalActivePowRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.totalActivePowRead_o;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 2376;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/totalActivePow' */
        windEmulatorStep4_WECSim_B.totalActivePow_g = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->totalActivePow_Gain_m) *
          windEmulatorStep4_WECSim_B.totalActivePowRead_o;

        /* S-Function (slecatpdorx): '<S11>/L1L2VoltageRead' */
        {
          /*------------ S-Function Block: <S11>/L1L2VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L1L2VoltageRead_m;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 2600;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/L1L2Voltage' */
        windEmulatorStep4_WECSim_B.L1L2Voltage_n = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L1L2Voltage_Gain_j) *
          windEmulatorStep4_WECSim_B.L1L2VoltageRead_m;

        /* S-Function (slecatpdorx): '<S11>/L2L3VoltageRead' */
        {
          /*------------ S-Function Block: <S11>/L2L3VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L2L3VoltageRead_e;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 2632;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/L2L3Voltage' */
        windEmulatorStep4_WECSim_B.L2L3Voltage_h = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L2L3Voltage_Gain_m) *
          windEmulatorStep4_WECSim_B.L2L3VoltageRead_e;

        /* S-Function (slecatpdorx): '<S11>/L3L1VoltageRead' */
        {
          /*------------ S-Function Block: <S11>/L3L1VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.L3L1VoltageRead_n;
          uint8_T *ecatRxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatRxBufPtr = (uint8_T *)xpcEtherCATgetPDin( 0 );
          bitOffset = 2664;
          for (sigIdx=0; sigIdx < 1; sigIdx++) {
            switch ( 1 ) {
             case SS_DOUBLE:
              ((real_T *)sigOutputPtr)[sigIdx] = *((real_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_SINGLE:
              ((real32_T *)sigOutputPtr)[sigIdx] = *((real32_T *)(ecatRxBufPtr+
                bitOffset/8));
              break;

             case SS_INT8:
              if ((bitOffset % 8 == 0) && (32 == 8) && (4 == sizeof(int8_T))) {
                ((int8_T *)sigOutputPtr)[sigIdx] = *((int8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT8:
              if ((bitOffset%8 == 0) && (32 == 8) && (4 == sizeof(uint8_T))) {
                ((uint8_T *)sigOutputPtr)[sigIdx] = *((uint8_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_BOOLEAN:
              slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                 sigIdx*4, 4);
              break;

             case SS_INT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(int16_T))) {
                ((int16_T *)sigOutputPtr)[sigIdx] = *((int16_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT16:
              if ((bitOffset%16 == 0) && (32 == 16) && (4 == sizeof(uint16_T)))
              {
                ((uint16_T *)sigOutputPtr)[sigIdx] = *((uint16_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_INT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(int32_T))) {
                ((int32_T *)sigOutputPtr)[sigIdx] = *((int32_T *)(ecatRxBufPtr+
                  bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             case SS_UINT32:
              if ((bitOffset%32 == 0) && (32 == 32) && (4 == sizeof(uint32_T)))
              {
                ((uint32_T *)sigOutputPtr)[sigIdx] = *((uint32_T *)(ecatRxBufPtr
                  +bitOffset/8));
              } else {
                slrtEcatCopyBitsRx(ecatRxBufPtr, bitOffset, 32, sigOutputPtr+
                                   sigIdx*4, 4);
              }
              break;

             default:
              /* Fatal error, unsupported type. This is checked in getDataSizes, so it should never happen. */
              break;
            }

            bitOffset += 32;
          }
        }

        /* Gain: '<S11>/L3L1Voltage' */
        windEmulatorStep4_WECSim_B.L3L1Voltage_m = static_cast<real_T>
          (windEmulatorStep4_WECSim_cal->L3L1Voltage_Gain_k) *
          windEmulatorStep4_WECSim_B.L3L1VoltageRead_n;

        /* BusAssignment: '<S11>/Bus Assignment' */
        windEmulatorStep4_WECSim_B.BusAssignment_l.L1InaccurateU =
          windEmulatorStep4_WECSim_B.L1InaccurateU_f;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L1InaccurateI =
          windEmulatorStep4_WECSim_B.L1InaccurateI_f;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L1Voltage =
          windEmulatorStep4_WECSim_B.L1Voltage_f;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L1Current =
          windEmulatorStep4_WECSim_B.L1Current_a;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L1PowFactor =
          windEmulatorStep4_WECSim_B.L1PowFactor_l;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L1ActivePow =
          windEmulatorStep4_WECSim_B.L1ActivePow_m;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L1THDu =
          windEmulatorStep4_WECSim_B.L1THDu_l;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L1THDi =
          windEmulatorStep4_WECSim_B.L1THDi_e;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L2InaccurateU =
          windEmulatorStep4_WECSim_B.L2InaccurateU_p;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L2InaccurateI =
          windEmulatorStep4_WECSim_B.L2InaccurateI_p;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L2Voltage =
          windEmulatorStep4_WECSim_B.L2Voltage_j;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L2Current =
          windEmulatorStep4_WECSim_B.L2Current_h;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L2PowFactor =
          windEmulatorStep4_WECSim_B.L2PowFactor_p;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L2ActivePow =
          windEmulatorStep4_WECSim_B.L2ActivePow_b;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L2THDu =
          windEmulatorStep4_WECSim_B.L2THDu_i;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L2THDi =
          windEmulatorStep4_WECSim_B.L2THDi_m;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L3InaccurateU =
          windEmulatorStep4_WECSim_B.L3InaccurateU_n;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L3InaccurateI =
          windEmulatorStep4_WECSim_B.L3InaccurateI_m;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L3Voltage =
          windEmulatorStep4_WECSim_B.L3Voltage_n;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L3Current =
          windEmulatorStep4_WECSim_B.L3Current_h;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L3PowFactor =
          windEmulatorStep4_WECSim_B.L3PowFactor_o;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L3ActivePow =
          windEmulatorStep4_WECSim_B.L3ActivePow_l;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L3THDu =
          windEmulatorStep4_WECSim_B.L3THDu_i;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L3THDi =
          windEmulatorStep4_WECSim_B.L3THDi_j;
        windEmulatorStep4_WECSim_B.BusAssignment_l.Frequency =
          windEmulatorStep4_WECSim_B.totalFrequency_g;
        windEmulatorStep4_WECSim_B.BusAssignment_l.TotalPowFactor =
          windEmulatorStep4_WECSim_B.totalPowFactor_e;
        windEmulatorStep4_WECSim_B.BusAssignment_l.TotalActivePow =
          windEmulatorStep4_WECSim_B.totalActivePow_g;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L1L2Voltage =
          windEmulatorStep4_WECSim_B.L1L2Voltage_n;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L2L3Voltage =
          windEmulatorStep4_WECSim_B.L2L3Voltage_h;
        windEmulatorStep4_WECSim_B.BusAssignment_l.L3L1Volage =
          windEmulatorStep4_WECSim_B.L3L1Voltage_m;

        /* ToAsyncQueueBlock generated from: '<S37>/invPowerAcs880' */
        slrtLogSignal
          (windEmulatorStep4_WECSim_DW.TAQSigLogging_InsertedFor_inv_p.SLRTSigHandles,
           (((windEmulatorStep4_WECSim_M->Timing.clockTick1+
              windEmulatorStep4_WECSim_M->Timing.clockTickH1* 4294967296.0)) *
            0.004));

        /* DataTypeConversion: '<S555>/Cast To uint32' incorporates:
         *  Constant: '<S555>/sidType'
         */
        windEmulatorStep4_WECSim_B.CastTouint32 =
          windEmulatorStep4_WECSim_cal->sidType_Value;
      }

      /* Switch: '<S609>/Switch3' */
      if (windEmulatorStep4_WECSim_B.Outofbounds) {
        /* Switch: '<S609>/Switch3' incorporates:
         *  Constant: '<S609>/Set bound'
         */
        windEmulatorStep4_WECSim_B.Switch3 =
          windEmulatorStep4_WECSim_cal->Setbound_Value;
      } else {
        /* Switch: '<S609>/Switch3' incorporates:
         *  Inport: '<Root>/inportCaseCounter'
         */
        windEmulatorStep4_WECSim_B.Switch3 =
          windEmulatorStep4_WECSim_U.inportCaseCounter;
      }

      /* End of Switch: '<S609>/Switch3' */

      /* Gain: '<S609>/caseCounterSignalsNow' */
      windEmulatorStep4_WECSim_B.caseCounterSignalsNow =
        windEmulatorStep4_WECSim_cal->caseCounterSignalsNow_Gain *
        windEmulatorStep4_WECSim_B.Switch3;

      /* BusAssignment: '<S555>/Bus Assignment' */
      windEmulatorStep4_WECSim_B.BusAssignment_j.acs800TorqueSetpoint_Nm =
        windEmulatorStep4_WECSim_B.Product_g;
      windEmulatorStep4_WECSim_B.BusAssignment_j.acs880SpeedSetpoint_rpm =
        windEmulatorStep4_WECSim_B.Product1;
      windEmulatorStep4_WECSim_B.BusAssignment_j.sidType =
        windEmulatorStep4_WECSim_B.CastTouint32;
      windEmulatorStep4_WECSim_B.BusAssignment_j.fromFileCaseCounter =
        windEmulatorStep4_WECSim_B.caseCounterSignalsNow;
      if (tmp_g) {
        /* ToAsyncQueueBlock generated from: '<S40>/sidInfoSignals' */
        slrtLogSignal
          (windEmulatorStep4_WECSim_DW.TAQSigLogging_InsertedFor_sidIn.SLRTSigHandles,
           (((windEmulatorStep4_WECSim_M->Timing.clockTick1+
              windEmulatorStep4_WECSim_M->Timing.clockTickH1* 4294967296.0)) *
            0.004));

        /* Constant: '<S5>/Constant' */
        windEmulatorStep4_WECSim_B.Constant =
          windEmulatorStep4_WECSim_cal->Constant_Value_m3;

        /* S-Function (slrealtimeenablelogging): '<S5>/Enable File Log' */

        /* Level2 S-Function Block: '<S5>/Enable File Log' (slrealtimeenablelogging) */
        {
          SimStruct *rts = windEmulatorStep4_WECSim_M->childSfunctions[0];
          sfcnOutputs(rts,0);
        }

        /* Memory: '<S29>/Memory' */
        windEmulatorStep4_WECSim_B.Memory =
          windEmulatorStep4_WECSim_DW.Memory_PreviousInput;

        /* Sum: '<S29>/Sum' incorporates:
         *  Constant: '<S29>/loopAdd'
         */
        windEmulatorStep4_WECSim_B.Sum_e =
          windEmulatorStep4_WECSim_cal->loopAdd_Value +
          windEmulatorStep4_WECSim_B.Memory;

        /* Gain: '<S29>/loopCounter' */
        windEmulatorStep4_WECSim_B.loopCounter =
          windEmulatorStep4_WECSim_cal->loopCounter_Gain *
          windEmulatorStep4_WECSim_B.Sum_e;

        /* DataTypeConversion: '<S29>/Data Type Conversion' */
        windEmulatorStep4_WECSim_B.DataTypeConversion_f =
          windEmulatorStep4_WECSim_B.loopCounter;

        /* Gain: '<S29>/time_s' */
        windEmulatorStep4_WECSim_B.time_s = *get_Ts() *
          windEmulatorStep4_WECSim_B.DataTypeConversion_f;

        /* S-Function (slecatpdotx): '<S14>/ACS800CtrlWord' */
        {
          /*------------ S-Function Block: <S14>/ACS800CtrlWord PDO transmit block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          int_T i;
          uint8_T *sigInputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.BusAssignment_kc.ctrlWord;
          uint8_T *ecatTxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatTxBufPtr = (uint8_T *)xpcEtherCATgetPDout( 0 );
          bitOffset = 520;
          for (i = 0; i < 1; i++) {
            switch ( 5 ){
             case SS_DOUBLE:
              *((real_T *)(ecatTxBufPtr+bitOffset/8)) = ((real_T *)sigInputPtr)
                [i];
              break;

             case SS_SINGLE:
              *((real32_T *)(ecatTxBufPtr+bitOffset/8)) = ((real32_T *)
                sigInputPtr)[i];
              break;

             case SS_INT8:
              if ((16 == 8) && (bitOffset%8 == 0)) {
                *((int8_T *)(ecatTxBufPtr+bitOffset/8)) = ((int8_T *)sigInputPtr)
                  [i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int8_T *)sigInputPtr)));
              }
              break;

             case SS_UINT8:
              if ((16 == 8) && (bitOffset%8 == 0)) {
                *((uint8_T *)(ecatTxBufPtr+bitOffset/8)) = ((uint8_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((uint8_T *)sigInputPtr)));
              }
              break;

             case SS_BOOLEAN:
              if ((16 == 8) && (bitOffset%8 == 0)) {
                *((int8_T *)(ecatTxBufPtr+bitOffset/8)) = ((int8_T *)sigInputPtr)
                  [i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int8_T *)sigInputPtr)));
              }
              break;

             case SS_INT16:
              if ((16 == 16) && (bitOffset%16 == 0)) {
                *((int16_T *)(ecatTxBufPtr+bitOffset/8)) = ((int16_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int16_T *)sigInputPtr)));
              }
              break;

             case SS_UINT16:
              if ((16 == 16) && (bitOffset%16 == 0)) {
                *((uint16_T *)(ecatTxBufPtr+bitOffset/8)) =((uint16_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((uint16_T *)sigInputPtr)));
              }
              break;

             case SS_INT32:
              if ((16 == 32) && (bitOffset%32 == 0)) {
                *((int32_T *)(ecatTxBufPtr+bitOffset/8)) = ((int32_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int32_T *)sigInputPtr)));
              }
              break;

             case SS_UINT32:
              if ((16 == 32) && (bitOffset%32 == 0)) {
                *((uint32_T *)(ecatTxBufPtr+bitOffset/8)) = ((uint32_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((uint32_T *)sigInputPtr)));
              }
              break;

             default:
              /* Fatal error. Should never happen as this is checked in mdlStart. */
              break;
            }

            bitOffset += 16;
          }
        }
      }

      /* Product: '<S611>/Product' incorporates:
       *  Constant: '<S611>/Constant EngineeringValue'
       *  Constant: '<S611>/Constant FieldbusValue'
       */
      windEmulatorStep4_WECSim_B.Product_k =
        windEmulatorStep4_WECSim_B.BusAssignment_kc.torqueSetpoint_percent *
        *get_acs800TorqueNomFb() / *get_acs800TorqueNomEng();

      /* DataTypeConversion: '<S14>/ACS800Ref2Int16' */
      riseValLimit = std::floor(windEmulatorStep4_WECSim_B.Product_k);
      if (rtIsNaN(riseValLimit) || rtIsInf(riseValLimit)) {
        riseValLimit = 0.0;
      } else {
        riseValLimit = std::fmod(riseValLimit, 65536.0);
      }

      /* DataTypeConversion: '<S14>/ACS800Ref2Int16' */
      windEmulatorStep4_WECSim_B.ACS800Ref2Int16 = static_cast<int16_T>
        (riseValLimit < 0.0 ? static_cast<int32_T>(static_cast<int16_T>(-
           static_cast<int16_T>(static_cast<uint16_T>(-riseValLimit)))) :
         static_cast<int32_T>(static_cast<int16_T>(static_cast<uint16_T>
           (riseValLimit))));
      if (tmp_g) {
        /* S-Function (slecatpdotx): '<S14>/ACS800TorqueSetpoint' */
        {
          /*------------ S-Function Block: <S14>/ACS800TorqueSetpoint PDO transmit block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          int_T i;
          uint8_T *sigInputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.ACS800Ref2Int16;
          uint8_T *ecatTxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatTxBufPtr = (uint8_T *)xpcEtherCATgetPDout( 0 );
          bitOffset = 552;
          for (i = 0; i < 1; i++) {
            switch ( 4 ){
             case SS_DOUBLE:
              *((real_T *)(ecatTxBufPtr+bitOffset/8)) = ((real_T *)sigInputPtr)
                [i];
              break;

             case SS_SINGLE:
              *((real32_T *)(ecatTxBufPtr+bitOffset/8)) = ((real32_T *)
                sigInputPtr)[i];
              break;

             case SS_INT8:
              if ((16 == 8) && (bitOffset%8 == 0)) {
                *((int8_T *)(ecatTxBufPtr+bitOffset/8)) = ((int8_T *)sigInputPtr)
                  [i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int8_T *)sigInputPtr)));
              }
              break;

             case SS_UINT8:
              if ((16 == 8) && (bitOffset%8 == 0)) {
                *((uint8_T *)(ecatTxBufPtr+bitOffset/8)) = ((uint8_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((uint8_T *)sigInputPtr)));
              }
              break;

             case SS_BOOLEAN:
              if ((16 == 8) && (bitOffset%8 == 0)) {
                *((int8_T *)(ecatTxBufPtr+bitOffset/8)) = ((int8_T *)sigInputPtr)
                  [i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int8_T *)sigInputPtr)));
              }
              break;

             case SS_INT16:
              if ((16 == 16) && (bitOffset%16 == 0)) {
                *((int16_T *)(ecatTxBufPtr+bitOffset/8)) = ((int16_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int16_T *)sigInputPtr)));
              }
              break;

             case SS_UINT16:
              if ((16 == 16) && (bitOffset%16 == 0)) {
                *((uint16_T *)(ecatTxBufPtr+bitOffset/8)) =((uint16_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((uint16_T *)sigInputPtr)));
              }
              break;

             case SS_INT32:
              if ((16 == 32) && (bitOffset%32 == 0)) {
                *((int32_T *)(ecatTxBufPtr+bitOffset/8)) = ((int32_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int32_T *)sigInputPtr)));
              }
              break;

             case SS_UINT32:
              if ((16 == 32) && (bitOffset%32 == 0)) {
                *((uint32_T *)(ecatTxBufPtr+bitOffset/8)) = ((uint32_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((uint32_T *)sigInputPtr)));
              }
              break;

             default:
              /* Fatal error. Should never happen as this is checked in mdlStart. */
              break;
            }

            bitOffset += 16;
          }
        }

        /* Product: '<S610>/Product' incorporates:
         *  Constant: '<S14>/Constant'
         *  Constant: '<S610>/Constant EngineeringValue'
         *  Constant: '<S610>/Constant FieldbusValue'
         */
        windEmulatorStep4_WECSim_B.Product_d =
          windEmulatorStep4_WECSim_cal->Constant_Value_m * *get_acs800SpeedNomFb
          () / *get_acs800SpeedNomEng();

        /* DataTypeConversion: '<S14>/ACS800Ref1Int16' */
        riseValLimit = std::floor(windEmulatorStep4_WECSim_B.Product_d);
        if (rtIsNaN(riseValLimit) || rtIsInf(riseValLimit)) {
          riseValLimit = 0.0;
        } else {
          riseValLimit = std::fmod(riseValLimit, 65536.0);
        }

        /* DataTypeConversion: '<S14>/ACS800Ref1Int16' */
        windEmulatorStep4_WECSim_B.ACS800Ref1Int16 = static_cast<int16_T>
          (riseValLimit < 0.0 ? static_cast<int32_T>(static_cast<int16_T>(-
             static_cast<int16_T>(static_cast<uint16_T>(-riseValLimit)))) :
           static_cast<int32_T>(static_cast<int16_T>(static_cast<uint16_T>
             (riseValLimit))));

        /* S-Function (slecatpdotx): '<S14>/ACS800SpeedSetpoint' */
        {
          /*------------ S-Function Block: <S14>/ACS800SpeedSetpoint PDO transmit block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          int_T i;
          uint8_T *sigInputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.ACS800Ref1Int16;
          uint8_T *ecatTxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatTxBufPtr = (uint8_T *)xpcEtherCATgetPDout( 0 );
          bitOffset = 536;
          for (i = 0; i < 1; i++) {
            switch ( 4 ){
             case SS_DOUBLE:
              *((real_T *)(ecatTxBufPtr+bitOffset/8)) = ((real_T *)sigInputPtr)
                [i];
              break;

             case SS_SINGLE:
              *((real32_T *)(ecatTxBufPtr+bitOffset/8)) = ((real32_T *)
                sigInputPtr)[i];
              break;

             case SS_INT8:
              if ((16 == 8) && (bitOffset%8 == 0)) {
                *((int8_T *)(ecatTxBufPtr+bitOffset/8)) = ((int8_T *)sigInputPtr)
                  [i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int8_T *)sigInputPtr)));
              }
              break;

             case SS_UINT8:
              if ((16 == 8) && (bitOffset%8 == 0)) {
                *((uint8_T *)(ecatTxBufPtr+bitOffset/8)) = ((uint8_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((uint8_T *)sigInputPtr)));
              }
              break;

             case SS_BOOLEAN:
              if ((16 == 8) && (bitOffset%8 == 0)) {
                *((int8_T *)(ecatTxBufPtr+bitOffset/8)) = ((int8_T *)sigInputPtr)
                  [i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int8_T *)sigInputPtr)));
              }
              break;

             case SS_INT16:
              if ((16 == 16) && (bitOffset%16 == 0)) {
                *((int16_T *)(ecatTxBufPtr+bitOffset/8)) = ((int16_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int16_T *)sigInputPtr)));
              }
              break;

             case SS_UINT16:
              if ((16 == 16) && (bitOffset%16 == 0)) {
                *((uint16_T *)(ecatTxBufPtr+bitOffset/8)) =((uint16_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((uint16_T *)sigInputPtr)));
              }
              break;

             case SS_INT32:
              if ((16 == 32) && (bitOffset%32 == 0)) {
                *((int32_T *)(ecatTxBufPtr+bitOffset/8)) = ((int32_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int32_T *)sigInputPtr)));
              }
              break;

             case SS_UINT32:
              if ((16 == 32) && (bitOffset%32 == 0)) {
                *((uint32_T *)(ecatTxBufPtr+bitOffset/8)) = ((uint32_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((uint32_T *)sigInputPtr)));
              }
              break;

             default:
              /* Fatal error. Should never happen as this is checked in mdlStart. */
              break;
            }

            bitOffset += 16;
          }
        }

        /* S-Function (slecatpdotx): '<S15>/ACS880CtrlWord' */
        {
          /*------------ S-Function Block: <S15>/ACS880CtrlWord PDO transmit block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          int_T i;
          uint8_T *sigInputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.BusAssignment_k.ctrlWord;
          uint8_T *ecatTxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatTxBufPtr = (uint8_T *)xpcEtherCATgetPDout( 0 );
          bitOffset = 344;
          for (i = 0; i < 1; i++) {
            switch ( 5 ){
             case SS_DOUBLE:
              *((real_T *)(ecatTxBufPtr+bitOffset/8)) = ((real_T *)sigInputPtr)
                [i];
              break;

             case SS_SINGLE:
              *((real32_T *)(ecatTxBufPtr+bitOffset/8)) = ((real32_T *)
                sigInputPtr)[i];
              break;

             case SS_INT8:
              if ((16 == 8) && (bitOffset%8 == 0)) {
                *((int8_T *)(ecatTxBufPtr+bitOffset/8)) = ((int8_T *)sigInputPtr)
                  [i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int8_T *)sigInputPtr)));
              }
              break;

             case SS_UINT8:
              if ((16 == 8) && (bitOffset%8 == 0)) {
                *((uint8_T *)(ecatTxBufPtr+bitOffset/8)) = ((uint8_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((uint8_T *)sigInputPtr)));
              }
              break;

             case SS_BOOLEAN:
              if ((16 == 8) && (bitOffset%8 == 0)) {
                *((int8_T *)(ecatTxBufPtr+bitOffset/8)) = ((int8_T *)sigInputPtr)
                  [i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int8_T *)sigInputPtr)));
              }
              break;

             case SS_INT16:
              if ((16 == 16) && (bitOffset%16 == 0)) {
                *((int16_T *)(ecatTxBufPtr+bitOffset/8)) = ((int16_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int16_T *)sigInputPtr)));
              }
              break;

             case SS_UINT16:
              if ((16 == 16) && (bitOffset%16 == 0)) {
                *((uint16_T *)(ecatTxBufPtr+bitOffset/8)) =((uint16_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((uint16_T *)sigInputPtr)));
              }
              break;

             case SS_INT32:
              if ((16 == 32) && (bitOffset%32 == 0)) {
                *((int32_T *)(ecatTxBufPtr+bitOffset/8)) = ((int32_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int32_T *)sigInputPtr)));
              }
              break;

             case SS_UINT32:
              if ((16 == 32) && (bitOffset%32 == 0)) {
                *((uint32_T *)(ecatTxBufPtr+bitOffset/8)) = ((uint32_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((uint32_T *)sigInputPtr)));
              }
              break;

             default:
              /* Fatal error. Should never happen as this is checked in mdlStart. */
              break;
            }

            bitOffset += 16;
          }
        }
      }

      /* Product: '<S612>/Product' incorporates:
       *  Constant: '<S612>/Constant EngineeringValue'
       *  Constant: '<S612>/Constant FieldbusValue'
       */
      windEmulatorStep4_WECSim_B.Product_l =
        windEmulatorStep4_WECSim_B.BusAssignment_k.torqueSetpoint_percent *
        *get_acs880TorqueFieldbusScale() / *get_acs880TorqueSetpointScaling();

      /* DataTypeConversion: '<S15>/ACS880Ref2Int16' */
      riseValLimit = std::floor(windEmulatorStep4_WECSim_B.Product_l);
      if (rtIsNaN(riseValLimit) || rtIsInf(riseValLimit)) {
        riseValLimit = 0.0;
      } else {
        riseValLimit = std::fmod(riseValLimit, 65536.0);
      }

      /* DataTypeConversion: '<S15>/ACS880Ref2Int16' */
      windEmulatorStep4_WECSim_B.ACS880Ref2Int16 = static_cast<int16_T>
        (riseValLimit < 0.0 ? static_cast<int32_T>(static_cast<int16_T>(-
           static_cast<int16_T>(static_cast<uint16_T>(-riseValLimit)))) :
         static_cast<int32_T>(static_cast<int16_T>(static_cast<uint16_T>
           (riseValLimit))));
      if (tmp_g) {
        /* S-Function (slecatpdotx): '<S15>/ACS880TorqueSetpoint' */
        {
          /*------------ S-Function Block: <S15>/ACS880TorqueSetpoint PDO transmit block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          int_T i;
          uint8_T *sigInputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.ACS880Ref2Int16;
          uint8_T *ecatTxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatTxBufPtr = (uint8_T *)xpcEtherCATgetPDout( 0 );
          bitOffset = 376;
          for (i = 0; i < 1; i++) {
            switch ( 4 ){
             case SS_DOUBLE:
              *((real_T *)(ecatTxBufPtr+bitOffset/8)) = ((real_T *)sigInputPtr)
                [i];
              break;

             case SS_SINGLE:
              *((real32_T *)(ecatTxBufPtr+bitOffset/8)) = ((real32_T *)
                sigInputPtr)[i];
              break;

             case SS_INT8:
              if ((16 == 8) && (bitOffset%8 == 0)) {
                *((int8_T *)(ecatTxBufPtr+bitOffset/8)) = ((int8_T *)sigInputPtr)
                  [i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int8_T *)sigInputPtr)));
              }
              break;

             case SS_UINT8:
              if ((16 == 8) && (bitOffset%8 == 0)) {
                *((uint8_T *)(ecatTxBufPtr+bitOffset/8)) = ((uint8_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((uint8_T *)sigInputPtr)));
              }
              break;

             case SS_BOOLEAN:
              if ((16 == 8) && (bitOffset%8 == 0)) {
                *((int8_T *)(ecatTxBufPtr+bitOffset/8)) = ((int8_T *)sigInputPtr)
                  [i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int8_T *)sigInputPtr)));
              }
              break;

             case SS_INT16:
              if ((16 == 16) && (bitOffset%16 == 0)) {
                *((int16_T *)(ecatTxBufPtr+bitOffset/8)) = ((int16_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int16_T *)sigInputPtr)));
              }
              break;

             case SS_UINT16:
              if ((16 == 16) && (bitOffset%16 == 0)) {
                *((uint16_T *)(ecatTxBufPtr+bitOffset/8)) =((uint16_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((uint16_T *)sigInputPtr)));
              }
              break;

             case SS_INT32:
              if ((16 == 32) && (bitOffset%32 == 0)) {
                *((int32_T *)(ecatTxBufPtr+bitOffset/8)) = ((int32_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int32_T *)sigInputPtr)));
              }
              break;

             case SS_UINT32:
              if ((16 == 32) && (bitOffset%32 == 0)) {
                *((uint32_T *)(ecatTxBufPtr+bitOffset/8)) = ((uint32_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((uint32_T *)sigInputPtr)));
              }
              break;

             default:
              /* Fatal error. Should never happen as this is checked in mdlStart. */
              break;
            }

            bitOffset += 16;
          }
        }

        /* DataTypeConversion: '<S15>/ACS880Ref1Int16' incorporates:
         *  Constant: '<S15>/Constant2'
         */
        riseValLimit = std::floor(windEmulatorStep4_WECSim_cal->Constant2_Value);
        if (rtIsNaN(riseValLimit) || rtIsInf(riseValLimit)) {
          riseValLimit = 0.0;
        } else {
          riseValLimit = std::fmod(riseValLimit, 65536.0);
        }

        /* DataTypeConversion: '<S15>/ACS880Ref1Int16' */
        windEmulatorStep4_WECSim_B.ACS880Ref1Int16 = static_cast<int16_T>
          (riseValLimit < 0.0 ? static_cast<int32_T>(static_cast<int16_T>(-
             static_cast<int16_T>(static_cast<uint16_T>(-riseValLimit)))) :
           static_cast<int32_T>(static_cast<int16_T>(static_cast<uint16_T>
             (riseValLimit))));

        /* S-Function (slecatpdotx): '<S15>/ACS880SpeedSetpoint' */
        {
          /*------------ S-Function Block: <S15>/ACS880SpeedSetpoint PDO transmit block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          int_T i;
          uint8_T *sigInputPtr = (uint8_T *)
            &windEmulatorStep4_WECSim_B.ACS880Ref1Int16;
          uint8_T *ecatTxBufPtr;       // Pointer to the stack PDO rx buffer
          ecatTxBufPtr = (uint8_T *)xpcEtherCATgetPDout( 0 );
          bitOffset = 360;
          for (i = 0; i < 1; i++) {
            switch ( 4 ){
             case SS_DOUBLE:
              *((real_T *)(ecatTxBufPtr+bitOffset/8)) = ((real_T *)sigInputPtr)
                [i];
              break;

             case SS_SINGLE:
              *((real32_T *)(ecatTxBufPtr+bitOffset/8)) = ((real32_T *)
                sigInputPtr)[i];
              break;

             case SS_INT8:
              if ((16 == 8) && (bitOffset%8 == 0)) {
                *((int8_T *)(ecatTxBufPtr+bitOffset/8)) = ((int8_T *)sigInputPtr)
                  [i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int8_T *)sigInputPtr)));
              }
              break;

             case SS_UINT8:
              if ((16 == 8) && (bitOffset%8 == 0)) {
                *((uint8_T *)(ecatTxBufPtr+bitOffset/8)) = ((uint8_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((uint8_T *)sigInputPtr)));
              }
              break;

             case SS_BOOLEAN:
              if ((16 == 8) && (bitOffset%8 == 0)) {
                *((int8_T *)(ecatTxBufPtr+bitOffset/8)) = ((int8_T *)sigInputPtr)
                  [i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int8_T *)sigInputPtr)));
              }
              break;

             case SS_INT16:
              if ((16 == 16) && (bitOffset%16 == 0)) {
                *((int16_T *)(ecatTxBufPtr+bitOffset/8)) = ((int16_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int16_T *)sigInputPtr)));
              }
              break;

             case SS_UINT16:
              if ((16 == 16) && (bitOffset%16 == 0)) {
                *((uint16_T *)(ecatTxBufPtr+bitOffset/8)) =((uint16_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((uint16_T *)sigInputPtr)));
              }
              break;

             case SS_INT32:
              if ((16 == 32) && (bitOffset%32 == 0)) {
                *((int32_T *)(ecatTxBufPtr+bitOffset/8)) = ((int32_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((int32_T *)sigInputPtr)));
              }
              break;

             case SS_UINT32:
              if ((16 == 32) && (bitOffset%32 == 0)) {
                *((uint32_T *)(ecatTxBufPtr+bitOffset/8)) = ((uint32_T *)
                  sigInputPtr)[i];
              } else {
                slrtEcatCopyBitsTx((uint8_T *)ecatTxBufPtr, bitOffset, 16,
                                   (uint32_T)(*((uint32_T *)sigInputPtr)));
              }
              break;

             default:
              /* Fatal error. Should never happen as this is checked in mdlStart. */
              break;
            }

            bitOffset += 16;
          }
        }
      }

      /* SimscapeExecutionBlock: '<S216>/STATE_1' incorporates:
       *  SimscapeExecutionBlock: '<S216>/OUTPUT_1_0'
       *  SimscapeExecutionBlock: '<S216>/OUTPUT_1_1'
       *  SimscapeExecutionBlock: '<S332>/OUTPUT_1_0'
       *  SimscapeExecutionBlock: '<S332>/STATE_1'
       *  SimscapeInputBlock: '<S332>/INPUT_1_1_1'
       */
      simulationData = static_cast<NeslSimulationData *>
        (windEmulatorStep4_WECSim_DW.STATE_1_SimData_h);
      u1 = Clock_tmp;
      time_3 = u1;
      simulationData->mData->mTime.mN = 1;
      simulationData->mData->mTime.mX = &time_3;
      simulationData->mData->mContStates.mN = 2;
      simulationData->mData->mContStates.mX =
        &windEmulatorStep4_WECSim_X.windEmulatorStep4_WECSimhptoSim[0];
      simulationData->mData->mDiscStates.mN = 0;
      simulationData->mData->mDiscStates.mX =
        &windEmulatorStep4_WECSim_DW.STATE_1_Discrete;
      simulationData->mData->mModeVector.mN = 0;
      simulationData->mData->mModeVector.mX =
        &windEmulatorStep4_WECSim_DW.STATE_1_Modes_j;
      f = false;
      simulationData->mData->mFoundZcEvents = f;
      simulationData->mData->mHadEvents = false;
      f = rtmIsMajorTimeStep(windEmulatorStep4_WECSim_M);
      simulationData->mData->mIsMajorTimeStep = f;
      tmp_11 = false;
      simulationData->mData->mIsSolverAssertCheck = tmp_11;
      simulationData->mData->mIsSolverCheckingCIC = false;
      tmp_11 = rtsiIsSolverComputingJacobian
        (&windEmulatorStep4_WECSim_M->solverInfo);
      simulationData->mData->mIsComputingJacobian = tmp_11;
      simulationData->mData->mIsEvaluatingF0 = false;
      simulationData->mData->mIsSolverRequestingReset = false;
      simulationData->mData->mIsModeUpdateTimeStep = tmp_h;
      tmp_5[0] = 0;
      tmp_4[0] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[0];
      tmp_4[1] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[1];
      tmp_4[2] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[2];
      tmp_4[3] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[3];
      tmp_5[1] = 4;
      tmp_4[4] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[0];
      tmp_4[5] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[1];
      tmp_4[6] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[2];
      tmp_4[7] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[3];
      tmp_5[2] = 8;
      tmp_4[8] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[0];
      tmp_4[9] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[1];
      tmp_4[10] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[2];
      tmp_4[11] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[3];
      tmp_5[3] = 12;
      tmp_4[12] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[0];
      tmp_4[13] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[1];
      tmp_4[14] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[2];
      tmp_4[15] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[3];
      tmp_5[4] = 16;
      tmp_4[16] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[0];
      tmp_4[17] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[1];
      tmp_4[18] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[2];
      tmp_4[19] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[3];
      tmp_5[5] = 20;
      tmp_4[20] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[0];
      tmp_4[21] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[1];
      tmp_4[22] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[2];
      tmp_4[23] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[3];
      tmp_5[6] = 24;
      tmp_4[24] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[0];
      tmp_4[25] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[1];
      tmp_4[26] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[2];
      tmp_4[27] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[3];
      tmp_5[7] = 28;
      tmp_4[28] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[0];
      tmp_4[29] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[1];
      tmp_4[30] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[2];
      tmp_4[31] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[3];
      tmp_5[8] = 32;
      tmp_4[32] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[0];
      tmp_4[33] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[1];
      tmp_4[34] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[2];
      tmp_4[35] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[3];
      tmp_5[9] = 36;
      tmp_4[36] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[0];
      tmp_4[37] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[1];
      tmp_4[38] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[2];
      tmp_4[39] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[3];
      tmp_5[10] = 40;
      tmp_4[40] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[0];
      tmp_4[41] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[1];
      tmp_4[42] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[2];
      tmp_4[43] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[3];
      tmp_5[11] = 44;
      tmp_4[44] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[0];
      tmp_4[45] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[1];
      tmp_4[46] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[2];
      tmp_4[47] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[3];
      tmp_5[12] = 48;
      tmp_4[48] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[0];
      tmp_4[49] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[1];
      tmp_4[50] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[2];
      tmp_4[51] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[3];
      tmp_5[13] = 52;
      simulationData->mData->mInputValues.mN = 52;
      simulationData->mData->mInputValues.mX = &tmp_4[0];
      simulationData->mData->mInputOffsets.mN = 14;
      simulationData->mData->mInputOffsets.mX = &tmp_5[0];
      simulationData->mData->mOutputs.mN = 2;
      simulationData->mData->mOutputs.mX =
        &windEmulatorStep4_WECSim_B.STATE_1_p[0];
      simulationData->mData->mTolerances.mN = 0;
      simulationData->mData->mTolerances.mX = NULL;
      simulationData->mData->mCstateHasChanged = false;
      simulationData->mData->mDstateHasChanged = false;
      time_4 = deltaT_tmp;
      simulationData->mData->mTime.mN = 1;
      simulationData->mData->mTime.mX = &time_4;
      simulationData->mData->mSampleHits.mN = 0;
      simulationData->mData->mSampleHits.mX = NULL;
      simulationData->mData->mIsFundamentalSampleHit = false;
      simulationData->mData->mHadEvents = false;
      simulator = static_cast<NeslSimulator *>
        (windEmulatorStep4_WECSim_DW.STATE_1_Simulator_f);
      diag = static_cast<NeuDiagnosticManager *>
        (windEmulatorStep4_WECSim_DW.STATE_1_DiagMgr_o);
      diagTree = neu_diagnostic_manager_get_initial_tree(diag);
      i = ne_simulator_method(simulator, NESL_SIM_OUTPUTS, simulationData, diag);
      if (i != 0) {
        tmp_11 = error_buffer_is_empty(rtmGetErrorStatus
          (windEmulatorStep4_WECSim_M));
        if (tmp_11) {
          msg = rtw_diagnostics_msg(diagTree);
          rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
        }
      }

      /* SimscapeExecutionBlock: '<S216>/OUTPUT_1_1' */
      simulationData = static_cast<NeslSimulationData *>
        (windEmulatorStep4_WECSim_DW.OUTPUT_1_1_SimData);
      time_5 = u1;
      simulationData->mData->mTime.mN = 1;
      simulationData->mData->mTime.mX = &time_5;
      simulationData->mData->mContStates.mN = 0;
      simulationData->mData->mContStates.mX = NULL;
      simulationData->mData->mDiscStates.mN = 0;
      simulationData->mData->mDiscStates.mX =
        &windEmulatorStep4_WECSim_DW.OUTPUT_1_1_Discrete;
      simulationData->mData->mModeVector.mN = 0;
      simulationData->mData->mModeVector.mX =
        &windEmulatorStep4_WECSim_DW.OUTPUT_1_1_Modes;
      tmp_11 = false;
      simulationData->mData->mFoundZcEvents = tmp_11;
      simulationData->mData->mHadEvents = false;
      simulationData->mData->mIsMajorTimeStep = f;
      tmp_11 = false;
      simulationData->mData->mIsSolverAssertCheck = tmp_11;
      simulationData->mData->mIsSolverCheckingCIC = false;
      simulationData->mData->mIsComputingJacobian = false;
      simulationData->mData->mIsEvaluatingF0 = false;
      simulationData->mData->mIsSolverRequestingReset = false;
      simulationData->mData->mIsModeUpdateTimeStep = tmp_h;
      tmp_7[0] = 0;
      tmp_6[0] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[0];
      tmp_6[1] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[1];
      tmp_6[2] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[2];
      tmp_6[3] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[3];
      tmp_7[1] = 4;
      tmp_6[4] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[0];
      tmp_6[5] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[1];
      tmp_6[6] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[2];
      tmp_6[7] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[3];
      tmp_7[2] = 8;
      tmp_6[8] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[0];
      tmp_6[9] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[1];
      tmp_6[10] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[2];
      tmp_6[11] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[3];
      tmp_7[3] = 12;
      tmp_6[12] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[0];
      tmp_6[13] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[1];
      tmp_6[14] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[2];
      tmp_6[15] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[3];
      tmp_7[4] = 16;
      tmp_6[16] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[0];
      tmp_6[17] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[1];
      tmp_6[18] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[2];
      tmp_6[19] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[3];
      tmp_7[5] = 20;
      tmp_6[20] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[0];
      tmp_6[21] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[1];
      tmp_6[22] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[2];
      tmp_6[23] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[3];
      tmp_7[6] = 24;
      tmp_6[24] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[0];
      tmp_6[25] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[1];
      tmp_6[26] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[2];
      tmp_6[27] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[3];
      tmp_7[7] = 28;
      tmp_6[28] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[0];
      tmp_6[29] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[1];
      tmp_6[30] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[2];
      tmp_6[31] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[3];
      tmp_7[8] = 32;
      tmp_6[32] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[0];
      tmp_6[33] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[1];
      tmp_6[34] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[2];
      tmp_6[35] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[3];
      tmp_7[9] = 36;
      tmp_6[36] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[0];
      tmp_6[37] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[1];
      tmp_6[38] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[2];
      tmp_6[39] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[3];
      tmp_7[10] = 40;
      tmp_6[40] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[0];
      tmp_6[41] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[1];
      tmp_6[42] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[2];
      tmp_6[43] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[3];
      tmp_7[11] = 44;
      tmp_6[44] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[0];
      tmp_6[45] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[1];
      tmp_6[46] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[2];
      tmp_6[47] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[3];
      tmp_7[12] = 48;
      tmp_6[48] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[0];
      tmp_6[49] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[1];
      tmp_6[50] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[2];
      tmp_6[51] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[3];
      tmp_7[13] = 52;
      tmp_6[52] = windEmulatorStep4_WECSim_B.STATE_1_p[0];
      tmp_6[53] = windEmulatorStep4_WECSim_B.STATE_1_p[1];
      tmp_7[14] = 54;
      simulationData->mData->mInputValues.mN = 54;
      simulationData->mData->mInputValues.mX = &tmp_6[0];
      simulationData->mData->mInputOffsets.mN = 15;
      simulationData->mData->mInputOffsets.mX = &tmp_7[0];
      simulationData->mData->mOutputs.mN = 28;
      simulationData->mData->mOutputs.mX =
        &windEmulatorStep4_WECSim_B.OUTPUT_1_1[0];
      simulationData->mData->mTolerances.mN = 0;
      simulationData->mData->mTolerances.mX = NULL;
      simulationData->mData->mCstateHasChanged = false;
      simulationData->mData->mDstateHasChanged = false;
      time_6 = deltaT_tmp;
      simulationData->mData->mTime.mN = 1;
      simulationData->mData->mTime.mX = &time_6;
      simulationData->mData->mSampleHits.mN = 0;
      simulationData->mData->mSampleHits.mX = NULL;
      simulationData->mData->mIsFundamentalSampleHit = false;
      simulationData->mData->mHadEvents = false;
      simulator = static_cast<NeslSimulator *>
        (windEmulatorStep4_WECSim_DW.OUTPUT_1_1_Simulator);
      diag = static_cast<NeuDiagnosticManager *>
        (windEmulatorStep4_WECSim_DW.OUTPUT_1_1_DiagMgr);
      diagTree = neu_diagnostic_manager_get_initial_tree(diag);
      i = ne_simulator_method(simulator, NESL_SIM_OUTPUTS, simulationData, diag);
      if (i != 0) {
        tmp_11 = error_buffer_is_empty(rtmGetErrorStatus
          (windEmulatorStep4_WECSim_M));
        if (tmp_11) {
          msg = rtw_diagnostics_msg(diagTree);
          rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
        }
      }

      for (i = 0; i < 6; i++) {
        /* Assignment: '<S59>/Assignment1' incorporates:
         *  Constant: '<S59>/Constant1'
         */
        windEmulatorStep4_WECSim_B.velocity[i] =
          windEmulatorStep4_WECSim_cal->Constant1_Value_b2[i];
      }

      /* Assignment: '<S59>/Assignment1' */
      windEmulatorStep4_WECSim_B.velocity[4] =
        windEmulatorStep4_WECSim_B.OUTPUT_1_1[1];
      if (tmp_g) {
      }

      /* Switch: '<S58>/Switch' */
      if (windEmulatorStep4_WECSim_B.velocity[4] >
          windEmulatorStep4_WECSim_cal->Switch_Threshold_l) {
        /* Switch: '<S58>/Switch' incorporates:
         *  Constant: '<S58>/Constant3'
         */
        windEmulatorStep4_WECSim_B.Switch_l =
          windEmulatorStep4_WECSim_cal->Constant3_Value_h;
      } else {
        /* Switch: '<S58>/Switch' incorporates:
         *  Constant: '<S58>/Constant4'
         */
        windEmulatorStep4_WECSim_B.Switch_l =
          windEmulatorStep4_WECSim_cal->Constant4_Value;
      }

      /* End of Switch: '<S58>/Switch' */

      /* Step: '<S58>/Step' */
      if (deltaT_tmp < windEmulatorStep4_WECSim_cal->Step_Time) {
        /* Step: '<S58>/Step' */
        windEmulatorStep4_WECSim_B.Step_g =
          windEmulatorStep4_WECSim_cal->Step_Y0_c;
      } else {
        /* Step: '<S58>/Step' */
        windEmulatorStep4_WECSim_B.Step_g =
          windEmulatorStep4_WECSim_cal->Step_YFinal;
      }

      /* SimscapeInputBlock: '<S332>/INPUT_2_1_1' */
      windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[0] =
        windEmulatorStep4_WECSim_B.Step_g;
      windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[1] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[2] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[3] = 0.0;

      /* SimscapeInputBlock: '<S332>/INPUT_3_1_1' */
      if (windEmulatorStep4_WECSim_DW.INPUT_3_1_1_FirstOutput_4203252 == 0.0) {
        windEmulatorStep4_WECSim_DW.INPUT_3_1_1_FirstOutput_4203252 = 1.0;
        windEmulatorStep4_WECSim_X.windEmulatorStep4_WECSimhptoS_n =
          windEmulatorStep4_WECSim_B.velocity[4];
      }

      windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[0] =
        windEmulatorStep4_WECSim_X.windEmulatorStep4_WECSimhptoS_n;
      windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[1] =
        (windEmulatorStep4_WECSim_B.velocity[4] -
         windEmulatorStep4_WECSim_X.windEmulatorStep4_WECSimhptoS_n) * 1000.0;
      windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[2] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[3] = 0.0;

      /* End of SimscapeInputBlock: '<S332>/INPUT_3_1_1' */
      if (tmp_g) {
        /* Delay: '<S58>/Delay One Step' */
        windEmulatorStep4_WECSim_B.DelayOneStep =
          windEmulatorStep4_WECSim_DW.DelayOneStep_DSTATE;

        /* Sum: '<S58>/Sum1' incorporates:
         *  Constant: '<S58>/Constant'
         */
        windEmulatorStep4_WECSim_B.Sum1 = *get_genShaftSpeedRef() -
          windEmulatorStep4_WECSim_B.DelayOneStep;

        /* Gain: '<S287>/Proportional Gain' */
        windEmulatorStep4_WECSim_B.ProportionalGain = *get_pGainGen() *
          windEmulatorStep4_WECSim_B.Sum1;

        /* DiscreteIntegrator: '<S282>/Integrator' */
        windEmulatorStep4_WECSim_B.Integrator_f =
          windEmulatorStep4_WECSim_DW.Integrator_DSTATE;

        /* Gain: '<S273>/Derivative Gain' */
        windEmulatorStep4_WECSim_B.DerivativeGain =
          windEmulatorStep4_WECSim_cal->DiscretePIDController_D *
          windEmulatorStep4_WECSim_B.Sum1;

        /* SampleTimeMath: '<S277>/Tsamp'
         *
         * About '<S277>/Tsamp':
         *  y = u * K where K = 1 / ( w * Ts )
         *   */
        windEmulatorStep4_WECSim_B.Tsamp =
          windEmulatorStep4_WECSim_B.DerivativeGain *
          windEmulatorStep4_WECSim_cal->Tsamp_WtEt;

        /* Delay: '<S275>/UD' */
        windEmulatorStep4_WECSim_B.UD = windEmulatorStep4_WECSim_DW.UD_DSTATE_j;

        /* Sum: '<S275>/Diff' */
        windEmulatorStep4_WECSim_B.Diff_o = windEmulatorStep4_WECSim_B.Tsamp -
          windEmulatorStep4_WECSim_B.UD;

        /* Sum: '<S291>/Sum' */
        windEmulatorStep4_WECSim_B.Sum_g =
          (windEmulatorStep4_WECSim_B.ProportionalGain +
           windEmulatorStep4_WECSim_B.Integrator_f) +
          windEmulatorStep4_WECSim_B.Diff_o;
      }

      /* SimscapeInputBlock: '<S332>/INPUT_1_1_1' */
      windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[0] =
        windEmulatorStep4_WECSim_B.Sum_g;
      windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[1] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[2] = 0.0;
      if (f) {
        windEmulatorStep4_WECSim_DW.INPUT_1_1_1_Discrete_2152258201[0] =
          !(windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[0] ==
            windEmulatorStep4_WECSim_DW.INPUT_1_1_1_Discrete_2152258201[1]);
        windEmulatorStep4_WECSim_DW.INPUT_1_1_1_Discrete_2152258201[1] =
          windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[0];
      }

      windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[0] =
        windEmulatorStep4_WECSim_DW.INPUT_1_1_1_Discrete_2152258201[1];
      windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[3] =
        windEmulatorStep4_WECSim_DW.INPUT_1_1_1_Discrete_2152258201[0];

      /* SimscapeExecutionBlock: '<S332>/STATE_1' */
      simulationData = static_cast<NeslSimulationData *>
        (windEmulatorStep4_WECSim_DW.STATE_1_SimData_a);
      time_7 = u1;
      simulationData->mData->mTime.mN = 1;
      simulationData->mData->mTime.mX = &time_7;
      simulationData->mData->mContStates.mN = 35;
      simulationData->mData->mContStates.mX =
        &windEmulatorStep4_WECSim_X.windEmulatorStep4_WECSimhptoS_h[0];
      simulationData->mData->mDiscStates.mN = 6;
      simulationData->mData->mDiscStates.mX =
        &windEmulatorStep4_WECSim_DW.STATE_1_Discrete_208214823[0];
      simulationData->mData->mModeVector.mN = 21;
      simulationData->mData->mModeVector.mX =
        &windEmulatorStep4_WECSim_DW.STATE_1_Modes_i[0];
      tmp_11 = false;
      simulationData->mData->mFoundZcEvents = tmp_11;
      simulationData->mData->mHadEvents = false;
      simulationData->mData->mIsMajorTimeStep = f;
      tmp_11 = false;
      simulationData->mData->mIsSolverAssertCheck = tmp_11;
      simulationData->mData->mIsSolverCheckingCIC = false;
      tmp_11 = rtsiIsSolverComputingJacobian
        (&windEmulatorStep4_WECSim_M->solverInfo);
      simulationData->mData->mIsComputingJacobian = tmp_11;
      simulationData->mData->mIsEvaluatingF0 = false;
      simulationData->mData->mIsSolverRequestingReset = false;
      simulationData->mData->mIsModeUpdateTimeStep = tmp_h;
      tmp_9[0] = 0;
      tmp_8[0] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[0];
      tmp_8[1] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[1];
      tmp_8[2] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[2];
      tmp_8[3] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[3];
      tmp_9[1] = 4;
      tmp_8[4] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[0];
      tmp_8[5] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[1];
      tmp_8[6] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[2];
      tmp_8[7] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[3];
      tmp_9[2] = 8;
      tmp_8[8] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[0];
      tmp_8[9] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[1];
      tmp_8[10] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[2];
      tmp_8[11] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[3];
      tmp_9[3] = 12;
      simulationData->mData->mInputValues.mN = 12;
      simulationData->mData->mInputValues.mX = &tmp_8[0];
      simulationData->mData->mInputOffsets.mN = 4;
      simulationData->mData->mInputOffsets.mX = &tmp_9[0];
      simulationData->mData->mOutputs.mN = 62;
      simulationData->mData->mOutputs.mX =
        &windEmulatorStep4_WECSim_B.STATE_1_d[0];
      simulationData->mData->mTolerances.mN = 0;
      simulationData->mData->mTolerances.mX = NULL;
      simulationData->mData->mCstateHasChanged = false;
      simulationData->mData->mDstateHasChanged = false;
      time_8 = deltaT_tmp;
      simulationData->mData->mTime.mN = 1;
      simulationData->mData->mTime.mX = &time_8;
      simulationData->mData->mSampleHits.mN = 0;
      simulationData->mData->mSampleHits.mX = NULL;
      simulationData->mData->mIsFundamentalSampleHit = false;
      simulationData->mData->mHadEvents = false;
      simulator = static_cast<NeslSimulator *>
        (windEmulatorStep4_WECSim_DW.STATE_1_Simulator_i);
      diag = static_cast<NeuDiagnosticManager *>
        (windEmulatorStep4_WECSim_DW.STATE_1_DiagMgr_g);
      diagTree = neu_diagnostic_manager_get_initial_tree(diag);
      i = ne_simulator_method(simulator, NESL_SIM_OUTPUTS, simulationData, diag);
      if (i != 0) {
        tmp_11 = error_buffer_is_empty(rtmGetErrorStatus
          (windEmulatorStep4_WECSim_M));
        if (tmp_11) {
          msg = rtw_diagnostics_msg(diagTree);
          rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
        }
      }

      /* SimscapeExecutionBlock: '<S332>/OUTPUT_1_0' */
      simulationData = static_cast<NeslSimulationData *>
        (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_SimData_l);
      time_9 = u1;
      simulationData->mData->mTime.mN = 1;
      simulationData->mData->mTime.mX = &time_9;
      simulationData->mData->mContStates.mN = 0;
      simulationData->mData->mContStates.mX = NULL;
      simulationData->mData->mDiscStates.mN = 0;
      simulationData->mData->mDiscStates.mX =
        &windEmulatorStep4_WECSim_DW.OUTPUT_1_0_Discrete_k;
      simulationData->mData->mModeVector.mN = 0;
      simulationData->mData->mModeVector.mX =
        &windEmulatorStep4_WECSim_DW.OUTPUT_1_0_Modes_b;
      tmp_11 = false;
      simulationData->mData->mFoundZcEvents = tmp_11;
      simulationData->mData->mHadEvents = false;
      simulationData->mData->mIsMajorTimeStep = f;
      tmp_11 = false;
      simulationData->mData->mIsSolverAssertCheck = tmp_11;
      simulationData->mData->mIsSolverCheckingCIC = false;
      simulationData->mData->mIsComputingJacobian = false;
      simulationData->mData->mIsEvaluatingF0 = false;
      simulationData->mData->mIsSolverRequestingReset = false;
      simulationData->mData->mIsModeUpdateTimeStep = tmp_h;
      tmp_b[0] = 0;
      tmp_a[0] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[0];
      tmp_a[1] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[1];
      tmp_a[2] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[2];
      tmp_a[3] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[3];
      tmp_b[1] = 4;
      tmp_a[4] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[0];
      tmp_a[5] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[1];
      tmp_a[6] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[2];
      tmp_a[7] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[3];
      tmp_b[2] = 8;
      tmp_a[8] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[0];
      tmp_a[9] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[1];
      tmp_a[10] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[2];
      tmp_a[11] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[3];
      tmp_b[3] = 12;
      std::memcpy(&tmp_a[12], &windEmulatorStep4_WECSim_B.STATE_1_d[0], 62U *
                  sizeof(real_T));
      tmp_b[4] = 74;
      simulationData->mData->mInputValues.mN = 74;
      simulationData->mData->mInputValues.mX = &tmp_a[0];
      simulationData->mData->mInputOffsets.mN = 5;
      simulationData->mData->mInputOffsets.mX = &tmp_b[0];
      simulationData->mData->mOutputs.mN = 24;
      simulationData->mData->mOutputs.mX =
        &windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[0];
      simulationData->mData->mTolerances.mN = 0;
      simulationData->mData->mTolerances.mX = NULL;
      simulationData->mData->mCstateHasChanged = false;
      simulationData->mData->mDstateHasChanged = false;
      time_a = deltaT_tmp;
      simulationData->mData->mTime.mN = 1;
      simulationData->mData->mTime.mX = &time_a;
      simulationData->mData->mSampleHits.mN = 0;
      simulationData->mData->mSampleHits.mX = NULL;
      simulationData->mData->mIsFundamentalSampleHit = false;
      simulationData->mData->mHadEvents = false;
      simulator = static_cast<NeslSimulator *>
        (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_Simulator_d);
      diag = static_cast<NeuDiagnosticManager *>
        (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_DiagMgr_i);
      diagTree = neu_diagnostic_manager_get_initial_tree(diag);
      i = ne_simulator_method(simulator, NESL_SIM_OUTPUTS, simulationData, diag);
      if (i != 0) {
        tmp_11 = error_buffer_is_empty(rtmGetErrorStatus
          (windEmulatorStep4_WECSim_M));
        if (tmp_11) {
          msg = rtw_diagnostics_msg(diagTree);
          rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
        }
      }

      /* Abs: '<S58>/Abs' */
      windEmulatorStep4_WECSim_B.Abs_p = std::abs
        (windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[15]);

      /* Gain: '<S58>/Gain' */
      windEmulatorStep4_WECSim_B.Gain_f =
        windEmulatorStep4_WECSim_cal->Gain_Gain_b *
        windEmulatorStep4_WECSim_B.Abs_p;

      /* Product: '<S58>/Product5' */
      windEmulatorStep4_WECSim_B.ptoTorqueHydraulic =
        windEmulatorStep4_WECSim_B.Switch_l * windEmulatorStep4_WECSim_B.Gain_f;
      if (tmp_g) {
      }

      /* Product: '<S58>/Product4' */
      windEmulatorStep4_WECSim_B.ptoPowerMech =
        windEmulatorStep4_WECSim_B.velocity[4] *
        windEmulatorStep4_WECSim_B.ptoTorqueHydraulic;
      if (tmp_g) {
      }

      /* Switch: '<S58>/Switch1' */
      if (windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[18] >
          windEmulatorStep4_WECSim_cal->Switch1_Threshold_i) {
        /* Switch: '<S58>/Switch1' incorporates:
         *  Constant: '<S58>/Constant2'
         */
        windEmulatorStep4_WECSim_B.Switch1_k =
          windEmulatorStep4_WECSim_cal->Constant2_Value_i;
      } else {
        /* Switch: '<S58>/Switch1' incorporates:
         *  Constant: '<S58>/Constant5'
         */
        windEmulatorStep4_WECSim_B.Switch1_k =
          windEmulatorStep4_WECSim_cal->Constant5_Value;
      }

      /* End of Switch: '<S58>/Switch1' */

      /* Abs: '<S58>/Abs1' */
      windEmulatorStep4_WECSim_B.Abs1 = std::abs
        (windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[12]);

      /* Product: '<S58>/Product1' */
      windEmulatorStep4_WECSim_B.pistonPowerMech =
        windEmulatorStep4_WECSim_B.Switch1_k * windEmulatorStep4_WECSim_B.Abs1;
      if (tmp_g) {
      }

      /* Product: '<S58>/Product' incorporates:
       *  Constant: '<S58>/Constant1'
       */
      windEmulatorStep4_WECSim_B.shaftSpeed =
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[14] * *get_radsec2rpm();
      if (tmp_g) {
      }

      /* Product: '<S58>/Product2' */
      windEmulatorStep4_WECSim_B.shaftPower = windEmulatorStep4_WECSim_B.Sum_g *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[14];
      if (tmp_g) {
      }

      /* Product: '<S58>/Product3' */
      windEmulatorStep4_WECSim_B.powerHM =
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[14] *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[16];
      if (tmp_g) {
        /* Gain: '<S279>/Integral Gain' */
        windEmulatorStep4_WECSim_B.IntegralGain = *get_iGainGen() *
          windEmulatorStep4_WECSim_B.Sum1;
      }

      /* Gain: '<S219>/Gain' */
      windEmulatorStep4_WECSim_B.flowRateHCB = tmp_z *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[0];
      if (tmp_g) {
      }

      /* Gain: '<S220>/Gain' */
      windEmulatorStep4_WECSim_B.flowRateCV4 = tmp_z *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[1];
      if (tmp_g) {
      }

      /* Gain: '<S221>/Gain' */
      windEmulatorStep4_WECSim_B.flowRateCV1 = tmp_z *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[2];
      if (tmp_g) {
      }

      /* Gain: '<S222>/Gain' */
      windEmulatorStep4_WECSim_B.flowRateCV3 = tmp_z *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[3];
      if (tmp_g) {
      }

      /* Gain: '<S223>/Gain' */
      windEmulatorStep4_WECSim_B.flowRateHCA = tmp_z *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[4];
      if (tmp_g) {
      }

      /* Gain: '<S224>/Gain' */
      windEmulatorStep4_WECSim_B.flowRateC = tmp_z *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[5];
      if (tmp_g) {
      }

      /* Gain: '<S225>/Gain' */
      windEmulatorStep4_WECSim_B.flowRateAccHP = tmp_z *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[6];
      if (tmp_g) {
      }

      /* Gain: '<S226>/Gain' */
      windEmulatorStep4_WECSim_B.flowRateHMin = tmp_z *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[7];
      if (tmp_g) {
      }

      /* Gain: '<S227>/Gain' */
      windEmulatorStep4_WECSim_B.flowRateHMout = tmp_z *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[8];
      if (tmp_g) {
      }

      /* Gain: '<S228>/Gain' */
      windEmulatorStep4_WECSim_B.flowRateAccLP = tmp_z *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[9];
      if (tmp_g) {
      }

      /* Gain: '<S229>/Gain' */
      windEmulatorStep4_WECSim_B.flowRateD = tmp_z *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[10];
      if (tmp_g) {
      }

      /* Gain: '<S230>/Gain' */
      windEmulatorStep4_WECSim_B.flowRateCV2 = tmp_z *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[11];
      if (tmp_g) {
      }

      /* Gain: '<S241>/Gain' */
      windEmulatorStep4_WECSim_B.pressureA = tmp_y *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[19];
      if (tmp_g) {
      }

      /* Gain: '<S242>/Gain' */
      windEmulatorStep4_WECSim_B.pressureB = tmp_y *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[20];
      if (tmp_g) {
      }

      /* Gain: '<S243>/Gain' */
      windEmulatorStep4_WECSim_B.pressureC = tmp_y *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[21];
      if (tmp_g) {
      }

      /* Gain: '<S244>/Gain' */
      windEmulatorStep4_WECSim_B.pressureHM = tmp_y *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[22];
      if (tmp_g) {
      }

      /* Gain: '<S245>/Gain' */
      windEmulatorStep4_WECSim_B.pressureD = tmp_y *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_m[23];
      if (tmp_g) {
      }

      for (i = 0; i < 6; i++) {
        /* Assignment: '<S59>/Assignment ' incorporates:
         *  Constant: '<S59>/Constant'
         */
        windEmulatorStep4_WECSim_B.position[i] =
          windEmulatorStep4_WECSim_cal->Constant_Value_c[i];
      }

      /* Assignment: '<S59>/Assignment ' */
      windEmulatorStep4_WECSim_B.position[4] =
        windEmulatorStep4_WECSim_B.OUTPUT_1_1[0];

      /* SimscapeInputBlock: '<S216>/INPUT_5_1_1' */
      windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[0] =
        windEmulatorStep4_WECSim_B.ptoTorqueHydraulic;
      windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[1] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[2] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[3] = 0.0;

      /* MATLAB Function: '<S81>/quaternion2EulXYZ' incorporates:
       *  SimscapeExecutionBlock: '<S216>/OUTPUT_1_1'
       */
      windEmulatorS_quaternion2EulXYZ(&windEmulatorStep4_WECSim_B.OUTPUT_1_1[2],
        &windEmulatorStep4_WECSim_B.sf_quaternion2EulXYZ,
        &windEmulatorStep4_WECSim_DW.sf_quaternion2EulXYZ);

      /* TransportDelay: '<S60>/Transport Delay' */
      for (i = 0; i < 6; i++) {
        riseValLimit = rt_TDelayInterpolate(windEmulatorStep4_WECSim_M->
          Timing.t[0] - windEmulatorStep4_WECSim_cal->TransportDelay_Delay,
          windEmulatorStep4_WECSim_DW.TransportDelay_RWORK[0],static_cast<real_T
          *>(windEmulatorStep4_WECSim_DW.TransportDelay_PWORK[i]),
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK[i + 18],
          &windEmulatorStep4_WECSim_DW.TransportDelay_IWORK[i + 12],
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK[i],
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK[i + 6],
          windEmulatorStep4_WECSim_cal->TransportDelay_InitOutput,false,
          rtmIsMinorTimeStep(windEmulatorStep4_WECSim_M) && ((static_cast<real_T
          *>(windEmulatorStep4_WECSim_DW.TransportDelay_PWORK[0]))
          [windEmulatorStep4_WECSim_DW.TransportDelay_IWORK[6] +
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK[18]] ==
          windEmulatorStep4_WECSim_M->Timing.t[0]));
        windEmulatorStep4_WECSim_B.TransportDelay[i] = riseValLimit;
      }

      /* End of TransportDelay: '<S60>/Transport Delay' */

      /* MATLAB Function: '<S133>/Yaw Kinematic Transforms' incorporates:
       *  Constant: '<S133>/Constant'
       *  SimscapeExecutionBlock: '<S216>/OUTPUT_1_1'
       */
      windEmul_YawKinematicTransforms
        (windEmulatorStep4_WECSim_cal->Constant_Value_h1,
         &windEmulatorStep4_WECSim_B.OUTPUT_1_1[9],
         windEmulatorStep4_WECSim_B.sf_quaternion2EulXYZ.E,
         &windEmulatorStep4_WECSim_B.OUTPUT_1_1[12],
         &windEmulatorStep4_WECSim_B.OUTPUT_1_1[6],
         windEmulatorStep4_WECSim_B.TransportDelay,
         &windEmulatorStep4_WECSim_B.sf_YawKinematicTransforms,
         &windEmulatorStep4_WECSim_DW.sf_YawKinematicTransforms);

      /* Product: '<S131>/Product' incorporates:
       *  Constant: '<S67>/Constant'
       *  Product: '<S140>/Product1'
       *  Product: '<S143>/Product1'
       *  Product: '<S148>/Product1'
       *  Product: '<S153>/Product2'
       *  Product: '<S210>/Product'
       *  Product: '<S61>/Product1'
       *  Product: '<S64>/Product1'
       *  Product: '<S69>/Product1'
       *  Product: '<S74>/Product2'
       */
      tmp_i = &windEmulatorStep4_WECSim_cal->Constant_Value_l.hf1.fDamping[0];
      for (i = 0; i < 6; i++) {
        tmp_e[i] = windEmulatorStep4_WECSim_B.sf_YawKinematicTransforms.velLoc[i];
        tmp_f[i] = 0.0;
      }

      for (i = 0; i < 6; i++) {
        riseValLimit = tmp_e[i];
        for (i_0 = 0; i_0 < 6; i_0++) {
          tmp_f[i_0] += tmp_i[6 * i + i_0] * riseValLimit;
        }
      }

      for (i = 0; i < 6; i++) {
        /* Product: '<S131>/Product' */
        windEmulatorStep4_WECSim_B.F_SingleFrequency[i] = tmp_f[i];
      }

      /* Product: '<S69>/Product1' incorporates:
       *  Constant: '<S67>/Constant'
       *  Product: '<S131>/Product'
       *  Product: '<S140>/Product1'
       *  Product: '<S143>/Product1'
       *  Product: '<S148>/Product1'
       *  Product: '<S153>/Product2'
       *  Product: '<S210>/Product'
       *  Product: '<S61>/Product1'
       *  Product: '<S64>/Product1'
       *  Product: '<S74>/Product2'
       */
      tmp_i = &windEmulatorStep4_WECSim_cal->Constant_Value_l.hf1.fAddedMass[0];
      for (i = 0; i < 6; i++) {
        tmp_e[i] = windEmulatorStep4_WECSim_B.sf_YawKinematicTransforms.accLoc[i];
        tmp_f[i] = 0.0;
      }

      for (i = 0; i < 6; i++) {
        riseValLimit = tmp_e[i];
        for (i_0 = 0; i_0 < 6; i_0++) {
          tmp_f[i_0] += tmp_i[6 * i + i_0] * riseValLimit;
        }
      }

      for (i = 0; i < 6; i++) {
        /* Product: '<S69>/Product1' */
        windEmulatorStep4_WECSim_B.F_AddedMass[i] = tmp_f[i];
      }

      /* Clock: '<S57>/Clock' */
      windEmulatorStep4_WECSim_B.Clock_n = Clock_tmp;
      if (tmp_g) {
        /* Product: '<S124>/Divide' incorporates:
         *  Constant: '<S124>/Constant4'
         *  Constant: '<S124>/Ramp Time'
         */
        windEmulatorStep4_WECSim_B.frequency =
          windEmulatorStep4_WECSim_cal->Constant4_Value_d /
          windEmulatorStep4_WECSim_cal->RampTime_Value;
      }

      /* RelationalOperator: '<S68>/Less Than' incorporates:
       *  Constant: '<S68>/Ramp Function Time'
       */
      windEmulatorStep4_WECSim_B.LessThan = (windEmulatorStep4_WECSim_B.Clock_n <
        windEmulatorStep4_WECSim_cal->RampFunctionTime_Value);

      /* Switch: '<S68>/Switch' */
      if (windEmulatorStep4_WECSim_B.LessThan) {
        /* Product: '<S124>/Product2' */
        windEmulatorStep4_WECSim_B.Product2_cu =
          windEmulatorStep4_WECSim_B.Clock_n *
          windEmulatorStep4_WECSim_B.frequency;

        /* Sum: '<S124>/Add2' incorporates:
         *  Constant: '<S124>/Constant3'
         */
        windEmulatorStep4_WECSim_B.Add2_e =
          windEmulatorStep4_WECSim_B.Product2_cu +
          windEmulatorStep4_WECSim_cal->Constant3_Value_d;

        /* Sin: '<S124>/Sine Wave Function' */
        windEmulatorStep4_WECSim_B.SineWaveFunction_g = std::sin
          (windEmulatorStep4_WECSim_cal->SineWaveFunction_Freq *
           windEmulatorStep4_WECSim_B.Add2_e +
           windEmulatorStep4_WECSim_cal->SineWaveFunction_Phase) *
          windEmulatorStep4_WECSim_cal->SineWaveFunction_Amp +
          windEmulatorStep4_WECSim_cal->SineWaveFunction_Bias;

        /* Sum: '<S124>/Add1' incorporates:
         *  Constant: '<S124>/Constant2'
         */
        windEmulatorStep4_WECSim_B.Add1_h =
          windEmulatorStep4_WECSim_cal->Constant2_Value_d +
          windEmulatorStep4_WECSim_B.SineWaveFunction_g;

        /* Product: '<S124>/Product1' incorporates:
         *  Constant: '<S124>/Constant1'
         */
        windEmulatorStep4_WECSim_B.Product1_k =
          windEmulatorStep4_WECSim_B.Add1_h *
          windEmulatorStep4_WECSim_cal->Constant1_Value_g;

        /* Switch: '<S68>/Switch' */
        windEmulatorStep4_WECSim_B.R = windEmulatorStep4_WECSim_B.Product1_k;
      } else {
        /* Switch: '<S68>/Switch' incorporates:
         *  Constant: '<S68>/Constant'
         */
        windEmulatorStep4_WECSim_B.R =
          windEmulatorStep4_WECSim_cal->Constant_Value_d;
      }

      /* End of Switch: '<S68>/Switch' */
      if (tmp_g) {
        for (i = 0; i < 6; i++) {
          /* Product: '<S126>/Product3' incorporates:
           *  Constant: '<S126>/Wave Amplitude1'
           *  Constant: '<S67>/Constant'
           */
          windEmulatorStep4_WECSim_B.Product3[i] =
            windEmulatorStep4_WECSim_cal->WaveAmplitude1_Value *
            windEmulatorStep4_WECSim_cal->Constant_Value_l.hf1.fExt.md[i];
        }
      }

      /* Product: '<S126>/Product' incorporates:
       *  Constant: '<S126>/Wave Frequency'
       */
      windEmulatorStep4_WECSim_B.Product_p = windEmulatorStep4_WECSim_B.Clock_n *
        windEmulatorStep4_WECSim_cal->WaveFrequency_Value;

      /* Sum: '<S126>/Add3' incorporates:
       *  Constant: '<S126>/Center of Gravity'
       *  Constant: '<S126>/Constant1'
       */
      windEmulatorStep4_WECSim_B.x_cg[0] =
        windEmulatorStep4_WECSim_B.OUTPUT_1_1[9] -
        windEmulatorStep4_WECSim_cal->CenterofGravity_Value[0];
      windEmulatorStep4_WECSim_B.x_cg[3] =
        windEmulatorStep4_WECSim_B.sf_quaternion2EulXYZ.E[0] -
        windEmulatorStep4_WECSim_cal->Constant1_Value_d[0];
      windEmulatorStep4_WECSim_B.x_cg[1] =
        windEmulatorStep4_WECSim_B.OUTPUT_1_1[10] -
        windEmulatorStep4_WECSim_cal->CenterofGravity_Value[1];
      windEmulatorStep4_WECSim_B.x_cg[4] =
        windEmulatorStep4_WECSim_B.sf_quaternion2EulXYZ.E[1] -
        windEmulatorStep4_WECSim_cal->Constant1_Value_d[1];
      windEmulatorStep4_WECSim_B.x_cg[2] =
        windEmulatorStep4_WECSim_B.OUTPUT_1_1[11] -
        windEmulatorStep4_WECSim_cal->CenterofGravity_Value[2];
      windEmulatorStep4_WECSim_B.x_cg[5] =
        windEmulatorStep4_WECSim_B.sf_quaternion2EulXYZ.E[2] -
        windEmulatorStep4_WECSim_cal->Constant1_Value_d[2];

      /* MATLAB Function: '<S126>/MATLAB Function1' incorporates:
       *  Constant: '<S126>/Displacement Phase Enable1'
       *  Constant: '<S126>/Wave Number'
       *  Constant: '<S126>/Wave direction'
       *  Sum: '<S126>/Add3'
       */
      windEmulatorSte_MATLABFunction1(&windEmulatorStep4_WECSim_B.x_cg[0],
        windEmulatorStep4_WECSim_cal->DisplacementPhaseEnable1_Value,
        windEmulatorStep4_WECSim_cal->WaveNumber_Value,
        windEmulatorStep4_WECSim_cal->Wavedirection_Value,
        &windEmulatorStep4_WECSim_B.sf_MATLABFunction1,
        &windEmulatorStep4_WECSim_DW.sf_MATLABFunction1);

      /* Sum: '<S126>/Add' incorporates:
       *  Constant: '<S126>/Constant'
       */
      windEmulatorStep4_WECSim_B.Add = (windEmulatorStep4_WECSim_B.Product_p +
        windEmulatorStep4_WECSim_cal->Constant_Value_n) +
        windEmulatorStep4_WECSim_B.sf_MATLABFunction1.dispPhase;

      /* Sin: '<S126>/Sine Wave Function1' */
      windEmulatorStep4_WECSim_B.coswt = std::sin
        (windEmulatorStep4_WECSim_cal->SineWaveFunction1_Freq *
         windEmulatorStep4_WECSim_B.Add +
         windEmulatorStep4_WECSim_cal->SineWaveFunction1_Phase) *
        windEmulatorStep4_WECSim_cal->SineWaveFunction1_Amp +
        windEmulatorStep4_WECSim_cal->SineWaveFunction1_Bias;

      /* Sum: '<S126>/Add2' */
      windEmulatorStep4_WECSim_B.Add2_i = windEmulatorStep4_WECSim_B.Product_p +
        windEmulatorStep4_WECSim_B.sf_MATLABFunction1.dispPhase;

      /* Sin: '<S126>/Sine Wave Function' */
      windEmulatorStep4_WECSim_B.sinwt = std::sin
        (windEmulatorStep4_WECSim_cal->SineWaveFunction_Freq_g *
         windEmulatorStep4_WECSim_B.Add2_i +
         windEmulatorStep4_WECSim_cal->SineWaveFunction_Phase_p) *
        windEmulatorStep4_WECSim_cal->SineWaveFunction_Amp_j +
        windEmulatorStep4_WECSim_cal->SineWaveFunction_Bias_l;
      for (i = 0; i < 6; i++) {
        /* Product: '<S126>/Product1' incorporates:
         *  Constant: '<S126>/Wave Amplitude'
         *  Constant: '<S67>/Constant'
         */
        Clock_tmp = windEmulatorStep4_WECSim_cal->WaveAmplitude_Value *
          windEmulatorStep4_WECSim_cal->Constant_Value_l.hf1.fExt.re[i] *
          windEmulatorStep4_WECSim_B.coswt;
        windEmulatorStep4_WECSim_B.Product1_b[i] = Clock_tmp;

        /* Product: '<S126>/Product2' incorporates:
         *  Constant: '<S126>/Wave Amplitude'
         *  Constant: '<S67>/Constant'
         */
        riseValLimit = windEmulatorStep4_WECSim_cal->WaveAmplitude_Value *
          windEmulatorStep4_WECSim_B.sinwt *
          windEmulatorStep4_WECSim_cal->Constant_Value_l.hf1.fExt.im[i];
        windEmulatorStep4_WECSim_B.Product2_o[i] = riseValLimit;

        /* Sum: '<S126>/Add1' incorporates:
         *  Product: '<S126>/Product1'
         *  Product: '<S126>/Product2'
         *  Product: '<S126>/Product3'
         */
        Clock_tmp = (windEmulatorStep4_WECSim_B.Product3[i] + Clock_tmp) -
          riseValLimit;
        windEmulatorStep4_WECSim_B.Add1[i] = Clock_tmp;

        /* Sum: '<S68>/Add' incorporates:
         *  Constant: '<S128>/Constant'
         *  Constant: '<S129>/Constant'
         *  Sum: '<S126>/Add1'
         */
        Clock_tmp = (Clock_tmp + windEmulatorStep4_WECSim_cal->
                     Constant_Value_a[i]) +
          windEmulatorStep4_WECSim_cal->Constant_Value_g[i];
        windEmulatorStep4_WECSim_B.F_wave[i] = Clock_tmp;

        /* Product: '<S68>/Product' incorporates:
         *  Sum: '<S68>/Add'
         */
        windEmulatorStep4_WECSim_B.F_Excitation[i] =
          windEmulatorStep4_WECSim_B.R * Clock_tmp;
      }

      if (tmp_g) {
        /* Product: '<S75>/Product' incorporates:
         *  Constant: '<S67>/Constant'
         *  Constant: '<S75>/Gravity'
         */
        windEmulatorStep4_WECSim_B.F_Gravity =
          windEmulatorStep4_WECSim_cal->Gravity_Value *
          windEmulatorStep4_WECSim_cal->Constant_Value_l.hf1.mass;

        /* Product: '<S75>/Product1' incorporates:
         *  Constant: '<S67>/Constant'
         *  Constant: '<S75>/Gravity'
         *  Constant: '<S75>/Water Density'
         */
        windEmulatorStep4_WECSim_B.F_Buoyancy =
          windEmulatorStep4_WECSim_cal->WaterDensity_Value *
          windEmulatorStep4_WECSim_cal->Gravity_Value *
          windEmulatorStep4_WECSim_cal->Constant_Value_l.hf1.volume;

        /* Sum: '<S75>/Add1' */
        windEmulatorStep4_WECSim_B.Add1_a = windEmulatorStep4_WECSim_B.F_Gravity
          - windEmulatorStep4_WECSim_B.F_Buoyancy;

        /* Assignment: '<S75>/Assignment (Add Net Bouyancy Force  to Z-Direction)1' incorporates:
         *  Constant: '<S75>/Constant2'
         */
        windEmulatorStep4_WECSim_B.AssignmentAddNetBouyancyForceto[0] =
          windEmulatorStep4_WECSim_cal->Constant2_Value_c[0];

        /* Sum: '<S75>/Add3' incorporates:
         *  Constant: '<S67>/Constant'
         *  Constant: '<S75>/Center of Gravity'
         */
        windEmulatorStep4_WECSim_B.Add3[0] =
          windEmulatorStep4_WECSim_cal->Constant_Value_l.hf1.centerBuoyancy[0] -
          windEmulatorStep4_WECSim_cal->CenterofGravity_Value_d[0];

        /* Assignment: '<S75>/Assignment (Add Net Bouyancy Force  to Z-Direction)1' incorporates:
         *  Constant: '<S75>/Constant2'
         */
        windEmulatorStep4_WECSim_B.AssignmentAddNetBouyancyForceto[1] =
          windEmulatorStep4_WECSim_cal->Constant2_Value_c[1];

        /* Sum: '<S75>/Add3' incorporates:
         *  Constant: '<S67>/Constant'
         *  Constant: '<S75>/Center of Gravity'
         */
        windEmulatorStep4_WECSim_B.Add3[1] =
          windEmulatorStep4_WECSim_cal->Constant_Value_l.hf1.centerBuoyancy[1] -
          windEmulatorStep4_WECSim_cal->CenterofGravity_Value_d[1];
        windEmulatorStep4_WECSim_B.Add3[2] =
          windEmulatorStep4_WECSim_cal->Constant_Value_l.hf1.centerBuoyancy[2] -
          windEmulatorStep4_WECSim_cal->CenterofGravity_Value_d[2];

        /* Assignment: '<S75>/Assignment (Add Net Bouyancy Force  to Z-Direction)1' */
        windEmulatorStep4_WECSim_B.AssignmentAddNetBouyancyForceto[2] =
          windEmulatorStep4_WECSim_B.F_Buoyancy;

        /* Product: '<S76>/Element product' */
        windEmulatorStep4_WECSim_B.Elementproduct[0] =
          windEmulatorStep4_WECSim_B.AssignmentAddNetBouyancyForceto[1] *
          windEmulatorStep4_WECSim_B.Add3[2];
        windEmulatorStep4_WECSim_B.Elementproduct[1] =
          windEmulatorStep4_WECSim_B.Add3[0] *
          windEmulatorStep4_WECSim_B.AssignmentAddNetBouyancyForceto[2];
        windEmulatorStep4_WECSim_B.Elementproduct[2] =
          windEmulatorStep4_WECSim_B.AssignmentAddNetBouyancyForceto[0] *
          windEmulatorStep4_WECSim_B.Add3[1];
        windEmulatorStep4_WECSim_B.Elementproduct[3] =
          windEmulatorStep4_WECSim_B.Add3[1] *
          windEmulatorStep4_WECSim_B.AssignmentAddNetBouyancyForceto[2];
        windEmulatorStep4_WECSim_B.Elementproduct[4] =
          windEmulatorStep4_WECSim_B.AssignmentAddNetBouyancyForceto[0] *
          windEmulatorStep4_WECSim_B.Add3[2];
        windEmulatorStep4_WECSim_B.Elementproduct[5] =
          windEmulatorStep4_WECSim_B.Add3[0] *
          windEmulatorStep4_WECSim_B.AssignmentAddNetBouyancyForceto[1];

        /* Sum: '<S76>/Add3' */
        windEmulatorStep4_WECSim_B.Add3_d[0] =
          windEmulatorStep4_WECSim_B.Elementproduct[0] -
          windEmulatorStep4_WECSim_B.Elementproduct[3];
        windEmulatorStep4_WECSim_B.Add3_d[1] =
          windEmulatorStep4_WECSim_B.Elementproduct[1] -
          windEmulatorStep4_WECSim_B.Elementproduct[4];
        windEmulatorStep4_WECSim_B.Add3_d[2] =
          windEmulatorStep4_WECSim_B.Elementproduct[2] -
          windEmulatorStep4_WECSim_B.Elementproduct[5];
        for (i = 0; i < 6; i++) {
          /* Assignment: '<S75>/Assignment (Add Net Bouyancy Force  to Z-Direction)' incorporates:
           *  Constant: '<S75>/Constant1'
           */
          windEmulatorStep4_WECSim_B.VerticalBuoyancyForce[i] =
            windEmulatorStep4_WECSim_cal->Constant1_Value_dp[i];

          /* Assignment: '<S75>/Assignment (Add Net Bouyancy Force  to Z-Direction)2' incorporates:
           *  Constant: '<S75>/Constant1'
           */
          windEmulatorStep4_WECSim_B.Rotationalbuoyancyforce[i] =
            windEmulatorStep4_WECSim_cal->Constant1_Value_dp[i];
        }

        /* Assignment: '<S75>/Assignment (Add Net Bouyancy Force  to Z-Direction)' */
        windEmulatorStep4_WECSim_B.VerticalBuoyancyForce[2] =
          windEmulatorStep4_WECSim_B.Add1_a;

        /* Assignment: '<S75>/Assignment (Add Net Bouyancy Force  to Z-Direction)2' incorporates:
         *  Sum: '<S76>/Add3'
         */
        windEmulatorStep4_WECSim_B.Rotationalbuoyancyforce[3] =
          windEmulatorStep4_WECSim_B.Add3_d[0];
        windEmulatorStep4_WECSim_B.Rotationalbuoyancyforce[4] =
          windEmulatorStep4_WECSim_B.Add3_d[1];
        windEmulatorStep4_WECSim_B.Rotationalbuoyancyforce[5] =
          windEmulatorStep4_WECSim_B.Add3_d[2];
        for (i = 0; i < 6; i++) {
          /* Sum: '<S75>/Add2' */
          windEmulatorStep4_WECSim_B.Netbuoyancyforce[i] =
            windEmulatorStep4_WECSim_B.VerticalBuoyancyForce[i] +
            windEmulatorStep4_WECSim_B.Rotationalbuoyancyforce[i];
        }
      }

      /* Sum: '<S63>/Add' incorporates:
       *  Constant: '<S63>/Constant'
       *  Constant: '<S67>/Constant'
       */
      windEmulatorStep4_WECSim_B.x_cg_b[0] =
        windEmulatorStep4_WECSim_B.sf_YawKinematicTransforms.dispLoc[0] -
        windEmulatorStep4_WECSim_cal->Constant_Value_l.hf1.centerGravity[0];
      windEmulatorStep4_WECSim_B.x_cg_b[3] =
        windEmulatorStep4_WECSim_B.sf_YawKinematicTransforms.dispLoc[3] -
        windEmulatorStep4_WECSim_cal->Constant_Value_i[0];
      windEmulatorStep4_WECSim_B.x_cg_b[1] =
        windEmulatorStep4_WECSim_B.sf_YawKinematicTransforms.dispLoc[1] -
        windEmulatorStep4_WECSim_cal->Constant_Value_l.hf1.centerGravity[1];
      windEmulatorStep4_WECSim_B.x_cg_b[4] =
        windEmulatorStep4_WECSim_B.sf_YawKinematicTransforms.dispLoc[4] -
        windEmulatorStep4_WECSim_cal->Constant_Value_i[1];
      windEmulatorStep4_WECSim_B.x_cg_b[2] =
        windEmulatorStep4_WECSim_B.sf_YawKinematicTransforms.dispLoc[2] -
        windEmulatorStep4_WECSim_cal->Constant_Value_l.hf1.centerGravity[2];
      windEmulatorStep4_WECSim_B.x_cg_b[5] =
        windEmulatorStep4_WECSim_B.sf_YawKinematicTransforms.dispLoc[5] -
        windEmulatorStep4_WECSim_cal->Constant_Value_i[2];

      /* Product: '<S74>/Product2' incorporates:
       *  Constant: '<S67>/Constant'
       *  Product: '<S131>/Product'
       *  Product: '<S140>/Product1'
       *  Product: '<S143>/Product1'
       *  Product: '<S148>/Product1'
       *  Product: '<S153>/Product2'
       *  Product: '<S210>/Product'
       *  Product: '<S61>/Product1'
       *  Product: '<S64>/Product1'
       *  Product: '<S69>/Product1'
       *  Sum: '<S63>/Add'
       */
      tmp_i =
        &windEmulatorStep4_WECSim_cal->Constant_Value_l.hf1.linearHydroRestCoef
        [0];
      for (i = 0; i < 6; i++) {
        tmp_e[i] = windEmulatorStep4_WECSim_B.x_cg_b[i];
        windEmulatorStep4_WECSim_B.LinearRestoringForce[i] = 0.0;
      }

      for (i = 0; i < 6; i++) {
        riseValLimit = tmp_e[i];
        for (i_0 = 0; i_0 < 6; i_0++) {
          windEmulatorStep4_WECSim_B.LinearRestoringForce[i_0] += tmp_i[6 * i +
            i_0] * riseValLimit;
        }
      }

      for (i = 0; i < 6; i++) {
        /* Sum: '<S74>/Add' incorporates:
         *  Product: '<S74>/Product2'
         */
        windEmulatorStep4_WECSim_B.F_Restoring[i] =
          windEmulatorStep4_WECSim_B.Netbuoyancyforce[i] +
          windEmulatorStep4_WECSim_B.LinearRestoringForce[i];
      }

      /* MATLAB Function: '<S70>/Yaw Force Transforms' incorporates:
       *  Constant: '<S70>/Constant'
       *  SimscapeExecutionBlock: '<S216>/OUTPUT_1_1'
       */
      windEmulator_YawForceTransforms
        (windEmulatorStep4_WECSim_cal->Constant_Value_h,
         &windEmulatorStep4_WECSim_B.OUTPUT_1_1[9],
         windEmulatorStep4_WECSim_B.sf_quaternion2EulXYZ.E,
         windEmulatorStep4_WECSim_B.F_SingleFrequency,
         windEmulatorStep4_WECSim_B.F_AddedMass,
         windEmulatorStep4_WECSim_B.F_Excitation,
         windEmulatorStep4_WECSim_B.F_Restoring,
         &windEmulatorStep4_WECSim_B.sf_YawForceTransforms,
         &windEmulatorStep4_WECSim_DW.sf_YawForceTransforms);

      /* SignalConversion generated from: '<S61>/Product1' */
      windEmulatorStep4_WECSim_B.v[0] = windEmulatorStep4_WECSim_B.OUTPUT_1_1[12];
      windEmulatorStep4_WECSim_B.v[3] = windEmulatorStep4_WECSim_B.OUTPUT_1_1[6];
      windEmulatorStep4_WECSim_B.v[1] = windEmulatorStep4_WECSim_B.OUTPUT_1_1[13];
      windEmulatorStep4_WECSim_B.v[4] = windEmulatorStep4_WECSim_B.OUTPUT_1_1[7];
      windEmulatorStep4_WECSim_B.v[2] = windEmulatorStep4_WECSim_B.OUTPUT_1_1[14];
      windEmulatorStep4_WECSim_B.v[5] = windEmulatorStep4_WECSim_B.OUTPUT_1_1[8];
      for (i = 0; i < 6; i++) {
        /* Abs: '<S64>/Abs1' */
        Clock_tmp = std::abs(windEmulatorStep4_WECSim_B.v[i]);
        windEmulatorStep4_WECSim_B.Abs1_e[i] = Clock_tmp;

        /* Product: '<S64>/Product' */
        windEmulatorStep4_WECSim_B.vv[i] = windEmulatorStep4_WECSim_B.v[i] *
          Clock_tmp;
      }

      /* Product: '<S64>/Product1' incorporates:
       *  Constant: '<S67>/Constant'
       *  Product: '<S131>/Product'
       *  Product: '<S140>/Product1'
       *  Product: '<S143>/Product1'
       *  Product: '<S148>/Product1'
       *  Product: '<S153>/Product2'
       *  Product: '<S210>/Product'
       *  Product: '<S61>/Product1'
       *  Product: '<S69>/Product1'
       *  Product: '<S74>/Product2'
       */
      tmp_i = &windEmulatorStep4_WECSim_cal->Constant_Value_l.hf1.quadDrag[0];
      for (i = 0; i < 6; i++) {
        tmp_e[i] = windEmulatorStep4_WECSim_B.vv[i];
        tmp_f[i] = 0.0;
      }

      for (i = 0; i < 6; i++) {
        riseValLimit = tmp_e[i];
        for (i_0 = 0; i_0 < 6; i_0++) {
          tmp_f[i_0] += tmp_i[6 * i + i_0] * riseValLimit;
        }
      }

      for (i = 0; i < 6; i++) {
        /* Product: '<S64>/Product1' */
        Clock_tmp = tmp_f[i];
        windEmulatorStep4_WECSim_B.F_quadraticViscous[i] = Clock_tmp;

        /* Sum: '<S64>/VisSum ' incorporates:
         *  Constant: '<S78>/Constant'
         */
        windEmulatorStep4_WECSim_B.F_MorisonAndViscous[i] = Clock_tmp -
          windEmulatorStep4_WECSim_cal->Constant_Value_b[i];
      }

      /* Product: '<S61>/Product1' incorporates:
       *  Constant: '<S67>/Constant'
       *  Product: '<S131>/Product'
       *  Product: '<S140>/Product1'
       *  Product: '<S143>/Product1'
       *  Product: '<S148>/Product1'
       *  Product: '<S153>/Product2'
       *  Product: '<S210>/Product'
       *  Product: '<S64>/Product1'
       *  Product: '<S69>/Product1'
       *  Product: '<S74>/Product2'
       */
      tmp_i = &windEmulatorStep4_WECSim_cal->Constant_Value_l.hf1.linearDamping
        [0];
      for (i = 0; i < 6; i++) {
        tmp_e[i] = windEmulatorStep4_WECSim_B.v[i];
        tmp_f[i] = 0.0;
      }

      for (i = 0; i < 6; i++) {
        riseValLimit = tmp_e[i];
        for (i_0 = 0; i_0 < 6; i_0++) {
          tmp_f[i_0] += tmp_i[6 * i + i_0] * riseValLimit;
        }
      }

      for (i = 0; i < 6; i++) {
        /* Product: '<S61>/Product1' */
        Clock_tmp = tmp_f[i];
        windEmulatorStep4_WECSim_B.F_LinearDamping[i] = Clock_tmp;

        /* Sum: '<S60>/Sum' */
        windEmulatorStep4_WECSim_B.F_Total[i] =
          ((((windEmulatorStep4_WECSim_B.sf_YawForceTransforms.F_Excitation[i] -
              windEmulatorStep4_WECSim_B.sf_YawForceTransforms.F_RadiationDamping
              [i]) -
             windEmulatorStep4_WECSim_B.sf_YawForceTransforms.F_AddedMass[i]) -
            windEmulatorStep4_WECSim_B.sf_YawForceTransforms.F_Restoring[i]) -
           windEmulatorStep4_WECSim_B.F_MorisonAndViscous[i]) - Clock_tmp;
      }

      /* SimscapeInputBlock: '<S216>/INPUT_1_1_1' */
      windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[0] =
        windEmulatorStep4_WECSim_B.F_Total[0];
      windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[1] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[2] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[3] = 0.0;

      /* SimscapeInputBlock: '<S216>/INPUT_1_1_2' */
      windEmulatorStep4_WECSim_B.INPUT_1_1_2[0] =
        windEmulatorStep4_WECSim_B.F_Total[1];
      windEmulatorStep4_WECSim_B.INPUT_1_1_2[1] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_1_1_2[2] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_1_1_2[3] = 0.0;

      /* SimscapeInputBlock: '<S216>/INPUT_1_1_3' */
      windEmulatorStep4_WECSim_B.INPUT_1_1_3[0] =
        windEmulatorStep4_WECSim_B.F_Total[2];
      windEmulatorStep4_WECSim_B.INPUT_1_1_3[1] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_1_1_3[2] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_1_1_3[3] = 0.0;

      /* SimscapeInputBlock: '<S216>/INPUT_2_1_1' */
      windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[0] =
        windEmulatorStep4_WECSim_B.F_Total[3];
      windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[1] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[2] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[3] = 0.0;

      /* SimscapeInputBlock: '<S216>/INPUT_2_1_2' */
      windEmulatorStep4_WECSim_B.INPUT_2_1_2[0] =
        windEmulatorStep4_WECSim_B.F_Total[4];
      windEmulatorStep4_WECSim_B.INPUT_2_1_2[1] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_2_1_2[2] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_2_1_2[3] = 0.0;

      /* SimscapeInputBlock: '<S216>/INPUT_2_1_3' */
      windEmulatorStep4_WECSim_B.INPUT_2_1_3[0] =
        windEmulatorStep4_WECSim_B.F_Total[5];
      windEmulatorStep4_WECSim_B.INPUT_2_1_3[1] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_2_1_3[2] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_2_1_3[3] = 0.0;

      /* MATLAB Function: '<S160>/quaternion2EulXYZ' incorporates:
       *  SimscapeExecutionBlock: '<S216>/OUTPUT_1_1'
       */
      windEmulatorS_quaternion2EulXYZ(&windEmulatorStep4_WECSim_B.OUTPUT_1_1[15],
        &windEmulatorStep4_WECSim_B.sf_quaternion2EulXYZ_c,
        &windEmulatorStep4_WECSim_DW.sf_quaternion2EulXYZ_c);

      /* TransportDelay: '<S139>/Transport Delay' */
      for (i = 0; i < 6; i++) {
        riseValLimit = rt_TDelayInterpolate(windEmulatorStep4_WECSim_M->
          Timing.t[0] - windEmulatorStep4_WECSim_cal->TransportDelay_Delay_h,
          windEmulatorStep4_WECSim_DW.TransportDelay_RWORK_k[0],
          static_cast<real_T *>
          (windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_f[i]),
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c[i + 18],
          &windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c[i + 12],
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c[i],
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c[i + 6],
          windEmulatorStep4_WECSim_cal->TransportDelay_InitOutput_m,false,
          rtmIsMinorTimeStep(windEmulatorStep4_WECSim_M) && ((static_cast<real_T
          *>(windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_f[0]))
          [windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c[6] +
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c[18]] ==
          windEmulatorStep4_WECSim_M->Timing.t[0]));
        windEmulatorStep4_WECSim_B.TransportDelay_b[i] = riseValLimit;
      }

      /* End of TransportDelay: '<S139>/Transport Delay' */

      /* MATLAB Function: '<S212>/Yaw Kinematic Transforms' incorporates:
       *  Constant: '<S212>/Constant'
       *  SimscapeExecutionBlock: '<S216>/OUTPUT_1_1'
       */
      windEmul_YawKinematicTransforms
        (windEmulatorStep4_WECSim_cal->Constant_Value_o,
         &windEmulatorStep4_WECSim_B.OUTPUT_1_1[22],
         windEmulatorStep4_WECSim_B.sf_quaternion2EulXYZ_c.E,
         &windEmulatorStep4_WECSim_B.OUTPUT_1_1[25],
         &windEmulatorStep4_WECSim_B.OUTPUT_1_1[19],
         windEmulatorStep4_WECSim_B.TransportDelay_b,
         &windEmulatorStep4_WECSim_B.sf_YawKinematicTransforms_l,
         &windEmulatorStep4_WECSim_DW.sf_YawKinematicTransforms_l);

      /* Product: '<S210>/Product' incorporates:
       *  Constant: '<S146>/Constant'
       *  Product: '<S131>/Product'
       *  Product: '<S140>/Product1'
       *  Product: '<S143>/Product1'
       *  Product: '<S148>/Product1'
       *  Product: '<S153>/Product2'
       *  Product: '<S61>/Product1'
       *  Product: '<S64>/Product1'
       *  Product: '<S69>/Product1'
       *  Product: '<S74>/Product2'
       */
      tmp_i = &windEmulatorStep4_WECSim_cal->Constant_Value.hf1.fDamping[0];
      for (i = 0; i < 6; i++) {
        tmp_e[i] =
          windEmulatorStep4_WECSim_B.sf_YawKinematicTransforms_l.velLoc[i];
        tmp_f[i] = 0.0;
      }

      for (i = 0; i < 6; i++) {
        riseValLimit = tmp_e[i];
        for (i_0 = 0; i_0 < 6; i_0++) {
          tmp_f[i_0] += tmp_i[6 * i + i_0] * riseValLimit;
        }
      }

      for (i = 0; i < 6; i++) {
        /* Product: '<S210>/Product' */
        windEmulatorStep4_WECSim_B.F_SingleFrequency_b[i] = tmp_f[i];
      }

      /* Product: '<S148>/Product1' incorporates:
       *  Constant: '<S146>/Constant'
       *  Product: '<S131>/Product'
       *  Product: '<S140>/Product1'
       *  Product: '<S143>/Product1'
       *  Product: '<S153>/Product2'
       *  Product: '<S210>/Product'
       *  Product: '<S61>/Product1'
       *  Product: '<S64>/Product1'
       *  Product: '<S69>/Product1'
       *  Product: '<S74>/Product2'
       */
      tmp_i = &windEmulatorStep4_WECSim_cal->Constant_Value.hf1.fAddedMass[0];
      for (i = 0; i < 6; i++) {
        tmp_e[i] =
          windEmulatorStep4_WECSim_B.sf_YawKinematicTransforms_l.accLoc[i];
        tmp_f[i] = 0.0;
      }

      for (i = 0; i < 6; i++) {
        riseValLimit = tmp_e[i];
        for (i_0 = 0; i_0 < 6; i_0++) {
          tmp_f[i_0] += tmp_i[6 * i + i_0] * riseValLimit;
        }
      }

      for (i = 0; i < 6; i++) {
        /* Product: '<S148>/Product1' */
        windEmulatorStep4_WECSim_B.F_AddedMass_k[i] = tmp_f[i];
      }

      if (tmp_g) {
        /* Product: '<S203>/Divide' incorporates:
         *  Constant: '<S203>/Constant4'
         *  Constant: '<S203>/Ramp Time'
         */
        windEmulatorStep4_WECSim_B.frequency_g =
          windEmulatorStep4_WECSim_cal->Constant4_Value_c /
          windEmulatorStep4_WECSim_cal->RampTime_Value_n;
      }

      /* RelationalOperator: '<S147>/Less Than' incorporates:
       *  Constant: '<S147>/Ramp Function Time'
       */
      windEmulatorStep4_WECSim_B.LessThan_k =
        (windEmulatorStep4_WECSim_B.Clock_n <
         windEmulatorStep4_WECSim_cal->RampFunctionTime_Value_d);

      /* Switch: '<S147>/Switch' */
      if (windEmulatorStep4_WECSim_B.LessThan_k) {
        /* Product: '<S203>/Product2' */
        windEmulatorStep4_WECSim_B.Product2_c =
          windEmulatorStep4_WECSim_B.Clock_n *
          windEmulatorStep4_WECSim_B.frequency_g;

        /* Sum: '<S203>/Add2' incorporates:
         *  Constant: '<S203>/Constant3'
         */
        windEmulatorStep4_WECSim_B.Add2_g =
          windEmulatorStep4_WECSim_B.Product2_c +
          windEmulatorStep4_WECSim_cal->Constant3_Value_p;

        /* Sin: '<S203>/Sine Wave Function' */
        windEmulatorStep4_WECSim_B.SineWaveFunction = std::sin
          (windEmulatorStep4_WECSim_cal->SineWaveFunction_Freq_e *
           windEmulatorStep4_WECSim_B.Add2_g +
           windEmulatorStep4_WECSim_cal->SineWaveFunction_Phase_i) *
          windEmulatorStep4_WECSim_cal->SineWaveFunction_Amp_a +
          windEmulatorStep4_WECSim_cal->SineWaveFunction_Bias_k;

        /* Sum: '<S203>/Add1' incorporates:
         *  Constant: '<S203>/Constant2'
         */
        windEmulatorStep4_WECSim_B.Add1_lhe =
          windEmulatorStep4_WECSim_cal->Constant2_Value_p +
          windEmulatorStep4_WECSim_B.SineWaveFunction;

        /* Product: '<S203>/Product1' incorporates:
         *  Constant: '<S203>/Constant1'
         */
        windEmulatorStep4_WECSim_B.Product1_d =
          windEmulatorStep4_WECSim_B.Add1_lhe *
          windEmulatorStep4_WECSim_cal->Constant1_Value_l;

        /* Switch: '<S147>/Switch' */
        windEmulatorStep4_WECSim_B.R_p = windEmulatorStep4_WECSim_B.Product1_d;
      } else {
        /* Switch: '<S147>/Switch' incorporates:
         *  Constant: '<S147>/Constant'
         */
        windEmulatorStep4_WECSim_B.R_p =
          windEmulatorStep4_WECSim_cal->Constant_Value_j;
      }

      /* End of Switch: '<S147>/Switch' */
      if (tmp_g) {
        for (i = 0; i < 6; i++) {
          /* Product: '<S205>/Product3' incorporates:
           *  Constant: '<S146>/Constant'
           *  Constant: '<S205>/Wave Amplitude1'
           */
          windEmulatorStep4_WECSim_B.Product3_p[i] =
            windEmulatorStep4_WECSim_cal->WaveAmplitude1_Value_d *
            windEmulatorStep4_WECSim_cal->Constant_Value.hf1.fExt.md[i];
        }
      }

      /* Product: '<S205>/Product' incorporates:
       *  Constant: '<S205>/Wave Frequency'
       */
      windEmulatorStep4_WECSim_B.Product_oc = windEmulatorStep4_WECSim_B.Clock_n
        * windEmulatorStep4_WECSim_cal->WaveFrequency_Value_i;

      /* Sum: '<S205>/Add3' incorporates:
       *  Constant: '<S205>/Center of Gravity'
       *  Constant: '<S205>/Constant1'
       */
      windEmulatorStep4_WECSim_B.x_cg_j[0] =
        windEmulatorStep4_WECSim_B.OUTPUT_1_1[22] -
        windEmulatorStep4_WECSim_cal->CenterofGravity_Value_h[0];
      windEmulatorStep4_WECSim_B.x_cg_j[3] =
        windEmulatorStep4_WECSim_B.sf_quaternion2EulXYZ_c.E[0] -
        windEmulatorStep4_WECSim_cal->Constant1_Value_gn[0];
      windEmulatorStep4_WECSim_B.x_cg_j[1] =
        windEmulatorStep4_WECSim_B.OUTPUT_1_1[23] -
        windEmulatorStep4_WECSim_cal->CenterofGravity_Value_h[1];
      windEmulatorStep4_WECSim_B.x_cg_j[4] =
        windEmulatorStep4_WECSim_B.sf_quaternion2EulXYZ_c.E[1] -
        windEmulatorStep4_WECSim_cal->Constant1_Value_gn[1];
      windEmulatorStep4_WECSim_B.x_cg_j[2] =
        windEmulatorStep4_WECSim_B.OUTPUT_1_1[24] -
        windEmulatorStep4_WECSim_cal->CenterofGravity_Value_h[2];
      windEmulatorStep4_WECSim_B.x_cg_j[5] =
        windEmulatorStep4_WECSim_B.sf_quaternion2EulXYZ_c.E[2] -
        windEmulatorStep4_WECSim_cal->Constant1_Value_gn[2];

      /* MATLAB Function: '<S205>/MATLAB Function1' incorporates:
       *  Constant: '<S205>/Displacement Phase Enable1'
       *  Constant: '<S205>/Wave Number'
       *  Constant: '<S205>/Wave direction'
       *  Sum: '<S205>/Add3'
       */
      windEmulatorSte_MATLABFunction1(&windEmulatorStep4_WECSim_B.x_cg_j[0],
        windEmulatorStep4_WECSim_cal->DisplacementPhaseEnable1_Valu_h,
        windEmulatorStep4_WECSim_cal->WaveNumber_Value_i,
        windEmulatorStep4_WECSim_cal->Wavedirection_Value_i,
        &windEmulatorStep4_WECSim_B.sf_MATLABFunction1_e,
        &windEmulatorStep4_WECSim_DW.sf_MATLABFunction1_e);

      /* Sum: '<S205>/Add' incorporates:
       *  Constant: '<S205>/Constant'
       */
      windEmulatorStep4_WECSim_B.Add_m = (windEmulatorStep4_WECSim_B.Product_oc
        + windEmulatorStep4_WECSim_cal->Constant_Value_hr) +
        windEmulatorStep4_WECSim_B.sf_MATLABFunction1_e.dispPhase;

      /* Sin: '<S205>/Sine Wave Function1' */
      windEmulatorStep4_WECSim_B.coswt_k = std::sin
        (windEmulatorStep4_WECSim_cal->SineWaveFunction1_Freq_i *
         windEmulatorStep4_WECSim_B.Add_m +
         windEmulatorStep4_WECSim_cal->SineWaveFunction1_Phase_h) *
        windEmulatorStep4_WECSim_cal->SineWaveFunction1_Amp_a +
        windEmulatorStep4_WECSim_cal->SineWaveFunction1_Bias_d;

      /* Sum: '<S205>/Add2' */
      windEmulatorStep4_WECSim_B.Add2_f = windEmulatorStep4_WECSim_B.Product_oc
        + windEmulatorStep4_WECSim_B.sf_MATLABFunction1_e.dispPhase;

      /* Sin: '<S205>/Sine Wave Function' */
      windEmulatorStep4_WECSim_B.sinwt_g = std::sin
        (windEmulatorStep4_WECSim_cal->SineWaveFunction_Freq_d *
         windEmulatorStep4_WECSim_B.Add2_f +
         windEmulatorStep4_WECSim_cal->SineWaveFunction_Phase_ih) *
        windEmulatorStep4_WECSim_cal->SineWaveFunction_Amp_b +
        windEmulatorStep4_WECSim_cal->SineWaveFunction_Bias_j;
      for (i = 0; i < 6; i++) {
        /* Product: '<S205>/Product1' incorporates:
         *  Constant: '<S146>/Constant'
         *  Constant: '<S205>/Wave Amplitude'
         */
        Clock_tmp = windEmulatorStep4_WECSim_cal->WaveAmplitude_Value_f *
          windEmulatorStep4_WECSim_cal->Constant_Value.hf1.fExt.re[i] *
          windEmulatorStep4_WECSim_B.coswt_k;
        windEmulatorStep4_WECSim_B.Product1_a[i] = Clock_tmp;

        /* Product: '<S205>/Product2' incorporates:
         *  Constant: '<S146>/Constant'
         *  Constant: '<S205>/Wave Amplitude'
         */
        riseValLimit = windEmulatorStep4_WECSim_cal->WaveAmplitude_Value_f *
          windEmulatorStep4_WECSim_B.sinwt_g *
          windEmulatorStep4_WECSim_cal->Constant_Value.hf1.fExt.im[i];
        windEmulatorStep4_WECSim_B.Product2_p[i] = riseValLimit;

        /* Sum: '<S205>/Add1' incorporates:
         *  Product: '<S205>/Product1'
         *  Product: '<S205>/Product2'
         *  Product: '<S205>/Product3'
         */
        Clock_tmp = (windEmulatorStep4_WECSim_B.Product3_p[i] + Clock_tmp) -
          riseValLimit;
        windEmulatorStep4_WECSim_B.Add1_e[i] = Clock_tmp;

        /* Sum: '<S147>/Add' incorporates:
         *  Constant: '<S207>/Constant'
         *  Constant: '<S208>/Constant'
         *  Sum: '<S205>/Add1'
         */
        Clock_tmp = (Clock_tmp + windEmulatorStep4_WECSim_cal->
                     Constant_Value_i1[i]) +
          windEmulatorStep4_WECSim_cal->Constant_Value_am[i];
        windEmulatorStep4_WECSim_B.F_wave_i[i] = Clock_tmp;

        /* Product: '<S147>/Product' incorporates:
         *  Sum: '<S147>/Add'
         */
        windEmulatorStep4_WECSim_B.F_Excitation_f[i] =
          windEmulatorStep4_WECSim_B.R_p * Clock_tmp;
      }

      if (tmp_g) {
        /* Product: '<S154>/Product' incorporates:
         *  Constant: '<S146>/Constant'
         *  Constant: '<S154>/Gravity'
         */
        windEmulatorStep4_WECSim_B.F_Gravity_h =
          windEmulatorStep4_WECSim_cal->Gravity_Value_a *
          windEmulatorStep4_WECSim_cal->Constant_Value.hf1.mass;

        /* Product: '<S154>/Product1' incorporates:
         *  Constant: '<S146>/Constant'
         *  Constant: '<S154>/Gravity'
         *  Constant: '<S154>/Water Density'
         */
        windEmulatorStep4_WECSim_B.F_Buoyancy_k =
          windEmulatorStep4_WECSim_cal->WaterDensity_Value_j *
          windEmulatorStep4_WECSim_cal->Gravity_Value_a *
          windEmulatorStep4_WECSim_cal->Constant_Value.hf1.volume;

        /* Sum: '<S154>/Add1' */
        windEmulatorStep4_WECSim_B.Add1_l =
          windEmulatorStep4_WECSim_B.F_Gravity_h -
          windEmulatorStep4_WECSim_B.F_Buoyancy_k;

        /* Assignment: '<S154>/Assignment (Add Net Bouyancy Force  to Z-Direction)1' incorporates:
         *  Constant: '<S154>/Constant2'
         */
        windEmulatorStep4_WECSim_B.AssignmentAddNetBouyancyForce_c[0] =
          windEmulatorStep4_WECSim_cal->Constant2_Value_m3[0];

        /* Sum: '<S154>/Add3' incorporates:
         *  Constant: '<S146>/Constant'
         *  Constant: '<S154>/Center of Gravity'
         */
        windEmulatorStep4_WECSim_B.Add3_m[0] =
          windEmulatorStep4_WECSim_cal->Constant_Value.hf1.centerBuoyancy[0] -
          windEmulatorStep4_WECSim_cal->CenterofGravity_Value_k[0];

        /* Assignment: '<S154>/Assignment (Add Net Bouyancy Force  to Z-Direction)1' incorporates:
         *  Constant: '<S154>/Constant2'
         */
        windEmulatorStep4_WECSim_B.AssignmentAddNetBouyancyForce_c[1] =
          windEmulatorStep4_WECSim_cal->Constant2_Value_m3[1];

        /* Sum: '<S154>/Add3' incorporates:
         *  Constant: '<S146>/Constant'
         *  Constant: '<S154>/Center of Gravity'
         */
        windEmulatorStep4_WECSim_B.Add3_m[1] =
          windEmulatorStep4_WECSim_cal->Constant_Value.hf1.centerBuoyancy[1] -
          windEmulatorStep4_WECSim_cal->CenterofGravity_Value_k[1];
        windEmulatorStep4_WECSim_B.Add3_m[2] =
          windEmulatorStep4_WECSim_cal->Constant_Value.hf1.centerBuoyancy[2] -
          windEmulatorStep4_WECSim_cal->CenterofGravity_Value_k[2];

        /* Assignment: '<S154>/Assignment (Add Net Bouyancy Force  to Z-Direction)1' */
        windEmulatorStep4_WECSim_B.AssignmentAddNetBouyancyForce_c[2] =
          windEmulatorStep4_WECSim_B.F_Buoyancy_k;

        /* Product: '<S155>/Element product' */
        windEmulatorStep4_WECSim_B.Elementproduct_p[0] =
          windEmulatorStep4_WECSim_B.AssignmentAddNetBouyancyForce_c[1] *
          windEmulatorStep4_WECSim_B.Add3_m[2];
        windEmulatorStep4_WECSim_B.Elementproduct_p[1] =
          windEmulatorStep4_WECSim_B.Add3_m[0] *
          windEmulatorStep4_WECSim_B.AssignmentAddNetBouyancyForce_c[2];
        windEmulatorStep4_WECSim_B.Elementproduct_p[2] =
          windEmulatorStep4_WECSim_B.AssignmentAddNetBouyancyForce_c[0] *
          windEmulatorStep4_WECSim_B.Add3_m[1];
        windEmulatorStep4_WECSim_B.Elementproduct_p[3] =
          windEmulatorStep4_WECSim_B.Add3_m[1] *
          windEmulatorStep4_WECSim_B.AssignmentAddNetBouyancyForce_c[2];
        windEmulatorStep4_WECSim_B.Elementproduct_p[4] =
          windEmulatorStep4_WECSim_B.AssignmentAddNetBouyancyForce_c[0] *
          windEmulatorStep4_WECSim_B.Add3_m[2];
        windEmulatorStep4_WECSim_B.Elementproduct_p[5] =
          windEmulatorStep4_WECSim_B.Add3_m[0] *
          windEmulatorStep4_WECSim_B.AssignmentAddNetBouyancyForce_c[1];

        /* Sum: '<S155>/Add3' */
        windEmulatorStep4_WECSim_B.Add3_h[0] =
          windEmulatorStep4_WECSim_B.Elementproduct_p[0] -
          windEmulatorStep4_WECSim_B.Elementproduct_p[3];
        windEmulatorStep4_WECSim_B.Add3_h[1] =
          windEmulatorStep4_WECSim_B.Elementproduct_p[1] -
          windEmulatorStep4_WECSim_B.Elementproduct_p[4];
        windEmulatorStep4_WECSim_B.Add3_h[2] =
          windEmulatorStep4_WECSim_B.Elementproduct_p[2] -
          windEmulatorStep4_WECSim_B.Elementproduct_p[5];
        for (i = 0; i < 6; i++) {
          /* Assignment: '<S154>/Assignment (Add Net Bouyancy Force  to Z-Direction)' incorporates:
           *  Constant: '<S154>/Constant1'
           */
          windEmulatorStep4_WECSim_B.VerticalBuoyancyForce_g[i] =
            windEmulatorStep4_WECSim_cal->Constant1_Value_a[i];

          /* Assignment: '<S154>/Assignment (Add Net Bouyancy Force  to Z-Direction)2' incorporates:
           *  Constant: '<S154>/Constant1'
           */
          windEmulatorStep4_WECSim_B.Rotationalbuoyancyforce_l[i] =
            windEmulatorStep4_WECSim_cal->Constant1_Value_a[i];
        }

        /* Assignment: '<S154>/Assignment (Add Net Bouyancy Force  to Z-Direction)' */
        windEmulatorStep4_WECSim_B.VerticalBuoyancyForce_g[2] =
          windEmulatorStep4_WECSim_B.Add1_l;

        /* Assignment: '<S154>/Assignment (Add Net Bouyancy Force  to Z-Direction)2' incorporates:
         *  Sum: '<S155>/Add3'
         */
        windEmulatorStep4_WECSim_B.Rotationalbuoyancyforce_l[3] =
          windEmulatorStep4_WECSim_B.Add3_h[0];
        windEmulatorStep4_WECSim_B.Rotationalbuoyancyforce_l[4] =
          windEmulatorStep4_WECSim_B.Add3_h[1];
        windEmulatorStep4_WECSim_B.Rotationalbuoyancyforce_l[5] =
          windEmulatorStep4_WECSim_B.Add3_h[2];
        for (i = 0; i < 6; i++) {
          /* Sum: '<S154>/Add2' */
          windEmulatorStep4_WECSim_B.Netbuoyancyforce_c[i] =
            windEmulatorStep4_WECSim_B.VerticalBuoyancyForce_g[i] +
            windEmulatorStep4_WECSim_B.Rotationalbuoyancyforce_l[i];
        }
      }

      /* Sum: '<S142>/Add' incorporates:
       *  Constant: '<S142>/Constant'
       *  Constant: '<S146>/Constant'
       */
      windEmulatorStep4_WECSim_B.x_cg_o[0] =
        windEmulatorStep4_WECSim_B.sf_YawKinematicTransforms_l.dispLoc[0] -
        windEmulatorStep4_WECSim_cal->Constant_Value.hf1.centerGravity[0];
      windEmulatorStep4_WECSim_B.x_cg_o[3] =
        windEmulatorStep4_WECSim_B.sf_YawKinematicTransforms_l.dispLoc[3] -
        windEmulatorStep4_WECSim_cal->Constant_Value_iv[0];
      windEmulatorStep4_WECSim_B.x_cg_o[1] =
        windEmulatorStep4_WECSim_B.sf_YawKinematicTransforms_l.dispLoc[1] -
        windEmulatorStep4_WECSim_cal->Constant_Value.hf1.centerGravity[1];
      windEmulatorStep4_WECSim_B.x_cg_o[4] =
        windEmulatorStep4_WECSim_B.sf_YawKinematicTransforms_l.dispLoc[4] -
        windEmulatorStep4_WECSim_cal->Constant_Value_iv[1];
      windEmulatorStep4_WECSim_B.x_cg_o[2] =
        windEmulatorStep4_WECSim_B.sf_YawKinematicTransforms_l.dispLoc[2] -
        windEmulatorStep4_WECSim_cal->Constant_Value.hf1.centerGravity[2];
      windEmulatorStep4_WECSim_B.x_cg_o[5] =
        windEmulatorStep4_WECSim_B.sf_YawKinematicTransforms_l.dispLoc[5] -
        windEmulatorStep4_WECSim_cal->Constant_Value_iv[2];

      /* Product: '<S153>/Product2' incorporates:
       *  Constant: '<S146>/Constant'
       *  Product: '<S131>/Product'
       *  Product: '<S140>/Product1'
       *  Product: '<S143>/Product1'
       *  Product: '<S148>/Product1'
       *  Product: '<S210>/Product'
       *  Product: '<S61>/Product1'
       *  Product: '<S64>/Product1'
       *  Product: '<S69>/Product1'
       *  Product: '<S74>/Product2'
       *  Sum: '<S142>/Add'
       */
      tmp_i =
        &windEmulatorStep4_WECSim_cal->Constant_Value.hf1.linearHydroRestCoef[0];
      for (i = 0; i < 6; i++) {
        tmp_e[i] = windEmulatorStep4_WECSim_B.x_cg_o[i];
        windEmulatorStep4_WECSim_B.LinearRestoringForce_l[i] = 0.0;
      }

      for (i = 0; i < 6; i++) {
        riseValLimit = tmp_e[i];
        for (i_0 = 0; i_0 < 6; i_0++) {
          windEmulatorStep4_WECSim_B.LinearRestoringForce_l[i_0] += tmp_i[6 * i
            + i_0] * riseValLimit;
        }
      }

      for (i = 0; i < 6; i++) {
        /* Sum: '<S153>/Add' incorporates:
         *  Product: '<S153>/Product2'
         */
        windEmulatorStep4_WECSim_B.F_Restoring_a[i] =
          windEmulatorStep4_WECSim_B.Netbuoyancyforce_c[i] +
          windEmulatorStep4_WECSim_B.LinearRestoringForce_l[i];
      }

      /* MATLAB Function: '<S149>/Yaw Force Transforms' incorporates:
       *  Constant: '<S149>/Constant'
       *  SimscapeExecutionBlock: '<S216>/OUTPUT_1_1'
       */
      windEmulator_YawForceTransforms
        (windEmulatorStep4_WECSim_cal->Constant_Value_p,
         &windEmulatorStep4_WECSim_B.OUTPUT_1_1[22],
         windEmulatorStep4_WECSim_B.sf_quaternion2EulXYZ_c.E,
         windEmulatorStep4_WECSim_B.F_SingleFrequency_b,
         windEmulatorStep4_WECSim_B.F_AddedMass_k,
         windEmulatorStep4_WECSim_B.F_Excitation_f,
         windEmulatorStep4_WECSim_B.F_Restoring_a,
         &windEmulatorStep4_WECSim_B.sf_YawForceTransforms_i,
         &windEmulatorStep4_WECSim_DW.sf_YawForceTransforms_i);

      /* SignalConversion generated from: '<S140>/Product1' */
      windEmulatorStep4_WECSim_B.v_j[0] = windEmulatorStep4_WECSim_B.OUTPUT_1_1
        [25];
      windEmulatorStep4_WECSim_B.v_j[3] = windEmulatorStep4_WECSim_B.OUTPUT_1_1
        [19];
      windEmulatorStep4_WECSim_B.v_j[1] = windEmulatorStep4_WECSim_B.OUTPUT_1_1
        [26];
      windEmulatorStep4_WECSim_B.v_j[4] = windEmulatorStep4_WECSim_B.OUTPUT_1_1
        [20];
      windEmulatorStep4_WECSim_B.v_j[2] = windEmulatorStep4_WECSim_B.OUTPUT_1_1
        [27];
      windEmulatorStep4_WECSim_B.v_j[5] = windEmulatorStep4_WECSim_B.OUTPUT_1_1
        [21];
      for (i = 0; i < 6; i++) {
        /* Abs: '<S143>/Abs1' */
        Clock_tmp = std::abs(windEmulatorStep4_WECSim_B.v_j[i]);
        windEmulatorStep4_WECSim_B.Abs1_b[i] = Clock_tmp;

        /* Product: '<S143>/Product' */
        windEmulatorStep4_WECSim_B.vv_j[i] = windEmulatorStep4_WECSim_B.v_j[i] *
          Clock_tmp;
      }

      /* Product: '<S143>/Product1' incorporates:
       *  Constant: '<S146>/Constant'
       *  Product: '<S131>/Product'
       *  Product: '<S140>/Product1'
       *  Product: '<S148>/Product1'
       *  Product: '<S153>/Product2'
       *  Product: '<S210>/Product'
       *  Product: '<S61>/Product1'
       *  Product: '<S64>/Product1'
       *  Product: '<S69>/Product1'
       *  Product: '<S74>/Product2'
       */
      tmp_i = &windEmulatorStep4_WECSim_cal->Constant_Value.hf1.quadDrag[0];
      for (i = 0; i < 6; i++) {
        tmp_e[i] = windEmulatorStep4_WECSim_B.vv_j[i];
        tmp_f[i] = 0.0;
      }

      for (i = 0; i < 6; i++) {
        riseValLimit = tmp_e[i];
        for (i_0 = 0; i_0 < 6; i_0++) {
          tmp_f[i_0] += tmp_i[6 * i + i_0] * riseValLimit;
        }
      }

      for (i = 0; i < 6; i++) {
        /* Product: '<S143>/Product1' */
        Clock_tmp = tmp_f[i];
        windEmulatorStep4_WECSim_B.F_quadraticViscous_f[i] = Clock_tmp;

        /* Sum: '<S143>/VisSum ' incorporates:
         *  Constant: '<S157>/Constant'
         */
        windEmulatorStep4_WECSim_B.F_MorisonAndViscous_d[i] = Clock_tmp -
          windEmulatorStep4_WECSim_cal->Constant_Value_f[i];
      }

      /* Product: '<S140>/Product1' incorporates:
       *  Constant: '<S146>/Constant'
       *  Product: '<S131>/Product'
       *  Product: '<S143>/Product1'
       *  Product: '<S148>/Product1'
       *  Product: '<S153>/Product2'
       *  Product: '<S210>/Product'
       *  Product: '<S61>/Product1'
       *  Product: '<S64>/Product1'
       *  Product: '<S69>/Product1'
       *  Product: '<S74>/Product2'
       */
      tmp_i = &windEmulatorStep4_WECSim_cal->Constant_Value.hf1.linearDamping[0];
      for (i = 0; i < 6; i++) {
        tmp_e[i] = windEmulatorStep4_WECSim_B.v_j[i];
        tmp_f[i] = 0.0;
      }

      for (i = 0; i < 6; i++) {
        riseValLimit = tmp_e[i];
        for (i_0 = 0; i_0 < 6; i_0++) {
          tmp_f[i_0] += tmp_i[6 * i + i_0] * riseValLimit;
        }
      }

      for (i = 0; i < 6; i++) {
        /* Product: '<S140>/Product1' */
        Clock_tmp = tmp_f[i];
        windEmulatorStep4_WECSim_B.F_LinearDamping_k[i] = Clock_tmp;

        /* Sum: '<S139>/Sum' */
        windEmulatorStep4_WECSim_B.F_Total_c[i] =
          ((((windEmulatorStep4_WECSim_B.sf_YawForceTransforms_i.F_Excitation[i]
              - windEmulatorStep4_WECSim_B.sf_YawForceTransforms_i.F_RadiationDamping
              [i]) -
             windEmulatorStep4_WECSim_B.sf_YawForceTransforms_i.F_AddedMass[i])
            - windEmulatorStep4_WECSim_B.sf_YawForceTransforms_i.F_Restoring[i])
           - windEmulatorStep4_WECSim_B.F_MorisonAndViscous_d[i]) - Clock_tmp;
      }

      /* SimscapeInputBlock: '<S216>/INPUT_3_1_1' */
      windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[0] =
        windEmulatorStep4_WECSim_B.F_Total_c[0];
      windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[1] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[2] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[3] = 0.0;

      /* SimscapeInputBlock: '<S216>/INPUT_3_1_2' */
      windEmulatorStep4_WECSim_B.INPUT_3_1_2[0] =
        windEmulatorStep4_WECSim_B.F_Total_c[1];
      windEmulatorStep4_WECSim_B.INPUT_3_1_2[1] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_3_1_2[2] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_3_1_2[3] = 0.0;

      /* SimscapeInputBlock: '<S216>/INPUT_3_1_3' */
      windEmulatorStep4_WECSim_B.INPUT_3_1_3[0] =
        windEmulatorStep4_WECSim_B.F_Total_c[2];
      windEmulatorStep4_WECSim_B.INPUT_3_1_3[1] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_3_1_3[2] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_3_1_3[3] = 0.0;

      /* SimscapeInputBlock: '<S216>/INPUT_4_1_1' */
      windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[0] =
        windEmulatorStep4_WECSim_B.F_Total_c[3];
      windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[1] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[2] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[3] = 0.0;

      /* SimscapeInputBlock: '<S216>/INPUT_4_1_2' */
      windEmulatorStep4_WECSim_B.INPUT_4_1_2[0] =
        windEmulatorStep4_WECSim_B.F_Total_c[4];
      windEmulatorStep4_WECSim_B.INPUT_4_1_2[1] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_4_1_2[2] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_4_1_2[3] = 0.0;

      /* SimscapeInputBlock: '<S216>/INPUT_4_1_3' */
      windEmulatorStep4_WECSim_B.INPUT_4_1_3[0] =
        windEmulatorStep4_WECSim_B.F_Total_c[5];
      windEmulatorStep4_WECSim_B.INPUT_4_1_3[1] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_4_1_3[2] = 0.0;
      windEmulatorStep4_WECSim_B.INPUT_4_1_3[3] = 0.0;

      /* SimscapeExecutionBlock: '<S216>/OUTPUT_1_0' */
      simulationData = static_cast<NeslSimulationData *>
        (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_SimData_i);
      time_b = u1;
      simulationData->mData->mTime.mN = 1;
      simulationData->mData->mTime.mX = &time_b;
      simulationData->mData->mContStates.mN = 0;
      simulationData->mData->mContStates.mX = NULL;
      simulationData->mData->mDiscStates.mN = 0;
      simulationData->mData->mDiscStates.mX =
        &windEmulatorStep4_WECSim_DW.OUTPUT_1_0_Discrete_p;
      simulationData->mData->mModeVector.mN = 0;
      simulationData->mData->mModeVector.mX =
        &windEmulatorStep4_WECSim_DW.OUTPUT_1_0_Modes_j;
      tmp_11 = false;
      simulationData->mData->mFoundZcEvents = tmp_11;
      simulationData->mData->mHadEvents = false;
      simulationData->mData->mIsMajorTimeStep = f;
      f = false;
      simulationData->mData->mIsSolverAssertCheck = f;
      simulationData->mData->mIsSolverCheckingCIC = false;
      simulationData->mData->mIsComputingJacobian = false;
      simulationData->mData->mIsEvaluatingF0 = false;
      simulationData->mData->mIsSolverRequestingReset = false;
      simulationData->mData->mIsModeUpdateTimeStep = tmp_h;
      tmp_d[0] = 0;
      tmp_c[0] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[0];
      tmp_c[1] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[1];
      tmp_c[2] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[2];
      tmp_c[3] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[3];
      tmp_d[1] = 4;
      tmp_c[4] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[0];
      tmp_c[5] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[1];
      tmp_c[6] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[2];
      tmp_c[7] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[3];
      tmp_d[2] = 8;
      tmp_c[8] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[0];
      tmp_c[9] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[1];
      tmp_c[10] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[2];
      tmp_c[11] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[3];
      tmp_d[3] = 12;
      tmp_c[12] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[0];
      tmp_c[13] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[1];
      tmp_c[14] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[2];
      tmp_c[15] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[3];
      tmp_d[4] = 16;
      tmp_c[16] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[0];
      tmp_c[17] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[1];
      tmp_c[18] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[2];
      tmp_c[19] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[3];
      tmp_d[5] = 20;
      tmp_c[20] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[0];
      tmp_c[21] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[1];
      tmp_c[22] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[2];
      tmp_c[23] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[3];
      tmp_d[6] = 24;
      tmp_c[24] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[0];
      tmp_c[25] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[1];
      tmp_c[26] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[2];
      tmp_c[27] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[3];
      tmp_d[7] = 28;
      tmp_c[28] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[0];
      tmp_c[29] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[1];
      tmp_c[30] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[2];
      tmp_c[31] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[3];
      tmp_d[8] = 32;
      tmp_c[32] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[0];
      tmp_c[33] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[1];
      tmp_c[34] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[2];
      tmp_c[35] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[3];
      tmp_d[9] = 36;
      tmp_c[36] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[0];
      tmp_c[37] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[1];
      tmp_c[38] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[2];
      tmp_c[39] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[3];
      tmp_d[10] = 40;
      tmp_c[40] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[0];
      tmp_c[41] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[1];
      tmp_c[42] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[2];
      tmp_c[43] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[3];
      tmp_d[11] = 44;
      tmp_c[44] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[0];
      tmp_c[45] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[1];
      tmp_c[46] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[2];
      tmp_c[47] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[3];
      tmp_d[12] = 48;
      tmp_c[48] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[0];
      tmp_c[49] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[1];
      tmp_c[50] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[2];
      tmp_c[51] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[3];
      tmp_d[13] = 52;
      tmp_c[52] = windEmulatorStep4_WECSim_B.STATE_1_p[0];
      tmp_c[53] = windEmulatorStep4_WECSim_B.STATE_1_p[1];
      tmp_d[14] = 54;
      simulationData->mData->mInputValues.mN = 54;
      simulationData->mData->mInputValues.mX = &tmp_c[0];
      simulationData->mData->mInputOffsets.mN = 15;
      simulationData->mData->mInputOffsets.mX = &tmp_d[0];
      simulationData->mData->mOutputs.mN = 32;
      simulationData->mData->mOutputs.mX =
        &windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[0];
      simulationData->mData->mTolerances.mN = 0;
      simulationData->mData->mTolerances.mX = NULL;
      simulationData->mData->mCstateHasChanged = false;
      simulationData->mData->mDstateHasChanged = false;
      time_c = deltaT_tmp;
      simulationData->mData->mTime.mN = 1;
      simulationData->mData->mTime.mX = &time_c;
      simulationData->mData->mSampleHits.mN = 0;
      simulationData->mData->mSampleHits.mX = NULL;
      simulationData->mData->mIsFundamentalSampleHit = false;
      simulationData->mData->mHadEvents = false;
      simulator = static_cast<NeslSimulator *>
        (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_Simulator_l);
      diag = static_cast<NeuDiagnosticManager *>
        (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_DiagMgr_a);
      diagTree = neu_diagnostic_manager_get_initial_tree(diag);
      i = ne_simulator_method(simulator, NESL_SIM_OUTPUTS, simulationData, diag);
      if (i != 0) {
        f = error_buffer_is_empty(rtmGetErrorStatus(windEmulatorStep4_WECSim_M));
        if (f) {
          msg = rtw_diagnostics_msg(diagTree);
          rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
        }
      }

      /* SignalConversion generated from: '<S59>/Assignment6' */
      windEmulatorStep4_WECSim_B.TmpSignalConversionAtAssignment[0] =
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[14];
      windEmulatorStep4_WECSim_B.TmpSignalConversionAtAssignment[1] =
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[16];
      windEmulatorStep4_WECSim_B.TmpSignalConversionAtAssignment[2] =
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[15];

      /* Gain: '<S59>/Gain6' */
      windEmulatorStep4_WECSim_B.Gain6 =
        windEmulatorStep4_WECSim_cal->Gain6_Gain *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[15];

      /* SignalConversion generated from: '<S59>/Assignment7' */
      windEmulatorStep4_WECSim_B.TmpSignalConversionAtAssignme_e[0] =
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[17];
      windEmulatorStep4_WECSim_B.TmpSignalConversionAtAssignme_e[1] =
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[19];
      windEmulatorStep4_WECSim_B.TmpSignalConversionAtAssignme_e[2] =
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[18];

      /* Gain: '<S59>/Gain7' */
      windEmulatorStep4_WECSim_B.Gain7 =
        windEmulatorStep4_WECSim_cal->Gain7_Gain *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[18];

      /* Assignment: '<S59>/Assignment6' */
      windEmulatorStep4_WECSim_B.Assignment6[0] =
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtAssignment[0];

      /* Assignment: '<S59>/Assignment7' */
      windEmulatorStep4_WECSim_B.Assignment7[0] =
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtAssignme_e[0];

      /* Assignment: '<S59>/Assignment6' */
      windEmulatorStep4_WECSim_B.Assignment6[1] =
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtAssignment[1];

      /* Assignment: '<S59>/Assignment7' */
      windEmulatorStep4_WECSim_B.Assignment7[1] =
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtAssignme_e[1];

      /* Assignment: '<S59>/Assignment6' */
      windEmulatorStep4_WECSim_B.Assignment6[2] =
        windEmulatorStep4_WECSim_B.Gain6;

      /* Assignment: '<S59>/Assignment7' */
      windEmulatorStep4_WECSim_B.Assignment7[2] =
        windEmulatorStep4_WECSim_B.Gain7;
      for (i = 0; i < 6; i++) {
        /* Assignment: '<S59>/Assignment2' incorporates:
         *  Constant: '<S59>/Constant2'
         */
        windEmulatorStep4_WECSim_B.acceleration[i] =
          windEmulatorStep4_WECSim_cal->Constant2_Value_m[i];

        /* Assignment: '<S59>/Assignment3' incorporates:
         *  Constant: '<S59>/Constant3'
         */
        windEmulatorStep4_WECSim_B.forceActuation[i] =
          windEmulatorStep4_WECSim_cal->Constant3_Value_m[i];
      }

      /* Assignment: '<S59>/Assignment2' */
      windEmulatorStep4_WECSim_B.acceleration[4] =
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[6];

      /* Assignment: '<S59>/Assignment3' */
      windEmulatorStep4_WECSim_B.forceActuation[4] =
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[7];

      /* SignalConversion generated from: '<S59>/Assignment4' */
      windEmulatorStep4_WECSim_B.TmpSignalConversionAtAssignme_p[0] =
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[8];
      windEmulatorStep4_WECSim_B.TmpSignalConversionAtAssignme_p[1] =
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[10];
      windEmulatorStep4_WECSim_B.TmpSignalConversionAtAssignme_p[2] =
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[9];

      /* Gain: '<S59>/Gain4' */
      windEmulatorStep4_WECSim_B.Gain4_j =
        windEmulatorStep4_WECSim_cal->Gain4_Gain *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[9];

      /* SignalConversion generated from: '<S59>/Assignment5' */
      windEmulatorStep4_WECSim_B.TmpSignalConversionAtAssignme_f[0] =
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[11];
      windEmulatorStep4_WECSim_B.TmpSignalConversionAtAssignme_f[1] =
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[13];
      windEmulatorStep4_WECSim_B.TmpSignalConversionAtAssignme_f[2] =
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[12];

      /* Gain: '<S59>/Gain5' */
      windEmulatorStep4_WECSim_B.Gain5_j =
        windEmulatorStep4_WECSim_cal->Gain5_Gain *
        windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[12];

      /* Assignment: '<S59>/Assignment4' */
      windEmulatorStep4_WECSim_B.Assignment4[0] =
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtAssignme_p[0];

      /* Assignment: '<S59>/Assignment5' */
      windEmulatorStep4_WECSim_B.Assignment5[0] =
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtAssignme_f[0];

      /* Assignment: '<S59>/Assignment4' */
      windEmulatorStep4_WECSim_B.Assignment4[1] =
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtAssignme_p[1];

      /* Assignment: '<S59>/Assignment5' */
      windEmulatorStep4_WECSim_B.Assignment5[1] =
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtAssignme_f[1];

      /* Assignment: '<S59>/Assignment4' */
      windEmulatorStep4_WECSim_B.Assignment4[2] =
        windEmulatorStep4_WECSim_B.Gain4_j;

      /* Assignment: '<S59>/Assignment5' */
      windEmulatorStep4_WECSim_B.Assignment5[2] =
        windEmulatorStep4_WECSim_B.Gain5_j;

      /* Sum: '<S59>/Add' */
      windEmulatorStep4_WECSim_B.forceInternalMechanics[0] =
        (windEmulatorStep4_WECSim_B.Assignment6[0] -
         windEmulatorStep4_WECSim_B.forceActuation[0]) -
        windEmulatorStep4_WECSim_B.Assignment4[0];
      windEmulatorStep4_WECSim_B.forceInternalMechanics[3] =
        (windEmulatorStep4_WECSim_B.Assignment7[0] -
         windEmulatorStep4_WECSim_B.forceActuation[3]) -
        windEmulatorStep4_WECSim_B.Assignment5[0];
      windEmulatorStep4_WECSim_B.forceInternalMechanics[1] =
        (windEmulatorStep4_WECSim_B.Assignment6[1] -
         windEmulatorStep4_WECSim_B.forceActuation[1]) -
        windEmulatorStep4_WECSim_B.Assignment4[1];
      windEmulatorStep4_WECSim_B.forceInternalMechanics[4] =
        (windEmulatorStep4_WECSim_B.Assignment7[1] -
         windEmulatorStep4_WECSim_B.forceActuation[4]) -
        windEmulatorStep4_WECSim_B.Assignment5[1];
      windEmulatorStep4_WECSim_B.forceInternalMechanics[2] =
        (windEmulatorStep4_WECSim_B.Assignment6[2] -
         windEmulatorStep4_WECSim_B.forceActuation[2]) -
        windEmulatorStep4_WECSim_B.Assignment4[2];
      windEmulatorStep4_WECSim_B.forceInternalMechanics[5] =
        (windEmulatorStep4_WECSim_B.Assignment7[2] -
         windEmulatorStep4_WECSim_B.forceActuation[5]) -
        windEmulatorStep4_WECSim_B.Assignment5[2];
      for (i = 0; i < 6; i++) {
        /* Product: '<S59>/Product' */
        windEmulatorStep4_WECSim_B.powerInternalMechanics[i] =
          windEmulatorStep4_WECSim_B.forceInternalMechanics[i] *
          windEmulatorStep4_WECSim_B.velocity[i];
      }

      if (tmp_g) {
        /* SignalConversion generated from: '<S59>/To Workspace' */
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorkspac[18] =
          windEmulatorStep4_WECSim_B.Assignment6[0];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorkspac[21] =
          windEmulatorStep4_WECSim_B.Assignment7[0];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorkspac[19] =
          windEmulatorStep4_WECSim_B.Assignment6[1];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorkspac[22] =
          windEmulatorStep4_WECSim_B.Assignment7[1];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorkspac[20] =
          windEmulatorStep4_WECSim_B.Assignment6[2];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorkspac[23] =
          windEmulatorStep4_WECSim_B.Assignment7[2];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorkspac[30] =
          windEmulatorStep4_WECSim_B.Assignment4[0];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorkspac[33] =
          windEmulatorStep4_WECSim_B.Assignment5[0];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorkspac[31] =
          windEmulatorStep4_WECSim_B.Assignment4[1];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorkspac[34] =
          windEmulatorStep4_WECSim_B.Assignment5[1];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorkspac[32] =
          windEmulatorStep4_WECSim_B.Assignment4[2];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorkspac[35] =
          windEmulatorStep4_WECSim_B.Assignment5[2];
        for (i = 0; i < 6; i++) {
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorkspac[i] =
            windEmulatorStep4_WECSim_B.position[i];
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorkspac[i + 6] =
            windEmulatorStep4_WECSim_B.velocity[i];
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorkspac[i + 12] =
            windEmulatorStep4_WECSim_B.acceleration[i];
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorkspac[i + 24] =
            windEmulatorStep4_WECSim_B.forceActuation[i];
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorkspac[i + 36] =
            windEmulatorStep4_WECSim_B.forceInternalMechanics[i];
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorkspac[i + 42] =
            windEmulatorStep4_WECSim_B.powerInternalMechanics[i];
        }

        /* End of SignalConversion generated from: '<S59>/To Workspace' */
        /* SignalConversion generated from: '<S60>/To Workspace' */
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[0] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_1[9];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[3] =
          windEmulatorStep4_WECSim_B.sf_quaternion2EulXYZ.E[0];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[6] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_1[12];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[9] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_1[6];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[12] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[23];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[15] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[20];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[1] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_1[10];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[4] =
          windEmulatorStep4_WECSim_B.sf_quaternion2EulXYZ.E[1];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[7] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_1[13];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[10] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_1[7];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[13] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[24];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[16] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[21];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[2] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_1[11];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[5] =
          windEmulatorStep4_WECSim_B.sf_quaternion2EulXYZ.E[2];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[8] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_1[14];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[11] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_1[8];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[14] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[25];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[17] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[22];
        for (i = 0; i < 6; i++) {
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[i + 18] =
            windEmulatorStep4_WECSim_B.F_Total[i];
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[i + 24] =
            windEmulatorStep4_WECSim_B.sf_YawForceTransforms.F_Excitation[i];
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[i + 30] =
            windEmulatorStep4_WECSim_B.sf_YawForceTransforms.F_RadiationDamping[i];
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[i + 36] =
            windEmulatorStep4_WECSim_B.sf_YawForceTransforms.F_AddedMass[i];
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[i + 42] =
            windEmulatorStep4_WECSim_B.sf_YawForceTransforms.F_Restoring[i];
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[i + 48] =
            windEmulatorStep4_WECSim_B.F_MorisonAndViscous[i];
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_c[i + 54] =
            windEmulatorStep4_WECSim_B.F_LinearDamping[i];
        }

        /* End of SignalConversion generated from: '<S60>/To Workspace' */
      }

      /* TransportDelay: '<S64>/Transport Delay' */
      for (i = 0; i < 6; i++) {
        riseValLimit = rt_TDelayInterpolate(windEmulatorStep4_WECSim_M->
          Timing.t[0] - windEmulatorStep4_WECSim_cal->TransportDelay_Delay_o,
          windEmulatorStep4_WECSim_DW.TransportDelay_RWORK_l[0],
          static_cast<real_T *>
          (windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_j[i]),
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_o[i + 18],
          &windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_o[i + 12],
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_o[i],
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_o[i + 6],
          windEmulatorStep4_WECSim_cal->TransportDelay_InitOutput_b,false,
          rtmIsMinorTimeStep(windEmulatorStep4_WECSim_M) && ((static_cast<real_T
          *>(windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_j[0]))
          [windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_o[6] +
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_o[18]] ==
          windEmulatorStep4_WECSim_M->Timing.t[0]));
        windEmulatorStep4_WECSim_B.TransportDelay_o[i] = riseValLimit;
      }

      /* End of TransportDelay: '<S64>/Transport Delay' */
      if (tmp_g) {
        /* Outputs for Atomic SubSystem: '<S60>/Nonlinear Wave Elevation' */
        /* SimscapeExecutionBlock: '<S216>/OUTPUT_1_1' */
        windEmul_NonlinearWaveElevation(&windEmulatorStep4_WECSim_B.OUTPUT_1_1[9],
          windEmulatorStep4_WECSim_B.sf_quaternion2EulXYZ.E,
          &windEmulatorStep4_WECSim_B.NonlinearWaveElevation,
          &windEmulatorStep4_WECSim_cal->wind_NonlinearWaveElevation_cal);

        /* End of Outputs for SubSystem: '<S60>/Nonlinear Wave Elevation' */
        for (i = 0; i < 6; i++) {
          /* Constant: '<S55>/Constant' */
          u1 = windEmulatorStep4_WECSim_cal->Constant_Value_gg[i];
          windEmulatorStep4_WECSim_B.position_d[i] = u1;

          /* Constant: '<S55>/Constant1' */
          Clock_tmp = windEmulatorStep4_WECSim_cal->Constant1_Value_h[i];
          windEmulatorStep4_WECSim_B.velocity_g[i] = Clock_tmp;

          /* Constant: '<S55>/Constant2' */
          riseValLimit = windEmulatorStep4_WECSim_cal->Constant2_Value_n[i];
          windEmulatorStep4_WECSim_B.acceleration_b[i] = riseValLimit;

          /* SignalConversion generated from: '<S55>/To Workspace' */
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_m[i] = u1;
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_m[i + 6] =
            Clock_tmp;
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_m[i + 12] =
            riseValLimit;
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorksp_m[i + 18] =
            windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[i];
        }

        /* SignalConversion generated from: '<S139>/To Workspace' */
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[0] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_1[22];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[3] =
          windEmulatorStep4_WECSim_B.sf_quaternion2EulXYZ_c.E[0];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[6] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_1[25];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[9] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_1[19];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[12] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[29];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[15] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[26];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[1] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_1[23];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[4] =
          windEmulatorStep4_WECSim_B.sf_quaternion2EulXYZ_c.E[1];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[7] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_1[26];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[10] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_1[20];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[13] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[30];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[16] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[27];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[2] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_1[24];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[5] =
          windEmulatorStep4_WECSim_B.sf_quaternion2EulXYZ_c.E[2];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[8] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_1[27];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[11] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_1[21];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[14] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[31];
        windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[17] =
          windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[28];
        for (i = 0; i < 6; i++) {
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[i + 18] =
            windEmulatorStep4_WECSim_B.F_Total_c[i];
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[i + 24] =
            windEmulatorStep4_WECSim_B.sf_YawForceTransforms_i.F_Excitation[i];
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[i + 30] =
            windEmulatorStep4_WECSim_B.sf_YawForceTransforms_i.F_RadiationDamping
            [i];
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[i + 36] =
            windEmulatorStep4_WECSim_B.sf_YawForceTransforms_i.F_AddedMass[i];
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[i + 42] =
            windEmulatorStep4_WECSim_B.sf_YawForceTransforms_i.F_Restoring[i];
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[i + 48] =
            windEmulatorStep4_WECSim_B.F_MorisonAndViscous_d[i];
          windEmulatorStep4_WECSim_B.TmpSignalConversionAtToWorks_ms[i + 54] =
            windEmulatorStep4_WECSim_B.F_LinearDamping_k[i];
        }

        /* End of SignalConversion generated from: '<S139>/To Workspace' */
      }

      /* TransportDelay: '<S143>/Transport Delay' */
      for (i = 0; i < 6; i++) {
        riseValLimit = rt_TDelayInterpolate(windEmulatorStep4_WECSim_M->
          Timing.t[0] - windEmulatorStep4_WECSim_cal->TransportDelay_Delay_n,
          windEmulatorStep4_WECSim_DW.TransportDelay_RWORK_f[0],
          static_cast<real_T *>
          (windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_n[i]),
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c1[i + 18],
          &windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c1[i + 12],
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c1[i],
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c1[i + 6],
          windEmulatorStep4_WECSim_cal->TransportDelay_InitOutput_g,false,
          rtmIsMinorTimeStep(windEmulatorStep4_WECSim_M) && ((static_cast<real_T
          *>(windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_n[0]))
          [windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c1[6] +
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c1[18]] ==
          windEmulatorStep4_WECSim_M->Timing.t[0]));
        windEmulatorStep4_WECSim_B.TransportDelay_e[i] = riseValLimit;
      }

      /* End of TransportDelay: '<S143>/Transport Delay' */
      if (tmp_g) {
        /* Outputs for Atomic SubSystem: '<S139>/Nonlinear Wave Elevation' */
        /* SimscapeExecutionBlock: '<S216>/OUTPUT_1_1' */
        windEmul_NonlinearWaveElevation(&windEmulatorStep4_WECSim_B.OUTPUT_1_1
          [22], windEmulatorStep4_WECSim_B.sf_quaternion2EulXYZ_c.E,
          &windEmulatorStep4_WECSim_B.NonlinearWaveElevation_j,
          &windEmulatorStep4_WECSim_cal->wi_NonlinearWaveElevation_j_cal);

        /* End of Outputs for SubSystem: '<S139>/Nonlinear Wave Elevation' */

        /* Memory: '<S433>/Memory' */
        windEmulatorStep4_WECSim_B.Memory_f =
          windEmulatorStep4_WECSim_DW.Memory_PreviousInput_kk;

        /* RateLimiter: '<S373>/Rate Limiter1' incorporates:
         *  Constant: '<S373>/shaftSpeedRefMin'
         */
        rateLimiterRate = windEmulatorStep4_WECSim_cal->shaftSpeedRefMin_Value -
          windEmulatorStep4_WECSim_DW.PrevY_m;
        if (rateLimiterRate >
            windEmulatorStep4_WECSim_cal->RateLimiter1_RisingLim *
            windEmulatorStep4_WECSim_period) {
          /* RateLimiter: '<S373>/Rate Limiter1' */
          windEmulatorStep4_WECSim_B.RateLimiter1 =
            windEmulatorStep4_WECSim_cal->RateLimiter1_RisingLim *
            windEmulatorStep4_WECSim_period +
            windEmulatorStep4_WECSim_DW.PrevY_m;
        } else if (rateLimiterRate <
                   windEmulatorStep4_WECSim_cal->RateLimiter1_FallingLim *
                   windEmulatorStep4_WECSim_period) {
          /* RateLimiter: '<S373>/Rate Limiter1' */
          windEmulatorStep4_WECSim_B.RateLimiter1 =
            windEmulatorStep4_WECSim_cal->RateLimiter1_FallingLim *
            windEmulatorStep4_WECSim_period +
            windEmulatorStep4_WECSim_DW.PrevY_m;
        } else {
          /* RateLimiter: '<S373>/Rate Limiter1' */
          windEmulatorStep4_WECSim_B.RateLimiter1 =
            windEmulatorStep4_WECSim_cal->shaftSpeedRefMin_Value;
        }

        windEmulatorStep4_WECSim_DW.PrevY_m =
          windEmulatorStep4_WECSim_B.RateLimiter1;

        /* End of RateLimiter: '<S373>/Rate Limiter1' */

        /* Memory: '<S434>/Memory' */
        windEmulatorStep4_WECSim_B.Memory_in =
          windEmulatorStep4_WECSim_DW.Memory_PreviousInput_h;
      }

      /* Sum: '<S376>/Add' */
      windEmulatorStep4_WECSim_B.Add_mw =
        windEmulatorStep4_WECSim_B.BusAssignment_c.genSpeedActual -
        windEmulatorStep4_WECSim_B.BusAssignment_c.speedRef_rpm;

      /* Product: '<S431>/Product' incorporates:
       *  Constant: '<S431>/Constant1'
       */
      riseValLimit = -tmp_m->PG;

      /* Product: '<S431>/Product' */
      windEmulatorStep4_WECSim_B.ControlSignal31 = riseValLimit *
        windEmulatorStep4_WECSim_B.Add_mw;

      /* RelationalOperator: '<S431>/Relational Operator' */
      windEmulatorStep4_WECSim_B.RelationalOperator_l =
        (windEmulatorStep4_WECSim_B.BusAssignment_c.speedRef_rpm <=
         windEmulatorStep4_WECSim_B.BusAssignment_c.genSpeedActual);

      /* CombinatorialLogic: '<S433>/Logic' incorporates:
       *  Constant: '<S431>/Constant'
       */
      f = windEmulatorStep4_WECSim_B.RelationalOperator_l;
      q0 = f;
      f = windEmulatorStep4_WECSim_cal->Constant_Value_bl;
      q0 = (q0 << 1) + f;
      f = windEmulatorStep4_WECSim_B.Memory_f;
      q0 = (q0 << 1) + f;
      windEmulatorStep4_WECSim_B.Logic[0U] =
        windEmulatorStep4_WECSim_cal->Logic_table[q0];
      windEmulatorStep4_WECSim_B.Logic[1U] =
        windEmulatorStep4_WECSim_cal->Logic_table[q0 + 8U];

      /* Sum: '<S376>/Add1' */
      windEmulatorStep4_WECSim_B.Add1_p =
        windEmulatorStep4_WECSim_B.BusAssignment_c.genSpeedActual -
        windEmulatorStep4_WECSim_B.RateLimiter1;

      /* Product: '<S432>/Product' incorporates:
       *  Constant: '<S432>/Constant1'
       */
      riseValLimit = -tmp_l;

      /* Product: '<S432>/Product' */
      windEmulatorStep4_WECSim_B.ControlSignal31_o = riseValLimit *
        windEmulatorStep4_WECSim_B.Add1_p;

      /* RelationalOperator: '<S432>/Relational Operator' */
      windEmulatorStep4_WECSim_B.RelationalOperator_k =
        (windEmulatorStep4_WECSim_B.RateLimiter1 >=
         windEmulatorStep4_WECSim_B.BusAssignment_c.genSpeedActual);

      /* CombinatorialLogic: '<S434>/Logic' incorporates:
       *  Constant: '<S432>/Constant'
       */
      f = windEmulatorStep4_WECSim_B.RelationalOperator_k;
      q0 = f;
      f = windEmulatorStep4_WECSim_cal->Constant_Value_k3;
      q0 = (q0 << 1) + f;
      f = windEmulatorStep4_WECSim_B.Memory_in;
      q0 = (q0 << 1) + f;
      windEmulatorStep4_WECSim_B.Logic_g[0U] =
        windEmulatorStep4_WECSim_cal->Logic_table_o[q0];
      windEmulatorStep4_WECSim_B.Logic_g[1U] =
        windEmulatorStep4_WECSim_cal->Logic_table_o[q0 + 8U];

      /* Switch: '<S376>/Switch' incorporates:
       *  Switch: '<S376>/Switch1'
       *  Switch: '<S432>/Switch'
       */
      if (windEmulatorStep4_WECSim_B.Add_mw >
          windEmulatorStep4_WECSim_cal->Switch_Threshold_j) {
        /* Switch: '<S431>/Switch' */
        if (windEmulatorStep4_WECSim_B.Logic[0]) {
          /* Saturate: '<S431>/Saturation' */
          riseValLimit = windEmulatorStep4_WECSim_B.ControlSignal31;
          u1 = windEmulatorStep4_WECSim_cal->Saturation_LowerSat;
          rateLimiterRate = windEmulatorStep4_WECSim_cal->Saturation_UpperSat;
          if (riseValLimit > rateLimiterRate) {
            /* Saturate: '<S431>/Saturation' */
            windEmulatorStep4_WECSim_B.Saturation_a3 = rateLimiterRate;
          } else if (riseValLimit < u1) {
            /* Saturate: '<S431>/Saturation' */
            windEmulatorStep4_WECSim_B.Saturation_a3 = u1;
          } else {
            /* Saturate: '<S431>/Saturation' */
            windEmulatorStep4_WECSim_B.Saturation_a3 = riseValLimit;
          }

          /* End of Saturate: '<S431>/Saturation' */

          /* Switch: '<S431>/Switch' */
          windEmulatorStep4_WECSim_B.ControlSignal3_f =
            windEmulatorStep4_WECSim_B.Saturation_a3;
        } else {
          /* Switch: '<S431>/Switch' */
          windEmulatorStep4_WECSim_B.ControlSignal3_f =
            windEmulatorStep4_WECSim_B.ControlSignal31;
        }

        /* End of Switch: '<S431>/Switch' */

        /* Switch: '<S376>/Switch' */
        windEmulatorStep4_WECSim_B.Switch_g =
          windEmulatorStep4_WECSim_B.ControlSignal3_f;
      } else {
        if (windEmulatorStep4_WECSim_B.Add1_p >
            windEmulatorStep4_WECSim_cal->Switch1_Threshold) {
          /* Switch: '<S376>/Switch1' incorporates:
           *  Constant: '<S376>/Constant1'
           */
          windEmulatorStep4_WECSim_B.Switch1_o =
            windEmulatorStep4_WECSim_cal->Constant1_Value_c;
        } else {
          if (windEmulatorStep4_WECSim_B.Logic_g[0]) {
            /* Saturate: '<S432>/Saturation' incorporates:
             *  Switch: '<S376>/Switch1'
             *  Switch: '<S432>/Switch'
             */
            riseValLimit = windEmulatorStep4_WECSim_B.ControlSignal31_o;
            u1 = windEmulatorStep4_WECSim_cal->Saturation_LowerSat_n;
            rateLimiterRate =
              windEmulatorStep4_WECSim_cal->Saturation_UpperSat_g;
            if (riseValLimit > rateLimiterRate) {
              /* Saturate: '<S432>/Saturation' */
              windEmulatorStep4_WECSim_B.Saturation_k = rateLimiterRate;
            } else if (riseValLimit < u1) {
              /* Saturate: '<S432>/Saturation' */
              windEmulatorStep4_WECSim_B.Saturation_k = u1;
            } else {
              /* Saturate: '<S432>/Saturation' */
              windEmulatorStep4_WECSim_B.Saturation_k = riseValLimit;
            }

            /* End of Saturate: '<S432>/Saturation' */

            /* Switch: '<S432>/Switch' incorporates:
             *  Switch: '<S376>/Switch1'
             */
            windEmulatorStep4_WECSim_B.ControlSignal3_e =
              windEmulatorStep4_WECSim_B.Saturation_k;
          } else {
            /* Switch: '<S432>/Switch' incorporates:
             *  Switch: '<S376>/Switch1'
             */
            windEmulatorStep4_WECSim_B.ControlSignal3_e =
              windEmulatorStep4_WECSim_B.ControlSignal31_o;
          }

          /* Switch: '<S376>/Switch1' */
          windEmulatorStep4_WECSim_B.Switch1_o =
            windEmulatorStep4_WECSim_B.ControlSignal3_e;
        }

        /* Switch: '<S376>/Switch' incorporates:
         *  Switch: '<S376>/Switch1'
         *  Switch: '<S432>/Switch'
         */
        windEmulatorStep4_WECSim_B.Switch_g =
          windEmulatorStep4_WECSim_B.Switch1_o;
      }

      /* End of Switch: '<S376>/Switch' */

      /* Gain: '<S376>/Gain' */
      windEmulatorStep4_WECSim_B.Gain_lr = tmp_x *
        windEmulatorStep4_WECSim_B.Switch_g;

      /* RateLimiter: '<S373>/Rate Limiter' */
      if (windEmulatorStep4_WECSim_DW.LastMajorTime_d == (rtInf)) {
        /* RateLimiter: '<S373>/Rate Limiter' */
        windEmulatorStep4_WECSim_B.RateLimiter_b =
          windEmulatorStep4_WECSim_B.Gain_lr;
      } else {
        u1 = deltaT_tmp - windEmulatorStep4_WECSim_DW.LastMajorTime_d;
        if (windEmulatorStep4_WECSim_DW.LastMajorTime_d == deltaT_tmp) {
          if (windEmulatorStep4_WECSim_DW.PrevLimited_g) {
            /* RateLimiter: '<S373>/Rate Limiter' */
            windEmulatorStep4_WECSim_B.RateLimiter_b =
              windEmulatorStep4_WECSim_DW.PrevY_l;
          } else {
            /* RateLimiter: '<S373>/Rate Limiter' */
            windEmulatorStep4_WECSim_B.RateLimiter_b =
              windEmulatorStep4_WECSim_B.Gain_lr;
          }
        } else {
          riseValLimit = u1 * tmp_k;
          rateLimiterRate = windEmulatorStep4_WECSim_B.Gain_lr -
            windEmulatorStep4_WECSim_DW.PrevY_l;
          if (rateLimiterRate > riseValLimit) {
            /* RateLimiter: '<S373>/Rate Limiter' */
            windEmulatorStep4_WECSim_B.RateLimiter_b =
              windEmulatorStep4_WECSim_DW.PrevY_l + riseValLimit;
            f = true;
          } else {
            riseValLimit = -tmp_k;
            u1 *= riseValLimit;
            if (rateLimiterRate < u1) {
              /* RateLimiter: '<S373>/Rate Limiter' */
              windEmulatorStep4_WECSim_B.RateLimiter_b =
                windEmulatorStep4_WECSim_DW.PrevY_l + u1;
              f = true;
            } else {
              /* RateLimiter: '<S373>/Rate Limiter' */
              windEmulatorStep4_WECSim_B.RateLimiter_b =
                windEmulatorStep4_WECSim_B.Gain_lr;
              f = false;
            }
          }

          if (tmp_h) {
            windEmulatorStep4_WECSim_DW.PrevLimited_g = f;
          }
        }
      }

      /* Sum: '<S375>/Sum' */
      windEmulatorStep4_WECSim_B.wError =
        windEmulatorStep4_WECSim_B.BusAssignment_c.genSpeedActual -
        windEmulatorStep4_WECSim_B.BusAssignment_c.speedRef_rpm;
      if (tmp_g) {
        /* Gain: '<S417>/Proportional Gain' */
        windEmulatorStep4_WECSim_B.ProportionalGain_p = tmp_m->PG *
          windEmulatorStep4_WECSim_B.wError;

        /* DiscreteIntegrator: '<S412>/Integrator' */
        if (windEmulatorStep4_WECSim_B.BusAssignment_c.speedCtrlReset ||
            (windEmulatorStep4_WECSim_DW.Integrator_PrevResetState != 0)) {
          windEmulatorStep4_WECSim_DW.Integrator_DSTATE_d =
            windEmulatorStep4_WECSim_cal->PIDController_InitialConditio_a;
        }

        /* DiscreteIntegrator: '<S412>/Integrator' */
        windEmulatorStep4_WECSim_B.Integrator_b =
          windEmulatorStep4_WECSim_DW.Integrator_DSTATE_d;

        /* Gain: '<S405>/Derivative Gain' */
        windEmulatorStep4_WECSim_B.DerivativeGain_a =
          windEmulatorStep4_WECSim_cal->PIDController_D *
          windEmulatorStep4_WECSim_B.wError;

        /* DiscreteIntegrator: '<S407>/Filter' */
        if (windEmulatorStep4_WECSim_B.BusAssignment_c.speedCtrlReset ||
            (windEmulatorStep4_WECSim_DW.Filter_PrevResetState != 0)) {
          windEmulatorStep4_WECSim_DW.Filter_DSTATE =
            windEmulatorStep4_WECSim_cal->PIDController_InitialConditionF;
        }

        /* DiscreteIntegrator: '<S407>/Filter' */
        windEmulatorStep4_WECSim_B.Filter =
          windEmulatorStep4_WECSim_DW.Filter_DSTATE;

        /* Sum: '<S407>/SumD' */
        windEmulatorStep4_WECSim_B.SumD =
          windEmulatorStep4_WECSim_B.DerivativeGain_a -
          windEmulatorStep4_WECSim_B.Filter;

        /* Gain: '<S415>/Filter Coefficient' */
        windEmulatorStep4_WECSim_B.FilterCoefficient =
          windEmulatorStep4_WECSim_cal->PIDController_N *
          windEmulatorStep4_WECSim_B.SumD;

        /* Sum: '<S422>/Sum' */
        windEmulatorStep4_WECSim_B.Sum_c =
          (windEmulatorStep4_WECSim_B.ProportionalGain_p +
           windEmulatorStep4_WECSim_B.Integrator_b) +
          windEmulatorStep4_WECSim_B.FilterCoefficient;

        /* RelationalOperator: '<S420>/LowerRelop1' incorporates:
         *  Constant: '<S375>/Constant'
         */
        windEmulatorStep4_WECSim_B.LowerRelop1_g =
          (windEmulatorStep4_WECSim_B.Sum_c > tmp_j);

        /* RelationalOperator: '<S420>/UpperRelop' incorporates:
         *  Constant: '<S375>/Constant1'
         */
        riseValLimit = -tmp_j;

        /* RelationalOperator: '<S420>/UpperRelop' */
        windEmulatorStep4_WECSim_B.UpperRelop_g =
          (windEmulatorStep4_WECSim_B.Sum_c < riseValLimit);

        /* Switch: '<S420>/Switch' */
        if (windEmulatorStep4_WECSim_B.UpperRelop_g) {
          /* Switch: '<S420>/Switch' incorporates:
           *  Constant: '<S375>/Constant1'
           */
          windEmulatorStep4_WECSim_B.Switch_i = -tmp_j;
        } else {
          /* Switch: '<S420>/Switch' */
          windEmulatorStep4_WECSim_B.Switch_i = windEmulatorStep4_WECSim_B.Sum_c;
        }

        /* End of Switch: '<S420>/Switch' */

        /* Switch: '<S420>/Switch2' */
        if (windEmulatorStep4_WECSim_B.LowerRelop1_g) {
          /* Switch: '<S420>/Switch2' incorporates:
           *  Constant: '<S375>/Constant'
           */
          windEmulatorStep4_WECSim_B.Switch2_k = tmp_j;
        } else {
          /* Switch: '<S420>/Switch2' */
          windEmulatorStep4_WECSim_B.Switch2_k =
            windEmulatorStep4_WECSim_B.Switch_i;
        }

        /* End of Switch: '<S420>/Switch2' */

        /* Gain: '<S375>/Gain2' */
        windEmulatorStep4_WECSim_B.ContolTorque =
          windEmulatorStep4_WECSim_cal->Gain2_Gain_o *
          windEmulatorStep4_WECSim_B.Switch2_k;

        /* Gain: '<S409>/Integral Gain' */
        windEmulatorStep4_WECSim_B.IntegralGain_h = tmp_m->IG *
          windEmulatorStep4_WECSim_B.wError;

        /* Memory: '<S496>/Memory' */
        windEmulatorStep4_WECSim_B.Memory_o =
          windEmulatorStep4_WECSim_DW.Memory_PreviousInput_g;

        /* RateLimiter: '<S436>/Rate Limiter1' incorporates:
         *  Constant: '<S436>/shaftSpeedRefMin'
         */
        rateLimiterRate = windEmulatorStep4_WECSim_cal->shaftSpeedRefMin_Value_k
          - windEmulatorStep4_WECSim_DW.PrevY_fq;
        if (rateLimiterRate >
            windEmulatorStep4_WECSim_cal->RateLimiter1_RisingLim_o *
            windEmulatorStep4_WECSim_period) {
          /* RateLimiter: '<S436>/Rate Limiter1' */
          windEmulatorStep4_WECSim_B.RateLimiter1_m =
            windEmulatorStep4_WECSim_cal->RateLimiter1_RisingLim_o *
            windEmulatorStep4_WECSim_period +
            windEmulatorStep4_WECSim_DW.PrevY_fq;
        } else if (rateLimiterRate <
                   windEmulatorStep4_WECSim_cal->RateLimiter1_FallingLim_e *
                   windEmulatorStep4_WECSim_period) {
          /* RateLimiter: '<S436>/Rate Limiter1' */
          windEmulatorStep4_WECSim_B.RateLimiter1_m =
            windEmulatorStep4_WECSim_cal->RateLimiter1_FallingLim_e *
            windEmulatorStep4_WECSim_period +
            windEmulatorStep4_WECSim_DW.PrevY_fq;
        } else {
          /* RateLimiter: '<S436>/Rate Limiter1' */
          windEmulatorStep4_WECSim_B.RateLimiter1_m =
            windEmulatorStep4_WECSim_cal->shaftSpeedRefMin_Value_k;
        }

        windEmulatorStep4_WECSim_DW.PrevY_fq =
          windEmulatorStep4_WECSim_B.RateLimiter1_m;

        /* End of RateLimiter: '<S436>/Rate Limiter1' */

        /* Memory: '<S497>/Memory' */
        windEmulatorStep4_WECSim_B.Memory_a =
          windEmulatorStep4_WECSim_DW.Memory_PreviousInput_n;
      }

      /* Sum: '<S439>/Add' */
      windEmulatorStep4_WECSim_B.Add_mo =
        windEmulatorStep4_WECSim_B.BusAssignment_c.genSpeedActual -
        windEmulatorStep4_WECSim_B.BusAssignment_c.speedRef_rpm;

      /* Product: '<S494>/Product' incorporates:
       *  Constant: '<S494>/Constant1'
       */
      riseValLimit = -tmp_m->PG;

      /* Product: '<S494>/Product' */
      windEmulatorStep4_WECSim_B.ControlSignal31_d = riseValLimit *
        windEmulatorStep4_WECSim_B.Add_mo;

      /* RelationalOperator: '<S494>/Relational Operator' */
      windEmulatorStep4_WECSim_B.RelationalOperator_e =
        (windEmulatorStep4_WECSim_B.BusAssignment_c.speedRef_rpm <=
         windEmulatorStep4_WECSim_B.BusAssignment_c.genSpeedActual);

      /* CombinatorialLogic: '<S496>/Logic' incorporates:
       *  Constant: '<S494>/Constant'
       */
      f = windEmulatorStep4_WECSim_B.RelationalOperator_e;
      q0 = f;
      f = windEmulatorStep4_WECSim_cal->Constant_Value_cl;
      q0 = (q0 << 1) + f;
      f = windEmulatorStep4_WECSim_B.Memory_o;
      q0 = (q0 << 1) + f;
      windEmulatorStep4_WECSim_B.Logic_c[0U] =
        windEmulatorStep4_WECSim_cal->Logic_table_h[q0];
      windEmulatorStep4_WECSim_B.Logic_c[1U] =
        windEmulatorStep4_WECSim_cal->Logic_table_h[q0 + 8U];

      /* Sum: '<S439>/Add1' */
      windEmulatorStep4_WECSim_B.Add1_f =
        windEmulatorStep4_WECSim_B.BusAssignment_c.genSpeedActual -
        windEmulatorStep4_WECSim_B.RateLimiter1_m;

      /* Product: '<S495>/Product' incorporates:
       *  Constant: '<S495>/Constant1'
       */
      riseValLimit = -tmp_l;

      /* Product: '<S495>/Product' */
      windEmulatorStep4_WECSim_B.ControlSignal31_m = riseValLimit *
        windEmulatorStep4_WECSim_B.Add1_f;

      /* RelationalOperator: '<S495>/Relational Operator' */
      windEmulatorStep4_WECSim_B.RelationalOperator_g =
        (windEmulatorStep4_WECSim_B.RateLimiter1_m >=
         windEmulatorStep4_WECSim_B.BusAssignment_c.genSpeedActual);

      /* CombinatorialLogic: '<S497>/Logic' incorporates:
       *  Constant: '<S495>/Constant'
       */
      f = windEmulatorStep4_WECSim_B.RelationalOperator_g;
      q0 = f;
      f = windEmulatorStep4_WECSim_cal->Constant_Value_ks;
      q0 = (q0 << 1) + f;
      f = windEmulatorStep4_WECSim_B.Memory_a;
      q0 = (q0 << 1) + f;
      windEmulatorStep4_WECSim_B.Logic_p[0U] =
        windEmulatorStep4_WECSim_cal->Logic_table_n[q0];
      windEmulatorStep4_WECSim_B.Logic_p[1U] =
        windEmulatorStep4_WECSim_cal->Logic_table_n[q0 + 8U];

      /* Switch: '<S439>/Switch' incorporates:
       *  Switch: '<S439>/Switch1'
       *  Switch: '<S495>/Switch'
       */
      if (windEmulatorStep4_WECSim_B.Add_mo >
          windEmulatorStep4_WECSim_cal->Switch_Threshold_k) {
        /* Switch: '<S494>/Switch' */
        if (windEmulatorStep4_WECSim_B.Logic_c[0]) {
          /* Saturate: '<S494>/Saturation' */
          riseValLimit = windEmulatorStep4_WECSim_B.ControlSignal31_d;
          u1 = windEmulatorStep4_WECSim_cal->Saturation_LowerSat_f;
          rateLimiterRate = windEmulatorStep4_WECSim_cal->Saturation_UpperSat_f;
          if (riseValLimit > rateLimiterRate) {
            /* Saturate: '<S494>/Saturation' */
            windEmulatorStep4_WECSim_B.Saturation_af = rateLimiterRate;
          } else if (riseValLimit < u1) {
            /* Saturate: '<S494>/Saturation' */
            windEmulatorStep4_WECSim_B.Saturation_af = u1;
          } else {
            /* Saturate: '<S494>/Saturation' */
            windEmulatorStep4_WECSim_B.Saturation_af = riseValLimit;
          }

          /* End of Saturate: '<S494>/Saturation' */

          /* Switch: '<S494>/Switch' */
          windEmulatorStep4_WECSim_B.ControlSignal3_h =
            windEmulatorStep4_WECSim_B.Saturation_af;
        } else {
          /* Switch: '<S494>/Switch' */
          windEmulatorStep4_WECSim_B.ControlSignal3_h =
            windEmulatorStep4_WECSim_B.ControlSignal31_d;
        }

        /* End of Switch: '<S494>/Switch' */

        /* Switch: '<S439>/Switch' */
        windEmulatorStep4_WECSim_B.Switch_n =
          windEmulatorStep4_WECSim_B.ControlSignal3_h;
      } else {
        if (windEmulatorStep4_WECSim_B.Add1_f >
            windEmulatorStep4_WECSim_cal->Switch1_Threshold_k) {
          /* Switch: '<S439>/Switch1' incorporates:
           *  Constant: '<S439>/Constant1'
           */
          windEmulatorStep4_WECSim_B.Switch1_g =
            windEmulatorStep4_WECSim_cal->Constant1_Value_j;
        } else {
          if (windEmulatorStep4_WECSim_B.Logic_p[0]) {
            /* Saturate: '<S495>/Saturation' incorporates:
             *  Switch: '<S439>/Switch1'
             *  Switch: '<S495>/Switch'
             */
            riseValLimit = windEmulatorStep4_WECSim_B.ControlSignal31_m;
            u1 = windEmulatorStep4_WECSim_cal->Saturation_LowerSat_h;
            rateLimiterRate =
              windEmulatorStep4_WECSim_cal->Saturation_UpperSat_m;
            if (riseValLimit > rateLimiterRate) {
              /* Saturate: '<S495>/Saturation' */
              windEmulatorStep4_WECSim_B.Saturation_p = rateLimiterRate;
            } else if (riseValLimit < u1) {
              /* Saturate: '<S495>/Saturation' */
              windEmulatorStep4_WECSim_B.Saturation_p = u1;
            } else {
              /* Saturate: '<S495>/Saturation' */
              windEmulatorStep4_WECSim_B.Saturation_p = riseValLimit;
            }

            /* End of Saturate: '<S495>/Saturation' */

            /* Switch: '<S495>/Switch' incorporates:
             *  Switch: '<S439>/Switch1'
             */
            windEmulatorStep4_WECSim_B.ControlSignal3 =
              windEmulatorStep4_WECSim_B.Saturation_p;
          } else {
            /* Switch: '<S495>/Switch' incorporates:
             *  Switch: '<S439>/Switch1'
             */
            windEmulatorStep4_WECSim_B.ControlSignal3 =
              windEmulatorStep4_WECSim_B.ControlSignal31_m;
          }

          /* Switch: '<S439>/Switch1' */
          windEmulatorStep4_WECSim_B.Switch1_g =
            windEmulatorStep4_WECSim_B.ControlSignal3;
        }

        /* Switch: '<S439>/Switch' incorporates:
         *  Switch: '<S439>/Switch1'
         *  Switch: '<S495>/Switch'
         */
        windEmulatorStep4_WECSim_B.Switch_n =
          windEmulatorStep4_WECSim_B.Switch1_g;
      }

      /* End of Switch: '<S439>/Switch' */

      /* Gain: '<S439>/Gain' */
      windEmulatorStep4_WECSim_B.Gain_fp = tmp_x *
        windEmulatorStep4_WECSim_B.Switch_n;

      /* RateLimiter: '<S436>/Rate Limiter' */
      if (windEmulatorStep4_WECSim_DW.LastMajorTime_k == (rtInf)) {
        /* RateLimiter: '<S436>/Rate Limiter' */
        windEmulatorStep4_WECSim_B.RateLimiter_a =
          windEmulatorStep4_WECSim_B.Gain_fp;
      } else {
        u1 = deltaT_tmp - windEmulatorStep4_WECSim_DW.LastMajorTime_k;
        if (windEmulatorStep4_WECSim_DW.LastMajorTime_k == deltaT_tmp) {
          if (windEmulatorStep4_WECSim_DW.PrevLimited_dz) {
            /* RateLimiter: '<S436>/Rate Limiter' */
            windEmulatorStep4_WECSim_B.RateLimiter_a =
              windEmulatorStep4_WECSim_DW.PrevY_e;
          } else {
            /* RateLimiter: '<S436>/Rate Limiter' */
            windEmulatorStep4_WECSim_B.RateLimiter_a =
              windEmulatorStep4_WECSim_B.Gain_fp;
          }
        } else {
          riseValLimit = u1 * tmp_k;
          rateLimiterRate = windEmulatorStep4_WECSim_B.Gain_fp -
            windEmulatorStep4_WECSim_DW.PrevY_e;
          if (rateLimiterRate > riseValLimit) {
            /* RateLimiter: '<S436>/Rate Limiter' */
            windEmulatorStep4_WECSim_B.RateLimiter_a =
              windEmulatorStep4_WECSim_DW.PrevY_e + riseValLimit;
            f = true;
          } else {
            riseValLimit = -tmp_k;
            u1 *= riseValLimit;
            if (rateLimiterRate < u1) {
              /* RateLimiter: '<S436>/Rate Limiter' */
              windEmulatorStep4_WECSim_B.RateLimiter_a =
                windEmulatorStep4_WECSim_DW.PrevY_e + u1;
              f = true;
            } else {
              /* RateLimiter: '<S436>/Rate Limiter' */
              windEmulatorStep4_WECSim_B.RateLimiter_a =
                windEmulatorStep4_WECSim_B.Gain_fp;
              f = false;
            }
          }

          if (tmp_h) {
            windEmulatorStep4_WECSim_DW.PrevLimited_dz = f;
          }
        }
      }

      /* Sum: '<S438>/Sum' */
      windEmulatorStep4_WECSim_B.wError_c =
        windEmulatorStep4_WECSim_B.BusAssignment_c.genSpeedActual -
        windEmulatorStep4_WECSim_B.BusAssignment_c.speedRef_rpm;
      if (tmp_g) {
        /* Gain: '<S480>/Proportional Gain' */
        windEmulatorStep4_WECSim_B.ProportionalGain_h = tmp_m->PG *
          windEmulatorStep4_WECSim_B.wError_c;

        /* DiscreteIntegrator: '<S475>/Integrator' */
        if (windEmulatorStep4_WECSim_B.BusAssignment_c.speedCtrlReset ||
            (windEmulatorStep4_WECSim_DW.Integrator_PrevResetState_g != 0)) {
          windEmulatorStep4_WECSim_DW.Integrator_DSTATE_e =
            windEmulatorStep4_WECSim_cal->PIDController_InitialConditio_c;
        }

        /* DiscreteIntegrator: '<S475>/Integrator' */
        windEmulatorStep4_WECSim_B.Integrator_l =
          windEmulatorStep4_WECSim_DW.Integrator_DSTATE_e;

        /* Gain: '<S468>/Derivative Gain' */
        windEmulatorStep4_WECSim_B.DerivativeGain_g =
          windEmulatorStep4_WECSim_cal->PIDController_D_d *
          windEmulatorStep4_WECSim_B.wError_c;

        /* DiscreteIntegrator: '<S470>/Filter' */
        if (windEmulatorStep4_WECSim_B.BusAssignment_c.speedCtrlReset ||
            (windEmulatorStep4_WECSim_DW.Filter_PrevResetState_g != 0)) {
          windEmulatorStep4_WECSim_DW.Filter_DSTATE_b =
            windEmulatorStep4_WECSim_cal->PIDController_InitialConditio_k;
        }

        /* DiscreteIntegrator: '<S470>/Filter' */
        windEmulatorStep4_WECSim_B.Filter_j =
          windEmulatorStep4_WECSim_DW.Filter_DSTATE_b;

        /* Sum: '<S470>/SumD' */
        windEmulatorStep4_WECSim_B.SumD_c =
          windEmulatorStep4_WECSim_B.DerivativeGain_g -
          windEmulatorStep4_WECSim_B.Filter_j;

        /* Gain: '<S478>/Filter Coefficient' */
        windEmulatorStep4_WECSim_B.FilterCoefficient_g =
          windEmulatorStep4_WECSim_cal->PIDController_N_p *
          windEmulatorStep4_WECSim_B.SumD_c;

        /* Sum: '<S485>/Sum' */
        windEmulatorStep4_WECSim_B.Sum_n =
          (windEmulatorStep4_WECSim_B.ProportionalGain_h +
           windEmulatorStep4_WECSim_B.Integrator_l) +
          windEmulatorStep4_WECSim_B.FilterCoefficient_g;

        /* RelationalOperator: '<S483>/LowerRelop1' incorporates:
         *  Constant: '<S438>/Constant'
         */
        windEmulatorStep4_WECSim_B.LowerRelop1_h =
          (windEmulatorStep4_WECSim_B.Sum_n > tmp_j);

        /* RelationalOperator: '<S483>/UpperRelop' incorporates:
         *  Constant: '<S438>/Constant1'
         */
        riseValLimit = -tmp_j;

        /* RelationalOperator: '<S483>/UpperRelop' */
        windEmulatorStep4_WECSim_B.UpperRelop_m =
          (windEmulatorStep4_WECSim_B.Sum_n < riseValLimit);

        /* Switch: '<S483>/Switch' */
        if (windEmulatorStep4_WECSim_B.UpperRelop_m) {
          /* Switch: '<S483>/Switch' incorporates:
           *  Constant: '<S438>/Constant1'
           */
          windEmulatorStep4_WECSim_B.Switch_c = -tmp_j;
        } else {
          /* Switch: '<S483>/Switch' */
          windEmulatorStep4_WECSim_B.Switch_c = windEmulatorStep4_WECSim_B.Sum_n;
        }

        /* Switch: '<S483>/Switch2' */
        if (windEmulatorStep4_WECSim_B.LowerRelop1_h) {
          /* Switch: '<S483>/Switch2' incorporates:
           *  Constant: '<S438>/Constant'
           */
          windEmulatorStep4_WECSim_B.Switch2_m = tmp_j;
        } else {
          /* Switch: '<S483>/Switch2' */
          windEmulatorStep4_WECSim_B.Switch2_m =
            windEmulatorStep4_WECSim_B.Switch_c;
        }

        /* End of Switch: '<S483>/Switch2' */

        /* Gain: '<S438>/Gain2' */
        windEmulatorStep4_WECSim_B.ContolTorque_f =
          windEmulatorStep4_WECSim_cal->Gain2_Gain_e *
          windEmulatorStep4_WECSim_B.Switch2_m;

        /* Gain: '<S472>/Integral Gain' */
        windEmulatorStep4_WECSim_B.IntegralGain_hk = tmp_m->IG *
          windEmulatorStep4_WECSim_B.wError_c;
      }

      /* Switch generated from: '<S365>/Switch' incorporates:
       *  Constant: '<S436>/DeadBandController'
       *  Switch: '<S365>/Switch1'
       *  Switch: '<S436>/Switch'
       */
      if (windEmulatorStep4_WECSim_B.Abs >= tmp_v) {
        /* Switch: '<S373>/Switch' incorporates:
         *  Constant: '<S373>/DeadBandController'
         */
        if (tmp_u) {
          /* Switch: '<S373>/Switch' */
          windEmulatorStep4_WECSim_B.Switch_cn =
            windEmulatorStep4_WECSim_B.RateLimiter_b;
        } else {
          /* Switch: '<S373>/Switch' */
          windEmulatorStep4_WECSim_B.Switch_cn =
            windEmulatorStep4_WECSim_B.ContolTorque;
        }

        /* Switch generated from: '<S365>/Switch' */
        windEmulatorStep4_WECSim_B.ControlTorqueLoad =
          windEmulatorStep4_WECSim_B.Switch_cn;

        /* Switch: '<S365>/Switch1' incorporates:
         *  Constant: '<S365>/Constant'
         */
        windEmulatorStep4_WECSim_B.SwitchLogic =
          windEmulatorStep4_WECSim_cal->Constant_Value_ir;
      } else {
        if (tmp_u) {
          /* Switch: '<S436>/Switch' */
          windEmulatorStep4_WECSim_B.Switch_f =
            windEmulatorStep4_WECSim_B.RateLimiter_a;
        } else {
          /* Switch: '<S436>/Switch' */
          windEmulatorStep4_WECSim_B.Switch_f =
            windEmulatorStep4_WECSim_B.ContolTorque_f;
        }

        /* Switch generated from: '<S365>/Switch' */
        windEmulatorStep4_WECSim_B.ControlTorqueLoad =
          windEmulatorStep4_WECSim_B.Switch_f;

        /* Switch: '<S365>/Switch1' incorporates:
         *  Constant: '<S365>/Constant1'
         */
        windEmulatorStep4_WECSim_B.SwitchLogic =
          windEmulatorStep4_WECSim_cal->Constant1_Value_p;
      }

      if (tmp_g) {
        /* Gain: '<S499>/m3toL' */
        windEmulatorStep4_WECSim_B.FlowPump1 =
          windEmulatorStep4_WECSim_cal->m3toL_Gain *
          windEmulatorStep4_WECSim_B.OUTPUT_1_0[0];

        /* Sum: '<S522>/Sum' */
        windEmulatorStep4_WECSim_B.Sum_a =
          windEmulatorStep4_WECSim_B.BusAssignment_c.excForce_N +
          windEmulatorStep4_WECSim_B.OUTPUT_1_0[4];

        /* Gain: '<S500>/Gain' */
        windEmulatorStep4_WECSim_B.FlowAccumulator =
          windEmulatorStep4_WECSim_cal->Gain_Gain_p *
          windEmulatorStep4_WECSim_B.OUTPUT_1_0[1];

        /* DataTypeConversion: '<S609>/vecIndex' */
        windEmulatorStep4_WECSim_B.vecIndex = windEmulatorStep4_WECSim_B.Mod;

        /* DataTypeConversion: '<S609>/Cast To Double' */
        windEmulatorStep4_WECSim_B.CastToDouble_g =
          windEmulatorStep4_WECSim_B.vecIndex;

        /* DataTypeConversion: '<S609>/Cast To Double1' incorporates:
         *  Constant: '<S609>/Length of input'
         */
        windEmulatorStep4_WECSim_B.CastToDouble1_i =
          windEmulatorStep4_WECSim_cal->Lengthofinput_Value;

        /* Product: '<S609>/Divide1' */
        windEmulatorStep4_WECSim_B.Divide1 =
          windEmulatorStep4_WECSim_B.CastToDouble_g /
          windEmulatorStep4_WECSim_B.CastToDouble1_i;

        /* Bias: '<S609>/fileSamples' incorporates:
         *  Constant: '<S609>/Length of input'
         */
        windEmulatorStep4_WECSim_B.fileSamples =
          windEmulatorStep4_WECSim_cal->Lengthofinput_Value +
          windEmulatorStep4_WECSim_cal->fileSamples_Bias;

        /* Gain: '<S609>/vecPercent' */
        windEmulatorStep4_WECSim_B.vecPercent =
          windEmulatorStep4_WECSim_cal->vecPercent_Gain *
          windEmulatorStep4_WECSim_B.Divide1;
      }

      /* Product: '<S587>/IProd Out' incorporates:
       *  Constant: '<S13>/acs880SpeedIGain'
       */
      windEmulatorStep4_WECSim_B.IProdOut = windEmulatorStep4_WECSim_B.Sum *
        windEmulatorStep4_WECSim_cal->acs880SpeedIGain_Value;

      /* user code (Output function Trailer) */
      {
        /*------------ S-Function Block: <Root>/EtherCAT Init Write Process Data ,Run Admin Tasks and then Write Acyclic Data------------*/
        xpcEtherCATWriteProcessData(0,NULL);
        xpcEtherCATExecAdminJobs(0);
        xpcEtherCATWriteAcyclicData(0);
      }
    }
  }

  if (rtmIsMajorTimeStep(windEmulatorStep4_WECSim_M)) {
    NeslSimulationData *simulationData;
    NeslSimulator *simulator;
    NeuDiagnosticManager *diagnosticManager;
    NeuDiagnosticTree *diagnosticTree;
    char *msg;
    real_T tmp_3[52];
    real_T tmp_0[20];
    real_T tmp_6[12];
    real_T tmp_9[6];
    real_T LastMajorTime_tmp;
    real_T OUTPUT_1_0_k;
    real_T OUTPUT_1_0_k_0;
    real_T OUTPUT_1_0_k_1;
    real_T OUTPUT_1_0_k_2;
    real_T OUTPUT_1_0_k_3;
    real_T OUTPUT_1_0_k_4;
    real_T time;
    real_T time_0;
    real_T time_1;
    real_T time_tmp;
    int_T tmp_4[14];
    int_T tmp_1[6];
    int_T tmp_7[4];
    int_T TransportDelay_IWORK;
    int_T TransportDelay_IWORK_0;
    int_T TransportDelay_IWORK_1;
    int_T i;
    boolean_T tmp;
    boolean_T tmp_2;
    boolean_T tmp_5;
    boolean_T tmp_8;
    tmp_8 = rtmIsMajorTimeStep(windEmulatorStep4_WECSim_M);

    /* Update for RateLimiter: '<S609>/torqueSlewRate' incorporates:
     *  RateLimiter: '<S2>/acs880RateLim'
     *  RateLimiter: '<S373>/Rate Limiter'
     *  RateLimiter: '<S436>/Rate Limiter'
     *  RateLimiter: '<S609>/speedSlewRate'
     */
    windEmulatorStep4_WECSim_DW.PrevY =
      windEmulatorStep4_WECSim_B.torqueSlewRate;
    LastMajorTime_tmp = windEmulatorStep4_WECSim_M->Timing.t[0];
    windEmulatorStep4_WECSim_DW.LastMajorTime = LastMajorTime_tmp;

    /* Update for RateLimiter: '<S609>/speedSlewRate' */
    windEmulatorStep4_WECSim_DW.PrevY_a =
      windEmulatorStep4_WECSim_B.speedSlewRate;
    windEmulatorStep4_WECSim_DW.LastMajorTime_j = LastMajorTime_tmp;
    if (tmp_8) {
      /* Update for Memory: '<S2>/Memory' incorporates:
       *  Constant: '<S2>/powerUpButton'
       */
      windEmulatorStep4_WECSim_DW.Memory_PreviousInput_i =
        windEmulatorStep4_WECSim_cal->powerUpButton_Value;

      /* Update for Memory: '<S2>/Memory1' incorporates:
       *  Constant: '<S2>/powerDownButton'
       */
      windEmulatorStep4_WECSim_DW.Memory1_PreviousInput =
        windEmulatorStep4_WECSim_cal->powerDownButton_Value;

      /* Update for Memory: '<S2>/Memory2' incorporates:
       *  Constant: '<S2>/resetFaultButton'
       */
      windEmulatorStep4_WECSim_DW.Memory2_PreviousInput =
        windEmulatorStep4_WECSim_cal->resetFaultButton_Value;

      /* Update for Memory: '<S4>/Memory' incorporates:
       *  Constant: '<S4>/eStopButton'
       */
      windEmulatorStep4_WECSim_DW.Memory_PreviousInput_k =
        windEmulatorStep4_WECSim_cal->eStopButton_Value;

      /* Update for Memory: '<S4>/Memory1' incorporates:
       *  Constant: '<S4>/startButton'
       */
      windEmulatorStep4_WECSim_DW.Memory1_PreviousInput_d =
        windEmulatorStep4_WECSim_cal->startButton_Value;

      /* Update for Memory: '<S4>/Memory2' incorporates:
       *  Constant: '<S4>/stopButton'
       */
      windEmulatorStep4_WECSim_DW.Memory2_PreviousInput_l =
        windEmulatorStep4_WECSim_cal->stopButton_Value;

      /* Update for SimscapeExecutionBlock: '<S541>/STATE_1' */
      simulationData = static_cast<NeslSimulationData *>
        (windEmulatorStep4_WECSim_DW.STATE_1_SimData);
      time = windEmulatorStep4_WECSim_M->Timing.t[0];
      simulationData->mData->mTime.mN = 1;
      simulationData->mData->mTime.mX = &time;
      simulationData->mData->mContStates.mN = 0;
      simulationData->mData->mContStates.mX = NULL;
      simulationData->mData->mDiscStates.mN = 22;
      simulationData->mData->mDiscStates.mX =
        &windEmulatorStep4_WECSim_DW.STATE_1_Discrete_1041191992[0];
      simulationData->mData->mModeVector.mN = 15;
      simulationData->mData->mModeVector.mX =
        &windEmulatorStep4_WECSim_DW.STATE_1_Modes[0];
      tmp = false;
      simulationData->mData->mFoundZcEvents = tmp;
      simulationData->mData->mHadEvents = false;
      simulationData->mData->mIsMajorTimeStep = true;
      tmp = false;
      simulationData->mData->mIsSolverAssertCheck = tmp;
      simulationData->mData->mIsSolverCheckingCIC = false;
      simulationData->mData->mIsComputingJacobian = false;
      simulationData->mData->mIsEvaluatingF0 = false;
      simulationData->mData->mIsSolverRequestingReset = false;
      simulationData->mData->mIsModeUpdateTimeStep = true;
      tmp_1[0] = 0;
      tmp_0[0] = windEmulatorStep4_WECSim_B.INPUT_1_1_1[0];
      tmp_0[1] = windEmulatorStep4_WECSim_B.INPUT_1_1_1[1];
      tmp_0[2] = windEmulatorStep4_WECSim_B.INPUT_1_1_1[2];
      tmp_0[3] = windEmulatorStep4_WECSim_B.INPUT_1_1_1[3];
      tmp_1[1] = 4;
      tmp_0[4] = windEmulatorStep4_WECSim_B.INPUT_2_1_1[0];
      tmp_0[5] = windEmulatorStep4_WECSim_B.INPUT_2_1_1[1];
      tmp_0[6] = windEmulatorStep4_WECSim_B.INPUT_2_1_1[2];
      tmp_0[7] = windEmulatorStep4_WECSim_B.INPUT_2_1_1[3];
      tmp_1[2] = 8;
      tmp_0[8] = windEmulatorStep4_WECSim_B.INPUT_4_1_1[0];
      tmp_0[9] = windEmulatorStep4_WECSim_B.INPUT_4_1_1[1];
      tmp_0[10] = windEmulatorStep4_WECSim_B.INPUT_4_1_1[2];
      tmp_0[11] = windEmulatorStep4_WECSim_B.INPUT_4_1_1[3];
      tmp_1[3] = 12;
      tmp_0[12] = windEmulatorStep4_WECSim_B.INPUT_5_1_1[0];
      tmp_0[13] = windEmulatorStep4_WECSim_B.INPUT_5_1_1[1];
      tmp_0[14] = windEmulatorStep4_WECSim_B.INPUT_5_1_1[2];
      tmp_0[15] = windEmulatorStep4_WECSim_B.INPUT_5_1_1[3];
      tmp_1[4] = 16;
      tmp_0[16] = windEmulatorStep4_WECSim_B.INPUT_3_1_1[0];
      tmp_0[17] = windEmulatorStep4_WECSim_B.INPUT_3_1_1[1];
      tmp_0[18] = windEmulatorStep4_WECSim_B.INPUT_3_1_1[2];
      tmp_0[19] = windEmulatorStep4_WECSim_B.INPUT_3_1_1[3];
      tmp_1[5] = 20;
      simulationData->mData->mInputValues.mN = 20;
      simulationData->mData->mInputValues.mX = &tmp_0[0];
      simulationData->mData->mInputOffsets.mN = 6;
      simulationData->mData->mInputOffsets.mX = &tmp_1[0];
      simulator = static_cast<NeslSimulator *>
        (windEmulatorStep4_WECSim_DW.STATE_1_Simulator);
      diagnosticManager = static_cast<NeuDiagnosticManager *>
        (windEmulatorStep4_WECSim_DW.STATE_1_DiagMgr);
      diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
      i = ne_simulator_method(simulator, NESL_SIM_UPDATE, simulationData,
        diagnosticManager);
      if (i != 0) {
        tmp = error_buffer_is_empty(rtmGetErrorStatus(windEmulatorStep4_WECSim_M));
        if (tmp) {
          msg = rtw_diagnostics_msg(diagnosticTree);
          rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
        }
      }

      /* End of Update for SimscapeExecutionBlock: '<S541>/STATE_1' */

      /* Update for DiscreteIntegrator: '<S437>/Discrete-Time Integrator' */
      windEmulatorStep4_WECSim_DW.DiscreteTimeIntegrator_DSTATE +=
        windEmulatorStep4_WECSim_cal->DiscreteTimeIntegrator_gainval *
        windEmulatorStep4_WECSim_B.Sum_f;

      /* Update for DiscreteIntegrator: '<S372>/Discrete-Time Integrator' */
      windEmulatorStep4_WECSim_DW.DiscreteTimeIntegrator_DSTATE_n +=
        windEmulatorStep4_WECSim_cal->DiscreteTimeIntegrator_gainva_l *
        windEmulatorStep4_WECSim_B.Sum_me;

      /* Update for DiscreteIntegrator: '<S435>/Discrete-Time Integrator' */
      windEmulatorStep4_WECSim_DW.DiscreteTimeIntegrator_DSTATE_l +=
        windEmulatorStep4_WECSim_cal->DiscreteTimeIntegrator_gainva_b *
        windEmulatorStep4_WECSim_B.Sum_d;

      /* Update for Memory: '<S1>/Memory' incorporates:
       *  Constant: '<S1>/powerUpButton'
       */
      windEmulatorStep4_WECSim_DW.Memory_PreviousInput_d =
        windEmulatorStep4_WECSim_cal->powerUpButton_Value_b;

      /* Update for Memory: '<S1>/Memory1' incorporates:
       *  Constant: '<S1>/powerDownButton'
       */
      windEmulatorStep4_WECSim_DW.Memory1_PreviousInput_p =
        windEmulatorStep4_WECSim_cal->powerDownButton_Value_k;

      /* Update for Memory: '<S1>/Memory2' incorporates:
       *  Constant: '<S1>/resetFaultButton'
       */
      windEmulatorStep4_WECSim_DW.Memory2_PreviousInput_h =
        windEmulatorStep4_WECSim_cal->resetFaultButton_Value_n;

      /* Update for Memory: '<S552>/lastRawCounts' */
      windEmulatorStep4_WECSim_DW.lastRawCounts_PreviousInput =
        windEmulatorStep4_WECSim_B.CastToDouble_p;

      /* Update for Memory: '<S552>/lastTurn' */
      windEmulatorStep4_WECSim_DW.lastTurn_PreviousInput =
        windEmulatorStep4_WECSim_B.Add1_p2;

      /* Update for UnitDelay: '<S553>/UD' */
      windEmulatorStep4_WECSim_DW.UD_DSTATE = windEmulatorStep4_WECSim_B.TSamp;

      /* Update for Memory: '<S29>/Memory' */
      windEmulatorStep4_WECSim_DW.Memory_PreviousInput =
        windEmulatorStep4_WECSim_B.Sum_e;
    }

    /* Update for RateLimiter: '<S2>/acs880RateLim' */
    windEmulatorStep4_WECSim_DW.PrevY_b =
      windEmulatorStep4_WECSim_B.acs880RateLim;
    windEmulatorStep4_WECSim_DW.LastMajorTime_a = LastMajorTime_tmp;

    /* Update for SimscapeExecutionBlock: '<S216>/STATE_1' incorporates:
     *  SimscapeExecutionBlock: '<S332>/STATE_1'
     *  TransportDelay: '<S139>/Transport Delay'
     *  TransportDelay: '<S143>/Transport Delay'
     *  TransportDelay: '<S60>/Transport Delay'
     *  TransportDelay: '<S64>/Transport Delay'
     */
    simulationData = static_cast<NeslSimulationData *>
      (windEmulatorStep4_WECSim_DW.STATE_1_SimData_h);
    time_tmp = windEmulatorStep4_WECSim_M->Timing.t[0];
    time_0 = time_tmp;
    simulationData->mData->mTime.mN = 1;
    simulationData->mData->mTime.mX = &time_0;
    simulationData->mData->mContStates.mN = 2;
    simulationData->mData->mContStates.mX =
      &windEmulatorStep4_WECSim_X.windEmulatorStep4_WECSimhptoSim[0];
    simulationData->mData->mDiscStates.mN = 0;
    simulationData->mData->mDiscStates.mX =
      &windEmulatorStep4_WECSim_DW.STATE_1_Discrete;
    simulationData->mData->mModeVector.mN = 0;
    simulationData->mData->mModeVector.mX =
      &windEmulatorStep4_WECSim_DW.STATE_1_Modes_j;
    tmp = false;
    simulationData->mData->mFoundZcEvents = tmp;
    simulationData->mData->mHadEvents = false;
    tmp = rtmIsMajorTimeStep(windEmulatorStep4_WECSim_M);
    simulationData->mData->mIsMajorTimeStep = tmp;
    tmp_2 = false;
    simulationData->mData->mIsSolverAssertCheck = tmp_2;
    simulationData->mData->mIsSolverCheckingCIC = false;
    tmp_2 = rtsiIsSolverComputingJacobian
      (&windEmulatorStep4_WECSim_M->solverInfo);
    simulationData->mData->mIsComputingJacobian = tmp_2;
    simulationData->mData->mIsEvaluatingF0 = false;
    simulationData->mData->mIsSolverRequestingReset = false;
    tmp_2 = rtsiIsModeUpdateTimeStep(&windEmulatorStep4_WECSim_M->solverInfo);
    simulationData->mData->mIsModeUpdateTimeStep = tmp_2;
    tmp_4[0] = 0;
    tmp_3[0] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[0];
    tmp_3[1] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[1];
    tmp_3[2] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[2];
    tmp_3[3] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[3];
    tmp_4[1] = 4;
    tmp_3[4] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[0];
    tmp_3[5] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[1];
    tmp_3[6] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[2];
    tmp_3[7] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[3];
    tmp_4[2] = 8;
    tmp_3[8] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[0];
    tmp_3[9] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[1];
    tmp_3[10] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[2];
    tmp_3[11] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[3];
    tmp_4[3] = 12;
    tmp_3[12] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[0];
    tmp_3[13] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[1];
    tmp_3[14] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[2];
    tmp_3[15] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[3];
    tmp_4[4] = 16;
    tmp_3[16] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[0];
    tmp_3[17] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[1];
    tmp_3[18] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[2];
    tmp_3[19] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[3];
    tmp_4[5] = 20;
    tmp_3[20] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[0];
    tmp_3[21] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[1];
    tmp_3[22] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[2];
    tmp_3[23] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[3];
    tmp_4[6] = 24;
    tmp_3[24] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[0];
    tmp_3[25] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[1];
    tmp_3[26] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[2];
    tmp_3[27] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[3];
    tmp_4[7] = 28;
    tmp_3[28] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[0];
    tmp_3[29] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[1];
    tmp_3[30] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[2];
    tmp_3[31] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[3];
    tmp_4[8] = 32;
    tmp_3[32] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[0];
    tmp_3[33] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[1];
    tmp_3[34] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[2];
    tmp_3[35] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[3];
    tmp_4[9] = 36;
    tmp_3[36] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[0];
    tmp_3[37] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[1];
    tmp_3[38] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[2];
    tmp_3[39] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[3];
    tmp_4[10] = 40;
    tmp_3[40] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[0];
    tmp_3[41] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[1];
    tmp_3[42] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[2];
    tmp_3[43] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[3];
    tmp_4[11] = 44;
    tmp_3[44] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[0];
    tmp_3[45] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[1];
    tmp_3[46] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[2];
    tmp_3[47] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[3];
    tmp_4[12] = 48;
    tmp_3[48] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[0];
    tmp_3[49] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[1];
    tmp_3[50] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[2];
    tmp_3[51] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[3];
    tmp_4[13] = 52;
    simulationData->mData->mInputValues.mN = 52;
    simulationData->mData->mInputValues.mX = &tmp_3[0];
    simulationData->mData->mInputOffsets.mN = 14;
    simulationData->mData->mInputOffsets.mX = &tmp_4[0];
    simulator = static_cast<NeslSimulator *>
      (windEmulatorStep4_WECSim_DW.STATE_1_Simulator_f);
    diagnosticManager = static_cast<NeuDiagnosticManager *>
      (windEmulatorStep4_WECSim_DW.STATE_1_DiagMgr_o);
    diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
    i = ne_simulator_method(simulator, NESL_SIM_UPDATE, simulationData,
      diagnosticManager);
    if (i != 0) {
      tmp_5 = error_buffer_is_empty(rtmGetErrorStatus(windEmulatorStep4_WECSim_M));
      if (tmp_5) {
        msg = rtw_diagnostics_msg(diagnosticTree);
        rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
      }
    }

    /* End of Update for SimscapeExecutionBlock: '<S216>/STATE_1' */
    if (tmp_8) {
      /* Update for Delay: '<S58>/Delay One Step' */
      windEmulatorStep4_WECSim_DW.DelayOneStep_DSTATE =
        windEmulatorStep4_WECSim_B.shaftSpeed;

      /* Update for DiscreteIntegrator: '<S282>/Integrator' */
      windEmulatorStep4_WECSim_DW.Integrator_DSTATE +=
        windEmulatorStep4_WECSim_cal->Integrator_gainval *
        windEmulatorStep4_WECSim_B.IntegralGain;

      /* Update for Delay: '<S275>/UD' */
      windEmulatorStep4_WECSim_DW.UD_DSTATE_j = windEmulatorStep4_WECSim_B.Tsamp;
    }

    /* Update for SimscapeExecutionBlock: '<S332>/STATE_1' */
    simulationData = static_cast<NeslSimulationData *>
      (windEmulatorStep4_WECSim_DW.STATE_1_SimData_a);
    time_1 = time_tmp;
    simulationData->mData->mTime.mN = 1;
    simulationData->mData->mTime.mX = &time_1;
    simulationData->mData->mContStates.mN = 35;
    simulationData->mData->mContStates.mX =
      &windEmulatorStep4_WECSim_X.windEmulatorStep4_WECSimhptoS_h[0];
    simulationData->mData->mDiscStates.mN = 6;
    simulationData->mData->mDiscStates.mX =
      &windEmulatorStep4_WECSim_DW.STATE_1_Discrete_208214823[0];
    simulationData->mData->mModeVector.mN = 21;
    simulationData->mData->mModeVector.mX =
      &windEmulatorStep4_WECSim_DW.STATE_1_Modes_i[0];
    tmp_5 = false;
    simulationData->mData->mFoundZcEvents = tmp_5;
    simulationData->mData->mHadEvents = false;
    simulationData->mData->mIsMajorTimeStep = tmp;
    tmp = false;
    simulationData->mData->mIsSolverAssertCheck = tmp;
    simulationData->mData->mIsSolverCheckingCIC = false;
    tmp = rtsiIsSolverComputingJacobian(&windEmulatorStep4_WECSim_M->solverInfo);
    simulationData->mData->mIsComputingJacobian = tmp;
    simulationData->mData->mIsEvaluatingF0 = false;
    simulationData->mData->mIsSolverRequestingReset = false;
    simulationData->mData->mIsModeUpdateTimeStep = tmp_2;
    tmp_7[0] = 0;
    tmp_6[0] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[0];
    tmp_6[1] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[1];
    tmp_6[2] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[2];
    tmp_6[3] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[3];
    tmp_7[1] = 4;
    tmp_6[4] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[0];
    tmp_6[5] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[1];
    tmp_6[6] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[2];
    tmp_6[7] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[3];
    tmp_7[2] = 8;
    tmp_6[8] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[0];
    tmp_6[9] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[1];
    tmp_6[10] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[2];
    tmp_6[11] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[3];
    tmp_7[3] = 12;
    simulationData->mData->mInputValues.mN = 12;
    simulationData->mData->mInputValues.mX = &tmp_6[0];
    simulationData->mData->mInputOffsets.mN = 4;
    simulationData->mData->mInputOffsets.mX = &tmp_7[0];
    simulator = static_cast<NeslSimulator *>
      (windEmulatorStep4_WECSim_DW.STATE_1_Simulator_i);
    diagnosticManager = static_cast<NeuDiagnosticManager *>
      (windEmulatorStep4_WECSim_DW.STATE_1_DiagMgr_g);
    diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
    i = ne_simulator_method(simulator, NESL_SIM_UPDATE, simulationData,
      diagnosticManager);
    if (i != 0) {
      tmp = error_buffer_is_empty(rtmGetErrorStatus(windEmulatorStep4_WECSim_M));
      if (tmp) {
        msg = rtw_diagnostics_msg(diagnosticTree);
        rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
      }
    }

    /* Update for TransportDelay: '<S60>/Transport Delay' */
    OUTPUT_1_0_k = windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[23];
    OUTPUT_1_0_k_0 = windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[20];
    OUTPUT_1_0_k_1 = windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[24];
    OUTPUT_1_0_k_2 = windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[21];
    OUTPUT_1_0_k_3 = windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[25];
    OUTPUT_1_0_k_4 = windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[22];
    for (i = 0; i < 6; i++) {
      TransportDelay_IWORK = windEmulatorStep4_WECSim_DW.TransportDelay_IWORK[i
        + 6];
      TransportDelay_IWORK_0 =
        windEmulatorStep4_WECSim_DW.TransportDelay_IWORK[i + 18];
      if (TransportDelay_IWORK < TransportDelay_IWORK_0 - 1) {
        TransportDelay_IWORK++;
        windEmulatorStep4_WECSim_DW.TransportDelay_IWORK[i + 6] =
          TransportDelay_IWORK;
      } else {
        TransportDelay_IWORK = 0;
        windEmulatorStep4_WECSim_DW.TransportDelay_IWORK[i + 6] = 0;
      }

      TransportDelay_IWORK_1 =
        windEmulatorStep4_WECSim_DW.TransportDelay_IWORK[i];
      if (TransportDelay_IWORK == TransportDelay_IWORK_1) {
        if (TransportDelay_IWORK_1 < TransportDelay_IWORK_0 - 1) {
          TransportDelay_IWORK_1++;
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK[i] =
            TransportDelay_IWORK_1;
        } else {
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK[i] = 0;
        }
      }

      tmp_9[0] = OUTPUT_1_0_k;
      tmp_9[3] = OUTPUT_1_0_k_0;
      tmp_9[1] = OUTPUT_1_0_k_1;
      tmp_9[4] = OUTPUT_1_0_k_2;
      tmp_9[2] = OUTPUT_1_0_k_3;
      tmp_9[5] = OUTPUT_1_0_k_4;
      (static_cast<real_T *>(windEmulatorStep4_WECSim_DW.TransportDelay_PWORK[i]))
        [TransportDelay_IWORK] = tmp_9[i];
      (static_cast<real_T *>(windEmulatorStep4_WECSim_DW.TransportDelay_PWORK[i]))
        [TransportDelay_IWORK_0 + TransportDelay_IWORK] = time_tmp;
    }

    /* Update for TransportDelay: '<S139>/Transport Delay' */
    OUTPUT_1_0_k = windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[29];
    OUTPUT_1_0_k_0 = windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[26];
    OUTPUT_1_0_k_1 = windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[30];
    OUTPUT_1_0_k_2 = windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[27];
    OUTPUT_1_0_k_3 = windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[31];
    OUTPUT_1_0_k_4 = windEmulatorStep4_WECSim_B.OUTPUT_1_0_k[28];
    for (i = 0; i < 6; i++) {
      TransportDelay_IWORK =
        windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c[i + 6];
      TransportDelay_IWORK_0 =
        windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c[i + 18];
      if (TransportDelay_IWORK < TransportDelay_IWORK_0 - 1) {
        TransportDelay_IWORK++;
        windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c[i + 6] =
          TransportDelay_IWORK;
      } else {
        TransportDelay_IWORK = 0;
        windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c[i + 6] = 0;
      }

      TransportDelay_IWORK_1 =
        windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c[i];
      if (TransportDelay_IWORK == TransportDelay_IWORK_1) {
        if (TransportDelay_IWORK_1 < TransportDelay_IWORK_0 - 1) {
          TransportDelay_IWORK_1++;
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c[i] =
            TransportDelay_IWORK_1;
        } else {
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c[i] = 0;
        }
      }

      tmp_9[0] = OUTPUT_1_0_k;
      tmp_9[3] = OUTPUT_1_0_k_0;
      tmp_9[1] = OUTPUT_1_0_k_1;
      tmp_9[4] = OUTPUT_1_0_k_2;
      tmp_9[2] = OUTPUT_1_0_k_3;
      tmp_9[5] = OUTPUT_1_0_k_4;
      (static_cast<real_T *>
        (windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_f[i]))
        [TransportDelay_IWORK] = tmp_9[i];
      (static_cast<real_T *>
        (windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_f[i]))
        [TransportDelay_IWORK_0 + TransportDelay_IWORK] = time_tmp;
    }

    /* Update for TransportDelay: '<S64>/Transport Delay' */
    for (i = 0; i < 6; i++) {
      TransportDelay_IWORK =
        windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_o[i + 6];
      TransportDelay_IWORK_0 =
        windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_o[i + 18];
      if (TransportDelay_IWORK < TransportDelay_IWORK_0 - 1) {
        TransportDelay_IWORK++;
        windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_o[i + 6] =
          TransportDelay_IWORK;
      } else {
        TransportDelay_IWORK = 0;
        windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_o[i + 6] = 0;
      }

      TransportDelay_IWORK_1 =
        windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_o[i];
      if (TransportDelay_IWORK == TransportDelay_IWORK_1) {
        if (TransportDelay_IWORK_1 < TransportDelay_IWORK_0 - 1) {
          TransportDelay_IWORK_1++;
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_o[i] =
            TransportDelay_IWORK_1;
        } else {
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_o[i] = 0;
        }
      }

      (static_cast<real_T *>
        (windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_j[i]))
        [TransportDelay_IWORK] = windEmulatorStep4_WECSim_B.TransportDelay[i];
      (static_cast<real_T *>
        (windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_j[i]))
        [TransportDelay_IWORK_0 + TransportDelay_IWORK] = time_tmp;
    }

    /* Update for TransportDelay: '<S143>/Transport Delay' */
    for (i = 0; i < 6; i++) {
      TransportDelay_IWORK =
        windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c1[i + 6];
      TransportDelay_IWORK_0 =
        windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c1[i + 18];
      if (TransportDelay_IWORK < TransportDelay_IWORK_0 - 1) {
        TransportDelay_IWORK++;
        windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c1[i + 6] =
          TransportDelay_IWORK;
      } else {
        TransportDelay_IWORK = 0;
        windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c1[i + 6] = 0;
      }

      TransportDelay_IWORK_1 =
        windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c1[i];
      if (TransportDelay_IWORK == TransportDelay_IWORK_1) {
        if (TransportDelay_IWORK_1 < TransportDelay_IWORK_0 - 1) {
          TransportDelay_IWORK_1++;
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c1[i] =
            TransportDelay_IWORK_1;
        } else {
          windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c1[i] = 0;
        }
      }

      (static_cast<real_T *>
        (windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_n[i]))
        [TransportDelay_IWORK] = windEmulatorStep4_WECSim_B.TransportDelay_b[i];
      (static_cast<real_T *>
        (windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_n[i]))
        [TransportDelay_IWORK_0 + TransportDelay_IWORK] = time_tmp;
    }

    if (tmp_8) {
      /* Update for Memory: '<S433>/Memory' */
      windEmulatorStep4_WECSim_DW.Memory_PreviousInput_kk =
        windEmulatorStep4_WECSim_B.Logic[0];

      /* Update for Memory: '<S434>/Memory' */
      windEmulatorStep4_WECSim_DW.Memory_PreviousInput_h =
        windEmulatorStep4_WECSim_B.Logic_g[0];

      /* Update for DiscreteIntegrator: '<S412>/Integrator' */
      windEmulatorStep4_WECSim_DW.Integrator_DSTATE_d +=
        windEmulatorStep4_WECSim_cal->Integrator_gainval_d *
        windEmulatorStep4_WECSim_B.IntegralGain_h;
      windEmulatorStep4_WECSim_DW.Integrator_PrevResetState = static_cast<int8_T>
        (windEmulatorStep4_WECSim_B.BusAssignment_c.speedCtrlReset);

      /* Update for DiscreteIntegrator: '<S407>/Filter' */
      windEmulatorStep4_WECSim_DW.Filter_DSTATE +=
        windEmulatorStep4_WECSim_cal->Filter_gainval *
        windEmulatorStep4_WECSim_B.FilterCoefficient;
      windEmulatorStep4_WECSim_DW.Filter_PrevResetState = static_cast<int8_T>
        (windEmulatorStep4_WECSim_B.BusAssignment_c.speedCtrlReset);

      /* Update for Memory: '<S496>/Memory' */
      windEmulatorStep4_WECSim_DW.Memory_PreviousInput_g =
        windEmulatorStep4_WECSim_B.Logic_c[0];

      /* Update for Memory: '<S497>/Memory' */
      windEmulatorStep4_WECSim_DW.Memory_PreviousInput_n =
        windEmulatorStep4_WECSim_B.Logic_p[0];

      /* Update for DiscreteIntegrator: '<S475>/Integrator' */
      windEmulatorStep4_WECSim_DW.Integrator_DSTATE_e +=
        windEmulatorStep4_WECSim_cal->Integrator_gainval_b *
        windEmulatorStep4_WECSim_B.IntegralGain_hk;
      windEmulatorStep4_WECSim_DW.Integrator_PrevResetState_g =
        static_cast<int8_T>
        (windEmulatorStep4_WECSim_B.BusAssignment_c.speedCtrlReset);

      /* Update for DiscreteIntegrator: '<S470>/Filter' */
      windEmulatorStep4_WECSim_DW.Filter_DSTATE_b +=
        windEmulatorStep4_WECSim_cal->Filter_gainval_d *
        windEmulatorStep4_WECSim_B.FilterCoefficient_g;
      windEmulatorStep4_WECSim_DW.Filter_PrevResetState_g = static_cast<int8_T>
        (windEmulatorStep4_WECSim_B.BusAssignment_c.speedCtrlReset);
    }

    /* Update for RateLimiter: '<S373>/Rate Limiter' */
    windEmulatorStep4_WECSim_DW.PrevY_l =
      windEmulatorStep4_WECSim_B.RateLimiter_b;
    windEmulatorStep4_WECSim_DW.LastMajorTime_d = LastMajorTime_tmp;

    /* Update for RateLimiter: '<S436>/Rate Limiter' */
    windEmulatorStep4_WECSim_DW.PrevY_e =
      windEmulatorStep4_WECSim_B.RateLimiter_a;
    windEmulatorStep4_WECSim_DW.LastMajorTime_k = LastMajorTime_tmp;
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(windEmulatorStep4_WECSim_M)) {
    rt_ertODEUpdateContinuousStates(&windEmulatorStep4_WECSim_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++windEmulatorStep4_WECSim_M->Timing.clockTick0)) {
      ++windEmulatorStep4_WECSim_M->Timing.clockTickH0;
    }

    windEmulatorStep4_WECSim_M->Timing.t[0] = rtsiGetSolverStopTime
      (&windEmulatorStep4_WECSim_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.004s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The absolute time is the multiplication of "clockTick1"
       * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
       * overflow during the application lifespan selected.
       * Timer of this task consists of two 32 bit unsigned integers.
       * The two integers represent the low bits Timing.clockTick1 and the high bits
       * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
       */
      if (!(++windEmulatorStep4_WECSim_M->Timing.clockTick1)) {
        ++windEmulatorStep4_WECSim_M->Timing.clockTickH1;
      }

      windEmulatorStep4_WECSim_M->Timing.t[1] =
        windEmulatorStep4_WECSim_M->Timing.clockTick1 *
        windEmulatorStep4_WECSim_M->Timing.stepSize1 +
        windEmulatorStep4_WECSim_M->Timing.clockTickH1 *
        windEmulatorStep4_WECSim_M->Timing.stepSize1 * 4294967296.0;
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void windEmulatorStep4_WECSim_derivatives(void)
{
  NeslSimulationData *simulationData;
  NeslSimulator *simulator;
  NeuDiagnosticManager *diagnosticManager;
  NeuDiagnosticTree *diagnosticTree;
  XDot_windEmulatorStep4_WECSim_T *_rtXdot;
  char *msg;
  real_T tmp_1[52];
  real_T tmp_4[12];
  real_T time;
  real_T time_0;
  real_T time_tmp;
  int_T tmp_2[14];
  int_T tmp_5[4];
  int_T is;
  uint32_T ri;
  boolean_T tmp;
  boolean_T tmp_0;
  boolean_T tmp_3;
  _rtXdot = ((XDot_windEmulatorStep4_WECSim_T *)
             windEmulatorStep4_WECSim_M->derivs);

  /* Derivatives for Integrator: '<S590>/Integrator' */
  if (!windEmulatorStep4_WECSim_B.BusAssignment_b.resetSidIntegrator) {
    _rtXdot->Integrator_CSTATE = windEmulatorStep4_WECSim_B.IProdOut;
  } else {
    /* level reset is active */
    _rtXdot->Integrator_CSTATE = 0.0;
  }

  /* End of Derivatives for Integrator: '<S590>/Integrator' */

  /* Derivatives for StateSpace: '<S531>/Internal' */
  _rtXdot->Internal_CSTATE[0] = 0.0;
  _rtXdot->Internal_CSTATE[1] = 0.0;
  _rtXdot->Internal_CSTATE[2] = 0.0;
  for (ri = windEmulatorStep4_WECSim_cal->Internal_A_jc[0U]; ri <
       windEmulatorStep4_WECSim_cal->Internal_A_jc[1U]; ri++) {
    _rtXdot->Internal_CSTATE[windEmulatorStep4_WECSim_cal->Internal_A_ir[ri]] +=
      windEmulatorStep4_WECSim_cal->Internal_A_pr[ri] *
      windEmulatorStep4_WECSim_X.Internal_CSTATE[0U];
  }

  for (ri = windEmulatorStep4_WECSim_cal->Internal_A_jc[1U]; ri <
       windEmulatorStep4_WECSim_cal->Internal_A_jc[2U]; ri++) {
    _rtXdot->Internal_CSTATE[windEmulatorStep4_WECSim_cal->Internal_A_ir[ri]] +=
      windEmulatorStep4_WECSim_cal->Internal_A_pr[ri] *
      windEmulatorStep4_WECSim_X.Internal_CSTATE[1U];
  }

  for (ri = windEmulatorStep4_WECSim_cal->Internal_A_jc[2U]; ri <
       windEmulatorStep4_WECSim_cal->Internal_A_jc[3U]; ri++) {
    _rtXdot->Internal_CSTATE[windEmulatorStep4_WECSim_cal->Internal_A_ir[ri]] +=
      windEmulatorStep4_WECSim_cal->Internal_A_pr[ri] *
      windEmulatorStep4_WECSim_X.Internal_CSTATE[2U];
  }

  for (ri = windEmulatorStep4_WECSim_cal->Internal_B_jc[0U]; ri <
       windEmulatorStep4_WECSim_cal->Internal_B_jc[1U]; ri++) {
    _rtXdot->Internal_CSTATE[windEmulatorStep4_WECSim_cal->Internal_B_ir] +=
      windEmulatorStep4_WECSim_cal->Internal_B_pr *
      windEmulatorStep4_WECSim_B.Sum_a;
  }

  /* End of Derivatives for StateSpace: '<S531>/Internal' */

  /* Derivatives for StateSpace: '<S545>/Internal' */
  _rtXdot->Internal_CSTATE_j = 0.0;
  for (ri = windEmulatorStep4_WECSim_cal->Internal_A_jc_i[0U]; ri <
       windEmulatorStep4_WECSim_cal->Internal_A_jc_i[1U]; ri++) {
    _rtXdot->Internal_CSTATE_j += windEmulatorStep4_WECSim_cal->Internal_A_pr_j *
      windEmulatorStep4_WECSim_X.Internal_CSTATE_j;
  }

  for (ri = windEmulatorStep4_WECSim_cal->Internal_B_jc_i[0U]; ri <
       windEmulatorStep4_WECSim_cal->Internal_B_jc_i[1U]; ri++) {
    _rtXdot->Internal_CSTATE_j += windEmulatorStep4_WECSim_cal->Internal_B_pr_g *
      windEmulatorStep4_WECSim_B.ControlSignal2;
  }

  /* End of Derivatives for StateSpace: '<S545>/Internal' */

  /* Derivatives for StateSpace: '<S542>/Internal' */
  _rtXdot->Internal_CSTATE_a = 0.0;
  for (ri = windEmulatorStep4_WECSim_cal->Internal_A_jc_e[0U]; ri <
       windEmulatorStep4_WECSim_cal->Internal_A_jc_e[1U]; ri++) {
    _rtXdot->Internal_CSTATE_a += windEmulatorStep4_WECSim_cal->Internal_A_pr_e *
      windEmulatorStep4_WECSim_X.Internal_CSTATE_a;
  }

  for (ri = windEmulatorStep4_WECSim_cal->Internal_B_jc_h[0U]; ri <
       windEmulatorStep4_WECSim_cal->Internal_B_jc_h[1U]; ri++) {
    _rtXdot->Internal_CSTATE_a += windEmulatorStep4_WECSim_cal->Internal_B_pr_k *
      windEmulatorStep4_WECSim_B.ControlSignal1;
  }

  /* End of Derivatives for StateSpace: '<S542>/Internal' */

  /* Derivatives for SimscapeExecutionBlock: '<S216>/STATE_1' incorporates:
   *  SimscapeExecutionBlock: '<S332>/STATE_1'
   */
  simulationData = static_cast<NeslSimulationData *>
    (windEmulatorStep4_WECSim_DW.STATE_1_SimData_h);
  time_tmp = windEmulatorStep4_WECSim_M->Timing.t[0];
  time = time_tmp;
  simulationData->mData->mTime.mN = 1;
  simulationData->mData->mTime.mX = &time;
  simulationData->mData->mContStates.mN = 2;
  simulationData->mData->mContStates.mX =
    &windEmulatorStep4_WECSim_X.windEmulatorStep4_WECSimhptoSim[0];
  simulationData->mData->mDiscStates.mN = 0;
  simulationData->mData->mDiscStates.mX =
    &windEmulatorStep4_WECSim_DW.STATE_1_Discrete;
  simulationData->mData->mModeVector.mN = 0;
  simulationData->mData->mModeVector.mX =
    &windEmulatorStep4_WECSim_DW.STATE_1_Modes_j;
  tmp = false;
  simulationData->mData->mFoundZcEvents = tmp;
  simulationData->mData->mHadEvents = false;
  tmp = rtmIsMajorTimeStep(windEmulatorStep4_WECSim_M);
  simulationData->mData->mIsMajorTimeStep = tmp;
  tmp_0 = false;
  simulationData->mData->mIsSolverAssertCheck = tmp_0;
  simulationData->mData->mIsSolverCheckingCIC = false;
  tmp_0 = rtsiIsSolverComputingJacobian(&windEmulatorStep4_WECSim_M->solverInfo);
  simulationData->mData->mIsComputingJacobian = tmp_0;
  simulationData->mData->mIsEvaluatingF0 = false;
  simulationData->mData->mIsSolverRequestingReset = false;
  tmp_0 = rtsiIsModeUpdateTimeStep(&windEmulatorStep4_WECSim_M->solverInfo);
  simulationData->mData->mIsModeUpdateTimeStep = tmp_0;
  tmp_2[0] = 0;
  tmp_1[0] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[0];
  tmp_1[1] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[1];
  tmp_1[2] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[2];
  tmp_1[3] = windEmulatorStep4_WECSim_B.INPUT_5_1_1_c[3];
  tmp_2[1] = 4;
  tmp_1[4] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[0];
  tmp_1[5] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[1];
  tmp_1[6] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[2];
  tmp_1[7] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_p[3];
  tmp_2[2] = 8;
  tmp_1[8] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[0];
  tmp_1[9] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[1];
  tmp_1[10] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[2];
  tmp_1[11] = windEmulatorStep4_WECSim_B.INPUT_1_1_2[3];
  tmp_2[3] = 12;
  tmp_1[12] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[0];
  tmp_1[13] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[1];
  tmp_1[14] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[2];
  tmp_1[15] = windEmulatorStep4_WECSim_B.INPUT_1_1_3[3];
  tmp_2[4] = 16;
  tmp_1[16] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[0];
  tmp_1[17] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[1];
  tmp_1[18] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[2];
  tmp_1[19] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_l[3];
  tmp_2[5] = 20;
  tmp_1[20] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[0];
  tmp_1[21] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[1];
  tmp_1[22] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[2];
  tmp_1[23] = windEmulatorStep4_WECSim_B.INPUT_2_1_2[3];
  tmp_2[6] = 24;
  tmp_1[24] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[0];
  tmp_1[25] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[1];
  tmp_1[26] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[2];
  tmp_1[27] = windEmulatorStep4_WECSim_B.INPUT_2_1_3[3];
  tmp_2[7] = 28;
  tmp_1[28] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[0];
  tmp_1[29] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[1];
  tmp_1[30] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[2];
  tmp_1[31] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_k[3];
  tmp_2[8] = 32;
  tmp_1[32] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[0];
  tmp_1[33] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[1];
  tmp_1[34] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[2];
  tmp_1[35] = windEmulatorStep4_WECSim_B.INPUT_3_1_2[3];
  tmp_2[9] = 36;
  tmp_1[36] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[0];
  tmp_1[37] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[1];
  tmp_1[38] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[2];
  tmp_1[39] = windEmulatorStep4_WECSim_B.INPUT_3_1_3[3];
  tmp_2[10] = 40;
  tmp_1[40] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[0];
  tmp_1[41] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[1];
  tmp_1[42] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[2];
  tmp_1[43] = windEmulatorStep4_WECSim_B.INPUT_4_1_1_k[3];
  tmp_2[11] = 44;
  tmp_1[44] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[0];
  tmp_1[45] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[1];
  tmp_1[46] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[2];
  tmp_1[47] = windEmulatorStep4_WECSim_B.INPUT_4_1_2[3];
  tmp_2[12] = 48;
  tmp_1[48] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[0];
  tmp_1[49] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[1];
  tmp_1[50] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[2];
  tmp_1[51] = windEmulatorStep4_WECSim_B.INPUT_4_1_3[3];
  tmp_2[13] = 52;
  simulationData->mData->mInputValues.mN = 52;
  simulationData->mData->mInputValues.mX = &tmp_1[0];
  simulationData->mData->mInputOffsets.mN = 14;
  simulationData->mData->mInputOffsets.mX = &tmp_2[0];
  simulationData->mData->mDx.mN = 2;
  simulationData->mData->mDx.mX = &_rtXdot->windEmulatorStep4_WECSimhptoSim[0];
  simulator = static_cast<NeslSimulator *>
    (windEmulatorStep4_WECSim_DW.STATE_1_Simulator_f);
  diagnosticManager = static_cast<NeuDiagnosticManager *>
    (windEmulatorStep4_WECSim_DW.STATE_1_DiagMgr_o);
  diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
  is = ne_simulator_method(simulator, NESL_SIM_DERIVATIVES, simulationData,
    diagnosticManager);
  if (is != 0) {
    tmp_3 = error_buffer_is_empty(rtmGetErrorStatus(windEmulatorStep4_WECSim_M));
    if (tmp_3) {
      msg = rtw_diagnostics_msg(diagnosticTree);
      rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
    }
  }

  /* End of Derivatives for SimscapeExecutionBlock: '<S216>/STATE_1' */

  /* Derivatives for SimscapeInputBlock: '<S332>/INPUT_3_1_1' */
  _rtXdot->windEmulatorStep4_WECSimhptoS_n =
    (windEmulatorStep4_WECSim_B.velocity[4] -
     windEmulatorStep4_WECSim_X.windEmulatorStep4_WECSimhptoS_n) * 1000.0;

  /* Derivatives for SimscapeExecutionBlock: '<S332>/STATE_1' */
  simulationData = static_cast<NeslSimulationData *>
    (windEmulatorStep4_WECSim_DW.STATE_1_SimData_a);
  time_0 = time_tmp;
  simulationData->mData->mTime.mN = 1;
  simulationData->mData->mTime.mX = &time_0;
  simulationData->mData->mContStates.mN = 35;
  simulationData->mData->mContStates.mX =
    &windEmulatorStep4_WECSim_X.windEmulatorStep4_WECSimhptoS_h[0];
  simulationData->mData->mDiscStates.mN = 6;
  simulationData->mData->mDiscStates.mX =
    &windEmulatorStep4_WECSim_DW.STATE_1_Discrete_208214823[0];
  simulationData->mData->mModeVector.mN = 21;
  simulationData->mData->mModeVector.mX =
    &windEmulatorStep4_WECSim_DW.STATE_1_Modes_i[0];
  tmp_3 = false;
  simulationData->mData->mFoundZcEvents = tmp_3;
  simulationData->mData->mHadEvents = false;
  simulationData->mData->mIsMajorTimeStep = tmp;
  tmp = false;
  simulationData->mData->mIsSolverAssertCheck = tmp;
  simulationData->mData->mIsSolverCheckingCIC = false;
  tmp = rtsiIsSolverComputingJacobian(&windEmulatorStep4_WECSim_M->solverInfo);
  simulationData->mData->mIsComputingJacobian = tmp;
  simulationData->mData->mIsEvaluatingF0 = false;
  simulationData->mData->mIsSolverRequestingReset = false;
  simulationData->mData->mIsModeUpdateTimeStep = tmp_0;
  tmp_5[0] = 0;
  tmp_4[0] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[0];
  tmp_4[1] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[1];
  tmp_4[2] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[2];
  tmp_4[3] = windEmulatorStep4_WECSim_B.INPUT_2_1_1_c[3];
  tmp_5[1] = 4;
  tmp_4[4] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[0];
  tmp_4[5] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[1];
  tmp_4[6] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[2];
  tmp_4[7] = windEmulatorStep4_WECSim_B.INPUT_3_1_1_d[3];
  tmp_5[2] = 8;
  tmp_4[8] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[0];
  tmp_4[9] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[1];
  tmp_4[10] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[2];
  tmp_4[11] = windEmulatorStep4_WECSim_B.INPUT_1_1_1_n[3];
  tmp_5[3] = 12;
  simulationData->mData->mInputValues.mN = 12;
  simulationData->mData->mInputValues.mX = &tmp_4[0];
  simulationData->mData->mInputOffsets.mN = 4;
  simulationData->mData->mInputOffsets.mX = &tmp_5[0];
  simulationData->mData->mDx.mN = 35;
  simulationData->mData->mDx.mX = &_rtXdot->windEmulatorStep4_WECSimhptoS_h[0];
  simulator = static_cast<NeslSimulator *>
    (windEmulatorStep4_WECSim_DW.STATE_1_Simulator_i);
  diagnosticManager = static_cast<NeuDiagnosticManager *>
    (windEmulatorStep4_WECSim_DW.STATE_1_DiagMgr_g);
  diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
  is = ne_simulator_method(simulator, NESL_SIM_DERIVATIVES, simulationData,
    diagnosticManager);
  if (is != 0) {
    tmp = error_buffer_is_empty(rtmGetErrorStatus(windEmulatorStep4_WECSim_M));
    if (tmp) {
      msg = rtw_diagnostics_msg(diagnosticTree);
      rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
    }
  }
}

/* Model initialize function */
void windEmulatorStep4_WECSim_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&windEmulatorStep4_WECSim_M->solverInfo,
                          &windEmulatorStep4_WECSim_M->Timing.simTimeStep);
    rtsiSetTPtr(&windEmulatorStep4_WECSim_M->solverInfo, &rtmGetTPtr
                (windEmulatorStep4_WECSim_M));
    rtsiSetStepSizePtr(&windEmulatorStep4_WECSim_M->solverInfo,
                       &windEmulatorStep4_WECSim_M->Timing.stepSize0);
    rtsiSetdXPtr(&windEmulatorStep4_WECSim_M->solverInfo,
                 &windEmulatorStep4_WECSim_M->derivs);
    rtsiSetContStatesPtr(&windEmulatorStep4_WECSim_M->solverInfo, (real_T **)
                         &windEmulatorStep4_WECSim_M->contStates);
    rtsiSetNumContStatesPtr(&windEmulatorStep4_WECSim_M->solverInfo,
      &windEmulatorStep4_WECSim_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&windEmulatorStep4_WECSim_M->solverInfo,
      &windEmulatorStep4_WECSim_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&windEmulatorStep4_WECSim_M->solverInfo,
      &windEmulatorStep4_WECSim_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&windEmulatorStep4_WECSim_M->solverInfo,
      &windEmulatorStep4_WECSim_M->periodicContStateRanges);
    rtsiSetContStateDisabledPtr(&windEmulatorStep4_WECSim_M->solverInfo,
      (boolean_T**) &windEmulatorStep4_WECSim_M->contStateDisabled);
    rtsiSetErrorStatusPtr(&windEmulatorStep4_WECSim_M->solverInfo,
                          (&rtmGetErrorStatus(windEmulatorStep4_WECSim_M)));
    rtsiSetSolverMassMatrixIr(&windEmulatorStep4_WECSim_M->solverInfo,
      windEmulatorStep4_WECSim_MassMatrix.ir);
    rtsiSetSolverMassMatrixJc(&windEmulatorStep4_WECSim_M->solverInfo,
      windEmulatorStep4_WECSim_MassMatrix.jc);
    rtsiSetSolverMassMatrixPr(&windEmulatorStep4_WECSim_M->solverInfo,
      windEmulatorStep4_WECSim_MassMatrix.pr);
    rtsiSetRTModelPtr(&windEmulatorStep4_WECSim_M->solverInfo,
                      windEmulatorStep4_WECSim_M);
  }

  rtsiSetSimTimeStep(&windEmulatorStep4_WECSim_M->solverInfo, MAJOR_TIME_STEP);
  rtsiSetIsMinorTimeStepWithModeChange(&windEmulatorStep4_WECSim_M->solverInfo,
    false);
  rtsiSetIsContModeFrozen(&windEmulatorStep4_WECSim_M->solverInfo, false);
  windEmulatorStep4_WECSim_M->intgData.x0 = windEmulatorStep4_WECSim_M->odeX0;
  windEmulatorStep4_WECSim_M->intgData.f0 = windEmulatorStep4_WECSim_M->odeF0;
  windEmulatorStep4_WECSim_M->intgData.x1start =
    windEmulatorStep4_WECSim_M->odeX1START;
  windEmulatorStep4_WECSim_M->intgData.f1 = windEmulatorStep4_WECSim_M->odeF1;
  windEmulatorStep4_WECSim_M->intgData.Delta =
    windEmulatorStep4_WECSim_M->odeDELTA;
  windEmulatorStep4_WECSim_M->intgData.E = windEmulatorStep4_WECSim_M->odeE;
  windEmulatorStep4_WECSim_M->intgData.fac = windEmulatorStep4_WECSim_M->odeFAC;

  /* initialize */
  {
    int_T i;
    real_T *f = windEmulatorStep4_WECSim_M->intgData.fac;
    for (i = 0; i < static_cast<int_T>(sizeof(windEmulatorStep4_WECSim_M->odeFAC)/
          sizeof(real_T)); i++) {
      f[i] = 1.5e-8;
    }
  }

  windEmulatorStep4_WECSim_M->intgData.DFDX =
    windEmulatorStep4_WECSim_M->odeDFDX;
  windEmulatorStep4_WECSim_M->intgData.W = windEmulatorStep4_WECSim_M->odeW;
  windEmulatorStep4_WECSim_M->intgData.pivots =
    windEmulatorStep4_WECSim_M->odePIVOTS;
  windEmulatorStep4_WECSim_M->intgData.xtmp =
    windEmulatorStep4_WECSim_M->odeXTMP;
  windEmulatorStep4_WECSim_M->intgData.ztmp =
    windEmulatorStep4_WECSim_M->odeZTMP;
  windEmulatorStep4_WECSim_M->intgData.M =
    windEmulatorStep4_WECSim_M->odeMASSMATRIX_M;
  windEmulatorStep4_WECSim_M->intgData.M1 =
    windEmulatorStep4_WECSim_M->odeMASSMATRIX_M1;
  windEmulatorStep4_WECSim_M->intgData.xdot =
    windEmulatorStep4_WECSim_M->odeXDOT;
  windEmulatorStep4_WECSim_M->intgData.Edot =
    windEmulatorStep4_WECSim_M->odeEDOT;
  windEmulatorStep4_WECSim_M->intgData.fminusMxdot =
    windEmulatorStep4_WECSim_M->odeFMXDOT;
  windEmulatorStep4_WECSim_M->intgData.isFirstStep = true;
  rtsiSetSolverExtrapolationOrder(&windEmulatorStep4_WECSim_M->solverInfo, 4);
  rtsiSetSolverNumberNewtonIterations(&windEmulatorStep4_WECSim_M->solverInfo, 1);
  windEmulatorStep4_WECSim_M->contStates = ((X_windEmulatorStep4_WECSim_T *)
    &windEmulatorStep4_WECSim_X);
  windEmulatorStep4_WECSim_M->contStateDisabled =
    ((XDis_windEmulatorStep4_WECSim_T *) &windEmulatorStep4_WECSim_XDis);
  windEmulatorStep4_WECSim_M->Timing.tStart = (0.0);
  windEmulatorStep4_WECSim_M->massMatrixType = ((ssMatrixType)3);
  windEmulatorStep4_WECSim_M->massMatrixNzMax = (26);
  windEmulatorStep4_WECSim_M->massMatrixIr =
    (windEmulatorStep4_WECSim_MassMatrix.ir);
  windEmulatorStep4_WECSim_M->massMatrixJc =
    (windEmulatorStep4_WECSim_MassMatrix.jc);
  windEmulatorStep4_WECSim_M->massMatrixPr =
    (windEmulatorStep4_WECSim_MassMatrix.pr);
  rtsiSetSolverMassMatrixType(&windEmulatorStep4_WECSim_M->solverInfo,
    (ssMatrixType)3);
  rtsiSetSolverMassMatrixNzMax(&windEmulatorStep4_WECSim_M->solverInfo, 26);
  rtsiSetSolverData(&windEmulatorStep4_WECSim_M->solverInfo, static_cast<void *>
                    (&windEmulatorStep4_WECSim_M->intgData));
  rtsiSetSolverName(&windEmulatorStep4_WECSim_M->solverInfo,"ode14x");
  windEmulatorStep4_WECSim_M->solverInfoPtr =
    (&windEmulatorStep4_WECSim_M->solverInfo);

  /* Initialize timing info */
  {
    int_T *mdlTsMap = windEmulatorStep4_WECSim_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    windEmulatorStep4_WECSim_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    windEmulatorStep4_WECSim_M->Timing.sampleTimes =
      (&windEmulatorStep4_WECSim_M->Timing.sampleTimesArray[0]);
    windEmulatorStep4_WECSim_M->Timing.offsetTimes =
      (&windEmulatorStep4_WECSim_M->Timing.offsetTimesArray[0]);

    /* task periods */
    windEmulatorStep4_WECSim_M->Timing.sampleTimes[0] = (0.0);
    windEmulatorStep4_WECSim_M->Timing.sampleTimes[1] = (0.004);

    /* task offsets */
    windEmulatorStep4_WECSim_M->Timing.offsetTimes[0] = (0.0);
    windEmulatorStep4_WECSim_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(windEmulatorStep4_WECSim_M,
             &windEmulatorStep4_WECSim_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = windEmulatorStep4_WECSim_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    windEmulatorStep4_WECSim_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(windEmulatorStep4_WECSim_M, -1);
  windEmulatorStep4_WECSim_M->Timing.stepSize0 = 0.004;
  windEmulatorStep4_WECSim_M->Timing.stepSize1 = 0.004;
  windEmulatorStep4_WECSim_M->solverInfoPtr =
    (&windEmulatorStep4_WECSim_M->solverInfo);
  windEmulatorStep4_WECSim_M->Timing.stepSize = (0.004);
  rtsiSetFixedStepSize(&windEmulatorStep4_WECSim_M->solverInfo, 0.004);
  rtsiSetSolverMode(&windEmulatorStep4_WECSim_M->solverInfo,
                    SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  (void) std::memset((static_cast<void *>(&windEmulatorStep4_WECSim_B)), 0,
                     sizeof(B_windEmulatorStep4_WECSim_T));

  {
    windEmulatorStep4_WECSim_B.state_e = abbStateEnum_undefined;
    windEmulatorStep4_WECSim_B.state_ed = abbStateEnum_undefined;
    windEmulatorStep4_WECSim_B.expType_a = expTypeEnum_off;
    windEmulatorStep4_WECSim_B.toExpTypeEnum = expTypeEnum_off;
  }

  /* states (continuous) */
  {
    (void) std::memset(static_cast<void *>(&windEmulatorStep4_WECSim_X), 0,
                       sizeof(X_windEmulatorStep4_WECSim_T));
  }

  /* disabled states */
  {
    (void) std::memset(static_cast<void *>(&windEmulatorStep4_WECSim_XDis), 0,
                       sizeof(XDis_windEmulatorStep4_WECSim_T));
  }

  /* global mass matrix */
  {
    int_T *ir = windEmulatorStep4_WECSim_MassMatrix.ir;
    int_T *jc = windEmulatorStep4_WECSim_MassMatrix.jc;
    real_T *pr = windEmulatorStep4_WECSim_MassMatrix.pr;
    (void) std::memset(static_cast<void *>(ir), 0,
                       26*sizeof(int_T));
    (void) std::memset(static_cast<void *>(jc), 0,
                       (44+1)*sizeof(int_T));
    (void) std::memset(static_cast<void *>(pr), 0,
                       26*sizeof(real_T));
  }

  /* states (dwork) */
  (void) std::memset(static_cast<void *>(&windEmulatorStep4_WECSim_DW), 0,
                     sizeof(DW_windEmulatorStep4_WECSim_T));

  /* external inputs */
  (void)std::memset(&windEmulatorStep4_WECSim_U, 0, sizeof
                    (ExtU_windEmulatorStep4_WECSim_T));

  /* Root-level init GlobalMassMatrixPr offset */
  {
    windEmulatorStep4_WECSim_DW.STATE_1_MASS_MATRIX_PR = 9;/* '<S332>/STATE_1' */
  }

  /* child S-Function registration */
  {
    RTWSfcnInfo *sfcnInfo =
      &windEmulatorStep4_WECSim_M->NonInlinedSFcns.sfcnInfo;
    windEmulatorStep4_WECSim_M->sfcnInfo = (sfcnInfo);
    rtssSetErrorStatusPtr(sfcnInfo, (&rtmGetErrorStatus
      (windEmulatorStep4_WECSim_M)));
    windEmulatorStep4_WECSim_M->Sizes.numSampTimes = (2);
    rtssSetNumRootSampTimesPtr(sfcnInfo,
      &windEmulatorStep4_WECSim_M->Sizes.numSampTimes);
    windEmulatorStep4_WECSim_M->NonInlinedSFcns.taskTimePtrs[0] = (&rtmGetTPtr
      (windEmulatorStep4_WECSim_M)[0]);
    windEmulatorStep4_WECSim_M->NonInlinedSFcns.taskTimePtrs[1] = (&rtmGetTPtr
      (windEmulatorStep4_WECSim_M)[1]);
    rtssSetTPtrPtr(sfcnInfo,
                   windEmulatorStep4_WECSim_M->NonInlinedSFcns.taskTimePtrs);
    rtssSetTStartPtr(sfcnInfo, &rtmGetTStart(windEmulatorStep4_WECSim_M));
    rtssSetTFinalPtr(sfcnInfo, &rtmGetTFinal(windEmulatorStep4_WECSim_M));
    rtssSetTimeOfLastOutputPtr(sfcnInfo, &rtmGetTimeOfLastOutput
      (windEmulatorStep4_WECSim_M));
    rtssSetStepSizePtr(sfcnInfo, &windEmulatorStep4_WECSim_M->Timing.stepSize);
    rtssSetStopRequestedPtr(sfcnInfo, &rtmGetStopRequested
      (windEmulatorStep4_WECSim_M));
    rtssSetDerivCacheNeedsResetPtr(sfcnInfo,
      &windEmulatorStep4_WECSim_M->derivCacheNeedsReset);
    rtssSetZCCacheNeedsResetPtr(sfcnInfo,
      &windEmulatorStep4_WECSim_M->zCCacheNeedsReset);
    rtssSetContTimeOutputInconsistentWithStateAtMajorStepPtr(sfcnInfo,
      &windEmulatorStep4_WECSim_M->CTOutputIncnstWithState);
    rtssSetSampleHitsPtr(sfcnInfo,
                         &windEmulatorStep4_WECSim_M->Timing.sampleHits);
    rtssSetPerTaskSampleHitsPtr(sfcnInfo,
      &windEmulatorStep4_WECSim_M->Timing.perTaskSampleHits);
    rtssSetSimModePtr(sfcnInfo, &windEmulatorStep4_WECSim_M->simMode);
    rtssSetSolverInfoPtr(sfcnInfo, &windEmulatorStep4_WECSim_M->solverInfoPtr);
  }

  windEmulatorStep4_WECSim_M->Sizes.numSFcns = (1);

  /* register each child */
  {
    (void) std::memset(static_cast<void *>
                       (&windEmulatorStep4_WECSim_M->NonInlinedSFcns.childSFunctions
                        [0]), 0,
                       1*sizeof(SimStruct));
    windEmulatorStep4_WECSim_M->childSfunctions =
      (&windEmulatorStep4_WECSim_M->NonInlinedSFcns.childSFunctionPtrs[0]);
    windEmulatorStep4_WECSim_M->childSfunctions[0] =
      (&windEmulatorStep4_WECSim_M->NonInlinedSFcns.childSFunctions[0]);

    /* Level2 S-Function Block: windEmulatorStep4_WECSim/<S5>/Enable File Log (slrealtimeenablelogging) */
    {
      SimStruct *rts = windEmulatorStep4_WECSim_M->childSfunctions[0];

      /* timing info */
      time_T *sfcnPeriod =
        windEmulatorStep4_WECSim_M->NonInlinedSFcns.Sfcn0.sfcnPeriod;
      time_T *sfcnOffset =
        windEmulatorStep4_WECSim_M->NonInlinedSFcns.Sfcn0.sfcnOffset;
      int_T *sfcnTsMap =
        windEmulatorStep4_WECSim_M->NonInlinedSFcns.Sfcn0.sfcnTsMap;
      (void) std::memset(static_cast<void*>(sfcnPeriod), 0,
                         sizeof(time_T)*1);
      (void) std::memset(static_cast<void*>(sfcnOffset), 0,
                         sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      {
        ssSetBlkInfo2Ptr(rts,
                         &windEmulatorStep4_WECSim_M->NonInlinedSFcns.blkInfo2[0]);
      }

      _ssSetBlkInfo2PortInfo2Ptr(rts,
        &windEmulatorStep4_WECSim_M->NonInlinedSFcns.inputOutputPortInfo2[0]);

      /* Set up the mdlInfo pointer */
      ssSetRTWSfcnInfo(rts, windEmulatorStep4_WECSim_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &windEmulatorStep4_WECSim_M->
                           NonInlinedSFcns.methods2[0]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &windEmulatorStep4_WECSim_M->
                           NonInlinedSFcns.methods3[0]);
      }

      /* Allocate memory of model methods 4 */
      {
        ssSetModelMethods4(rts,
                           &windEmulatorStep4_WECSim_M->
                           NonInlinedSFcns.methods4[0]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &windEmulatorStep4_WECSim_M->NonInlinedSFcns.statesInfo2
                         [0]);
        ssSetPeriodicStatesInfo(rts,
          &windEmulatorStep4_WECSim_M->NonInlinedSFcns.periodicStatesInfo[0]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &windEmulatorStep4_WECSim_M->NonInlinedSFcns.Sfcn0.inputPortInfo[0]);
        ssSetPortInfoForInputs(rts,
          &windEmulatorStep4_WECSim_M->NonInlinedSFcns.Sfcn0.inputPortInfo[0]);
        _ssSetPortInfo2ForInputUnits(rts,
          &windEmulatorStep4_WECSim_M->NonInlinedSFcns.Sfcn0.inputPortUnits[0]);
        ssSetInputPortUnit(rts, 0, 0);
        _ssSetPortInfo2ForInputCoSimAttribute(rts,
          &windEmulatorStep4_WECSim_M->NonInlinedSFcns.Sfcn0.inputPortCoSimAttribute
          [0]);
        ssSetInputPortIsContinuousQuantity(rts, 0, 0);

        /* port 0 */
        {
          ssSetInputPortRequiredContiguous(rts, 0, 1);
          ssSetInputPortSignal(rts, 0, &windEmulatorStep4_WECSim_B.Constant);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidthAsInt(rts, 0, 1);
        }
      }

      /* path info */
      ssSetModelName(rts, "Enable File Log");
      ssSetPath(rts, "windEmulatorStep4_WECSim/fileAndUI/Enable File Log");
      ssSetRTModel(rts,windEmulatorStep4_WECSim_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* registration */
      slrealtimeenablelogging(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 0.004);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCsAsInt(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }
  }

  {
    NeModelParameters modelParameters;
    NeModelParameters modelParameters_0;
    NeModelParameters modelParameters_1;
    NeModelParameters modelParameters_2;
    NeModelParameters modelParameters_3;
    NeModelParameters modelParameters_4;
    NeModelParameters modelParameters_5;
    NeslRtpManager *manager;
    NeslSimulationData *tmp;
    NeslSimulator *simulator;
    NeuDiagnosticManager *diagnosticManager;
    NeuDiagnosticTree *diagnosticTree;
    char *msg;
    real_T tmp_0;
    int_T i;
    int_T startIdx;
    boolean_T tmp_1;
    boolean_T zcDisabled;

    /* Start for Constant: '<S2>/ACS880CtrlMode' */
    tmp_1 = *get_ctrlModeTorque();

    /* Start for S-Function (slecatinit): '<Root>/EtherCAT Init' */
    slrealtime::StartCallbackService::registerCB( std::bind
      ( Root_EtherCATInit_callback, nullptr ), 10 );

    /* Start for ToAsyncQueueBlock generated from: '<S27>/acs880Signals' */
    windEmulatorStep4_WECSim_DW.TAQSigLogging_InsertedFor_acs88.SLRTSigHandles =
      slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
      (&windEmulatorStep4_WECSim_B.BusAssignment_a));

    /* Start for Constant: '<S2>/ACS880CtrlMode' */
    windEmulatorStep4_WECSim_B.ACS880CtrlMode = tmp_1;

    /* Start for Constant: '<S4>/expType' */
    windEmulatorStep4_WECSim_B.expType_a =
      windEmulatorStep4_WECSim_cal->expType_Value;

    /* Start for Constant: '<S4>/expRunTime' */
    windEmulatorStep4_WECSim_B.expRunTime =
      windEmulatorStep4_WECSim_cal->expRunTime_Value;

    /* Start for SimscapeRtp: '<S508>/RTP_1' */
    manager = nesl_lease_rtp_manager(
      "windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Solver Configuration_1",
      0);
    zcDisabled = pointer_is_null(manager);
    if (zcDisabled) {
      windEmulatorStep4_WECSim_1e9c788f_1_gateway();
      manager = nesl_lease_rtp_manager(
        "windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Solver Configuration_1",
        0);
    }

    windEmulatorStep4_WECSim_DW.RTP_1_RtpManager = (void *)manager;
    windEmulatorStep4_WECSim_DW.RTP_1_SetParametersNeeded = true;

    /* End of Start for SimscapeRtp: '<S508>/RTP_1' */

    /* Start for SimscapeExecutionBlock: '<S541>/STATE_1' */
    simulator = nesl_lease_simulator(
      "windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Solver Configuration_1",
      0, 0);
    windEmulatorStep4_WECSim_DW.STATE_1_Simulator = (void *)simulator;
    zcDisabled = pointer_is_null(windEmulatorStep4_WECSim_DW.STATE_1_Simulator);
    if (zcDisabled) {
      windEmulatorStep4_WECSim_1e9c788f_1_gateway();
      simulator = nesl_lease_simulator(
        "windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Solver Configuration_1",
        0, 0);
      windEmulatorStep4_WECSim_DW.STATE_1_Simulator = (void *)simulator;
    }

    tmp = nesl_create_simulation_data();
    windEmulatorStep4_WECSim_DW.STATE_1_SimData = (void *)tmp;
    diagnosticManager = rtw_create_diagnostics();
    windEmulatorStep4_WECSim_DW.STATE_1_DiagMgr = (void *)diagnosticManager;
    modelParameters.mSolverType = NE_SOLVER_TYPE_DAE;
    modelParameters.mSolverTolerance = 0.001;
    modelParameters.mSolverAbsTol = 0.001;
    modelParameters.mSolverRelTol = 0.001;
    modelParameters.mVariableStepSolver = false;
    modelParameters.mIsUsingODEN = false;
    modelParameters.mSolverModifyAbsTol = NE_MODIFY_ABS_TOL_NO;
    modelParameters.mFixedStepSize = 0.004;
    modelParameters.mStartTime = 0.0;
    modelParameters.mLoadInitialState = false;
    modelParameters.mUseSimState = false;
    modelParameters.mLinTrimCompile = false;
    modelParameters.mLoggingMode = SSC_LOGGING_OFF;
    modelParameters.mRTWModifiedTimeStamp = 7.06642264E+8;
    modelParameters.mZcDisabled = true;
    modelParameters.mUseModelRefSolver = false;
    modelParameters.mTargetFPGAHIL = false;
    tmp_0 = 0.001;
    modelParameters.mSolverTolerance = tmp_0;
    tmp_0 = 0.004;
    modelParameters.mFixedStepSize = tmp_0;
    zcDisabled = false;
    modelParameters.mVariableStepSolver = zcDisabled;
    zcDisabled = false;
    modelParameters.mIsUsingODEN = zcDisabled;
    modelParameters.mZcDisabled = true;
    simulator = static_cast<NeslSimulator *>
      (windEmulatorStep4_WECSim_DW.STATE_1_Simulator);
    diagnosticManager = static_cast<NeuDiagnosticManager *>
      (windEmulatorStep4_WECSim_DW.STATE_1_DiagMgr);
    diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
    startIdx = nesl_initialize_simulator(simulator, &modelParameters,
      diagnosticManager);
    if (startIdx != 0) {
      zcDisabled = error_buffer_is_empty(rtmGetErrorStatus
        (windEmulatorStep4_WECSim_M));
      if (zcDisabled) {
        msg = rtw_diagnostics_msg(diagnosticTree);
        rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
      }
    }

    /* End of Start for SimscapeExecutionBlock: '<S541>/STATE_1' */

    /* Start for SimscapeExecutionBlock: '<S541>/OUTPUT_1_0' */
    simulator = nesl_lease_simulator(
      "windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Solver Configuration_1",
      1, 0);
    windEmulatorStep4_WECSim_DW.OUTPUT_1_0_Simulator = (void *)simulator;
    zcDisabled = pointer_is_null
      (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_Simulator);
    if (zcDisabled) {
      windEmulatorStep4_WECSim_1e9c788f_1_gateway();
      simulator = nesl_lease_simulator(
        "windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Solver Configuration_1",
        1, 0);
      windEmulatorStep4_WECSim_DW.OUTPUT_1_0_Simulator = (void *)simulator;
    }

    tmp = nesl_create_simulation_data();
    windEmulatorStep4_WECSim_DW.OUTPUT_1_0_SimData = (void *)tmp;
    diagnosticManager = rtw_create_diagnostics();
    windEmulatorStep4_WECSim_DW.OUTPUT_1_0_DiagMgr = (void *)diagnosticManager;
    modelParameters_0.mSolverType = NE_SOLVER_TYPE_DAE;
    modelParameters_0.mSolverTolerance = 0.001;
    modelParameters_0.mSolverAbsTol = 0.001;
    modelParameters_0.mSolverRelTol = 0.001;
    modelParameters_0.mVariableStepSolver = false;
    modelParameters_0.mIsUsingODEN = false;
    modelParameters_0.mSolverModifyAbsTol = NE_MODIFY_ABS_TOL_NO;
    modelParameters_0.mFixedStepSize = 0.004;
    modelParameters_0.mStartTime = 0.0;
    modelParameters_0.mLoadInitialState = false;
    modelParameters_0.mUseSimState = false;
    modelParameters_0.mLinTrimCompile = false;
    modelParameters_0.mLoggingMode = SSC_LOGGING_OFF;
    modelParameters_0.mRTWModifiedTimeStamp = 7.06642264E+8;
    modelParameters_0.mZcDisabled = true;
    modelParameters_0.mUseModelRefSolver = false;
    modelParameters_0.mTargetFPGAHIL = false;
    tmp_0 = 0.001;
    modelParameters_0.mSolverTolerance = tmp_0;
    tmp_0 = 0.004;
    modelParameters_0.mFixedStepSize = tmp_0;
    zcDisabled = false;
    modelParameters_0.mVariableStepSolver = zcDisabled;
    zcDisabled = false;
    modelParameters_0.mIsUsingODEN = zcDisabled;
    modelParameters_0.mZcDisabled = true;
    simulator = static_cast<NeslSimulator *>
      (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_Simulator);
    diagnosticManager = static_cast<NeuDiagnosticManager *>
      (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_DiagMgr);
    diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
    startIdx = nesl_initialize_simulator(simulator, &modelParameters_0,
      diagnosticManager);
    if (startIdx != 0) {
      zcDisabled = error_buffer_is_empty(rtmGetErrorStatus
        (windEmulatorStep4_WECSim_M));
      if (zcDisabled) {
        msg = rtw_diagnostics_msg(diagnosticTree);
        rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
      }
    }

    /* End of Start for SimscapeExecutionBlock: '<S541>/OUTPUT_1_0' */

    /* Start for ToAsyncQueueBlock generated from: '<S25>/acs880CtrlSignals' */
    windEmulatorStep4_WECSim_DW.TAQSigLogging_InsertedFor_acs_a.SLRTSigHandles =
      slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
      (&windEmulatorStep4_WECSim_B.BusAssignment_k));

    /* Start for ToAsyncQueueBlock generated from: '<S23>/acs800Signals' */
    windEmulatorStep4_WECSim_DW.TAQSigLogging_InsertedFor_acs80.SLRTSigHandles =
      slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
      (&windEmulatorStep4_WECSim_B.BusAssignment_h));

    /* Start for Constant: '<S1>/ACS800CtrlMode' */
    windEmulatorStep4_WECSim_B.ACS800CtrlMode = tmp_1;

    /* Start for ToAsyncQueueBlock generated from: '<S21>/acs800CtrlSignals' */
    windEmulatorStep4_WECSim_DW.TAQSigLogging_InsertedFor_acs_l.SLRTSigHandles =
      slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
      (&windEmulatorStep4_WECSim_B.BusAssignment_kc));

    /* Start for ToAsyncQueueBlock generated from: '<S34>/hptoSignals' */
    windEmulatorStep4_WECSim_DW.TAQSigLogging_InsertedFor_hptoS.SLRTSigHandles =
      slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
      (&windEmulatorStep4_WECSim_B.BusAssignment_n));

    /* Start for ToAsyncQueueBlock generated from: '<S32>/hptoCtrl' */
    windEmulatorStep4_WECSim_DW.TAQSigLogging_InsertedFor_hptoC.SLRTSigHandles =
      slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
      (&windEmulatorStep4_WECSim_B.BusAssignment_c));

    /* Start for ToAsyncQueueBlock generated from: '<S30>/expCtrlSignals' */
    windEmulatorStep4_WECSim_DW.TAQSigLogging_InsertedFor_expCt.SLRTSigHandles =
      slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
      (&windEmulatorStep4_WECSim_B.BusAssignment_b));

    /* Start for ToAsyncQueueBlock generated from: '<S38>/shaftSignals' */
    windEmulatorStep4_WECSim_DW.TAQSigLogging_InsertedFor_shaft.SLRTSigHandles =
      slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
      (&windEmulatorStep4_WECSim_B.BusAssignment_g));

    /* Start for ToAsyncQueueBlock generated from: '<S36>/invPowerAcs800' */
    windEmulatorStep4_WECSim_DW.TAQSigLogging_InsertedFor_invPo.SLRTSigHandles =
      slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
      (&windEmulatorStep4_WECSim_B.BusAssignment));

    /* Start for ToAsyncQueueBlock generated from: '<S37>/invPowerAcs880' */
    windEmulatorStep4_WECSim_DW.TAQSigLogging_InsertedFor_inv_p.SLRTSigHandles =
      slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
      (&windEmulatorStep4_WECSim_B.BusAssignment_l));

    /* Start for ToAsyncQueueBlock generated from: '<S40>/sidInfoSignals' */
    windEmulatorStep4_WECSim_DW.TAQSigLogging_InsertedFor_sidIn.SLRTSigHandles =
      slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
      (&windEmulatorStep4_WECSim_B.BusAssignment_j));

    /* Start for Constant: '<S5>/Constant' */
    windEmulatorStep4_WECSim_B.Constant =
      windEmulatorStep4_WECSim_cal->Constant_Value_m3;

    /* Start for S-Function (slrealtimeenablelogging): '<S5>/Enable File Log' */
    /* Level2 S-Function Block: '<S5>/Enable File Log' (slrealtimeenablelogging) */
    {
      SimStruct *rts = windEmulatorStep4_WECSim_M->childSfunctions[0];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Start for SimscapeExecutionBlock: '<S216>/STATE_1' */
    simulator = nesl_lease_simulator(
      "windEmulatorStep4_WECSim/hptoSim/WECSimModel/Global Reference Frame/Solver Configuration_1",
      0, 0);
    windEmulatorStep4_WECSim_DW.STATE_1_Simulator_f = (void *)simulator;
    zcDisabled = pointer_is_null(windEmulatorStep4_WECSim_DW.STATE_1_Simulator_f);
    if (zcDisabled) {
      windEmulatorStep4_WECSim_dfbb7ac7_1_gateway();
      simulator = nesl_lease_simulator(
        "windEmulatorStep4_WECSim/hptoSim/WECSimModel/Global Reference Frame/Solver Configuration_1",
        0, 0);
      windEmulatorStep4_WECSim_DW.STATE_1_Simulator_f = (void *)simulator;
    }

    tmp = nesl_create_simulation_data();
    windEmulatorStep4_WECSim_DW.STATE_1_SimData_h = (void *)tmp;
    diagnosticManager = rtw_create_diagnostics();
    windEmulatorStep4_WECSim_DW.STATE_1_DiagMgr_o = (void *)diagnosticManager;
    modelParameters_1.mSolverType = NE_SOLVER_TYPE_DAE;
    modelParameters_1.mSolverTolerance = 0.001;
    modelParameters_1.mSolverAbsTol = 0.001;
    modelParameters_1.mSolverRelTol = 0.001;
    modelParameters_1.mVariableStepSolver = false;
    modelParameters_1.mIsUsingODEN = false;
    modelParameters_1.mSolverModifyAbsTol = NE_MODIFY_ABS_TOL_NO;
    modelParameters_1.mFixedStepSize = 0.004;
    modelParameters_1.mStartTime = 0.0;
    modelParameters_1.mLoadInitialState = false;
    modelParameters_1.mUseSimState = false;
    modelParameters_1.mLinTrimCompile = false;
    modelParameters_1.mLoggingMode = SSC_LOGGING_OFF;
    modelParameters_1.mRTWModifiedTimeStamp = 7.06642264E+8;
    modelParameters_1.mZcDisabled = true;
    modelParameters_1.mUseModelRefSolver = false;
    modelParameters_1.mTargetFPGAHIL = false;
    tmp_0 = 0.001;
    modelParameters_1.mSolverTolerance = tmp_0;
    tmp_0 = 0.004;
    modelParameters_1.mFixedStepSize = tmp_0;
    zcDisabled = false;
    modelParameters_1.mVariableStepSolver = zcDisabled;
    zcDisabled = false;
    modelParameters_1.mIsUsingODEN = zcDisabled;
    modelParameters_1.mZcDisabled = true;
    simulator = static_cast<NeslSimulator *>
      (windEmulatorStep4_WECSim_DW.STATE_1_Simulator_f);
    diagnosticManager = static_cast<NeuDiagnosticManager *>
      (windEmulatorStep4_WECSim_DW.STATE_1_DiagMgr_o);
    diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
    startIdx = nesl_initialize_simulator(simulator, &modelParameters_1,
      diagnosticManager);
    if (startIdx != 0) {
      zcDisabled = error_buffer_is_empty(rtmGetErrorStatus
        (windEmulatorStep4_WECSim_M));
      if (zcDisabled) {
        msg = rtw_diagnostics_msg(diagnosticTree);
        rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
      }
    }

    /* End of Start for SimscapeExecutionBlock: '<S216>/STATE_1' */

    /* Start for SimscapeExecutionBlock: '<S216>/OUTPUT_1_1' */
    simulator = nesl_lease_simulator(
      "windEmulatorStep4_WECSim/hptoSim/WECSimModel/Global Reference Frame/Solver Configuration_1",
      1, 1);
    windEmulatorStep4_WECSim_DW.OUTPUT_1_1_Simulator = (void *)simulator;
    zcDisabled = pointer_is_null
      (windEmulatorStep4_WECSim_DW.OUTPUT_1_1_Simulator);
    if (zcDisabled) {
      windEmulatorStep4_WECSim_dfbb7ac7_1_gateway();
      simulator = nesl_lease_simulator(
        "windEmulatorStep4_WECSim/hptoSim/WECSimModel/Global Reference Frame/Solver Configuration_1",
        1, 1);
      windEmulatorStep4_WECSim_DW.OUTPUT_1_1_Simulator = (void *)simulator;
    }

    tmp = nesl_create_simulation_data();
    windEmulatorStep4_WECSim_DW.OUTPUT_1_1_SimData = (void *)tmp;
    diagnosticManager = rtw_create_diagnostics();
    windEmulatorStep4_WECSim_DW.OUTPUT_1_1_DiagMgr = (void *)diagnosticManager;
    modelParameters_2.mSolverType = NE_SOLVER_TYPE_DAE;
    modelParameters_2.mSolverTolerance = 0.001;
    modelParameters_2.mSolverAbsTol = 0.001;
    modelParameters_2.mSolverRelTol = 0.001;
    modelParameters_2.mVariableStepSolver = false;
    modelParameters_2.mIsUsingODEN = false;
    modelParameters_2.mSolverModifyAbsTol = NE_MODIFY_ABS_TOL_NO;
    modelParameters_2.mFixedStepSize = 0.004;
    modelParameters_2.mStartTime = 0.0;
    modelParameters_2.mLoadInitialState = false;
    modelParameters_2.mUseSimState = false;
    modelParameters_2.mLinTrimCompile = false;
    modelParameters_2.mLoggingMode = SSC_LOGGING_OFF;
    modelParameters_2.mRTWModifiedTimeStamp = 7.06642264E+8;
    modelParameters_2.mZcDisabled = true;
    modelParameters_2.mUseModelRefSolver = false;
    modelParameters_2.mTargetFPGAHIL = false;
    tmp_0 = 0.001;
    modelParameters_2.mSolverTolerance = tmp_0;
    tmp_0 = 0.004;
    modelParameters_2.mFixedStepSize = tmp_0;
    zcDisabled = false;
    modelParameters_2.mVariableStepSolver = zcDisabled;
    zcDisabled = false;
    modelParameters_2.mIsUsingODEN = zcDisabled;
    modelParameters_2.mZcDisabled = true;
    simulator = static_cast<NeslSimulator *>
      (windEmulatorStep4_WECSim_DW.OUTPUT_1_1_Simulator);
    diagnosticManager = static_cast<NeuDiagnosticManager *>
      (windEmulatorStep4_WECSim_DW.OUTPUT_1_1_DiagMgr);
    diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
    startIdx = nesl_initialize_simulator(simulator, &modelParameters_2,
      diagnosticManager);
    if (startIdx != 0) {
      zcDisabled = error_buffer_is_empty(rtmGetErrorStatus
        (windEmulatorStep4_WECSim_M));
      if (zcDisabled) {
        msg = rtw_diagnostics_msg(diagnosticTree);
        rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
      }
    }

    /* End of Start for SimscapeExecutionBlock: '<S216>/OUTPUT_1_1' */

    /* Start for SimscapeExecutionBlock: '<S332>/STATE_1' */
    simulator = nesl_lease_simulator(
      "windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Solver Configuration_1",
      0, 0);
    windEmulatorStep4_WECSim_DW.STATE_1_Simulator_i = (void *)simulator;
    zcDisabled = pointer_is_null(windEmulatorStep4_WECSim_DW.STATE_1_Simulator_i);
    if (zcDisabled) {
      windEmulatorStep4_WECSim_5bdcd402_1_gateway();
      simulator = nesl_lease_simulator(
        "windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Solver Configuration_1",
        0, 0);
      windEmulatorStep4_WECSim_DW.STATE_1_Simulator_i = (void *)simulator;
    }

    tmp = nesl_create_simulation_data();
    windEmulatorStep4_WECSim_DW.STATE_1_SimData_a = (void *)tmp;
    diagnosticManager = rtw_create_diagnostics();
    windEmulatorStep4_WECSim_DW.STATE_1_DiagMgr_g = (void *)diagnosticManager;
    modelParameters_3.mSolverType = NE_SOLVER_TYPE_DAE;
    modelParameters_3.mSolverTolerance = 0.001;
    modelParameters_3.mSolverAbsTol = 0.001;
    modelParameters_3.mSolverRelTol = 0.001;
    modelParameters_3.mVariableStepSolver = false;
    modelParameters_3.mIsUsingODEN = false;
    modelParameters_3.mSolverModifyAbsTol = NE_MODIFY_ABS_TOL_NO;
    modelParameters_3.mFixedStepSize = 0.004;
    modelParameters_3.mStartTime = 0.0;
    modelParameters_3.mLoadInitialState = false;
    modelParameters_3.mUseSimState = false;
    modelParameters_3.mLinTrimCompile = false;
    modelParameters_3.mLoggingMode = SSC_LOGGING_OFF;
    modelParameters_3.mRTWModifiedTimeStamp = 7.06642264E+8;
    modelParameters_3.mZcDisabled = true;
    modelParameters_3.mUseModelRefSolver = false;
    modelParameters_3.mTargetFPGAHIL = false;
    tmp_0 = 0.001;
    modelParameters_3.mSolverTolerance = tmp_0;
    tmp_0 = 0.004;
    modelParameters_3.mFixedStepSize = tmp_0;
    zcDisabled = false;
    modelParameters_3.mVariableStepSolver = zcDisabled;
    zcDisabled = false;
    modelParameters_3.mIsUsingODEN = zcDisabled;
    modelParameters_3.mZcDisabled = true;
    simulator = static_cast<NeslSimulator *>
      (windEmulatorStep4_WECSim_DW.STATE_1_Simulator_i);
    diagnosticManager = static_cast<NeuDiagnosticManager *>
      (windEmulatorStep4_WECSim_DW.STATE_1_DiagMgr_g);
    diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
    startIdx = nesl_initialize_simulator(simulator, &modelParameters_3,
      diagnosticManager);
    if (startIdx != 0) {
      zcDisabled = error_buffer_is_empty(rtmGetErrorStatus
        (windEmulatorStep4_WECSim_M));
      if (zcDisabled) {
        msg = rtw_diagnostics_msg(diagnosticTree);
        rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
      }
    }

    /* End of Start for SimscapeExecutionBlock: '<S332>/STATE_1' */

    /* Start for SimscapeExecutionBlock: '<S332>/OUTPUT_1_0' */
    simulator = nesl_lease_simulator(
      "windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Solver Configuration_1",
      1, 0);
    windEmulatorStep4_WECSim_DW.OUTPUT_1_0_Simulator_d = (void *)simulator;
    zcDisabled = pointer_is_null
      (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_Simulator_d);
    if (zcDisabled) {
      windEmulatorStep4_WECSim_5bdcd402_1_gateway();
      simulator = nesl_lease_simulator(
        "windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Solver Configuration_1",
        1, 0);
      windEmulatorStep4_WECSim_DW.OUTPUT_1_0_Simulator_d = (void *)simulator;
    }

    tmp = nesl_create_simulation_data();
    windEmulatorStep4_WECSim_DW.OUTPUT_1_0_SimData_l = (void *)tmp;
    diagnosticManager = rtw_create_diagnostics();
    windEmulatorStep4_WECSim_DW.OUTPUT_1_0_DiagMgr_i = (void *)diagnosticManager;
    modelParameters_4.mSolverType = NE_SOLVER_TYPE_DAE;
    modelParameters_4.mSolverTolerance = 0.001;
    modelParameters_4.mSolverAbsTol = 0.001;
    modelParameters_4.mSolverRelTol = 0.001;
    modelParameters_4.mVariableStepSolver = false;
    modelParameters_4.mIsUsingODEN = false;
    modelParameters_4.mSolverModifyAbsTol = NE_MODIFY_ABS_TOL_NO;
    modelParameters_4.mFixedStepSize = 0.004;
    modelParameters_4.mStartTime = 0.0;
    modelParameters_4.mLoadInitialState = false;
    modelParameters_4.mUseSimState = false;
    modelParameters_4.mLinTrimCompile = false;
    modelParameters_4.mLoggingMode = SSC_LOGGING_OFF;
    modelParameters_4.mRTWModifiedTimeStamp = 7.06642264E+8;
    modelParameters_4.mZcDisabled = true;
    modelParameters_4.mUseModelRefSolver = false;
    modelParameters_4.mTargetFPGAHIL = false;
    tmp_0 = 0.001;
    modelParameters_4.mSolverTolerance = tmp_0;
    tmp_0 = 0.004;
    modelParameters_4.mFixedStepSize = tmp_0;
    zcDisabled = false;
    modelParameters_4.mVariableStepSolver = zcDisabled;
    zcDisabled = false;
    modelParameters_4.mIsUsingODEN = zcDisabled;
    modelParameters_4.mZcDisabled = true;
    simulator = static_cast<NeslSimulator *>
      (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_Simulator_d);
    diagnosticManager = static_cast<NeuDiagnosticManager *>
      (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_DiagMgr_i);
    diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
    startIdx = nesl_initialize_simulator(simulator, &modelParameters_4,
      diagnosticManager);
    if (startIdx != 0) {
      zcDisabled = error_buffer_is_empty(rtmGetErrorStatus
        (windEmulatorStep4_WECSim_M));
      if (zcDisabled) {
        msg = rtw_diagnostics_msg(diagnosticTree);
        rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
      }
    }

    /* End of Start for SimscapeExecutionBlock: '<S332>/OUTPUT_1_0' */
    /* Start for TransportDelay: '<S60>/Transport Delay' */
    windEmulatorStep4_WECSim_DW.TransportDelay_RWORK[0] = 0.0;
    startIdx = 1;
    for (i = 0; i < 6; i++) {
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK[i] = 0;
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK[i + 6] = 0;
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK[i + 12] = 0;
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK[i + 18] = 1024;
      windEmulatorStep4_WECSim_DW.TransportDelay_PWORK[i] =
        &windEmulatorStep4_WECSim_DW.TransportDelay_RWORK[startIdx];
      startIdx += 2048;
      (static_cast<real_T *>(windEmulatorStep4_WECSim_DW.TransportDelay_PWORK[i]))
        [0] = windEmulatorStep4_WECSim_cal->TransportDelay_InitOutput;
      (static_cast<real_T *>(windEmulatorStep4_WECSim_DW.TransportDelay_PWORK[i]))
        [1024] = windEmulatorStep4_WECSim_M->Timing.t[0];
    }

    /* End of Start for TransportDelay: '<S60>/Transport Delay' */

    /* Start for TransportDelay: '<S139>/Transport Delay' */
    windEmulatorStep4_WECSim_DW.TransportDelay_RWORK_k[0] = 0.0;
    startIdx = 1;
    for (i = 0; i < 6; i++) {
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c[i] = 0;
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c[i + 6] = 0;
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c[i + 12] = 0;
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c[i + 18] = 1024;
      windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_f[i] =
        &windEmulatorStep4_WECSim_DW.TransportDelay_RWORK_k[startIdx];
      startIdx += 2048;
      (static_cast<real_T *>
        (windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_f[i]))[0] =
        windEmulatorStep4_WECSim_cal->TransportDelay_InitOutput_m;
      (static_cast<real_T *>
        (windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_f[i]))[1024] =
        windEmulatorStep4_WECSim_M->Timing.t[0];
    }

    /* End of Start for TransportDelay: '<S139>/Transport Delay' */

    /* Start for SimscapeExecutionBlock: '<S216>/OUTPUT_1_0' */
    simulator = nesl_lease_simulator(
      "windEmulatorStep4_WECSim/hptoSim/WECSimModel/Global Reference Frame/Solver Configuration_1",
      1, 0);
    windEmulatorStep4_WECSim_DW.OUTPUT_1_0_Simulator_l = (void *)simulator;
    zcDisabled = pointer_is_null
      (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_Simulator_l);
    if (zcDisabled) {
      windEmulatorStep4_WECSim_dfbb7ac7_1_gateway();
      simulator = nesl_lease_simulator(
        "windEmulatorStep4_WECSim/hptoSim/WECSimModel/Global Reference Frame/Solver Configuration_1",
        1, 0);
      windEmulatorStep4_WECSim_DW.OUTPUT_1_0_Simulator_l = (void *)simulator;
    }

    tmp = nesl_create_simulation_data();
    windEmulatorStep4_WECSim_DW.OUTPUT_1_0_SimData_i = (void *)tmp;
    diagnosticManager = rtw_create_diagnostics();
    windEmulatorStep4_WECSim_DW.OUTPUT_1_0_DiagMgr_a = (void *)diagnosticManager;
    modelParameters_5.mSolverType = NE_SOLVER_TYPE_DAE;
    modelParameters_5.mSolverTolerance = 0.001;
    modelParameters_5.mSolverAbsTol = 0.001;
    modelParameters_5.mSolverRelTol = 0.001;
    modelParameters_5.mVariableStepSolver = false;
    modelParameters_5.mIsUsingODEN = false;
    modelParameters_5.mSolverModifyAbsTol = NE_MODIFY_ABS_TOL_NO;
    modelParameters_5.mFixedStepSize = 0.004;
    modelParameters_5.mStartTime = 0.0;
    modelParameters_5.mLoadInitialState = false;
    modelParameters_5.mUseSimState = false;
    modelParameters_5.mLinTrimCompile = false;
    modelParameters_5.mLoggingMode = SSC_LOGGING_OFF;
    modelParameters_5.mRTWModifiedTimeStamp = 7.06642264E+8;
    modelParameters_5.mZcDisabled = true;
    modelParameters_5.mUseModelRefSolver = false;
    modelParameters_5.mTargetFPGAHIL = false;
    tmp_0 = 0.001;
    modelParameters_5.mSolverTolerance = tmp_0;
    tmp_0 = 0.004;
    modelParameters_5.mFixedStepSize = tmp_0;
    zcDisabled = false;
    modelParameters_5.mVariableStepSolver = zcDisabled;
    zcDisabled = false;
    modelParameters_5.mIsUsingODEN = zcDisabled;
    modelParameters_5.mZcDisabled = true;
    simulator = static_cast<NeslSimulator *>
      (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_Simulator_l);
    diagnosticManager = static_cast<NeuDiagnosticManager *>
      (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_DiagMgr_a);
    diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
    startIdx = nesl_initialize_simulator(simulator, &modelParameters_5,
      diagnosticManager);
    if (startIdx != 0) {
      zcDisabled = error_buffer_is_empty(rtmGetErrorStatus
        (windEmulatorStep4_WECSim_M));
      if (zcDisabled) {
        msg = rtw_diagnostics_msg(diagnosticTree);
        rtmSetErrorStatus(windEmulatorStep4_WECSim_M, msg);
      }
    }

    /* End of Start for SimscapeExecutionBlock: '<S216>/OUTPUT_1_0' */

    /* Start for TransportDelay: '<S64>/Transport Delay' */
    windEmulatorStep4_WECSim_DW.TransportDelay_RWORK_l[0] = 0.0;
    startIdx = 1;
    for (i = 0; i < 6; i++) {
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_o[i] = 0;
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_o[i + 6] = 0;
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_o[i + 12] = 0;
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_o[i + 18] = 1024;
      windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_j[i] =
        &windEmulatorStep4_WECSim_DW.TransportDelay_RWORK_l[startIdx];
      startIdx += 2048;
      (static_cast<real_T *>
        (windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_j[i]))[0] =
        windEmulatorStep4_WECSim_cal->TransportDelay_InitOutput_b;
      (static_cast<real_T *>
        (windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_j[i]))[1024] =
        windEmulatorStep4_WECSim_M->Timing.t[0];
    }

    /* End of Start for TransportDelay: '<S64>/Transport Delay' */

    /* Start for TransportDelay: '<S143>/Transport Delay' */
    windEmulatorStep4_WECSim_DW.TransportDelay_RWORK_f[0] = 0.0;
    startIdx = 1;
    for (i = 0; i < 6; i++) {
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c1[i] = 0;
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c1[i + 6] = 0;
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c1[i + 12] = 0;
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c1[i + 18] = 1024;
      windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_n[i] =
        &windEmulatorStep4_WECSim_DW.TransportDelay_RWORK_f[startIdx];
      startIdx += 2048;
      (static_cast<real_T *>
        (windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_n[i]))[0] =
        windEmulatorStep4_WECSim_cal->TransportDelay_InitOutput_g;
      (static_cast<real_T *>
        (windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_n[i]))[1024] =
        windEmulatorStep4_WECSim_M->Timing.t[0];
    }

    /* End of Start for TransportDelay: '<S143>/Transport Delay' */
  }

  windEmulatorStep4_WECSim_PrevZCX.Integrator_Reset_ZCE = UNINITIALIZED_ZCSIG;

  {
    int_T is;
    uint32_T is_UpdateStateMachine_tmp;
    boolean_T tmp;
    boolean_T tmp_0;

    /* InitializeConditions for Memory: '<S2>/Memory' */
    windEmulatorStep4_WECSim_DW.Memory_PreviousInput_i =
      windEmulatorStep4_WECSim_cal->Memory_InitialCondition_a;

    /* InitializeConditions for Memory: '<S2>/Memory1' */
    windEmulatorStep4_WECSim_DW.Memory1_PreviousInput =
      windEmulatorStep4_WECSim_cal->Memory1_InitialCondition;

    /* InitializeConditions for Memory: '<S2>/Memory2' */
    windEmulatorStep4_WECSim_DW.Memory2_PreviousInput =
      windEmulatorStep4_WECSim_cal->Memory2_InitialCondition;

    /* InitializeConditions for Memory: '<S4>/Memory' */
    windEmulatorStep4_WECSim_DW.Memory_PreviousInput_k =
      windEmulatorStep4_WECSim_cal->Memory_InitialCondition_d;

    /* InitializeConditions for Memory: '<S4>/Memory1' */
    windEmulatorStep4_WECSim_DW.Memory1_PreviousInput_d =
      windEmulatorStep4_WECSim_cal->Memory1_InitialCondition_h;

    /* InitializeConditions for Memory: '<S4>/Memory2' */
    windEmulatorStep4_WECSim_DW.Memory2_PreviousInput_l =
      windEmulatorStep4_WECSim_cal->Memory2_InitialCondition_i;

    /* InitializeConditions for RateLimiter: '<S609>/torqueSlewRate' */
    windEmulatorStep4_WECSim_DW.LastMajorTime = (rtInf);

    /* InitializeConditions for RateLimiter: '<S609>/speedSlewRate' */
    windEmulatorStep4_WECSim_DW.LastMajorTime_j = (rtInf);

    /* InitializeConditions for Integrator: '<S590>/Integrator' */
    windEmulatorStep4_WECSim_X.Integrator_CSTATE =
      windEmulatorStep4_WECSim_cal->PIDController_InitialConditio_o;

    /* InitializeConditions for RateLimiter: '<S498>/Rate Limiter' */
    windEmulatorStep4_WECSim_DW.PrevY_k =
      windEmulatorStep4_WECSim_cal->RateLimiter_IC;

    /* InitializeConditions for StateSpace: '<S531>/Internal' */
    windEmulatorStep4_WECSim_X.Internal_CSTATE[0] =
      windEmulatorStep4_WECSim_cal->Internal_InitialCondition;
    windEmulatorStep4_WECSim_X.Internal_CSTATE[1] =
      windEmulatorStep4_WECSim_cal->Internal_InitialCondition;
    windEmulatorStep4_WECSim_X.Internal_CSTATE[2] =
      windEmulatorStep4_WECSim_cal->Internal_InitialCondition;

    /* InitializeConditions for StateSpace: '<S545>/Internal' */
    windEmulatorStep4_WECSim_X.Internal_CSTATE_j =
      windEmulatorStep4_WECSim_cal->Internal_InitialCondition_p;

    /* InitializeConditions for RateLimiter: '<S7>/Rate Limiter' */
    windEmulatorStep4_WECSim_DW.PrevY_f =
      windEmulatorStep4_WECSim_cal->RateLimiter_IC_d;

    /* InitializeConditions for StateSpace: '<S542>/Internal' */
    windEmulatorStep4_WECSim_X.Internal_CSTATE_a =
      windEmulatorStep4_WECSim_cal->Internal_InitialCondition_j;

    /* InitializeConditions for DiscreteIntegrator: '<S437>/Discrete-Time Integrator' */
    windEmulatorStep4_WECSim_DW.DiscreteTimeIntegrator_DSTATE =
      windEmulatorStep4_WECSim_cal->DiscreteTimeIntegrator_IC;

    /* InitializeConditions for DiscreteIntegrator: '<S372>/Discrete-Time Integrator' */
    windEmulatorStep4_WECSim_DW.DiscreteTimeIntegrator_DSTATE_n =
      windEmulatorStep4_WECSim_cal->DiscreteTimeIntegrator_IC_e;

    /* InitializeConditions for DiscreteIntegrator: '<S435>/Discrete-Time Integrator' */
    windEmulatorStep4_WECSim_DW.DiscreteTimeIntegrator_DSTATE_l =
      windEmulatorStep4_WECSim_cal->DiscreteTimeIntegrator_IC_l;

    /* InitializeConditions for RateLimiter: '<S2>/acs880RateLim' */
    windEmulatorStep4_WECSim_DW.LastMajorTime_a = (rtInf);

    /* InitializeConditions for Memory: '<S1>/Memory' */
    windEmulatorStep4_WECSim_DW.Memory_PreviousInput_d =
      windEmulatorStep4_WECSim_cal->Memory_InitialCondition_l;

    /* InitializeConditions for Memory: '<S1>/Memory1' */
    windEmulatorStep4_WECSim_DW.Memory1_PreviousInput_p =
      windEmulatorStep4_WECSim_cal->Memory1_InitialCondition_a;

    /* InitializeConditions for Memory: '<S1>/Memory2' */
    windEmulatorStep4_WECSim_DW.Memory2_PreviousInput_h =
      windEmulatorStep4_WECSim_cal->Memory2_InitialCondition_p;

    /* InitializeConditions for Memory: '<S552>/lastRawCounts' */
    windEmulatorStep4_WECSim_DW.lastRawCounts_PreviousInput =
      windEmulatorStep4_WECSim_cal->lastRawCounts_InitialCondition;

    /* InitializeConditions for Memory: '<S552>/lastTurn' */
    windEmulatorStep4_WECSim_DW.lastTurn_PreviousInput =
      windEmulatorStep4_WECSim_cal->lastTurn_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S553>/UD' */
    windEmulatorStep4_WECSim_DW.UD_DSTATE =
      windEmulatorStep4_WECSim_cal->posToVel_ICPrevScaledInput;

    /* InitializeConditions for Memory: '<S29>/Memory' */
    windEmulatorStep4_WECSim_DW.Memory_PreviousInput =
      windEmulatorStep4_WECSim_cal->Memory_InitialCondition;

    /* InitializeConditions for SimscapeInputBlock: '<S332>/INPUT_3_1_1' */
    if (rtmIsMajorTimeStep(windEmulatorStep4_WECSim_M)) {
      windEmulatorStep4_WECSim_DW.INPUT_3_1_1_FirstOutput_4203252 = 0.0;
    }

    /* End of InitializeConditions for SimscapeInputBlock: '<S332>/INPUT_3_1_1' */

    /* InitializeConditions for Delay: '<S58>/Delay One Step' */
    windEmulatorStep4_WECSim_DW.DelayOneStep_DSTATE =
      windEmulatorStep4_WECSim_cal->DelayOneStep_InitialCondition;

    /* InitializeConditions for DiscreteIntegrator: '<S282>/Integrator' */
    windEmulatorStep4_WECSim_DW.Integrator_DSTATE =
      windEmulatorStep4_WECSim_cal->DiscretePIDController_InitialCo;

    /* InitializeConditions for Delay: '<S275>/UD' */
    windEmulatorStep4_WECSim_DW.UD_DSTATE_j =
      windEmulatorStep4_WECSim_cal->DiscretePIDController_Different;

    /* InitializeConditions for SimscapeExecutionBlock: '<S332>/STATE_1' */
    tmp = false;
    tmp_0 = false;
    if (tmp_0 || tmp) {
      is = strcmp(rtsiGetSolverName(&windEmulatorStep4_WECSim_M->solverInfo),
                  "daessc");
      tmp = (is == 0);
      is = strcmp(rtsiGetSolverName(&windEmulatorStep4_WECSim_M->solverInfo),
                  "ode14x");
      tmp = (is == 0) | tmp;
      is = strcmp(rtsiGetSolverName(&windEmulatorStep4_WECSim_M->solverInfo),
                  "ode15s");
      tmp = (is == 0) | tmp;
      is = strcmp(rtsiGetSolverName(&windEmulatorStep4_WECSim_M->solverInfo),
                  "ode1be");
      tmp = (is == 0) | tmp;
      is = strcmp(rtsiGetSolverName(&windEmulatorStep4_WECSim_M->solverInfo),
                  "ode23t");
      tmp = (is == 0) | tmp;
      is = strcmp(rtsiGetSolverName(&windEmulatorStep4_WECSim_M->solverInfo),
                  "odeN");
      tmp = (is == 0) | tmp;
      if (!tmp) {
        rtmSetErrorStatus(windEmulatorStep4_WECSim_M,
                          "Detected inconsistent solvers in the model reference hierarchy. Model built with ode14x requires one of {daessc, ode14x, ode15s, ode1be, ode23t, odeN} solvers to run. Use one of the required solvers in the top model.");
      }
    }

    /* End of InitializeConditions for SimscapeExecutionBlock: '<S332>/STATE_1' */

    /* InitializeConditions for TransportDelay: '<S60>/Transport Delay' */
    for (is = 0; is < 6; is++) {
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK[is] = 0;
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK[is + 6] = 0;
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK[is + 12] = 0;
      (static_cast<real_T *>(windEmulatorStep4_WECSim_DW.TransportDelay_PWORK[is]))
        [0] = 0.0;
      (static_cast<real_T *>(windEmulatorStep4_WECSim_DW.TransportDelay_PWORK[is]))
        [windEmulatorStep4_WECSim_DW.TransportDelay_IWORK[is + 18]] =
        windEmulatorStep4_WECSim_M->Timing.t[0];
    }

    /* End of InitializeConditions for TransportDelay: '<S60>/Transport Delay' */

    /* InitializeConditions for TransportDelay: '<S139>/Transport Delay' */
    for (is = 0; is < 6; is++) {
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c[is] = 0;
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c[is + 6] = 0;
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c[is + 12] = 0;
      (static_cast<real_T *>
        (windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_f[is]))[0] = 0.0;
      (static_cast<real_T *>
        (windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_f[is]))
        [windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c[is + 18]] =
        windEmulatorStep4_WECSim_M->Timing.t[0];
    }

    /* End of InitializeConditions for TransportDelay: '<S139>/Transport Delay' */

    /* InitializeConditions for TransportDelay: '<S64>/Transport Delay' */
    for (is = 0; is < 6; is++) {
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_o[is] = 0;
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_o[is + 6] = 0;
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_o[is + 12] = 0;
      (static_cast<real_T *>
        (windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_j[is]))[0] = 0.0;
      (static_cast<real_T *>
        (windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_j[is]))
        [windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_o[is + 18]] =
        windEmulatorStep4_WECSim_M->Timing.t[0];
    }

    /* End of InitializeConditions for TransportDelay: '<S64>/Transport Delay' */

    /* InitializeConditions for TransportDelay: '<S143>/Transport Delay' */
    for (is = 0; is < 6; is++) {
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c1[is] = 0;
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c1[is + 6] = 0;
      windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c1[is + 12] = 0;
      (static_cast<real_T *>
        (windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_n[is]))[0] = 0.0;
      (static_cast<real_T *>
        (windEmulatorStep4_WECSim_DW.TransportDelay_PWORK_n[is]))
        [windEmulatorStep4_WECSim_DW.TransportDelay_IWORK_c1[is + 18]] =
        windEmulatorStep4_WECSim_M->Timing.t[0];
    }

    /* End of InitializeConditions for TransportDelay: '<S143>/Transport Delay' */

    /* InitializeConditions for Memory: '<S433>/Memory' */
    windEmulatorStep4_WECSim_DW.Memory_PreviousInput_kk =
      windEmulatorStep4_WECSim_cal->SRFlipFlop_initial_condition;

    /* InitializeConditions for RateLimiter: '<S373>/Rate Limiter1' */
    windEmulatorStep4_WECSim_DW.PrevY_m =
      windEmulatorStep4_WECSim_cal->RateLimiter1_IC;

    /* InitializeConditions for Memory: '<S434>/Memory' */
    windEmulatorStep4_WECSim_DW.Memory_PreviousInput_h =
      windEmulatorStep4_WECSim_cal->SRFlipFlop_initial_condition_d;

    /* InitializeConditions for RateLimiter: '<S373>/Rate Limiter' */
    windEmulatorStep4_WECSim_DW.LastMajorTime_d = (rtInf);

    /* InitializeConditions for DiscreteIntegrator: '<S412>/Integrator' */
    windEmulatorStep4_WECSim_DW.Integrator_DSTATE_d =
      windEmulatorStep4_WECSim_cal->PIDController_InitialConditio_a;
    windEmulatorStep4_WECSim_DW.Integrator_PrevResetState = 0;

    /* InitializeConditions for DiscreteIntegrator: '<S407>/Filter' */
    windEmulatorStep4_WECSim_DW.Filter_DSTATE =
      windEmulatorStep4_WECSim_cal->PIDController_InitialConditionF;
    windEmulatorStep4_WECSim_DW.Filter_PrevResetState = 0;

    /* InitializeConditions for Memory: '<S496>/Memory' */
    windEmulatorStep4_WECSim_DW.Memory_PreviousInput_g =
      windEmulatorStep4_WECSim_cal->SRFlipFlop_initial_condition_k;

    /* InitializeConditions for RateLimiter: '<S436>/Rate Limiter1' */
    windEmulatorStep4_WECSim_DW.PrevY_fq =
      windEmulatorStep4_WECSim_cal->RateLimiter1_IC_g;

    /* InitializeConditions for Memory: '<S497>/Memory' */
    windEmulatorStep4_WECSim_DW.Memory_PreviousInput_n =
      windEmulatorStep4_WECSim_cal->SRFlipFlop_initial_condition_j;

    /* InitializeConditions for RateLimiter: '<S436>/Rate Limiter' */
    windEmulatorStep4_WECSim_DW.LastMajorTime_k = (rtInf);

    /* InitializeConditions for DiscreteIntegrator: '<S475>/Integrator' */
    windEmulatorStep4_WECSim_DW.Integrator_DSTATE_e =
      windEmulatorStep4_WECSim_cal->PIDController_InitialConditio_c;
    windEmulatorStep4_WECSim_DW.Integrator_PrevResetState_g = 0;

    /* InitializeConditions for DiscreteIntegrator: '<S470>/Filter' */
    windEmulatorStep4_WECSim_DW.Filter_DSTATE_b =
      windEmulatorStep4_WECSim_cal->PIDController_InitialConditio_k;
    windEmulatorStep4_WECSim_DW.Filter_PrevResetState_g = 0;

    /* SystemInitialize for MATLAB Function: '<S49>/Parse Status Word' */
    windEmulat_ParseStatusWord_Init
      (&windEmulatorStep4_WECSim_DW.sf_ParseStatusWord_h);

    /* SystemInitialize for Chart: '<S18>/ABB Fieldbus Control' incorporates:
     *  Chart: '<S16>/ABB Fieldbus Control'
     *  Chart: '<S4>/FexcRamp'
     */
    windEmulatorStep4_WECSim_DW.temporalCounter_i1_l = 0U;
    windEmulatorStep4_WECSim_DW.sfEvent_l = windEmulatorStep4__CALL_EVENT_k;
    windEmulatorStep4_WECSim_DW.swRDY_ON = 0.0;
    windEmulatorStep4_WECSim_DW.swRDY_RUN = 0.0;
    windEmulatorStep4_WECSim_DW.swRDY_REF = 0.0;
    windEmulatorStep4_WECSim_DW.swTRIPPED = 0.0;
    windEmulatorStep4_WECSim_DW.swOFF_2_STA = 0.0;
    windEmulatorStep4_WECSim_DW.swOFF_3_STA = 0.0;
    windEmulatorStep4_WECSim_DW.swSWC_ON_INHIB = 0.0;
    windEmulatorStep4_WECSim_DW.swAT_SETPOINT = 0.0;
    windEmulatorStep4_WECSim_DW.swEXT_RUN_ENABLE = 0.0;
    windEmulatorStep4_WECSim_DW.cwOFF2_CONTROL = 0.0;
    windEmulatorStep4_WECSim_DW.cwOFF3_CONTROL = 0.0;
    windEmulatorStep4_WECSim_DW.cwENABLE_OPERATION = 0.0;
    windEmulatorStep4_WECSim_DW.cwRAMP_OUT_ZERO = 0.0;
    windEmulatorStep4_WECSim_DW.cwRAMP_HOLD = 0.0;
    windEmulatorStep4_WECSim_DW.cwRAMP_IN_ZERO = 0.0;
    windEmulatorStep4_WECSim_DW.cwRESET = 0.0;
    windEmulatorStep4_WECSim_DW.cwREMOTE_CMD = 0.0;
    windEmulatorStep4_WECSim_DW.cwOFF1_CONTROL = 0.0;
    windEmulatorStep4_WECSim_DW.swEXT_CTRL_LOC = 0.0;
    windEmulatorStep4_WECSim_DW.swREMOTE = 0.0;
    windEmulatorStep4_WECSim_B.ControlWord = 0U;
    windEmulatorStep4_WECSim_DW.swWARNING = 0.0;
    windEmulatorStep4_WECSim_DW.swABOVE_LIMIT = 0.0;
    windEmulatorStep4_WECSim_DW.swMSW_B13 = 0.0;
    windEmulatorStep4_WECSim_DW.swMSW_B14 = 0.0;
    windEmulatorStep4_WECSim_DW.swCOMM_ERR = 0.0;
    windEmulatorStep4_WECSim_B.state_e = abbStateEnum_undefined;
    windEmulatorStep4_WECSim_DW.is_active_c7_windEmulatorStep4_ = 0U;
    windEmulatorStep4_WECSim_DW.is_active_UpdateStateMachine = 0U;
    is_UpdateStateMachine_tmp = windEmulator_IN_NO_ACTIVE_CHILD;
    windEmulatorStep4_WECSim_DW.is_UpdateStateMachine =
      is_UpdateStateMachine_tmp;
    windEmulatorStep4_WECSim_DW.is_active_UpdateControlWord = 0U;

    /* SystemInitialize for Chart: '<S4>/FexcRamp' */
    windEmulatorStep4_WECSim_DW.temporalCounter_i1 = 0U;
    windEmulatorStep4_WECSim_DW.sfEvent = windEmulatorStep4__CALL_EVENT_k;
    windEmulatorStep4_WECSim_B.time_c = 0.0;
    windEmulatorStep4_WECSim_B.ramp_l = 0.0;
    windEmulatorStep4_WECSim_DW.rampLast = 0.0;
    windEmulatorStep4_WECSim_DW.rampUpTime = 0.0;
    windEmulatorStep4_WECSim_DW.runTime = 0.0;
    windEmulatorStep4_WECSim_B.runCounter_b = 0U;
    windEmulatorStep4_WECSim_B.stepCounter_n = 0U;
    windEmulatorStep4_WECSim_B.resetHilIntegrator_h = false;
    windEmulatorStep4_WECSim_B.resetSidIntegrator_h = false;
    windEmulatorStep4_WECSim_DW.is_active_c3_windEmulatorStep4_ = 0U;
    windEmulatorStep4_WECSim_DW.is_c3_windEmulatorStep4_WECSim =
      is_UpdateStateMachine_tmp;

    /* SystemInitialize for MATLAB Function: '<S46>/parseCtrlWord' */
    windEmulator_parseCtrlWord_Init
      (&windEmulatorStep4_WECSim_DW.sf_parseCtrlWord_h);

    /* SystemInitialize for MATLAB Function: '<S44>/Parse Status Word' */
    windEmulat_ParseStatusWord_Init
      (&windEmulatorStep4_WECSim_DW.sf_ParseStatusWord);

    /* SystemInitialize for Chart: '<S16>/ABB Fieldbus Control' */
    windEmulatorStep4_WECSim_DW.temporalCounter_i1_g = 0U;
    windEmulatorStep4_WECSim_DW.sfEvent_d = windEmulatorStep4__CALL_EVENT_k;
    windEmulatorStep4_WECSim_DW.swRDY_ON_f = 0.0;
    windEmulatorStep4_WECSim_DW.swRDY_RUN_f = 0.0;
    windEmulatorStep4_WECSim_DW.swRDY_REF_a = 0.0;
    windEmulatorStep4_WECSim_DW.swTRIPPED_e = 0.0;
    windEmulatorStep4_WECSim_DW.swOFF_2_STA_c = 0.0;
    windEmulatorStep4_WECSim_DW.swOFF_3_STA_i = 0.0;
    windEmulatorStep4_WECSim_DW.swSWC_ON_INHIB_l = 0.0;
    windEmulatorStep4_WECSim_DW.swAT_SETPOINT_b = 0.0;
    windEmulatorStep4_WECSim_DW.swEXT_RUN_ENABLE_h = 0.0;
    windEmulatorStep4_WECSim_DW.cwOFF2_CONTROL_b = 0.0;
    windEmulatorStep4_WECSim_DW.cwOFF3_CONTROL_n = 0.0;
    windEmulatorStep4_WECSim_DW.cwENABLE_OPERATION_l = 0.0;
    windEmulatorStep4_WECSim_DW.cwRAMP_OUT_ZERO_m = 0.0;
    windEmulatorStep4_WECSim_DW.cwRAMP_HOLD_j = 0.0;
    windEmulatorStep4_WECSim_DW.cwRAMP_IN_ZERO_m = 0.0;
    windEmulatorStep4_WECSim_DW.cwRESET_e = 0.0;
    windEmulatorStep4_WECSim_DW.cwREMOTE_CMD_c = 0.0;
    windEmulatorStep4_WECSim_DW.cwOFF1_CONTROL_i = 0.0;
    windEmulatorStep4_WECSim_DW.swEXT_CTRL_LOC_c = 0.0;
    windEmulatorStep4_WECSim_DW.swREMOTE_j = 0.0;
    windEmulatorStep4_WECSim_B.ControlWord_l = 0U;
    windEmulatorStep4_WECSim_DW.swWARNING_k = 0.0;
    windEmulatorStep4_WECSim_DW.swABOVE_LIMIT_a = 0.0;
    windEmulatorStep4_WECSim_DW.swMSW_B13_l = 0.0;
    windEmulatorStep4_WECSim_DW.swMSW_B14_i = 0.0;
    windEmulatorStep4_WECSim_DW.swCOMM_ERR_p = 0.0;
    windEmulatorStep4_WECSim_B.state_ed = abbStateEnum_undefined;
    windEmulatorStep4_WECSim_DW.is_active_c9_windEmulatorStep4_ = 0U;
    windEmulatorStep4_WECSim_DW.is_active_UpdateStateMachine_a = 0U;
    windEmulatorStep4_WECSim_DW.is_UpdateStateMachine_g =
      is_UpdateStateMachine_tmp;
    windEmulatorStep4_WECSim_DW.is_active_UpdateControlWord_p = 0U;

    /* SystemInitialize for MATLAB Function: '<S41>/parseCtrlWord' */
    windEmulator_parseCtrlWord_Init
      (&windEmulatorStep4_WECSim_DW.sf_parseCtrlWord);

    /* SystemInitialize for MATLAB Function: '<S81>/quaternion2EulXYZ' */
    windEmul_quaternion2EulXYZ_Init
      (&windEmulatorStep4_WECSim_DW.sf_quaternion2EulXYZ);

    /* SystemInitialize for MATLAB Function: '<S133>/Yaw Kinematic Transforms' */
    win_YawKinematicTransforms_Init
      (&windEmulatorStep4_WECSim_DW.sf_YawKinematicTransforms);

    /* SystemInitialize for MATLAB Function: '<S126>/MATLAB Function1' */
    windEmulat_MATLABFunction1_Init
      (&windEmulatorStep4_WECSim_DW.sf_MATLABFunction1);

    /* SystemInitialize for MATLAB Function: '<S70>/Yaw Force Transforms' */
    windEmu_YawForceTransforms_Init
      (&windEmulatorStep4_WECSim_DW.sf_YawForceTransforms);

    /* SystemInitialize for MATLAB Function: '<S160>/quaternion2EulXYZ' */
    windEmul_quaternion2EulXYZ_Init
      (&windEmulatorStep4_WECSim_DW.sf_quaternion2EulXYZ_c);

    /* SystemInitialize for MATLAB Function: '<S212>/Yaw Kinematic Transforms' */
    win_YawKinematicTransforms_Init
      (&windEmulatorStep4_WECSim_DW.sf_YawKinematicTransforms_l);

    /* SystemInitialize for MATLAB Function: '<S205>/MATLAB Function1' */
    windEmulat_MATLABFunction1_Init
      (&windEmulatorStep4_WECSim_DW.sf_MATLABFunction1_e);

    /* SystemInitialize for MATLAB Function: '<S149>/Yaw Force Transforms' */
    windEmu_YawForceTransforms_Init
      (&windEmulatorStep4_WECSim_DW.sf_YawForceTransforms_i);
    windEmulator_MovingAverage_Init(&windEmulatorStep4_WECSim_DW.MovingAverage_p);
    windEmulator_MovingAverage_Init(&windEmulatorStep4_WECSim_DW.MovingAverage);
    windEmulat_MovingAverage_e_Init
      (&windEmulatorStep4_WECSim_DW.MovingAverage_pn);
    windEmulat_MovingAverage_e_Init(&windEmulatorStep4_WECSim_DW.MovingAverage1);

    /* Root-level InitSystemMatrices */
    {
      static int_T modelMassMatrixIr[26] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 9,
        10, 12, 13, 14, 19, 33, 17, 15, 16, 18, 22, 21, 20, 23, 24 };

      static int_T modelMassMatrixJc[45] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
        11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 26, 26, 26,
        26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26 };

      static real_T modelMassMatrixPr[26] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0 };

      (void) std::memcpy(windEmulatorStep4_WECSim_MassMatrix.ir,
                         modelMassMatrixIr,
                         26*sizeof(int_T));
      (void) std::memcpy(windEmulatorStep4_WECSim_MassMatrix.jc,
                         modelMassMatrixJc,
                         45*sizeof(int_T));
      (void) std::memcpy(windEmulatorStep4_WECSim_MassMatrix.pr,
                         modelMassMatrixPr,
                         26*sizeof(real_T));
    }
  }
}

/* Model terminate function */
void windEmulatorStep4_WECSim_terminate(void)
{
  NeslSimulationData *simulationData;
  NeuDiagnosticManager *diagnosticManager;
  windEmulator_MovingAverage_Term(&windEmulatorStep4_WECSim_DW.MovingAverage_p);

  /* Terminate for SimscapeExecutionBlock: '<S541>/STATE_1' */
  diagnosticManager = static_cast<NeuDiagnosticManager *>
    (windEmulatorStep4_WECSim_DW.STATE_1_DiagMgr);
  neu_destroy_diagnostic_manager(diagnosticManager);
  simulationData = static_cast<NeslSimulationData *>
    (windEmulatorStep4_WECSim_DW.STATE_1_SimData);
  nesl_destroy_simulation_data(simulationData);
  nesl_erase_simulator("windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Solver Configuration_1");

  /* Terminate for SimscapeExecutionBlock: '<S541>/OUTPUT_1_0' */
  diagnosticManager = static_cast<NeuDiagnosticManager *>
    (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_DiagMgr);
  neu_destroy_diagnostic_manager(diagnosticManager);
  simulationData = static_cast<NeslSimulationData *>
    (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_SimData);
  nesl_destroy_simulation_data(simulationData);
  nesl_erase_simulator("windEmulatorStep4_WECSim/hptoSim/hptoModel/HPTO/Solver Configuration_1");
  windEmulator_MovingAverage_Term(&windEmulatorStep4_WECSim_DW.MovingAverage);
  windEmulat_MovingAverage_g_Term(&windEmulatorStep4_WECSim_DW.MovingAverage_pn);
  windEmulat_MovingAverage_g_Term(&windEmulatorStep4_WECSim_DW.MovingAverage1);

  /* Terminate for S-Function (slrealtimeenablelogging): '<S5>/Enable File Log' */
  /* Level2 S-Function Block: '<S5>/Enable File Log' (slrealtimeenablelogging) */
  {
    SimStruct *rts = windEmulatorStep4_WECSim_M->childSfunctions[0];
    sfcnTerminate(rts);
  }

  /* Terminate for SimscapeExecutionBlock: '<S216>/STATE_1' */
  diagnosticManager = static_cast<NeuDiagnosticManager *>
    (windEmulatorStep4_WECSim_DW.STATE_1_DiagMgr_o);
  neu_destroy_diagnostic_manager(diagnosticManager);
  simulationData = static_cast<NeslSimulationData *>
    (windEmulatorStep4_WECSim_DW.STATE_1_SimData_h);
  nesl_destroy_simulation_data(simulationData);
  nesl_erase_simulator("windEmulatorStep4_WECSim/hptoSim/WECSimModel/Global Reference Frame/Solver Configuration_1");

  /* Terminate for SimscapeExecutionBlock: '<S216>/OUTPUT_1_1' */
  diagnosticManager = static_cast<NeuDiagnosticManager *>
    (windEmulatorStep4_WECSim_DW.OUTPUT_1_1_DiagMgr);
  neu_destroy_diagnostic_manager(diagnosticManager);
  simulationData = static_cast<NeslSimulationData *>
    (windEmulatorStep4_WECSim_DW.OUTPUT_1_1_SimData);
  nesl_destroy_simulation_data(simulationData);
  nesl_erase_simulator("windEmulatorStep4_WECSim/hptoSim/WECSimModel/Global Reference Frame/Solver Configuration_1");

  /* Terminate for SimscapeExecutionBlock: '<S332>/STATE_1' */
  diagnosticManager = static_cast<NeuDiagnosticManager *>
    (windEmulatorStep4_WECSim_DW.STATE_1_DiagMgr_g);
  neu_destroy_diagnostic_manager(diagnosticManager);
  simulationData = static_cast<NeslSimulationData *>
    (windEmulatorStep4_WECSim_DW.STATE_1_SimData_a);
  nesl_destroy_simulation_data(simulationData);
  nesl_erase_simulator("windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Solver Configuration_1");

  /* Terminate for SimscapeExecutionBlock: '<S332>/OUTPUT_1_0' */
  diagnosticManager = static_cast<NeuDiagnosticManager *>
    (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_DiagMgr_i);
  neu_destroy_diagnostic_manager(diagnosticManager);
  simulationData = static_cast<NeslSimulationData *>
    (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_SimData_l);
  nesl_destroy_simulation_data(simulationData);
  nesl_erase_simulator("windEmulatorStep4_WECSim/hptoSim/WECSimModel/PTO_WECSim/Solver Configuration_1");

  /* Terminate for SimscapeExecutionBlock: '<S216>/OUTPUT_1_0' */
  diagnosticManager = static_cast<NeuDiagnosticManager *>
    (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_DiagMgr_a);
  neu_destroy_diagnostic_manager(diagnosticManager);
  simulationData = static_cast<NeslSimulationData *>
    (windEmulatorStep4_WECSim_DW.OUTPUT_1_0_SimData_i);
  nesl_destroy_simulation_data(simulationData);
  nesl_erase_simulator("windEmulatorStep4_WECSim/hptoSim/WECSimModel/Global Reference Frame/Solver Configuration_1");

  /* user code (Terminate function Trailer) */

  /*------------ S-Function Block: <Root>/EtherCAT Init Process Shutdown Network ------------*/
  {
    int_T status;
    status = xpcEtherCATstop(0, 1000 );/* 1 second timeout */
  }
}
