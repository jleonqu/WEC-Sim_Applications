/*
 * windEmulatorStep4.cpp
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

#include "windEmulatorStep4.h"
#include "rtwtypes.h"
#include "windEmulatorStep4_types.h"
#include "windEmulatorStep4_private.h"
#include "rte_windEmulatorStep4_parameters.h"
#include "windEmulatorStep4_cal.h"
#include "abbState.h"
#include <cmath>

extern "C"
{

#include "rt_nonfinite.h"

}

#include <cstring>
#include "sidType.h"
#include "expType.h"
#include <stddef.h>
#include "zero_crossing_types.h"

/* Named constants for Chart: '<S16>/ABB Fieldbus Control' */
const uint32_T windEmula_IN_notReadyToSwitchOn = 3U;
const uint32_T windEmulat_IN_operationDisabled = 4U;
const uint32_T windEmulato_IN_operationEnabled = 5U;
const int32_T windEmulatorStep4_CALL_EVENT = -1;
const uint32_T windEmulatorStep4_IN_DelayOFF1 = 1U;
const uint32_T windEmulatorStep4_IN_initialize = 2U;
const uint8_T windEmulator_IN_NO_ACTIVE_CHILD = 0U;
const uint32_T windEmulator_IN_readyToSwitchOn = 6U;

/* Named constants for Chart: '<S4>/FexcRamp' */
const uint32_T windEmulatorStep4_IN_idle = 1U;
const uint32_T windEmulatorStep4_IN_init = 2U;
const uint32_T windEmulatorStep4_IN_rampdown = 3U;
const uint32_T windEmulatorStep4_IN_rampup = 4U;
const uint32_T windEmulatorStep4_IN_runing = 5U;
const real_T windEmulatorStep4_period = 0.004;

/* Block signals (default storage) */
B_windEmulatorStep4_T windEmulatorStep4_B;

/* Continuous states */
X_windEmulatorStep4_T windEmulatorStep4_X;

/* Block states (default storage) */
DW_windEmulatorStep4_T windEmulatorStep4_DW;

/* Previous zero-crossings (trigger) states */
PrevZCX_windEmulatorStep4_T windEmulatorStep4_PrevZCX;

/* External inputs (root inport signals with default storage) */
ExtU_windEmulatorStep4_T windEmulatorStep4_U;

/* Real-time model */
RT_MODEL_windEmulatorStep4_T windEmulatorStep4_M_ = RT_MODEL_windEmulatorStep4_T
  ();
RT_MODEL_windEmulatorStep4_T *const windEmulatorStep4_M = &windEmulatorStep4_M_;

/* Forward declaration for local functions */
static void windEmulatorSt_SystemCore_setup(dsp_simulink_MovingAverage_wi_T *obj);

/* Forward declaration for local functions */
static void windEmulato_swParseStatusWord_a(void);
static void windEmulat_cwBuildControlWord_l(void);
static void windEmulatorStep_cwInitialize_d(void);
static void windEmulatorS_swParseStatusWord(void);
static void windEmulator_cwBuildControlWord(void);
static void windEmulatorStep4_cwInitialize(void);
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
                        1,
                        1,
                        (unsigned char *)xmlecatArr_0,
                        xmlecatArr_0_count,
                        0,
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

    rtmSetErrorStatus(windEmulatorStep4_M, errMsg);
    return;
  }
}

/*
 * This function updates continuous states using the ODE4 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE4_IntgData *id = static_cast<ODE4_IntgData *>(rtsiGetSolverData(si));
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T *f3 = id->f[3];
  real_T temp;
  int_T i;
  int_T nXc = 6;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) std::memcpy(y, x,
                     static_cast<uint_T>(nXc)*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  windEmulatorStep4_derivatives();

  /* f1 = f(t + (h/2), y + (h/2)*f0) */
  temp = 0.5 * h;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f0[i]);
  }

  rtsiSetT(si, t + temp);
  rtsiSetdX(si, f1);
  windEmulatorStep4_step();
  windEmulatorStep4_derivatives();

  /* f2 = f(t + (h/2), y + (h/2)*f1) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f1[i]);
  }

  rtsiSetdX(si, f2);
  windEmulatorStep4_step();
  windEmulatorStep4_derivatives();

  /* f3 = f(t + h, y + h*f2) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (h*f2[i]);
  }

  rtsiSetT(si, tnew);
  rtsiSetdX(si, f3);
  windEmulatorStep4_step();
  windEmulatorStep4_derivatives();

  /* tnew = t + h
     ynew = y + (h/6)*(f0 + 2*f1 + 2*f2 + 2*f3) */
  temp = h / 6.0;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + temp*(f0[i] + 2.0*f1[i] + 2.0*f2[i] + f3[i]);
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/*
 * Output and update for atomic system:
 *    '<S41>/parseCtrlWord'
 *    '<S46>/parseCtrlWord'
 */
void windEmulatorStep4_parseCtrlWord(uint16_T rtu_ctrlWord,
  B_parseCtrlWord_windEmulatorS_T *localB)
{
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
  dsp_simulink_MovingAverage_wi_T *obj_0;
  g_dsp_internal_SlidingWindowA_T *iobj_0;
  obj->isSetupComplete = false;
  obj->isInitialized = 1;
  obj_0 = obj;
  obj_0->NumChannels = 1;
  obj_0->FrameLength = 1;
  iobj_0 = &obj_0->_pobj0;
  iobj_0->isInitialized = 0;
  iobj_0->isInitialized = 0;
  obj_0->pStatistic = iobj_0;
  obj->isSetupComplete = true;
  obj->TunablePropsChanged = false;
}

/* System initialize for atomic system: */
void windEmulator_MovingAverage_Init(DW_MovingAverage_windEmulator_T *localDW)
{
  dsp_simulink_MovingAverage_wi_T *b_obj;
  g_dsp_internal_SlidingWindowA_T *obj;

  /* Start for MATLABSystem: '<S43>/Moving Average' */
  localDW->obj.matlabCodegenIsDeleted = true;
  b_obj = &localDW->obj;
  b_obj->isInitialized = 0;
  b_obj->NumChannels = -1;
  b_obj->FrameLength = -1;
  b_obj->matlabCodegenIsDeleted = false;
  localDW->objisempty = true;
  windEmulatorSt_SystemCore_setup(&localDW->obj);

  /* InitializeConditions for MATLABSystem: '<S43>/Moving Average' */
  b_obj = &localDW->obj;
  obj = b_obj->pStatistic;
  if (obj->isInitialized == 1) {
    obj->pCumSum = 0.0;
    for (int32_T i = 0; i < 2499; i++) {
      obj->pCumSumRev[i] = 0.0;
    }

    obj->pCumRevIndex = 1.0;
    obj->pModValueRev = 0.0;
  }

  /* End of InitializeConditions for MATLABSystem: '<S43>/Moving Average' */
}

/* Output and update for atomic system: */
void windEmulatorStep4_MovingAverage(real_T rtu_0,
  B_MovingAverage_windEmulatorS_T *localB, DW_MovingAverage_windEmulator_T
  *localDW)
{
  dsp_simulink_MovingAverage_wi_T *obj;
  dsp_simulink_MovingAverage_wi_T *obj_0;
  g_dsp_internal_SlidingWindowA_T *obj_1;
  g_dsp_internal_SlidingWindowA_T *obj_2;
  g_dsp_internal_SlidingWindowA_T *obj_3;
  g_dsp_internal_SlidingWindowA_T *obj_4;
  real_T csum;
  real_T cumRevIndex;
  real_T modValueRev;
  real_T tmp;
  real_T z;

  /* MATLABSystem: '<S43>/Moving Average' */
  obj = &localDW->obj;
  obj_0 = obj;
  if (obj_0->TunablePropsChanged) {
    obj_0->TunablePropsChanged = false;
  }

  obj_1 = obj->pStatistic;
  if (obj_1->isInitialized != 1) {
    obj_2 = obj_1;
    obj_3 = obj_2;
    obj_3->isSetupComplete = false;
    obj_3->isInitialized = 1;
    obj_4 = obj_3;
    obj_4->pCumSum = 0.0;
    for (int32_T i = 0; i < 2499; i++) {
      obj_4->pCumSumRev[i] = 0.0;
    }

    obj_4->pCumRevIndex = 1.0;
    obj_4->pModValueRev = 0.0;
    obj_3->isSetupComplete = true;
    obj_2->pCumSum = 0.0;
    for (int32_T i = 0; i < 2499; i++) {
      obj_2->pCumSumRev[i] = 0.0;
    }

    obj_2->pCumRevIndex = 1.0;
    obj_2->pModValueRev = 0.0;
  }

  cumRevIndex = obj_1->pCumRevIndex;
  csum = obj_1->pCumSum;
  for (int32_T i = 0; i < 2499; i++) {
    localB->csumrev[i] = obj_1->pCumSumRev[i];
  }

  modValueRev = obj_1->pModValueRev;
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

  obj_1->pCumSum = csum;
  for (int32_T i = 0; i < 2499; i++) {
    obj_1->pCumSumRev[i] = localB->csumrev[i];
  }

  obj_1->pCumRevIndex = cumRevIndex;
  obj_1->pModValueRev = modValueRev;

  /* MATLABSystem: '<S43>/Moving Average' */
  localB->MovingAverage = tmp;
}

/* Termination for atomic system: */
void windEmulator_MovingAverage_Term(DW_MovingAverage_windEmulator_T *localDW)
{
  dsp_simulink_MovingAverage_wi_T *obj;
  g_dsp_internal_SlidingWindowA_T *obj_0;

  /* Terminate for MATLABSystem: '<S43>/Moving Average' */
  obj = &localDW->obj;
  if (!obj->matlabCodegenIsDeleted) {
    obj->matlabCodegenIsDeleted = true;
    if ((obj->isInitialized == 1) && obj->isSetupComplete) {
      obj_0 = obj->pStatistic;
      if (obj_0->isInitialized == 1) {
        obj_0->isInitialized = 2;
      }

      obj->NumChannels = -1;
      obj->FrameLength = -1;
    }
  }

  /* End of Terminate for MATLABSystem: '<S43>/Moving Average' */
}

/*
 * Output and update for atomic system:
 *    '<S44>/Parse Status Word'
 *    '<S49>/Parse Status Word'
 */
void windEmulatorSte_ParseStatusWord(uint16_T rtu_StatusWord,
  B_ParseStatusWord_windEmulato_T *localB)
{
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

/* Function for Chart: '<S18>/ABB Fieldbus Control' */
static void windEmulato_swParseStatusWord_a(void)
{
  uint16_T sw;
  sw = windEmulatorStep4_B.BusAssignment_a.statusWord;
  windEmulatorStep4_DW.swRDY_ON = 0.0;
  windEmulatorStep4_DW.swRDY_RUN = 0.0;
  windEmulatorStep4_DW.swRDY_REF = 0.0;
  windEmulatorStep4_DW.swTRIPPED = 0.0;
  windEmulatorStep4_DW.swOFF_2_STA = 0.0;
  windEmulatorStep4_DW.swOFF_3_STA = 0.0;
  windEmulatorStep4_DW.swSWC_ON_INHIB = 0.0;
  windEmulatorStep4_DW.swWARNING = 0.0;
  windEmulatorStep4_DW.swAT_SETPOINT = 0.0;
  windEmulatorStep4_DW.swREMOTE = 0.0;
  windEmulatorStep4_DW.swABOVE_LIMIT = 0.0;
  windEmulatorStep4_DW.swEXT_CTRL_LOC = 0.0;
  windEmulatorStep4_DW.swEXT_RUN_ENABLE = 0.0;
  windEmulatorStep4_DW.swMSW_B13 = 0.0;
  windEmulatorStep4_DW.swMSW_B14 = 0.0;
  windEmulatorStep4_DW.swCOMM_ERR = 0.0;
  if (windEmulatorStep4_B.BusAssignment_a.statusWord >= 32768) {
    int32_T q0;
    uint32_T qY;
    windEmulatorStep4_DW.swCOMM_ERR = 1.0;
    q0 = windEmulatorStep4_B.BusAssignment_a.statusWord;
    qY = static_cast<uint32_T>(q0) - 32768U;
    if (qY > static_cast<uint32_T>(q0)) {
      qY = 0U;
    }

    q0 = static_cast<int32_T>(qY);
    sw = static_cast<uint16_T>(q0);
  }

  if (sw >= 16384) {
    windEmulatorStep4_DW.swMSW_B14 = 1.0;
    sw = static_cast<uint16_T>(sw - 16384);
  }

  if (sw >= 8192) {
    windEmulatorStep4_DW.swMSW_B13 = 1.0;
    sw = static_cast<uint16_T>(sw - 8192);
  }

  if (sw >= 4096) {
    windEmulatorStep4_DW.swEXT_RUN_ENABLE = 1.0;
    sw = static_cast<uint16_T>(sw - 4096);
  }

  if (sw >= 2048) {
    windEmulatorStep4_DW.swEXT_CTRL_LOC = 1.0;
    sw = static_cast<uint16_T>(sw - 2048);
  }

  if (sw >= 1024) {
    windEmulatorStep4_DW.swABOVE_LIMIT = 1.0;
    sw = static_cast<uint16_T>(sw - 1024);
  }

  if (sw >= 512) {
    windEmulatorStep4_DW.swREMOTE = 1.0;
    sw = static_cast<uint16_T>(sw - 512);
  }

  if (sw >= 256) {
    windEmulatorStep4_DW.swAT_SETPOINT = 1.0;
    sw = static_cast<uint16_T>(sw - 256);
  }

  if (sw >= 128) {
    windEmulatorStep4_DW.swWARNING = 1.0;
    sw = static_cast<uint16_T>(sw - 128);
  }

  if (sw >= 64) {
    windEmulatorStep4_DW.swSWC_ON_INHIB = 1.0;
    sw = static_cast<uint16_T>(sw - 64);
  }

  if (sw >= 32) {
    windEmulatorStep4_DW.swOFF_3_STA = 1.0;
    sw = static_cast<uint16_T>(sw - 32);
  }

  if (sw >= 16) {
    windEmulatorStep4_DW.swOFF_2_STA = 1.0;
    sw = static_cast<uint16_T>(sw - 16);
  }

  if (sw >= 8) {
    windEmulatorStep4_DW.swTRIPPED = 1.0;
    sw = static_cast<uint16_T>(sw - 8);
  }

  if (sw >= 4) {
    windEmulatorStep4_DW.swRDY_REF = 1.0;
    sw = static_cast<uint16_T>(sw - 4);
  }

  if (sw >= 2) {
    windEmulatorStep4_DW.swRDY_RUN = 1.0;
    sw = static_cast<uint16_T>(sw - 2);
  }

  if (sw >= 1) {
    windEmulatorStep4_DW.swRDY_ON = 1.0;
  }
}

/* Function for Chart: '<S18>/ABB Fieldbus Control' */
static void windEmulat_cwBuildControlWord_l(void)
{
  windEmulatorStep4_B.ControlWord = 0U;
  if (windEmulatorStep4_DW.cwOFF1_CONTROL != 0.0) {
    windEmulatorStep4_B.ControlWord = 1U;
  }

  if (windEmulatorStep4_DW.cwOFF2_CONTROL != 0.0) {
    windEmulatorStep4_B.ControlWord = static_cast<uint16_T>
      (windEmulatorStep4_B.ControlWord | 2);
  }

  if (windEmulatorStep4_DW.cwOFF3_CONTROL != 0.0) {
    windEmulatorStep4_B.ControlWord = static_cast<uint16_T>
      (windEmulatorStep4_B.ControlWord | 4);
  }

  if (windEmulatorStep4_DW.cwENABLE_OPERATION != 0.0) {
    windEmulatorStep4_B.ControlWord = static_cast<uint16_T>
      (windEmulatorStep4_B.ControlWord | 8);
  }

  if (windEmulatorStep4_DW.cwRAMP_OUT_ZERO != 0.0) {
    windEmulatorStep4_B.ControlWord = static_cast<uint16_T>
      (windEmulatorStep4_B.ControlWord | 16);
  }

  if (windEmulatorStep4_DW.cwRAMP_HOLD != 0.0) {
    windEmulatorStep4_B.ControlWord = static_cast<uint16_T>
      (windEmulatorStep4_B.ControlWord | 32);
  }

  if (windEmulatorStep4_DW.cwRAMP_IN_ZERO != 0.0) {
    windEmulatorStep4_B.ControlWord = static_cast<uint16_T>
      (windEmulatorStep4_B.ControlWord | 64);
  }

  if (windEmulatorStep4_DW.cwRESET != 0.0) {
    windEmulatorStep4_B.ControlWord = static_cast<uint16_T>
      (windEmulatorStep4_B.ControlWord | 128);
  }

  if (windEmulatorStep4_DW.cwREMOTE_CMD != 0.0) {
    windEmulatorStep4_B.ControlWord = static_cast<uint16_T>
      (windEmulatorStep4_B.ControlWord | 1024);
  }

  if (windEmulatorStep4_B.ACS880CtrlMode) {
    windEmulatorStep4_B.ControlWord = static_cast<uint16_T>
      (windEmulatorStep4_B.ControlWord | 2048);
  }
}

/* Function for Chart: '<S18>/ABB Fieldbus Control' */
static void windEmulatorStep_cwInitialize_d(void)
{
  windEmulatorStep4_DW.cwOFF1_CONTROL = 0.0;
  windEmulatorStep4_DW.cwOFF2_CONTROL = 1.0;
  windEmulatorStep4_DW.cwOFF3_CONTROL = 1.0;
  windEmulatorStep4_DW.cwENABLE_OPERATION = 0.0;
  windEmulatorStep4_DW.cwRAMP_OUT_ZERO = 1.0;
  windEmulatorStep4_DW.cwRAMP_HOLD = 1.0;
  windEmulatorStep4_DW.cwRAMP_IN_ZERO = 1.0;
  windEmulatorStep4_DW.cwRESET = 0.0;
  windEmulatorStep4_DW.cwREMOTE_CMD = 1.0;
}

/* Function for Chart: '<S16>/ABB Fieldbus Control' */
static void windEmulatorS_swParseStatusWord(void)
{
  uint16_T sw;
  sw = windEmulatorStep4_B.BusAssignment_h.statusWord;
  windEmulatorStep4_DW.swRDY_ON_k = 0.0;
  windEmulatorStep4_DW.swRDY_RUN_i = 0.0;
  windEmulatorStep4_DW.swRDY_REF_o = 0.0;
  windEmulatorStep4_DW.swTRIPPED_h = 0.0;
  windEmulatorStep4_DW.swOFF_2_STA_e = 0.0;
  windEmulatorStep4_DW.swOFF_3_STA_a = 0.0;
  windEmulatorStep4_DW.swSWC_ON_INHIB_a = 0.0;
  windEmulatorStep4_DW.swWARNING_j = 0.0;
  windEmulatorStep4_DW.swAT_SETPOINT_b = 0.0;
  windEmulatorStep4_DW.swREMOTE_l = 0.0;
  windEmulatorStep4_DW.swABOVE_LIMIT_c = 0.0;
  windEmulatorStep4_DW.swEXT_CTRL_LOC_m = 0.0;
  windEmulatorStep4_DW.swEXT_RUN_ENABLE_g = 0.0;
  windEmulatorStep4_DW.swMSW_B13_o = 0.0;
  windEmulatorStep4_DW.swMSW_B14_d = 0.0;
  windEmulatorStep4_DW.swCOMM_ERR_g = 0.0;
  if (windEmulatorStep4_B.BusAssignment_h.statusWord >= 32768) {
    int32_T q0;
    uint32_T qY;
    windEmulatorStep4_DW.swCOMM_ERR_g = 1.0;
    q0 = windEmulatorStep4_B.BusAssignment_h.statusWord;
    qY = static_cast<uint32_T>(q0) - 32768U;
    if (qY > static_cast<uint32_T>(q0)) {
      qY = 0U;
    }

    q0 = static_cast<int32_T>(qY);
    sw = static_cast<uint16_T>(q0);
  }

  if (sw >= 16384) {
    windEmulatorStep4_DW.swMSW_B14_d = 1.0;
    sw = static_cast<uint16_T>(sw - 16384);
  }

  if (sw >= 8192) {
    windEmulatorStep4_DW.swMSW_B13_o = 1.0;
    sw = static_cast<uint16_T>(sw - 8192);
  }

  if (sw >= 4096) {
    windEmulatorStep4_DW.swEXT_RUN_ENABLE_g = 1.0;
    sw = static_cast<uint16_T>(sw - 4096);
  }

  if (sw >= 2048) {
    windEmulatorStep4_DW.swEXT_CTRL_LOC_m = 1.0;
    sw = static_cast<uint16_T>(sw - 2048);
  }

  if (sw >= 1024) {
    windEmulatorStep4_DW.swABOVE_LIMIT_c = 1.0;
    sw = static_cast<uint16_T>(sw - 1024);
  }

  if (sw >= 512) {
    windEmulatorStep4_DW.swREMOTE_l = 1.0;
    sw = static_cast<uint16_T>(sw - 512);
  }

  if (sw >= 256) {
    windEmulatorStep4_DW.swAT_SETPOINT_b = 1.0;
    sw = static_cast<uint16_T>(sw - 256);
  }

  if (sw >= 128) {
    windEmulatorStep4_DW.swWARNING_j = 1.0;
    sw = static_cast<uint16_T>(sw - 128);
  }

  if (sw >= 64) {
    windEmulatorStep4_DW.swSWC_ON_INHIB_a = 1.0;
    sw = static_cast<uint16_T>(sw - 64);
  }

  if (sw >= 32) {
    windEmulatorStep4_DW.swOFF_3_STA_a = 1.0;
    sw = static_cast<uint16_T>(sw - 32);
  }

  if (sw >= 16) {
    windEmulatorStep4_DW.swOFF_2_STA_e = 1.0;
    sw = static_cast<uint16_T>(sw - 16);
  }

  if (sw >= 8) {
    windEmulatorStep4_DW.swTRIPPED_h = 1.0;
    sw = static_cast<uint16_T>(sw - 8);
  }

  if (sw >= 4) {
    windEmulatorStep4_DW.swRDY_REF_o = 1.0;
    sw = static_cast<uint16_T>(sw - 4);
  }

  if (sw >= 2) {
    windEmulatorStep4_DW.swRDY_RUN_i = 1.0;
    sw = static_cast<uint16_T>(sw - 2);
  }

  if (sw >= 1) {
    windEmulatorStep4_DW.swRDY_ON_k = 1.0;
  }
}

/* Function for Chart: '<S16>/ABB Fieldbus Control' */
static void windEmulator_cwBuildControlWord(void)
{
  windEmulatorStep4_B.ControlWord_c = 0U;
  if (windEmulatorStep4_DW.cwOFF1_CONTROL_a != 0.0) {
    windEmulatorStep4_B.ControlWord_c = 1U;
  }

  if (windEmulatorStep4_DW.cwOFF2_CONTROL_f != 0.0) {
    windEmulatorStep4_B.ControlWord_c = static_cast<uint16_T>
      (windEmulatorStep4_B.ControlWord_c | 2);
  }

  if (windEmulatorStep4_DW.cwOFF3_CONTROL_c != 0.0) {
    windEmulatorStep4_B.ControlWord_c = static_cast<uint16_T>
      (windEmulatorStep4_B.ControlWord_c | 4);
  }

  if (windEmulatorStep4_DW.cwENABLE_OPERATION_k != 0.0) {
    windEmulatorStep4_B.ControlWord_c = static_cast<uint16_T>
      (windEmulatorStep4_B.ControlWord_c | 8);
  }

  if (windEmulatorStep4_DW.cwRAMP_OUT_ZERO_a != 0.0) {
    windEmulatorStep4_B.ControlWord_c = static_cast<uint16_T>
      (windEmulatorStep4_B.ControlWord_c | 16);
  }

  if (windEmulatorStep4_DW.cwRAMP_HOLD_n != 0.0) {
    windEmulatorStep4_B.ControlWord_c = static_cast<uint16_T>
      (windEmulatorStep4_B.ControlWord_c | 32);
  }

  if (windEmulatorStep4_DW.cwRAMP_IN_ZERO_k != 0.0) {
    windEmulatorStep4_B.ControlWord_c = static_cast<uint16_T>
      (windEmulatorStep4_B.ControlWord_c | 64);
  }

  if (windEmulatorStep4_DW.cwRESET_b != 0.0) {
    windEmulatorStep4_B.ControlWord_c = static_cast<uint16_T>
      (windEmulatorStep4_B.ControlWord_c | 128);
  }

  if (windEmulatorStep4_DW.cwREMOTE_CMD_a != 0.0) {
    windEmulatorStep4_B.ControlWord_c = static_cast<uint16_T>
      (windEmulatorStep4_B.ControlWord_c | 1024);
  }

  if (windEmulatorStep4_B.ACS800CtrlMode) {
    windEmulatorStep4_B.ControlWord_c = static_cast<uint16_T>
      (windEmulatorStep4_B.ControlWord_c | 2048);
  }
}

/* Function for Chart: '<S16>/ABB Fieldbus Control' */
static void windEmulatorStep4_cwInitialize(void)
{
  windEmulatorStep4_DW.cwOFF1_CONTROL_a = 0.0;
  windEmulatorStep4_DW.cwOFF2_CONTROL_f = 1.0;
  windEmulatorStep4_DW.cwOFF3_CONTROL_c = 1.0;
  windEmulatorStep4_DW.cwENABLE_OPERATION_k = 0.0;
  windEmulatorStep4_DW.cwRAMP_OUT_ZERO_a = 1.0;
  windEmulatorStep4_DW.cwRAMP_HOLD_n = 1.0;
  windEmulatorStep4_DW.cwRAMP_IN_ZERO_k = 1.0;
  windEmulatorStep4_DW.cwRESET_b = 0.0;
  windEmulatorStep4_DW.cwREMOTE_CMD_a = 1.0;
}

/* Model step function */
void windEmulatorStep4_step(void)
{
  if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
    /* set solver stop time */
    if (!(windEmulatorStep4_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&windEmulatorStep4_M->solverInfo,
                            ((windEmulatorStep4_M->Timing.clockTickH0 + 1) *
        windEmulatorStep4_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&windEmulatorStep4_M->solverInfo,
                            ((windEmulatorStep4_M->Timing.clockTick0 + 1) *
        windEmulatorStep4_M->Timing.stepSize0 +
        windEmulatorStep4_M->Timing.clockTickH0 *
        windEmulatorStep4_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(windEmulatorStep4_M)) {
    windEmulatorStep4_M->Timing.t[0] = rtsiGetT(&windEmulatorStep4_M->solverInfo);
  }

  {
    NeParameterBundle expl_temp;
    NeslRtpManager *rtpManager;
    NeslSimulationData *simulationData;
    NeslSimulator *simulator;
    NeuDiagnosticManager *diag;
    invPowerBus *tmp_d;
    struct_KY1U3Kyrwv5e6VnUIBWG5G *tmp_9;
    struct_KY1U3Kyrwv5e6VnUIBWG5G *tmp_e;
    struct_KY1U3Kyrwv5e6VnUIBWG5G *tmp_f;
    real_T tmp_2[59];
    real_T tmp_0[20];
    real_T tmp[3];
    real_T rateLimiterRate;
    real_T riseValLimit;
    real_T time;
    real_T time_0;
    real_T time_1;
    real_T time_2;
    real_T tmp_4;
    real_T tmp_5;
    real_T tmp_6;
    real_T tmp_7;
    real_T tmp_8;
    real_T tmp_a;
    real_T tmp_b;
    real_T tmp_c;
    real_T tmp_h;
    real_T tmp_i;
    real_T tmp_j;
    real_T tmp_k;
    real_T tmp_m;
    real_T u1;
    real_T *parameterBundle_mRealParameters;
    int32_T isHit;
    int32_T isHit_0;
    int32_T rowIdx;
    int_T tmp_3[7];
    int_T tmp_1[6];
    uint32_T q0;
    boolean_T e_out;
    boolean_T tmp_g;
    boolean_T tmp_l;

    /* Gain: '<S123>/Gain' */
    tmp_m = *get_TorqueLoadMax();

    /* Constant: '<S1>/ACS800CtrlMode' */
    tmp_l = *get_ctrlModeTorque();

    /* Gain: '<S12>/rad//s->rpm' */
    tmp_k = *get_radps2rpm();

    /* Switch: '<S278>/Switch' incorporates:
     *  Constant: '<S13>/Constant2'
     */
    tmp_j = *get_acs880SpeedPILimLo();

    /* Switch: '<S278>/Switch2' incorporates:
     *  Constant: '<S13>/Constant1'
     */
    tmp_i = *get_acs880SpeedPILimUp();

    /* Switch: '<S53>/Switch1' */
    tmp_h = *get_minTorqueRef_Nm();

    /* Chart: '<S4>/FexcRamp' */
    u1 = *get_rampTime();

    /* Switch: '<S61>/Switch' incorporates:
     *  Constant: '<S61>/DeadBandController'
     *  Switch generated from: '<S53>/Switch'
     */
    tmp_g = *get_deadBandController();

    /* Gain: '<S121>/Gain' incorporates:
     *  Switch generated from: '<S53>/Switch'
     */
    tmp_f = get_TorqueInputControl();

    /* Gain: '<S60>/Gain' incorporates:
     *  Switch generated from: '<S53>/Switch'
     */
    tmp_e = get_PressureControl();

    /* BusAssignment: '<S11>/Bus Assignment' incorporates:
     *  Constant: '<S11>/Constant'
     */
    tmp_d = get_invPowerStruct();

    /* RateLimiter: '<S288>/torqueSlewRate' */
    tmp_c = *get_fromFileTorqueSlewRate();

    /* RateLimiter: '<S288>/speedSlewRate' */
    tmp_b = *get_fromFileSpeedSlewRate();

    /* RateLimiter: '<S178>/Rate Limiter' incorporates:
     *  Constant: '<S178>/Constant1'
     */
    tmp_a = *get_rpm2radps();

    /* Product: '<S174>/Product' incorporates:
     *  Constant: '<S174>/Constant1'
     */
    tmp_9 = get_SpeedControl();

    /* Product: '<S175>/Product' incorporates:
     *  Constant: '<S175>/Constant1'
     */
    tmp_8 = *get_belowMinPGain();

    /* RateLimiter: '<S120>/Rate Limiter' */
    tmp_7 = *get_deadbandTorqueSlewRate();

    /* Switch: '<S164>/Switch' incorporates:
     *  Constant: '<S122>/Constant1'
     */
    tmp_6 = *get_genMaxTorque();

    /* Gain: '<S53>/Gain' incorporates:
     *  Switch: '<S53>/Switch2'
     */
    tmp_5 = *get_Dm_max();

    /* Gain: '<S1>/Nm -> %' */
    tmp_4 = *get_acs880RatedTorque();

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
        memcpy(&windEmulatorStep4_B.EtherCATInit[0], data,6*sizeof(int32_T));
        mwErrorClear( (int_T)0 );

        // Clear all momentary triggered values
      }

      if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
        uint32_T qY;

        /* S-Function (slecatpdorx): '<S10>/EtherCAT PDO Receive1' */
        {
          /*------------ S-Function Block: <S10>/EtherCAT PDO Receive1 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.ACS880MotorVoltage;
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
        windEmulatorStep4_B.CastToDouble =
          windEmulatorStep4_B.ACS880MotorVoltage;

        /* Gain: '<S10>/Gain' */
        windEmulatorStep4_B.Gain = *get_acs880MotorVoltsScaling() *
          windEmulatorStep4_B.CastToDouble;

        /* S-Function (slecatpdorx): '<S10>/EtherCAT PDO Receive2' */
        {
          /*------------ S-Function Block: <S10>/EtherCAT PDO Receive2 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.ACS880MotorCurrent;
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
        windEmulatorStep4_B.CastToDouble1 =
          windEmulatorStep4_B.ACS880MotorCurrent;

        /* Gain: '<S10>/Gain2' */
        windEmulatorStep4_B.Gain2 = *get_acs880MotorCurrentScaling() *
          windEmulatorStep4_B.CastToDouble1;

        /* S-Function (slecatpdorx): '<S10>/EtherCAT PDO Receive3' */
        {
          /*------------ S-Function Block: <S10>/EtherCAT PDO Receive3 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.ACS880OutputFreq;
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
        windEmulatorStep4_B.CastToDouble2 = windEmulatorStep4_B.ACS880OutputFreq;

        /* Gain: '<S10>/Gain3' */
        windEmulatorStep4_B.Gain3 = *get_acs880FreqScaling() *
          windEmulatorStep4_B.CastToDouble2;

        /* S-Function (slecatpdorx): '<S10>/EtherCAT PDO Receive4' */
        {
          /*------------ S-Function Block: <S10>/EtherCAT PDO Receive4 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.ACS880MotorSpeed;
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
        windEmulatorStep4_B.CastToDouble3 = windEmulatorStep4_B.ACS880MotorSpeed;

        /* Gain: '<S10>/Gain4' */
        windEmulatorStep4_B.Gain4 = *get_acs880SpeedScaling() *
          windEmulatorStep4_B.CastToDouble3;

        /* S-Function (slecatpdorx): '<S10>/EtherCAT PDO Receive5' */
        {
          /*------------ S-Function Block: <S10>/EtherCAT PDO Receive5 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.ACS880MotorTorque;
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
        windEmulatorStep4_B.DataTypeConversion1 =
          windEmulatorStep4_B.ACS880MotorTorque;

        /* Gain: '<S10>/Gain1' */
        windEmulatorStep4_B.Gain1 = *get_acs880TorqueScaling() *
          windEmulatorStep4_B.DataTypeConversion1;

        /* S-Function (slecatpdorx): '<S10>/EtherCAT PDO Receive6' */
        {
          /*------------ S-Function Block: <S10>/EtherCAT PDO Receive6 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.ACS880MotorShaftPower;
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
        windEmulatorStep4_B.DataTypeConversion2 =
          windEmulatorStep4_B.ACS880MotorShaftPower;

        /* Gain: '<S10>/Gain5' */
        windEmulatorStep4_B.Gain5 = *get_acs880PowerScaling() *
          windEmulatorStep4_B.DataTypeConversion2;

        /* S-Function (slecatpdorx): '<S10>/EtherCAT PDO Receive' */
        {
          /*------------ S-Function Block: <S10>/EtherCAT PDO Receive PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.EtherCATPDOReceive;
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

        /* BusAssignment: '<S10>/Bus Assignment' incorporates:
         *  Constant: '<S10>/Constant'
         */
        windEmulatorStep4_B.BusAssignment_a = *get_acs880SignalStruct();

        /* BusAssignment: '<S10>/Bus Assignment' */
        windEmulatorStep4_B.BusAssignment_a.motorVoltage_V =
          windEmulatorStep4_B.Gain;
        windEmulatorStep4_B.BusAssignment_a.motorCurrent_A =
          windEmulatorStep4_B.Gain2;
        windEmulatorStep4_B.BusAssignment_a.frequency_Hz =
          windEmulatorStep4_B.Gain3;
        windEmulatorStep4_B.BusAssignment_a.motorSpeed_rpm =
          windEmulatorStep4_B.Gain4;
        windEmulatorStep4_B.BusAssignment_a.motorTorque_Nm =
          windEmulatorStep4_B.Gain1;
        windEmulatorStep4_B.BusAssignment_a.shaftPower_W =
          windEmulatorStep4_B.Gain5;
        windEmulatorStep4_B.BusAssignment_a.statusWord =
          windEmulatorStep4_B.EtherCATPDOReceive;

        /* ToAsyncQueueBlock generated from: '<S27>/acs880Signals' */
        slrtLogSignal
          (windEmulatorStep4_DW.TAQSigLogging_InsertedFor_acs88.SLRTSigHandles,
           (((windEmulatorStep4_M->Timing.clockTick1+
              windEmulatorStep4_M->Timing.clockTickH1* 4294967296.0)) * 0.004));

        /* Gain: '<S48>/rpm -> rad//s' */
        windEmulatorStep4_B.rpmrads = tmp_a *
          windEmulatorStep4_B.BusAssignment_a.motorSpeed_rpm;

        /* Product: '<S48>/Product' */
        windEmulatorStep4_B.Product = windEmulatorStep4_B.rpmrads *
          windEmulatorStep4_B.BusAssignment_a.motorTorque_Nm;
        windEmulatorStep4_MovingAverage(windEmulatorStep4_B.Product,
          &windEmulatorStep4_B.MovingAverage_p,
          &windEmulatorStep4_DW.MovingAverage_p);

        /* Gain: '<S48>/shaftPowerAverage_W' */
        windEmulatorStep4_B.shaftPowerAverage_W =
          windEmulatorStep4_cal->shaftPowerAverage_W_Gain *
          windEmulatorStep4_B.MovingAverage_p.MovingAverage;

        /* Gain: '<S48>/shaftPower_W' */
        windEmulatorStep4_B.shaftPower_W =
          windEmulatorStep4_cal->shaftPower_W_Gain * windEmulatorStep4_B.Product;

        /* MATLAB Function: '<S49>/Parse Status Word' */
        windEmulatorSte_ParseStatusWord
          (windEmulatorStep4_B.BusAssignment_a.statusWord,
           &windEmulatorStep4_B.sf_ParseStatusWord_h);

        /* Logic: '<S49>/aboveLimit' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_B.aboveLimit =
          (windEmulatorStep4_cal->Constant_Value_o &&
           (windEmulatorStep4_B.sf_ParseStatusWord_h.above_limit != 0));

        /* Logic: '<S49>/alarm' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_B.alarm = (windEmulatorStep4_cal->Constant_Value_o &&
          (windEmulatorStep4_B.sf_ParseStatusWord_h.alarm != 0));

        /* Logic: '<S49>/atSetpoint' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_B.atSetpoint =
          (windEmulatorStep4_cal->Constant_Value_o &&
           (windEmulatorStep4_B.sf_ParseStatusWord_h.at_setpoint != 0));

        /* Logic: '<S49>/commErr' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_B.commErr = (windEmulatorStep4_cal->Constant_Value_o &&
          (windEmulatorStep4_B.sf_ParseStatusWord_h.comm_err != 0));

        /* Logic: '<S49>/extCtrlLoc' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_B.extCtrlLoc =
          (windEmulatorStep4_cal->Constant_Value_o &&
           (windEmulatorStep4_B.sf_ParseStatusWord_h.ext_ctrl_loc != 0));

        /* Logic: '<S49>/extRunEnable' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_B.extRunEnable =
          (windEmulatorStep4_cal->Constant_Value_o &&
           (windEmulatorStep4_B.sf_ParseStatusWord_h.ext_run_enable != 0));

        /* Logic: '<S49>/mswB13' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_B.mswB13 = (windEmulatorStep4_cal->Constant_Value_o &&
          (windEmulatorStep4_B.sf_ParseStatusWord_h.msw_b13 != 0));

        /* Logic: '<S49>/mswB14' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_B.mswB14 = (windEmulatorStep4_cal->Constant_Value_o &&
          (windEmulatorStep4_B.sf_ParseStatusWord_h.msw_b14 != 0));

        /* Logic: '<S49>/off2' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_B.off2 = (windEmulatorStep4_cal->Constant_Value_o &&
          (windEmulatorStep4_B.sf_ParseStatusWord_h.off2 != 0));

        /* Logic: '<S49>/off3' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_B.off3 = (windEmulatorStep4_cal->Constant_Value_o &&
          (windEmulatorStep4_B.sf_ParseStatusWord_h.off3 != 0));

        /* Logic: '<S49>/rdyOn' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_B.rdyOn = (windEmulatorStep4_cal->Constant_Value_o &&
          (windEmulatorStep4_B.sf_ParseStatusWord_h.rdy_on != 0));

        /* Logic: '<S49>/rdyRef' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_B.rdyRef = (windEmulatorStep4_cal->Constant_Value_o &&
          (windEmulatorStep4_B.sf_ParseStatusWord_h.rdy_ref != 0));

        /* Logic: '<S49>/rdyRun' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_B.rdyRun = (windEmulatorStep4_cal->Constant_Value_o &&
          (windEmulatorStep4_B.sf_ParseStatusWord_h.rdy_run != 0));

        /* Logic: '<S49>/remote' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_B.remote = (windEmulatorStep4_cal->Constant_Value_o &&
          (windEmulatorStep4_B.sf_ParseStatusWord_h.remote != 0));

        /* Logic: '<S49>/switchOnInhibit' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_B.switchOnInhibit =
          (windEmulatorStep4_cal->Constant_Value_o &&
           (windEmulatorStep4_B.sf_ParseStatusWord_h.swc_on_inhib != 0));

        /* Logic: '<S49>/tripped' incorporates:
         *  Constant: '<S49>/Constant'
         */
        windEmulatorStep4_B.tripped = (windEmulatorStep4_cal->Constant_Value_o &&
          (windEmulatorStep4_B.sf_ParseStatusWord_h.tripped != 0));

        /* Gain: '<S28>/frequency_Hz' */
        windEmulatorStep4_B.frequency_Hz =
          windEmulatorStep4_cal->frequency_Hz_Gain *
          windEmulatorStep4_B.BusAssignment_a.frequency_Hz;

        /* Gain: '<S28>/motorCurrent_A' */
        windEmulatorStep4_B.motorCurrent_A =
          windEmulatorStep4_cal->motorCurrent_A_Gain *
          windEmulatorStep4_B.BusAssignment_a.motorCurrent_A;

        /* Gain: '<S28>/motorSpeed_rpm' */
        windEmulatorStep4_B.motorSpeed_rpm =
          windEmulatorStep4_cal->motorSpeed_rpm_Gain *
          windEmulatorStep4_B.BusAssignment_a.motorSpeed_rpm;

        /* Gain: '<S28>/motorTorque_Nm' */
        windEmulatorStep4_B.motorTorque_Nm =
          windEmulatorStep4_cal->motorTorque_Nm_Gain *
          windEmulatorStep4_B.BusAssignment_a.motorTorque_Nm;

        /* Gain: '<S28>/motorVoltage_V' */
        windEmulatorStep4_B.motorVoltage_V =
          windEmulatorStep4_cal->motorVoltage_V_Gain *
          windEmulatorStep4_B.BusAssignment_a.motorVoltage_V;

        /* Gain: '<S28>/shaftPower_W' */
        windEmulatorStep4_B.shaftPower_W_a =
          windEmulatorStep4_cal->shaftPower_W_Gain_e *
          windEmulatorStep4_B.BusAssignment_a.shaftPower_W;

        /* Bias: '<S28>/statusWord' */
        windEmulatorStep4_B.statusWord = static_cast<uint16_T>
          (static_cast<uint32_T>(windEmulatorStep4_B.BusAssignment_a.statusWord)
           + windEmulatorStep4_cal->statusWord_Bias);

        /* Memory: '<S2>/Memory' */
        windEmulatorStep4_B.Memory_b =
          windEmulatorStep4_DW.Memory_PreviousInput_i;

        /* RelationalOperator: '<S2>/NotEqual' incorporates:
         *  Constant: '<S2>/powerUpButton'
         */
        windEmulatorStep4_B.NotEqual =
          (windEmulatorStep4_cal->powerUpButton_Value !=
           windEmulatorStep4_B.Memory_b);

        /* Memory: '<S2>/Memory1' */
        windEmulatorStep4_B.Memory1 = windEmulatorStep4_DW.Memory1_PreviousInput;

        /* RelationalOperator: '<S2>/NotEqual1' incorporates:
         *  Constant: '<S2>/powerDownButton'
         */
        windEmulatorStep4_B.NotEqual1 =
          (windEmulatorStep4_cal->powerDownButton_Value !=
           windEmulatorStep4_B.Memory1);

        /* Memory: '<S2>/Memory2' */
        windEmulatorStep4_B.Memory2 = windEmulatorStep4_DW.Memory2_PreviousInput;

        /* RelationalOperator: '<S2>/NotEqual2' incorporates:
         *  Constant: '<S2>/resetFaultButton'
         */
        windEmulatorStep4_B.NotEqual2 =
          (windEmulatorStep4_cal->resetFaultButton_Value !=
           windEmulatorStep4_B.Memory2);

        /* DataTypeConversion: '<S2>/Cast To Double' */
        windEmulatorStep4_B.CastToDouble_o = windEmulatorStep4_B.NotEqual2;

        /* Memory: '<S4>/Memory' */
        windEmulatorStep4_B.Memory_i =
          windEmulatorStep4_DW.Memory_PreviousInput_k;

        /* RelationalOperator: '<S4>/NotEqual' incorporates:
         *  Constant: '<S4>/eStopButton'
         */
        windEmulatorStep4_B.NotEqual_i =
          (windEmulatorStep4_cal->eStopButton_Value !=
           windEmulatorStep4_B.Memory_i);

        /* Constant: '<S2>/ACS880CtrlMode' */
        windEmulatorStep4_B.ACS880CtrlMode = tmp_l;

        /* Chart: '<S18>/ABB Fieldbus Control' */
        if (windEmulatorStep4_DW.temporalCounter_i1_l < 31U) {
          windEmulatorStep4_DW.temporalCounter_i1_l = static_cast<uint8_T>
            (windEmulatorStep4_DW.temporalCounter_i1_l + 1U);
        }

        windEmulatorStep4_DW.sfEvent_g = windEmulatorStep4_CALL_EVENT;
        if (windEmulatorStep4_DW.is_active_c7_windEmulatorStep4 == 0U) {
          windEmulatorStep4_DW.is_active_c7_windEmulatorStep4 = 1U;
          windEmulatorStep4_DW.is_active_UpdateStateMachine = 1U;
          windEmulatorStep4_DW.is_UpdateStateMachine =
            windEmulatorStep4_IN_initialize;
          windEmulatorStep_cwInitialize_d();
          windEmulatorStep4_B.state_j = abbStateEnum_init;
          windEmulatorStep4_DW.is_active_UpdateControlWord = 1U;
        } else {
          windEmulato_swParseStatusWord_a();
          windEmulatorStep4_DW.cwRESET = windEmulatorStep4_B.CastToDouble_o;
          switch (windEmulatorStep4_DW.is_UpdateStateMachine) {
           case windEmulatorStep4_IN_DelayOFF1:
            if (windEmulatorStep4_DW.temporalCounter_i1_l >= 25U) {
              windEmulatorStep4_DW.is_UpdateStateMachine =
                windEmula_IN_notReadyToSwitchOn;
              windEmulatorStep4_DW.cwOFF1_CONTROL = 0.0;
            } else {
              windEmulatorStep4_B.state_j = abbStateEnum_delayOff1;
            }
            break;

           case windEmulatorStep4_IN_initialize:
            windEmulatorStep4_DW.is_UpdateStateMachine =
              windEmula_IN_notReadyToSwitchOn;
            windEmulatorStep4_DW.cwOFF1_CONTROL = 0.0;
            break;

           case windEmula_IN_notReadyToSwitchOn:
            e_out = ((windEmulatorStep4_DW.swRDY_ON == 1.0) &&
                     (windEmulatorStep4_DW.swWARNING == 0.0));
            if (e_out) {
              windEmulatorStep4_DW.is_UpdateStateMachine =
                windEmulator_IN_readyToSwitchOn;
              windEmulatorStep4_DW.cwOFF1_CONTROL = 0.0;
            } else {
              windEmulatorStep4_B.state_j = abbStateEnum_notReadyToSwitchOn;
            }
            break;

           case windEmulat_IN_operationDisabled:
            e_out = (windEmulatorStep4_B.NotEqual &&
                     (windEmulatorStep4_DW.swRDY_RUN == 1.0));
            if (e_out) {
              windEmulatorStep4_DW.is_UpdateStateMachine =
                windEmulato_IN_operationEnabled;
              windEmulatorStep4_DW.cwOFF1_CONTROL = 1.0;
              windEmulatorStep4_DW.cwENABLE_OPERATION = 1.0;
            } else {
              e_out = (windEmulatorStep4_B.NotEqual1 ||
                       windEmulatorStep4_B.NotEqual_i ||
                       (windEmulatorStep4_DW.swTRIPPED == 1.0) ||
                       (windEmulatorStep4_DW.swSWC_ON_INHIB == 1.0) ||
                       (windEmulatorStep4_DW.swWARNING == 1.0));
              if (e_out) {
                windEmulatorStep4_DW.is_UpdateStateMachine =
                  windEmulatorStep4_IN_DelayOFF1;
                windEmulatorStep4_DW.temporalCounter_i1_l = 0U;
              } else {
                windEmulatorStep4_B.state_j = abbStateEnum_operationDisabled;
              }
            }
            break;

           case windEmulato_IN_operationEnabled:
            if (windEmulatorStep4_B.NotEqual1) {
              windEmulatorStep4_DW.cwENABLE_OPERATION = 0.0;
              windEmulatorStep4_DW.is_UpdateStateMachine =
                windEmulat_IN_operationDisabled;
              windEmulatorStep4_DW.cwOFF1_CONTROL = 1.0;
            } else {
              e_out = (windEmulatorStep4_B.NotEqual_i ||
                       (windEmulatorStep4_DW.swTRIPPED == 1.0) ||
                       (windEmulatorStep4_DW.swSWC_ON_INHIB == 1.0) ||
                       (windEmulatorStep4_DW.swWARNING == 1.0));
              if (e_out) {
                windEmulatorStep4_DW.cwENABLE_OPERATION = 0.0;
                windEmulatorStep4_DW.is_UpdateStateMachine =
                  windEmulatorStep4_IN_DelayOFF1;
                windEmulatorStep4_DW.temporalCounter_i1_l = 0U;
              } else {
                windEmulatorStep4_B.state_j = abbStateEnum_operationEnabled;
              }
            }
            break;

           default:
            /* case IN_readyToSwitchOn: */
            e_out = (windEmulatorStep4_B.NotEqual &&
                     (windEmulatorStep4_DW.swREMOTE == 1.0));
            if (e_out) {
              windEmulatorStep4_DW.is_UpdateStateMachine =
                windEmulat_IN_operationDisabled;
              windEmulatorStep4_DW.cwOFF1_CONTROL = 1.0;
            } else {
              windEmulatorStep4_B.state_j = abbStateEnum_readyToSwitchOn;
            }
            break;
          }

          windEmulat_cwBuildControlWord_l();
        }

        /* End of Chart: '<S18>/ABB Fieldbus Control' */

        /* DataTypeConversion: '<S2>/Cast To Double1' */
        windEmulatorStep4_B.CastToDouble1_a = windEmulatorStep4_B.state_j;

        /* Constant: '<S4>/expType' */
        windEmulatorStep4_B.expType_a = windEmulatorStep4_cal->expType_Value;

        /* DataTypeConversion: '<S4>/ToUint16' incorporates:
         *  Constant: '<S4>/expType'
         */
        windEmulatorStep4_B.ToUint16 = windEmulatorStep4_B.expType_a;

        /* RelationalOperator: '<S4>/Equal' incorporates:
         *  Constant: '<S4>/expModeHil'
         *  Constant: '<S4>/expType'
         */
        windEmulatorStep4_B.Equal = (windEmulatorStep4_B.expType_a ==
          windEmulatorStep4_cal->expModeHil_Value);

        /* RelationalOperator: '<S4>/Equal1' incorporates:
         *  Constant: '<S4>/expModeSid'
         *  Constant: '<S4>/expType'
         */
        windEmulatorStep4_B.Equal1 = (windEmulatorStep4_B.expType_a ==
          windEmulatorStep4_cal->expModeSid_Value);

        /* Memory: '<S4>/Memory1' */
        windEmulatorStep4_B.Memory1_n =
          windEmulatorStep4_DW.Memory1_PreviousInput_d;

        /* RelationalOperator: '<S4>/NotEqual1' incorporates:
         *  Constant: '<S4>/startButton'
         */
        windEmulatorStep4_B.NotEqual1_a =
          (windEmulatorStep4_cal->startButton_Value !=
           windEmulatorStep4_B.Memory1_n);

        /* Memory: '<S4>/Memory2' */
        windEmulatorStep4_B.Memory2_c =
          windEmulatorStep4_DW.Memory2_PreviousInput_l;

        /* RelationalOperator: '<S4>/NotEqual2' incorporates:
         *  Constant: '<S4>/stopButton'
         */
        windEmulatorStep4_B.NotEqual2_c =
          (windEmulatorStep4_cal->stopButton_Value !=
           windEmulatorStep4_B.Memory2_c);

        /* Constant: '<S4>/expRunTime' */
        windEmulatorStep4_B.expRunTime = windEmulatorStep4_cal->expRunTime_Value;

        /* Chart: '<S4>/FexcRamp' incorporates:
         *  Constant: '<S4>/expType'
         */
        if (windEmulatorStep4_DW.temporalCounter_i1 < MAX_uint32_T) {
          windEmulatorStep4_DW.temporalCounter_i1++;
        }

        windEmulatorStep4_DW.sfEvent = windEmulatorStep4_CALL_EVENT;
        if (windEmulatorStep4_DW.is_active_c3_windEmulatorStep4 == 0U) {
          windEmulatorStep4_DW.is_active_c3_windEmulatorStep4 = 1U;
          windEmulatorStep4_DW.is_c3_windEmulatorStep4 =
            windEmulatorStep4_IN_init;
          windEmulatorStep4_B.runCounter_m = 0U;
          windEmulatorStep4_B.stepCounter_a = 1U;
          windEmulatorStep4_B.ramp_i = 0.0;
          windEmulatorStep4_DW.rampLast = 0.0;
          windEmulatorStep4_B.time_d = 0.0;
          windEmulatorStep4_DW.rampUpTime = 0.0;
          windEmulatorStep4_DW.runTime = 0.0;
          windEmulatorStep4_B.resetHilIntegrator_i = true;
          windEmulatorStep4_B.resetSidIntegrator_i = true;
        } else {
          switch (windEmulatorStep4_DW.is_c3_windEmulatorStep4) {
           case windEmulatorStep4_IN_idle:
            if (windEmulatorStep4_B.NotEqual1_a) {
              windEmulatorStep4_DW.is_c3_windEmulatorStep4 =
                windEmulatorStep4_IN_rampup;
              windEmulatorStep4_DW.temporalCounter_i1 = 0U;
              q0 = windEmulatorStep4_B.runCounter_m;
              qY = q0 + 1U;
              if (qY < q0) {
                qY = MAX_uint32_T;
              }

              windEmulatorStep4_B.runCounter_m = qY;
              windEmulatorStep4_B.resetHilIntegrator_i =
                (windEmulatorStep4_B.expType_a != expTypeEnum_hil);
              windEmulatorStep4_B.resetSidIntegrator_i =
                (windEmulatorStep4_B.expType_a != expTypeEnum_sid);
            } else {
              windEmulatorStep4_B.ramp_i = 0.0;
              windEmulatorStep4_DW.rampLast = 0.0;
              windEmulatorStep4_B.stepCounter_a = 1U;
              windEmulatorStep4_B.time_d = 0.0;
              windEmulatorStep4_DW.rampUpTime = 0.0;
              windEmulatorStep4_DW.runTime = 0.0;
              windEmulatorStep4_B.resetHilIntegrator_i = true;
              windEmulatorStep4_B.resetSidIntegrator_i = true;
            }
            break;

           case windEmulatorStep4_IN_init:
            windEmulatorStep4_DW.is_c3_windEmulatorStep4 =
              windEmulatorStep4_IN_idle;
            break;

           case windEmulatorStep4_IN_rampdown:
            if (windEmulatorStep4_B.ramp_i <= 0.0) {
              windEmulatorStep4_DW.is_c3_windEmulatorStep4 =
                windEmulatorStep4_IN_idle;
            } else {
              u1 = static_cast<real_T>(windEmulatorStep4_DW.temporalCounter_i1) *
                0.004 / u1;
              if ((u1 <= 0.0) || rtIsNaN(u1)) {
                u1 = 0.0;
              }

              windEmulatorStep4_B.ramp_i = windEmulatorStep4_DW.rampLast - u1;
              windEmulatorStep4_B.time_d = (windEmulatorStep4_DW.rampUpTime +
                windEmulatorStep4_DW.runTime) + static_cast<real_T>
                (windEmulatorStep4_DW.temporalCounter_i1) * 0.004;
              q0 = windEmulatorStep4_B.stepCounter_a;
              qY = q0 + 1U;
              if (qY < q0) {
                qY = MAX_uint32_T;
              }

              windEmulatorStep4_B.stepCounter_a = qY;
              windEmulatorStep4_B.resetHilIntegrator_i =
                (windEmulatorStep4_B.expType_a != expTypeEnum_hil);
              windEmulatorStep4_B.resetSidIntegrator_i =
                (windEmulatorStep4_B.expType_a != expTypeEnum_sid);
            }
            break;

           case windEmulatorStep4_IN_rampup:
            if (windEmulatorStep4_B.NotEqual2_c) {
              windEmulatorStep4_DW.is_c3_windEmulatorStep4 =
                windEmulatorStep4_IN_rampdown;
              windEmulatorStep4_DW.temporalCounter_i1 = 0U;
            } else if (windEmulatorStep4_B.ramp_i >= 1.0) {
              windEmulatorStep4_DW.is_c3_windEmulatorStep4 =
                windEmulatorStep4_IN_runing;
              windEmulatorStep4_DW.temporalCounter_i1 = 0U;
            } else {
              u1 = static_cast<real_T>(windEmulatorStep4_DW.temporalCounter_i1) *
                0.004 / u1;
              if ((u1 >= 1.0) || rtIsNaN(u1)) {
                windEmulatorStep4_B.ramp_i = 1.0;
              } else {
                windEmulatorStep4_B.ramp_i = u1;
              }

              q0 = windEmulatorStep4_B.stepCounter_a;
              qY = q0 + 1U;
              if (qY < q0) {
                qY = MAX_uint32_T;
              }

              windEmulatorStep4_B.stepCounter_a = qY;
              windEmulatorStep4_DW.rampLast = windEmulatorStep4_B.ramp_i;
              windEmulatorStep4_DW.rampUpTime = static_cast<real_T>
                (windEmulatorStep4_DW.temporalCounter_i1) * 0.004;
              windEmulatorStep4_B.time_d = windEmulatorStep4_DW.rampUpTime;
            }
            break;

           default:
            /* case IN_runing: */
            e_out = (windEmulatorStep4_B.NotEqual2_c ||
                     (windEmulatorStep4_B.time_d >=
                      windEmulatorStep4_B.expRunTime));
            if (e_out) {
              windEmulatorStep4_DW.is_c3_windEmulatorStep4 =
                windEmulatorStep4_IN_rampdown;
              windEmulatorStep4_DW.temporalCounter_i1 = 0U;
            } else {
              windEmulatorStep4_B.ramp_i = 1.0;
              windEmulatorStep4_DW.rampLast = windEmulatorStep4_B.ramp_i;
              windEmulatorStep4_DW.runTime = static_cast<real_T>
                (windEmulatorStep4_DW.temporalCounter_i1) * 0.004;
              windEmulatorStep4_B.time_d = windEmulatorStep4_DW.rampUpTime +
                windEmulatorStep4_DW.runTime;
              q0 = windEmulatorStep4_B.stepCounter_a;
              qY = q0 + 1U;
              if (qY < q0) {
                qY = MAX_uint32_T;
              }

              windEmulatorStep4_B.stepCounter_a = qY;
              windEmulatorStep4_B.resetHilIntegrator_i =
                (windEmulatorStep4_B.expType_a != expTypeEnum_hil);
              windEmulatorStep4_B.resetSidIntegrator_i =
                (windEmulatorStep4_B.expType_a != expTypeEnum_sid);
            }
            break;
          }
        }

        /* BusAssignment: '<S4>/Bus Assignment' incorporates:
         *  Constant: '<S4>/Constant'
         */
        windEmulatorStep4_B.BusAssignment_b = *get_expCtrlStruct();

        /* BusAssignment: '<S4>/Bus Assignment' */
        windEmulatorStep4_B.BusAssignment_b.expType =
          windEmulatorStep4_B.ToUint16;
        windEmulatorStep4_B.BusAssignment_b.runHil = windEmulatorStep4_B.Equal;
        windEmulatorStep4_B.BusAssignment_b.runSid = windEmulatorStep4_B.Equal1;
        windEmulatorStep4_B.BusAssignment_b.time = windEmulatorStep4_B.time_d;
        windEmulatorStep4_B.BusAssignment_b.ramp = windEmulatorStep4_B.ramp_i;
        windEmulatorStep4_B.BusAssignment_b.runCounter =
          windEmulatorStep4_B.runCounter_m;
        windEmulatorStep4_B.BusAssignment_b.stepCounter =
          windEmulatorStep4_B.stepCounter_a;
        windEmulatorStep4_B.BusAssignment_b.resetHilIntegrator =
          windEmulatorStep4_B.resetHilIntegrator_i;
        windEmulatorStep4_B.BusAssignment_b.resetSidIntegrator =
          windEmulatorStep4_B.resetSidIntegrator_i;

        /* DataTypeConversion: '<S3>/toExpTypeEnum' */
        windEmulatorStep4_B.toExpTypeEnum =
          windEmulatorStep4_B.BusAssignment_b.expType;

        /* Switch: '<S238>/Switch' */
        if (windEmulatorStep4_B.BusAssignment_b.runSid) {
          /* Switch: '<S238>/Switch' */
          windEmulatorStep4_B.rampValue =
            windEmulatorStep4_B.BusAssignment_b.ramp;
        } else {
          /* Switch: '<S238>/Switch' incorporates:
           *  Constant: '<S238>/Constant1'
           */
          windEmulatorStep4_B.rampValue =
            windEmulatorStep4_cal->Constant1_Value_e;
        }

        /* End of Switch: '<S238>/Switch' */

        /* MultiPortSwitch: '<S288>/Multiport Switch' incorporates:
         *  Constant: '<S238>/sidType'
         */
        switch (windEmulatorStep4_cal->sidType_Value) {
         case sidTypeEnum_off:
          /* MultiPortSwitch: '<S288>/Multiport Switch' incorporates:
           *  Constant: '<S288>/Constant'
           */
          windEmulatorStep4_B.MultiportSwitch_h =
            windEmulatorStep4_cal->Constant_Value_im;
          break;

         case sidTypeEnum_manual:
          /* MultiPortSwitch: '<S288>/Multiport Switch' incorporates:
           *  Constant: '<S288>/Constant'
           */
          windEmulatorStep4_B.MultiportSwitch_h =
            windEmulatorStep4_cal->Constant_Value_im;
          break;

         case sidTypeEnum_fromFile:
          /* MultiPortSwitch: '<S288>/Multiport Switch' */
          windEmulatorStep4_B.MultiportSwitch_h =
            windEmulatorStep4_B.BusAssignment_b.stepCounter;
          break;

         default:
          /* MultiPortSwitch: '<S288>/Multiport Switch' incorporates:
           *  Constant: '<S288>/Constant'
           */
          windEmulatorStep4_B.MultiportSwitch_h =
            windEmulatorStep4_cal->Constant_Value_im;
          break;
        }

        /* End of MultiPortSwitch: '<S288>/Multiport Switch' */

        /* Math: '<S288>/Mod' incorporates:
         *  Constant: '<S288>/Length of input'
         */
        q0 = windEmulatorStep4_B.MultiportSwitch_h;
        qY = windEmulatorStep4_cal->Lengthofinput_Value;
        if (qY == 0U) {
          /* Math: '<S288>/Mod' */
          windEmulatorStep4_B.Mod = q0;
        } else {
          /* Math: '<S288>/Mod' */
          windEmulatorStep4_B.Mod = q0 % qY;
        }

        /* End of Math: '<S288>/Mod' */

        /* RelationalOperator: '<S288>/Out of bounds' incorporates:
         *  Constant: '<S288>/Length of input'
         */
        windEmulatorStep4_B.Outofbounds = (windEmulatorStep4_B.Mod >
          windEmulatorStep4_cal->Lengthofinput_Value);
      }

      /* Switch: '<S288>/Switch2' */
      if (windEmulatorStep4_B.Outofbounds) {
        /* Switch: '<S288>/Switch2' incorporates:
         *  Constant: '<S288>/Set bound'
         */
        windEmulatorStep4_B.Switch2 = windEmulatorStep4_cal->Setbound_Value;
      } else {
        /* Switch: '<S288>/Switch2' incorporates:
         *  Inport: '<Root>/inportTorque_Nm'
         */
        windEmulatorStep4_B.Switch2 = windEmulatorStep4_U.inportTorque_Nm;
      }

      /* End of Switch: '<S288>/Switch2' */

      /* RateLimiter: '<S288>/torqueSlewRate' */
      if (windEmulatorStep4_DW.LastMajorTime == (rtInf)) {
        /* RateLimiter: '<S288>/torqueSlewRate' */
        windEmulatorStep4_B.torqueSlewRate = windEmulatorStep4_B.Switch2;
      } else {
        u1 = windEmulatorStep4_M->Timing.t[0] -
          windEmulatorStep4_DW.LastMajorTime;
        riseValLimit = u1 * tmp_c;
        rateLimiterRate = windEmulatorStep4_B.Switch2 -
          windEmulatorStep4_DW.PrevY;
        if (rateLimiterRate > riseValLimit) {
          /* RateLimiter: '<S288>/torqueSlewRate' */
          windEmulatorStep4_B.torqueSlewRate = windEmulatorStep4_DW.PrevY +
            riseValLimit;
        } else {
          riseValLimit = -tmp_c;
          u1 *= riseValLimit;
          if (rateLimiterRate < u1) {
            /* RateLimiter: '<S288>/torqueSlewRate' */
            windEmulatorStep4_B.torqueSlewRate = windEmulatorStep4_DW.PrevY + u1;
          } else {
            /* RateLimiter: '<S288>/torqueSlewRate' */
            windEmulatorStep4_B.torqueSlewRate = windEmulatorStep4_B.Switch2;
          }
        }
      }

      /* MultiPortSwitch: '<S238>/Multiport Switch' incorporates:
       *  Constant: '<S238>/sidType'
       */
      switch (windEmulatorStep4_cal->sidType_Value) {
       case sidTypeEnum_off:
        /* MultiPortSwitch: '<S238>/Multiport Switch' incorporates:
         *  Constant: '<S238>/Constant3'
         */
        windEmulatorStep4_B.MultiportSwitch =
          windEmulatorStep4_cal->Constant3_Value;
        break;

       case sidTypeEnum_manual:
        /* MultiPortSwitch: '<S238>/Multiport Switch' incorporates:
         *  Constant: '<S238>/manualTorqueSetpoint_Nm'
         */
        windEmulatorStep4_B.MultiportSwitch =
          windEmulatorStep4_cal->manualTorqueSetpoint_Nm_Value;
        break;

       case sidTypeEnum_fromFile:
        /* Gain: '<S288>/fromFileTorqueNow_Nm' */
        windEmulatorStep4_B.fromFileTorqueNow_Nm =
          windEmulatorStep4_cal->fromFileTorqueNow_Nm_Gain *
          windEmulatorStep4_B.torqueSlewRate;

        /* MultiPortSwitch: '<S238>/Multiport Switch' */
        windEmulatorStep4_B.MultiportSwitch =
          windEmulatorStep4_B.fromFileTorqueNow_Nm;
        break;

       default:
        /* MultiPortSwitch: '<S238>/Multiport Switch' incorporates:
         *  Constant: '<S238>/Constant3'
         */
        windEmulatorStep4_B.MultiportSwitch =
          windEmulatorStep4_cal->Constant3_Value;
        break;
      }

      /* End of MultiPortSwitch: '<S238>/Multiport Switch' */

      /* Product: '<S238>/Product' */
      windEmulatorStep4_B.Product_g = windEmulatorStep4_B.rampValue *
        windEmulatorStep4_B.MultiportSwitch;

      /* Switch: '<S288>/Switch1' */
      if (windEmulatorStep4_B.Outofbounds) {
        /* Switch: '<S288>/Switch1' incorporates:
         *  Constant: '<S288>/Set bound'
         */
        windEmulatorStep4_B.Switch1 = windEmulatorStep4_cal->Setbound_Value;
      } else {
        /* Switch: '<S288>/Switch1' incorporates:
         *  Inport: '<Root>/inportSpeed_rpm'
         */
        windEmulatorStep4_B.Switch1 = windEmulatorStep4_U.inportSpeed_rpm;
      }

      /* End of Switch: '<S288>/Switch1' */

      /* RateLimiter: '<S288>/speedSlewRate' */
      if (windEmulatorStep4_DW.LastMajorTime_p == (rtInf)) {
        /* RateLimiter: '<S288>/speedSlewRate' */
        windEmulatorStep4_B.speedSlewRate = windEmulatorStep4_B.Switch1;
      } else {
        u1 = windEmulatorStep4_M->Timing.t[0] -
          windEmulatorStep4_DW.LastMajorTime_p;
        riseValLimit = u1 * tmp_b;
        rateLimiterRate = windEmulatorStep4_B.Switch1 -
          windEmulatorStep4_DW.PrevY_a;
        if (rateLimiterRate > riseValLimit) {
          /* RateLimiter: '<S288>/speedSlewRate' */
          windEmulatorStep4_B.speedSlewRate = windEmulatorStep4_DW.PrevY_a +
            riseValLimit;
        } else {
          riseValLimit = -tmp_b;
          u1 *= riseValLimit;
          if (rateLimiterRate < u1) {
            /* RateLimiter: '<S288>/speedSlewRate' */
            windEmulatorStep4_B.speedSlewRate = windEmulatorStep4_DW.PrevY_a +
              u1;
          } else {
            /* RateLimiter: '<S288>/speedSlewRate' */
            windEmulatorStep4_B.speedSlewRate = windEmulatorStep4_B.Switch1;
          }
        }
      }

      /* MultiPortSwitch: '<S238>/Multiport Switch1' incorporates:
       *  Constant: '<S238>/sidType'
       */
      switch (windEmulatorStep4_cal->sidType_Value) {
       case sidTypeEnum_off:
        /* MultiPortSwitch: '<S238>/Multiport Switch1' incorporates:
         *  Constant: '<S238>/Constant3'
         */
        windEmulatorStep4_B.MultiportSwitch1 =
          windEmulatorStep4_cal->Constant3_Value;
        break;

       case sidTypeEnum_manual:
        /* MultiPortSwitch: '<S238>/Multiport Switch1' incorporates:
         *  Constant: '<S238>/manualSpeedSetpoint_rpm'
         */
        windEmulatorStep4_B.MultiportSwitch1 =
          windEmulatorStep4_cal->manualSpeedSetpoint_rpm_Value;
        break;

       case sidTypeEnum_fromFile:
        /* Gain: '<S288>/fromFileSpeedNow_rpm' */
        windEmulatorStep4_B.fromFileSpeedNow_rpm =
          windEmulatorStep4_cal->fromFileSpeedNow_rpm_Gain *
          windEmulatorStep4_B.speedSlewRate;

        /* MultiPortSwitch: '<S238>/Multiport Switch1' */
        windEmulatorStep4_B.MultiportSwitch1 =
          windEmulatorStep4_B.fromFileSpeedNow_rpm;
        break;

       default:
        /* MultiPortSwitch: '<S238>/Multiport Switch1' incorporates:
         *  Constant: '<S238>/Constant3'
         */
        windEmulatorStep4_B.MultiportSwitch1 =
          windEmulatorStep4_cal->Constant3_Value;
        break;
      }

      /* End of MultiPortSwitch: '<S238>/Multiport Switch1' */

      /* Product: '<S238>/Product1' */
      windEmulatorStep4_B.Product1 = windEmulatorStep4_B.rampValue *
        windEmulatorStep4_B.MultiportSwitch1;

      /* Sum: '<S13>/Sum' */
      windEmulatorStep4_B.Sum = windEmulatorStep4_B.Product1 -
        windEmulatorStep4_B.BusAssignment_a.motorSpeed_rpm;

      /* Product: '<S275>/PProd Out' incorporates:
       *  Constant: '<S13>/acs880SpeedPGain'
       */
      windEmulatorStep4_B.PProdOut = windEmulatorStep4_B.Sum *
        windEmulatorStep4_cal->acs880SpeedPGain_Value;

      /* Integrator: '<S270>/Integrator' */
      if (rtsiIsModeUpdateTimeStep(&windEmulatorStep4_M->solverInfo)) {
        e_out = (((windEmulatorStep4_PrevZCX.Integrator_Reset_ZCE == POS_ZCSIG)
                  != windEmulatorStep4_B.BusAssignment_b.resetSidIntegrator) &&
                 (windEmulatorStep4_PrevZCX.Integrator_Reset_ZCE !=
                  UNINITIALIZED_ZCSIG));
        windEmulatorStep4_PrevZCX.Integrator_Reset_ZCE =
          windEmulatorStep4_B.BusAssignment_b.resetSidIntegrator;

        /* evaluate zero-crossings and the level of the reset signal */
        if (e_out || windEmulatorStep4_B.BusAssignment_b.resetSidIntegrator) {
          windEmulatorStep4_X.Integrator_CSTATE =
            windEmulatorStep4_cal->PIDController_InitialConditio_o;
        }
      }

      /* Integrator: '<S270>/Integrator' */
      windEmulatorStep4_B.Integrator = windEmulatorStep4_X.Integrator_CSTATE;

      /* Sum: '<S280>/Sum' */
      windEmulatorStep4_B.Sum_k = windEmulatorStep4_B.PProdOut +
        windEmulatorStep4_B.Integrator;

      /* RelationalOperator: '<S278>/LowerRelop1' incorporates:
       *  Constant: '<S13>/Constant1'
       */
      windEmulatorStep4_B.LowerRelop1 = (windEmulatorStep4_B.Sum_k > tmp_i);

      /* RelationalOperator: '<S278>/UpperRelop' incorporates:
       *  Constant: '<S13>/Constant2'
       */
      windEmulatorStep4_B.UpperRelop = (windEmulatorStep4_B.Sum_k < tmp_j);

      /* Switch: '<S278>/Switch' */
      if (windEmulatorStep4_B.UpperRelop) {
        /* Switch: '<S278>/Switch' incorporates:
         *  Constant: '<S13>/Constant2'
         */
        windEmulatorStep4_B.Switch = tmp_j;
      } else {
        /* Switch: '<S278>/Switch' */
        windEmulatorStep4_B.Switch = windEmulatorStep4_B.Sum_k;
      }

      /* Switch: '<S278>/Switch2' */
      if (windEmulatorStep4_B.LowerRelop1) {
        /* Switch: '<S278>/Switch2' incorporates:
         *  Constant: '<S13>/Constant1'
         */
        windEmulatorStep4_B.Switch2_h = tmp_i;
      } else {
        /* Switch: '<S278>/Switch2' */
        windEmulatorStep4_B.Switch2_h = windEmulatorStep4_B.Switch;
      }

      /* BusAssignment: '<S13>/Bus Assignment' */
      windEmulatorStep4_B.BusAssignment_f.acs800Torque_Nm =
        windEmulatorStep4_B.Product_g;
      windEmulatorStep4_B.BusAssignment_f.acs880Torque_Nm =
        windEmulatorStep4_B.Switch2_h;
      if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
        /* RateLimiter: '<S7>/Rate Limiter' incorporates:
         *  Constant: '<S7>/speedReference'
         */
        rateLimiterRate = windEmulatorStep4_cal->speedReference_Value -
          windEmulatorStep4_DW.PrevY_f;
        if (rateLimiterRate > windEmulatorStep4_cal->RateLimiter_RisingLim *
            windEmulatorStep4_period) {
          /* RateLimiter: '<S7>/Rate Limiter' */
          windEmulatorStep4_B.RateLimiter =
            windEmulatorStep4_cal->RateLimiter_RisingLim *
            windEmulatorStep4_period + windEmulatorStep4_DW.PrevY_f;
        } else if (rateLimiterRate <
                   windEmulatorStep4_cal->RateLimiter_FallingLim *
                   windEmulatorStep4_period) {
          /* RateLimiter: '<S7>/Rate Limiter' */
          windEmulatorStep4_B.RateLimiter =
            windEmulatorStep4_cal->RateLimiter_FallingLim *
            windEmulatorStep4_period + windEmulatorStep4_DW.PrevY_f;
        } else {
          /* RateLimiter: '<S7>/Rate Limiter' */
          windEmulatorStep4_B.RateLimiter =
            windEmulatorStep4_cal->speedReference_Value;
        }

        windEmulatorStep4_DW.PrevY_f = windEmulatorStep4_B.RateLimiter;

        /* End of RateLimiter: '<S7>/Rate Limiter' */

        /* Gain: '<S7>/f->w' incorporates:
         *  Constant: '<S7>/excForceFreq_Hz'
         */
        windEmulatorStep4_B.fw = windEmulatorStep4_cal->fw_Gain *
          windEmulatorStep4_cal->excForceFreq_Hz_Value;

        /* Product: '<S7>/Product4' */
        windEmulatorStep4_B.Product4 = windEmulatorStep4_B.fw *
          windEmulatorStep4_B.BusAssignment_b.time;

        /* Trigonometry: '<S7>/Sin2' */
        windEmulatorStep4_B.Sin2 = std::sin(windEmulatorStep4_B.Product4);

        /* Gain: '<S7>/excForceAmpNow_N' incorporates:
         *  Constant: '<S7>/excForceAmp_N'
         */
        windEmulatorStep4_B.excForceAmpNow_N =
          windEmulatorStep4_cal->excForceAmpNow_N_Gain *
          windEmulatorStep4_cal->excForceAmp_N_Value;

        /* Product: '<S7>/ExcitationForce_N' */
        windEmulatorStep4_B.waveBotExcitationForce_N = windEmulatorStep4_B.Sin2 *
          windEmulatorStep4_B.excForceAmpNow_N;

        /* Switch: '<S7>/Switch' */
        if (windEmulatorStep4_B.BusAssignment_b.runHil) {
          /* Switch: '<S7>/Switch' */
          windEmulatorStep4_B.rampValue_a =
            windEmulatorStep4_B.BusAssignment_b.ramp;
        } else {
          /* Switch: '<S7>/Switch' incorporates:
           *  Constant: '<S7>/Constant1'
           */
          windEmulatorStep4_B.rampValue_a =
            windEmulatorStep4_cal->Constant1_Value;
        }

        /* End of Switch: '<S7>/Switch' */

        /* Product: '<S7>/Product' */
        windEmulatorStep4_B.Product_m =
          windEmulatorStep4_B.waveBotExcitationForce_N *
          windEmulatorStep4_B.rampValue_a;

        /* S-Function (slecatpdorx): '<S8>/EtherCAT PDO Receive8' */
        {
          /*------------ S-Function Block: <S8>/EtherCAT PDO Receive8 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.ACS800DcBusVoltage;
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
        windEmulatorStep4_B.DataTypeConversion4 =
          windEmulatorStep4_B.ACS800DcBusVoltage;

        /* Gain: '<S8>/Gain3' */
        windEmulatorStep4_B.Gain3_g = *get_acs800DcBusVoltsScaling() *
          windEmulatorStep4_B.DataTypeConversion4;

        /* S-Function (slecatpdorx): '<S8>/EtherCAT PDO Receive9' */
        {
          /*------------ S-Function Block: <S8>/EtherCAT PDO Receive9 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.ACS800Frequency;
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
        windEmulatorStep4_B.DataTypeConversion3 =
          windEmulatorStep4_B.ACS800Frequency;

        /* Gain: '<S8>/Gain4' */
        windEmulatorStep4_B.Gain4_p = *get_acs800FreqScaling() *
          windEmulatorStep4_B.DataTypeConversion3;

        /* S-Function (slecatpdorx): '<S8>/EtherCAT PDO Receive10' */
        {
          /*------------ S-Function Block: <S8>/EtherCAT PDO Receive10 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.ACS800Temperature;
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
        windEmulatorStep4_B.DataTypeConversion5 =
          windEmulatorStep4_B.ACS800Temperature;

        /* Gain: '<S8>/Gain5' */
        windEmulatorStep4_B.Gain5_g = *get_acs800TempScaling() *
          windEmulatorStep4_B.DataTypeConversion5;

        /* S-Function (slecatpdorx): '<S8>/EtherCAT PDO Receive11' */
        {
          /*------------ S-Function Block: <S8>/EtherCAT PDO Receive11 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.ACS800ActualFeedback;
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
        windEmulatorStep4_B.DataTypeConversion =
          windEmulatorStep4_B.ACS800ActualFeedback;

        /* Gain: '<S8>/Gain' */
        windEmulatorStep4_B.Gain_c = *get_acs800SpeedScaling() *
          windEmulatorStep4_B.DataTypeConversion;

        /* S-Function (slecatpdorx): '<S8>/EtherCAT PDO Receive12' */
        {
          /*------------ S-Function Block: <S8>/EtherCAT PDO Receive12 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)&windEmulatorStep4_B.ACS800Torque;
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
        windEmulatorStep4_B.DataTypeConversion1_e =
          windEmulatorStep4_B.ACS800Torque;

        /* Gain: '<S8>/Gain1' */
        windEmulatorStep4_B.Gain1_f = *get_acs800TorqueScaling() *
          windEmulatorStep4_B.DataTypeConversion1_e;

        /* S-Function (slecatpdorx): '<S8>/EtherCAT PDO Receive13' */
        {
          /*------------ S-Function Block: <S8>/EtherCAT PDO Receive13 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)&windEmulatorStep4_B.ACS800Power;
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
        windEmulatorStep4_B.DataTypeConversion2_i =
          windEmulatorStep4_B.ACS800Power;

        /* Gain: '<S8>/Gain2' */
        windEmulatorStep4_B.Gain2_d = *get_acs800PowerScaling() *
          windEmulatorStep4_B.DataTypeConversion2_i;

        /* S-Function (slecatpdorx): '<S8>/EtherCAT PDO Receive7' */
        {
          /*------------ S-Function Block: <S8>/EtherCAT PDO Receive7 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.ACS800StatusWord;
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

        /* BusAssignment: '<S8>/Bus Assignment' incorporates:
         *  Constant: '<S8>/Constant'
         */
        windEmulatorStep4_B.BusAssignment_h = *get_acs800SignalStruct();

        /* BusAssignment: '<S8>/Bus Assignment' */
        windEmulatorStep4_B.BusAssignment_h.dcBusVoltage_V =
          windEmulatorStep4_B.Gain3_g;
        windEmulatorStep4_B.BusAssignment_h.frequency_Hz =
          windEmulatorStep4_B.Gain4_p;
        windEmulatorStep4_B.BusAssignment_h.temperature =
          windEmulatorStep4_B.Gain5_g;
        windEmulatorStep4_B.BusAssignment_h.motorSpeed_rpm =
          windEmulatorStep4_B.Gain_c;
        windEmulatorStep4_B.BusAssignment_h.motorTorque_Nm =
          windEmulatorStep4_B.Gain1_f;
        windEmulatorStep4_B.BusAssignment_h.shaftPower_W =
          windEmulatorStep4_B.Gain2_d;
        windEmulatorStep4_B.BusAssignment_h.statusWord =
          windEmulatorStep4_B.ACS800StatusWord;

        /* BusAssignment: '<S7>/Bus Assignment' incorporates:
         *  Constant: '<S7>/Constant'
         */
        windEmulatorStep4_B.BusAssignment_c = *get_hptoCtrlStruct();

        /* BusAssignment: '<S7>/Bus Assignment' */
        windEmulatorStep4_B.BusAssignment_c.speedRef_rpm =
          windEmulatorStep4_B.RateLimiter;
        windEmulatorStep4_B.BusAssignment_c.excForce_N =
          windEmulatorStep4_B.Product_m;
        windEmulatorStep4_B.BusAssignment_c.genSpeedActual =
          windEmulatorStep4_B.BusAssignment_h.motorSpeed_rpm;
        windEmulatorStep4_B.BusAssignment_c.speedCtrlReset =
          windEmulatorStep4_B.BusAssignment_b.resetHilIntegrator;

        /* Memory: '<S117>/Memory' */
        windEmulatorStep4_B.Memory_f =
          windEmulatorStep4_DW.Memory_PreviousInput_kk;
      }

      /* Sum: '<S64>/Add' */
      windEmulatorStep4_B.Add =
        windEmulatorStep4_B.BusAssignment_c.genSpeedActual -
        windEmulatorStep4_B.BusAssignment_c.speedRef_rpm;

      /* Product: '<S115>/Product' incorporates:
       *  Constant: '<S115>/Constant1'
       */
      riseValLimit = -tmp_9->PG;

      /* Product: '<S115>/Product' */
      windEmulatorStep4_B.ControlSignal31 = riseValLimit *
        windEmulatorStep4_B.Add;

      /* RelationalOperator: '<S115>/Relational Operator' */
      windEmulatorStep4_B.RelationalOperator =
        (windEmulatorStep4_B.BusAssignment_c.speedRef_rpm <=
         windEmulatorStep4_B.BusAssignment_c.genSpeedActual);

      /* CombinatorialLogic: '<S117>/Logic' incorporates:
       *  Constant: '<S115>/Constant'
       */
      e_out = windEmulatorStep4_B.RelationalOperator;
      rowIdx = e_out;
      e_out = windEmulatorStep4_cal->Constant_Value_b;
      rowIdx = static_cast<int32_T>((static_cast<uint32_T>(rowIdx) << 1) + e_out);
      e_out = windEmulatorStep4_B.Memory_f;
      rowIdx = static_cast<int32_T>((static_cast<uint32_T>(rowIdx) << 1) + e_out);
      windEmulatorStep4_B.Logic[0U] = windEmulatorStep4_cal->Logic_table[
        static_cast<uint32_T>(rowIdx)];
      windEmulatorStep4_B.Logic[1U] = windEmulatorStep4_cal->Logic_table[
        static_cast<uint32_T>(rowIdx) + 8U];
      if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
        /* RateLimiter: '<S61>/Rate Limiter1' incorporates:
         *  Constant: '<S61>/shaftSpeedRefMin'
         */
        rateLimiterRate = windEmulatorStep4_cal->shaftSpeedRefMin_Value -
          windEmulatorStep4_DW.PrevY_m;
        if (rateLimiterRate > windEmulatorStep4_cal->RateLimiter1_RisingLim *
            windEmulatorStep4_period) {
          /* RateLimiter: '<S61>/Rate Limiter1' */
          windEmulatorStep4_B.RateLimiter1 =
            windEmulatorStep4_cal->RateLimiter1_RisingLim *
            windEmulatorStep4_period + windEmulatorStep4_DW.PrevY_m;
        } else if (rateLimiterRate <
                   windEmulatorStep4_cal->RateLimiter1_FallingLim *
                   windEmulatorStep4_period) {
          /* RateLimiter: '<S61>/Rate Limiter1' */
          windEmulatorStep4_B.RateLimiter1 =
            windEmulatorStep4_cal->RateLimiter1_FallingLim *
            windEmulatorStep4_period + windEmulatorStep4_DW.PrevY_m;
        } else {
          /* RateLimiter: '<S61>/Rate Limiter1' */
          windEmulatorStep4_B.RateLimiter1 =
            windEmulatorStep4_cal->shaftSpeedRefMin_Value;
        }

        windEmulatorStep4_DW.PrevY_m = windEmulatorStep4_B.RateLimiter1;

        /* End of RateLimiter: '<S61>/Rate Limiter1' */

        /* Memory: '<S118>/Memory' */
        windEmulatorStep4_B.Memory_in =
          windEmulatorStep4_DW.Memory_PreviousInput_h;
      }

      /* Sum: '<S64>/Add1' */
      windEmulatorStep4_B.Add1 =
        windEmulatorStep4_B.BusAssignment_c.genSpeedActual -
        windEmulatorStep4_B.RateLimiter1;

      /* Product: '<S116>/Product' incorporates:
       *  Constant: '<S116>/Constant1'
       */
      riseValLimit = -tmp_8;

      /* Product: '<S116>/Product' */
      windEmulatorStep4_B.ControlSignal31_o = riseValLimit *
        windEmulatorStep4_B.Add1;

      /* RelationalOperator: '<S116>/Relational Operator' */
      windEmulatorStep4_B.RelationalOperator_k =
        (windEmulatorStep4_B.RateLimiter1 >=
         windEmulatorStep4_B.BusAssignment_c.genSpeedActual);

      /* CombinatorialLogic: '<S118>/Logic' incorporates:
       *  Constant: '<S116>/Constant'
       */
      e_out = windEmulatorStep4_B.RelationalOperator_k;
      rowIdx = e_out;
      e_out = windEmulatorStep4_cal->Constant_Value_k;
      rowIdx = static_cast<int32_T>((static_cast<uint32_T>(rowIdx) << 1) + e_out);
      e_out = windEmulatorStep4_B.Memory_in;
      rowIdx = static_cast<int32_T>((static_cast<uint32_T>(rowIdx) << 1) + e_out);
      windEmulatorStep4_B.Logic_g[0U] = windEmulatorStep4_cal->Logic_table_o[
        static_cast<uint32_T>(rowIdx)];
      windEmulatorStep4_B.Logic_g[1U] = windEmulatorStep4_cal->Logic_table_o[
        static_cast<uint32_T>(rowIdx) + 8U];

      /* Switch: '<S64>/Switch' incorporates:
       *  Switch: '<S116>/Switch'
       *  Switch: '<S64>/Switch1'
       */
      if (windEmulatorStep4_B.Add > windEmulatorStep4_cal->Switch_Threshold_j) {
        /* Switch: '<S115>/Switch' */
        if (windEmulatorStep4_B.Logic[0]) {
          /* Saturate: '<S115>/Saturation' */
          riseValLimit = windEmulatorStep4_B.ControlSignal31;
          u1 = windEmulatorStep4_cal->Saturation_LowerSat;
          rateLimiterRate = windEmulatorStep4_cal->Saturation_UpperSat;
          if (riseValLimit > rateLimiterRate) {
            /* Saturate: '<S115>/Saturation' */
            windEmulatorStep4_B.Saturation_a3 = rateLimiterRate;
          } else if (riseValLimit < u1) {
            /* Saturate: '<S115>/Saturation' */
            windEmulatorStep4_B.Saturation_a3 = u1;
          } else {
            /* Saturate: '<S115>/Saturation' */
            windEmulatorStep4_B.Saturation_a3 = riseValLimit;
          }

          /* End of Saturate: '<S115>/Saturation' */

          /* Switch: '<S115>/Switch' */
          windEmulatorStep4_B.ControlSignal3_f =
            windEmulatorStep4_B.Saturation_a3;
        } else {
          /* Switch: '<S115>/Switch' */
          windEmulatorStep4_B.ControlSignal3_f =
            windEmulatorStep4_B.ControlSignal31;
        }

        /* End of Switch: '<S115>/Switch' */

        /* Switch: '<S64>/Switch' */
        windEmulatorStep4_B.Switch_g = windEmulatorStep4_B.ControlSignal3_f;
      } else {
        if (windEmulatorStep4_B.Add1 > windEmulatorStep4_cal->Switch1_Threshold)
        {
          /* Switch: '<S64>/Switch1' incorporates:
           *  Constant: '<S64>/Constant1'
           */
          windEmulatorStep4_B.Switch1_o =
            windEmulatorStep4_cal->Constant1_Value_c;
        } else {
          if (windEmulatorStep4_B.Logic_g[0]) {
            /* Saturate: '<S116>/Saturation' incorporates:
             *  Switch: '<S116>/Switch'
             *  Switch: '<S64>/Switch1'
             */
            riseValLimit = windEmulatorStep4_B.ControlSignal31_o;
            u1 = windEmulatorStep4_cal->Saturation_LowerSat_n;
            rateLimiterRate = windEmulatorStep4_cal->Saturation_UpperSat_g;
            if (riseValLimit > rateLimiterRate) {
              /* Saturate: '<S116>/Saturation' */
              windEmulatorStep4_B.Saturation_k = rateLimiterRate;
            } else if (riseValLimit < u1) {
              /* Saturate: '<S116>/Saturation' */
              windEmulatorStep4_B.Saturation_k = u1;
            } else {
              /* Saturate: '<S116>/Saturation' */
              windEmulatorStep4_B.Saturation_k = riseValLimit;
            }

            /* End of Saturate: '<S116>/Saturation' */

            /* Switch: '<S116>/Switch' incorporates:
             *  Switch: '<S64>/Switch1'
             */
            windEmulatorStep4_B.ControlSignal3_e =
              windEmulatorStep4_B.Saturation_k;
          } else {
            /* Switch: '<S116>/Switch' incorporates:
             *  Switch: '<S64>/Switch1'
             */
            windEmulatorStep4_B.ControlSignal3_e =
              windEmulatorStep4_B.ControlSignal31_o;
          }

          /* Switch: '<S64>/Switch1' */
          windEmulatorStep4_B.Switch1_o = windEmulatorStep4_B.ControlSignal3_e;
        }

        /* Switch: '<S64>/Switch' incorporates:
         *  Switch: '<S116>/Switch'
         *  Switch: '<S64>/Switch1'
         */
        windEmulatorStep4_B.Switch_g = windEmulatorStep4_B.Switch1_o;
      }

      /* End of Switch: '<S64>/Switch' */

      /* Gain: '<S64>/Gain' */
      windEmulatorStep4_B.Gain_l = tmp_m * windEmulatorStep4_B.Switch_g;

      /* RateLimiter: '<S61>/Rate Limiter' */
      if (windEmulatorStep4_DW.LastMajorTime_m == (rtInf)) {
        /* RateLimiter: '<S61>/Rate Limiter' */
        windEmulatorStep4_B.RateLimiter_b = windEmulatorStep4_B.Gain_l;
      } else {
        u1 = windEmulatorStep4_M->Timing.t[0] -
          windEmulatorStep4_DW.LastMajorTime_m;
        riseValLimit = u1 * tmp_7;
        rateLimiterRate = windEmulatorStep4_B.Gain_l -
          windEmulatorStep4_DW.PrevY_l;
        if (rateLimiterRate > riseValLimit) {
          /* RateLimiter: '<S61>/Rate Limiter' */
          windEmulatorStep4_B.RateLimiter_b = windEmulatorStep4_DW.PrevY_l +
            riseValLimit;
        } else {
          riseValLimit = -tmp_7;
          u1 *= riseValLimit;
          if (rateLimiterRate < u1) {
            /* RateLimiter: '<S61>/Rate Limiter' */
            windEmulatorStep4_B.RateLimiter_b = windEmulatorStep4_DW.PrevY_l +
              u1;
          } else {
            /* RateLimiter: '<S61>/Rate Limiter' */
            windEmulatorStep4_B.RateLimiter_b = windEmulatorStep4_B.Gain_l;
          }
        }
      }

      /* End of RateLimiter: '<S61>/Rate Limiter' */

      /* Sum: '<S63>/Sum' */
      windEmulatorStep4_B.wError =
        windEmulatorStep4_B.BusAssignment_c.genSpeedActual -
        windEmulatorStep4_B.BusAssignment_c.speedRef_rpm;
      if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
        /* Gain: '<S102>/Proportional Gain' */
        windEmulatorStep4_B.ProportionalGain = tmp_9->PG *
          windEmulatorStep4_B.wError;

        /* DiscreteIntegrator: '<S97>/Integrator' */
        if (windEmulatorStep4_B.BusAssignment_c.speedCtrlReset ||
            (windEmulatorStep4_DW.Integrator_PrevResetState != 0)) {
          windEmulatorStep4_DW.Integrator_DSTATE =
            windEmulatorStep4_cal->PIDController_InitialConditio_a;
        }

        /* DiscreteIntegrator: '<S97>/Integrator' */
        windEmulatorStep4_B.Integrator_b =
          windEmulatorStep4_DW.Integrator_DSTATE;

        /* Gain: '<S91>/Derivative Gain' */
        windEmulatorStep4_B.DerivativeGain =
          windEmulatorStep4_cal->PIDController_D * windEmulatorStep4_B.wError;

        /* DiscreteIntegrator: '<S92>/Filter' */
        if (windEmulatorStep4_B.BusAssignment_c.speedCtrlReset ||
            (windEmulatorStep4_DW.Filter_PrevResetState != 0)) {
          windEmulatorStep4_DW.Filter_DSTATE =
            windEmulatorStep4_cal->PIDController_InitialConditionF;
        }

        /* DiscreteIntegrator: '<S92>/Filter' */
        windEmulatorStep4_B.Filter = windEmulatorStep4_DW.Filter_DSTATE;

        /* Sum: '<S92>/SumD' */
        windEmulatorStep4_B.SumD = windEmulatorStep4_B.DerivativeGain -
          windEmulatorStep4_B.Filter;

        /* Gain: '<S100>/Filter Coefficient' */
        windEmulatorStep4_B.FilterCoefficient =
          windEmulatorStep4_cal->PIDController_N * windEmulatorStep4_B.SumD;

        /* Sum: '<S107>/Sum' */
        windEmulatorStep4_B.Sum_c = (windEmulatorStep4_B.ProportionalGain +
          windEmulatorStep4_B.Integrator_b) +
          windEmulatorStep4_B.FilterCoefficient;

        /* RelationalOperator: '<S105>/LowerRelop1' incorporates:
         *  Constant: '<S63>/Constant'
         */
        windEmulatorStep4_B.LowerRelop1_g = (windEmulatorStep4_B.Sum_c > tmp_6);

        /* RelationalOperator: '<S105>/UpperRelop' incorporates:
         *  Constant: '<S63>/Constant1'
         */
        riseValLimit = -tmp_6;

        /* RelationalOperator: '<S105>/UpperRelop' */
        windEmulatorStep4_B.UpperRelop_g = (windEmulatorStep4_B.Sum_c <
          riseValLimit);

        /* Switch: '<S105>/Switch' */
        if (windEmulatorStep4_B.UpperRelop_g) {
          /* Switch: '<S105>/Switch' incorporates:
           *  Constant: '<S63>/Constant1'
           */
          windEmulatorStep4_B.Switch_i = -tmp_6;
        } else {
          /* Switch: '<S105>/Switch' */
          windEmulatorStep4_B.Switch_i = windEmulatorStep4_B.Sum_c;
        }

        /* End of Switch: '<S105>/Switch' */

        /* Switch: '<S105>/Switch2' */
        if (windEmulatorStep4_B.LowerRelop1_g) {
          /* Switch: '<S105>/Switch2' incorporates:
           *  Constant: '<S63>/Constant'
           */
          windEmulatorStep4_B.Switch2_k = tmp_6;
        } else {
          /* Switch: '<S105>/Switch2' */
          windEmulatorStep4_B.Switch2_k = windEmulatorStep4_B.Switch_i;
        }

        /* End of Switch: '<S105>/Switch2' */

        /* Gain: '<S63>/Gain2' */
        windEmulatorStep4_B.ContolTorque = windEmulatorStep4_cal->Gain2_Gain_o *
          windEmulatorStep4_B.Switch2_k;
      }

      /* Step: '<S58>/Step' */
      rateLimiterRate = windEmulatorStep4_M->Timing.t[0];
      if (rateLimiterRate < windEmulatorStep4_cal->Ramp_start) {
        /* Step: '<S58>/Step' */
        windEmulatorStep4_B.Step = windEmulatorStep4_cal->Step_Y0;
      } else {
        /* Step: '<S58>/Step' */
        windEmulatorStep4_B.Step = windEmulatorStep4_cal->Ramp_slope;
      }

      /* End of Step: '<S58>/Step' */

      /* Clock: '<S58>/Clock' */
      windEmulatorStep4_B.Clock = windEmulatorStep4_M->Timing.t[0];

      /* Sum: '<S58>/Sum' incorporates:
       *  Constant: '<S58>/Constant'
       */
      windEmulatorStep4_B.Sum_m = windEmulatorStep4_B.Clock -
        windEmulatorStep4_cal->Ramp_start;

      /* Product: '<S58>/Product' */
      windEmulatorStep4_B.Product_a = windEmulatorStep4_B.Step *
        windEmulatorStep4_B.Sum_m;

      /* Sum: '<S58>/Output' incorporates:
       *  Constant: '<S58>/Constant1'
       */
      windEmulatorStep4_B.Output = windEmulatorStep4_B.Product_a +
        windEmulatorStep4_cal->Ramp_InitialOutput;

      /* Saturate: '<S53>/Saturation' */
      riseValLimit = windEmulatorStep4_B.Output;
      u1 = windEmulatorStep4_cal->Saturation_LowerSat_em;
      rateLimiterRate = windEmulatorStep4_cal->Saturation_UpperSat_da;
      if (riseValLimit > rateLimiterRate) {
        /* Saturate: '<S53>/Saturation' */
        windEmulatorStep4_B.Saturation = rateLimiterRate;
      } else if (riseValLimit < u1) {
        /* Saturate: '<S53>/Saturation' */
        windEmulatorStep4_B.Saturation = u1;
      } else {
        /* Saturate: '<S53>/Saturation' */
        windEmulatorStep4_B.Saturation = riseValLimit;
      }

      /* End of Saturate: '<S53>/Saturation' */
      if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
        /* Gain: '<S59>/kDampingNow' incorporates:
         *  Constant: '<S59>/kDamping'
         */
        windEmulatorStep4_B.kDampingNow =
          windEmulatorStep4_cal->kDampingNow_Gain *
          windEmulatorStep4_cal->kDamping_Value;

        /* RateLimiter: '<S178>/Rate Limiter' incorporates:
         *  Constant: '<S178>/Constant1'
         */
        riseValLimit = 2000.0 * tmp_a;
        rateLimiterRate = riseValLimit - windEmulatorStep4_DW.PrevY_k;
        if (rateLimiterRate > windEmulatorStep4_cal->RateLimiter_RisingLim_j *
            windEmulatorStep4_period) {
          /* RateLimiter: '<S178>/Rate Limiter' */
          windEmulatorStep4_B.RateLimiter_a =
            windEmulatorStep4_cal->RateLimiter_RisingLim_j *
            windEmulatorStep4_period + windEmulatorStep4_DW.PrevY_k;
        } else if (rateLimiterRate <
                   windEmulatorStep4_cal->RateLimiter_FallingLim_e *
                   windEmulatorStep4_period) {
          /* RateLimiter: '<S178>/Rate Limiter' */
          windEmulatorStep4_B.RateLimiter_a =
            windEmulatorStep4_cal->RateLimiter_FallingLim_e *
            windEmulatorStep4_period + windEmulatorStep4_DW.PrevY_k;
        } else {
          /* RateLimiter: '<S178>/Rate Limiter' */
          windEmulatorStep4_B.RateLimiter_a = 2000.0 * tmp_a;
        }

        windEmulatorStep4_DW.PrevY_k = windEmulatorStep4_B.RateLimiter_a;

        /* Gain: '<S54>/Gain' */
        windEmulatorStep4_B.Gain_d = tmp_a *
          windEmulatorStep4_B.BusAssignment_c.genSpeedActual;

        /* SimscapeInputBlock: '<S222>/INPUT_1_1_1' incorporates:
         *  SimscapeInputBlock: '<S222>/INPUT_5_1_1'
         */
        if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
          windEmulatorStep4_B.INPUT_1_1_1[0] = windEmulatorStep4_B.RateLimiter_a;
          windEmulatorStep4_B.INPUT_1_1_1[1] = 0.0;
          windEmulatorStep4_B.INPUT_1_1_1[2] = 0.0;
          windEmulatorStep4_DW.INPUT_1_1_1_Discrete[0] =
            !(windEmulatorStep4_B.INPUT_1_1_1[0] ==
              windEmulatorStep4_DW.INPUT_1_1_1_Discrete[1]);
          windEmulatorStep4_DW.INPUT_1_1_1_Discrete[1] =
            windEmulatorStep4_B.INPUT_1_1_1[0];
          windEmulatorStep4_B.INPUT_1_1_1[0] =
            windEmulatorStep4_DW.INPUT_1_1_1_Discrete[1];
          windEmulatorStep4_B.INPUT_1_1_1[3] =
            windEmulatorStep4_DW.INPUT_1_1_1_Discrete[0];
          windEmulatorStep4_B.INPUT_5_1_1[0] = windEmulatorStep4_B.Gain_d;
          windEmulatorStep4_B.INPUT_5_1_1[1] = 0.0;
          windEmulatorStep4_B.INPUT_5_1_1[2] = 0.0;
          windEmulatorStep4_DW.INPUT_5_1_1_Discrete[0] =
            !(windEmulatorStep4_B.INPUT_5_1_1[0] ==
              windEmulatorStep4_DW.INPUT_5_1_1_Discrete[1]);
          windEmulatorStep4_DW.INPUT_5_1_1_Discrete[1] =
            windEmulatorStep4_B.INPUT_5_1_1[0];
          windEmulatorStep4_B.INPUT_5_1_1[0] =
            windEmulatorStep4_DW.INPUT_5_1_1_Discrete[1];
          windEmulatorStep4_B.INPUT_5_1_1[3] =
            windEmulatorStep4_DW.INPUT_5_1_1_Discrete[0];
        }

        /* End of SimscapeInputBlock: '<S222>/INPUT_1_1_1' */
      }

      /* StateSpace: '<S211>/Internal' */
      windEmulatorStep4_B.Internal = 0.0;

      /* StateSpace: '<S211>/Internal' */
      for (q0 = windEmulatorStep4_cal->Internal_C_jc[0U]; q0 <
           windEmulatorStep4_cal->Internal_C_jc[1U]; q0++) {
        /* StateSpace: '<S211>/Internal' */
        windEmulatorStep4_B.Internal += windEmulatorStep4_cal->Internal_C_pr *
          windEmulatorStep4_X.Internal_CSTATE[0U];
      }

      for (q0 = windEmulatorStep4_cal->Internal_C_jc[1U]; q0 <
           windEmulatorStep4_cal->Internal_C_jc[2U]; q0++) {
        /* StateSpace: '<S211>/Internal' */
        windEmulatorStep4_B.Internal += windEmulatorStep4_cal->Internal_C_pr *
          windEmulatorStep4_X.Internal_CSTATE[1U];
      }

      for (q0 = windEmulatorStep4_cal->Internal_C_jc[2U]; q0 <
           windEmulatorStep4_cal->Internal_C_jc[3U]; q0++) {
        /* StateSpace: '<S211>/Internal' */
        windEmulatorStep4_B.Internal += windEmulatorStep4_cal->Internal_C_pr *
          windEmulatorStep4_X.Internal_CSTATE[2U];
      }

      /* SimscapeInputBlock: '<S222>/INPUT_2_1_1' */
      windEmulatorStep4_B.INPUT_2_1_1[0] = windEmulatorStep4_B.Internal;
      windEmulatorStep4_B.INPUT_2_1_1[1] = 0.0;
      windEmulatorStep4_B.INPUT_2_1_1[2] = 0.0;
      windEmulatorStep4_B.INPUT_2_1_1[3] = 0.0;

      /* StateSpace: '<S227>/Internal' */
      windEmulatorStep4_B.Internal_j = 0.0;

      /* StateSpace: '<S227>/Internal' */
      for (q0 = windEmulatorStep4_cal->Internal_C_jc_k[0U]; q0 <
           windEmulatorStep4_cal->Internal_C_jc_k[1U]; q0++) {
        /* StateSpace: '<S227>/Internal' */
        windEmulatorStep4_B.Internal_j += windEmulatorStep4_cal->Internal_C_pr_a
          * windEmulatorStep4_X.Internal_CSTATE_j;
      }

      /* Gain: '<S190>/Gain' */
      windEmulatorStep4_B.Gain_lp = windEmulatorStep4_cal->Gain_Gain *
        windEmulatorStep4_B.Internal_j;

      /* SimscapeInputBlock: '<S222>/INPUT_4_1_1' */
      windEmulatorStep4_B.INPUT_4_1_1[0] = windEmulatorStep4_B.Gain_lp;
      windEmulatorStep4_B.INPUT_4_1_1[1] = 0.0;
      windEmulatorStep4_B.INPUT_4_1_1[2] = 0.0;
      windEmulatorStep4_B.INPUT_4_1_1[3] = 0.0;

      /* StateSpace: '<S223>/Internal' */
      windEmulatorStep4_B.Internal_h = 0.0;

      /* StateSpace: '<S223>/Internal' */
      for (q0 = windEmulatorStep4_cal->Internal_C_jc_h[0U]; q0 <
           windEmulatorStep4_cal->Internal_C_jc_h[1U]; q0++) {
        /* StateSpace: '<S223>/Internal' */
        windEmulatorStep4_B.Internal_h += windEmulatorStep4_cal->Internal_C_pr_m
          * windEmulatorStep4_X.Internal_CSTATE_a;
      }

      /* Gain: '<S189>/Gain' */
      windEmulatorStep4_B.Gain_g = windEmulatorStep4_cal->Gain_Gain_j *
        windEmulatorStep4_B.Internal_h;

      /* SimscapeInputBlock: '<S222>/INPUT_3_1_1' */
      windEmulatorStep4_B.INPUT_3_1_1[0] = windEmulatorStep4_B.Gain_g;
      windEmulatorStep4_B.INPUT_3_1_1[1] = 0.0;
      windEmulatorStep4_B.INPUT_3_1_1[2] = 0.0;
      windEmulatorStep4_B.INPUT_3_1_1[3] = 0.0;
      if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
        NeuDiagnosticTree *diagTree;
        char *msg;

        /* SimscapeRtp: '<S188>/RTP_1' incorporates:
         *  Constant: '<S178>/Subsystem_around_RTP_D2E1D090_liquid_pressure'
         *  Constant: '<S178>/Subsystem_around_RTP_D2E1D090_liquid_volume'
         *  Constant: '<S54>/Subsystem_around_RTP_D290B913_fluid_volume'
         */
        if (windEmulatorStep4_DW.RTP_1_SetParametersNeeded) {
          tmp[0] = windEmulatorStep4_cal->RTP_D290B913_fluid_volume_Value;
          tmp[1] = windEmulatorStep4_cal->RTP_D2E1D090_liquid_pressure_Va;
          tmp[2] = windEmulatorStep4_cal->RTP_D2E1D090_liquid_volume_Valu;
          parameterBundle_mRealParameters = &tmp[0];
          rtpManager = static_cast<NeslRtpManager *>
            (windEmulatorStep4_DW.RTP_1_RtpManager);
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
          e_out = nesl_rtp_manager_set_rtps(rtpManager,
            windEmulatorStep4_M->Timing.t[0], expl_temp, diag);
          if (!e_out) {
            e_out = error_buffer_is_empty(rtmGetErrorStatus(windEmulatorStep4_M));
            if (e_out) {
              msg = rtw_diagnostics_msg(diagTree);
              rtmSetErrorStatus(windEmulatorStep4_M, msg);
            }
          }
        }

        windEmulatorStep4_DW.RTP_1_SetParametersNeeded = false;

        /* End of SimscapeRtp: '<S188>/RTP_1' */

        /* SimscapeExecutionBlock: '<S222>/STATE_1' */
        simulationData = static_cast<NeslSimulationData *>
          (windEmulatorStep4_DW.STATE_1_SimData);
        time = windEmulatorStep4_M->Timing.t[0];
        simulationData->mData->mTime.mN = 1;
        simulationData->mData->mTime.mX = &time;
        simulationData->mData->mContStates.mN = 0;
        simulationData->mData->mContStates.mX = NULL;
        simulationData->mData->mDiscStates.mN = 22;
        simulationData->mData->mDiscStates.mX =
          &windEmulatorStep4_DW.STATE_1_Discrete[0];
        simulationData->mData->mModeVector.mN = 17;
        simulationData->mData->mModeVector.mX =
          &windEmulatorStep4_DW.STATE_1_Modes[0];
        e_out = false;
        simulationData->mData->mFoundZcEvents = e_out;
        simulationData->mData->mIsMajorTimeStep = true;
        e_out = false;
        simulationData->mData->mIsSolverAssertCheck = e_out;
        simulationData->mData->mIsSolverCheckingCIC = false;
        simulationData->mData->mIsComputingJacobian = false;
        simulationData->mData->mIsEvaluatingF0 = false;
        simulationData->mData->mIsSolverRequestingReset = false;
        simulationData->mData->mIsModeUpdateTimeStep = true;
        tmp_1[0] = 0;
        tmp_0[0] = windEmulatorStep4_B.INPUT_1_1_1[0];
        tmp_0[1] = windEmulatorStep4_B.INPUT_1_1_1[1];
        tmp_0[2] = windEmulatorStep4_B.INPUT_1_1_1[2];
        tmp_0[3] = windEmulatorStep4_B.INPUT_1_1_1[3];
        tmp_1[1] = 4;
        tmp_0[4] = windEmulatorStep4_B.INPUT_5_1_1[0];
        tmp_0[5] = windEmulatorStep4_B.INPUT_5_1_1[1];
        tmp_0[6] = windEmulatorStep4_B.INPUT_5_1_1[2];
        tmp_0[7] = windEmulatorStep4_B.INPUT_5_1_1[3];
        tmp_1[2] = 8;
        tmp_0[8] = windEmulatorStep4_B.INPUT_2_1_1[0];
        tmp_0[9] = windEmulatorStep4_B.INPUT_2_1_1[1];
        tmp_0[10] = windEmulatorStep4_B.INPUT_2_1_1[2];
        tmp_0[11] = windEmulatorStep4_B.INPUT_2_1_1[3];
        tmp_1[3] = 12;
        tmp_0[12] = windEmulatorStep4_B.INPUT_4_1_1[0];
        tmp_0[13] = windEmulatorStep4_B.INPUT_4_1_1[1];
        tmp_0[14] = windEmulatorStep4_B.INPUT_4_1_1[2];
        tmp_0[15] = windEmulatorStep4_B.INPUT_4_1_1[3];
        tmp_1[4] = 16;
        tmp_0[16] = windEmulatorStep4_B.INPUT_3_1_1[0];
        tmp_0[17] = windEmulatorStep4_B.INPUT_3_1_1[1];
        tmp_0[18] = windEmulatorStep4_B.INPUT_3_1_1[2];
        tmp_0[19] = windEmulatorStep4_B.INPUT_3_1_1[3];
        tmp_1[5] = 20;
        simulationData->mData->mInputValues.mN = 20;
        simulationData->mData->mInputValues.mX = &tmp_0[0];
        simulationData->mData->mInputOffsets.mN = 6;
        simulationData->mData->mInputOffsets.mX = &tmp_1[0];
        simulationData->mData->mOutputs.mN = 39;
        simulationData->mData->mOutputs.mX = &windEmulatorStep4_B.STATE_1[0];
        simulationData->mData->mTolerances.mN = 0;
        simulationData->mData->mTolerances.mX = NULL;
        simulationData->mData->mCstateHasChanged = false;
        time_0 = windEmulatorStep4_M->Timing.t[1];
        simulationData->mData->mTime.mN = 1;
        simulationData->mData->mTime.mX = &time_0;
        isHit = 0;
        simulationData->mData->mSampleHits.mN = 1;
        simulationData->mData->mSampleHits.mX = &isHit;
        simulationData->mData->mIsFundamentalSampleHit = true;
        simulator = static_cast<NeslSimulator *>
          (windEmulatorStep4_DW.STATE_1_Simulator);
        diag = static_cast<NeuDiagnosticManager *>
          (windEmulatorStep4_DW.STATE_1_DiagMgr);
        diagTree = neu_diagnostic_manager_get_initial_tree(diag);
        rowIdx = ne_simulator_method(simulator, NESL_SIM_OUTPUTS, simulationData,
          diag);
        if (rowIdx != 0) {
          e_out = error_buffer_is_empty(rtmGetErrorStatus(windEmulatorStep4_M));
          if (e_out) {
            msg = rtw_diagnostics_msg(diagTree);
            rtmSetErrorStatus(windEmulatorStep4_M, msg);
          }
        }

        /* End of SimscapeExecutionBlock: '<S222>/STATE_1' */

        /* SimscapeExecutionBlock: '<S222>/OUTPUT_1_0' */
        simulationData = static_cast<NeslSimulationData *>
          (windEmulatorStep4_DW.OUTPUT_1_0_SimData);
        time_1 = windEmulatorStep4_M->Timing.t[0];
        simulationData->mData->mTime.mN = 1;
        simulationData->mData->mTime.mX = &time_1;
        simulationData->mData->mContStates.mN = 0;
        simulationData->mData->mContStates.mX = NULL;
        simulationData->mData->mDiscStates.mN = 0;
        simulationData->mData->mDiscStates.mX =
          &windEmulatorStep4_DW.OUTPUT_1_0_Discrete;
        simulationData->mData->mModeVector.mN = 0;
        simulationData->mData->mModeVector.mX =
          &windEmulatorStep4_DW.OUTPUT_1_0_Modes;
        e_out = false;
        simulationData->mData->mFoundZcEvents = e_out;
        simulationData->mData->mIsMajorTimeStep = true;
        e_out = false;
        simulationData->mData->mIsSolverAssertCheck = e_out;
        simulationData->mData->mIsSolverCheckingCIC = false;
        simulationData->mData->mIsComputingJacobian = false;
        simulationData->mData->mIsEvaluatingF0 = false;
        simulationData->mData->mIsSolverRequestingReset = false;
        simulationData->mData->mIsModeUpdateTimeStep = true;
        tmp_3[0] = 0;
        tmp_2[0] = windEmulatorStep4_B.INPUT_1_1_1[0];
        tmp_2[1] = windEmulatorStep4_B.INPUT_1_1_1[1];
        tmp_2[2] = windEmulatorStep4_B.INPUT_1_1_1[2];
        tmp_2[3] = windEmulatorStep4_B.INPUT_1_1_1[3];
        tmp_3[1] = 4;
        tmp_2[4] = windEmulatorStep4_B.INPUT_5_1_1[0];
        tmp_2[5] = windEmulatorStep4_B.INPUT_5_1_1[1];
        tmp_2[6] = windEmulatorStep4_B.INPUT_5_1_1[2];
        tmp_2[7] = windEmulatorStep4_B.INPUT_5_1_1[3];
        tmp_3[2] = 8;
        tmp_2[8] = windEmulatorStep4_B.INPUT_2_1_1[0];
        tmp_2[9] = windEmulatorStep4_B.INPUT_2_1_1[1];
        tmp_2[10] = windEmulatorStep4_B.INPUT_2_1_1[2];
        tmp_2[11] = windEmulatorStep4_B.INPUT_2_1_1[3];
        tmp_3[3] = 12;
        tmp_2[12] = windEmulatorStep4_B.INPUT_4_1_1[0];
        tmp_2[13] = windEmulatorStep4_B.INPUT_4_1_1[1];
        tmp_2[14] = windEmulatorStep4_B.INPUT_4_1_1[2];
        tmp_2[15] = windEmulatorStep4_B.INPUT_4_1_1[3];
        tmp_3[4] = 16;
        tmp_2[16] = windEmulatorStep4_B.INPUT_3_1_1[0];
        tmp_2[17] = windEmulatorStep4_B.INPUT_3_1_1[1];
        tmp_2[18] = windEmulatorStep4_B.INPUT_3_1_1[2];
        tmp_2[19] = windEmulatorStep4_B.INPUT_3_1_1[3];
        tmp_3[5] = 20;
        std::memcpy(&tmp_2[20], &windEmulatorStep4_B.STATE_1[0], 39U * sizeof
                    (real_T));
        tmp_3[6] = 59;
        simulationData->mData->mInputValues.mN = 59;
        simulationData->mData->mInputValues.mX = &tmp_2[0];
        simulationData->mData->mInputOffsets.mN = 7;
        simulationData->mData->mInputOffsets.mX = &tmp_3[0];
        simulationData->mData->mOutputs.mN = 11;
        simulationData->mData->mOutputs.mX = &windEmulatorStep4_B.OUTPUT_1_0[0];
        simulationData->mData->mTolerances.mN = 0;
        simulationData->mData->mTolerances.mX = NULL;
        simulationData->mData->mCstateHasChanged = false;
        time_2 = windEmulatorStep4_M->Timing.t[1];
        simulationData->mData->mTime.mN = 1;
        simulationData->mData->mTime.mX = &time_2;
        isHit_0 = 0;
        simulationData->mData->mSampleHits.mN = 1;
        simulationData->mData->mSampleHits.mX = &isHit_0;
        simulationData->mData->mIsFundamentalSampleHit = true;
        simulator = static_cast<NeslSimulator *>
          (windEmulatorStep4_DW.OUTPUT_1_0_Simulator);
        diag = static_cast<NeuDiagnosticManager *>
          (windEmulatorStep4_DW.OUTPUT_1_0_DiagMgr);
        diagTree = neu_diagnostic_manager_get_initial_tree(diag);
        rowIdx = ne_simulator_method(simulator, NESL_SIM_OUTPUTS, simulationData,
          diag);
        if (rowIdx != 0) {
          e_out = error_buffer_is_empty(rtmGetErrorStatus(windEmulatorStep4_M));
          if (e_out) {
            msg = rtw_diagnostics_msg(diagTree);
            rtmSetErrorStatus(windEmulatorStep4_M, msg);
          }
        }

        /* End of SimscapeExecutionBlock: '<S222>/OUTPUT_1_0' */

        /* Product: '<S59>/Product' */
        windEmulatorStep4_B.Product_mo = windEmulatorStep4_B.kDampingNow *
          windEmulatorStep4_B.OUTPUT_1_0[7];

        /* Gain: '<S59>/kSpringNow' incorporates:
         *  Constant: '<S59>/kSpring'
         */
        windEmulatorStep4_B.kSpringNow = windEmulatorStep4_cal->kSpringNow_Gain *
          windEmulatorStep4_cal->kSpring_Value;

        /* Product: '<S59>/Product1' */
        windEmulatorStep4_B.Product1_g = windEmulatorStep4_B.OUTPUT_1_0[6] *
          windEmulatorStep4_B.kSpringNow;

        /* Sum: '<S59>/Add' */
        windEmulatorStep4_B.ForceD = windEmulatorStep4_B.Product_mo +
          windEmulatorStep4_B.Product1_g;

        /* Product: '<S55>/Product2' incorporates:
         *  Constant: '<S55>/Constant2'
         */
        riseValLimit = *get_A() * 2.0 * 3.1415926535897931 / *get_T() /
          *get_targetOmega();

        /* Product: '<S55>/Product2' */
        windEmulatorStep4_B.Product2 = windEmulatorStep4_B.ForceD * riseValLimit;

        /* Gain: '<S55>/Gain' */
        windEmulatorStep4_B.Gain_n = windEmulatorStep4_cal->Gain_Gain_c *
          windEmulatorStep4_B.Product2;

        /* Memory: '<S176>/Memory' */
        windEmulatorStep4_B.Memory_o =
          windEmulatorStep4_DW.Memory_PreviousInput_g;
      }

      /* Product: '<S53>/Product' */
      windEmulatorStep4_B.TorqueInputRef = windEmulatorStep4_B.Saturation *
        windEmulatorStep4_B.Gain_n;

      /* Abs: '<S53>/Abs' */
      windEmulatorStep4_B.Abs = std::abs(windEmulatorStep4_B.TorqueInputRef);

      /* Sum: '<S123>/Add' */
      windEmulatorStep4_B.Add_m =
        windEmulatorStep4_B.BusAssignment_c.genSpeedActual -
        windEmulatorStep4_B.BusAssignment_c.speedRef_rpm;

      /* Product: '<S174>/Product' incorporates:
       *  Constant: '<S174>/Constant1'
       */
      riseValLimit = -tmp_9->PG;

      /* Product: '<S174>/Product' */
      windEmulatorStep4_B.ControlSignal31_d = riseValLimit *
        windEmulatorStep4_B.Add_m;

      /* RelationalOperator: '<S174>/Relational Operator' */
      windEmulatorStep4_B.RelationalOperator_e =
        (windEmulatorStep4_B.BusAssignment_c.speedRef_rpm <=
         windEmulatorStep4_B.BusAssignment_c.genSpeedActual);

      /* CombinatorialLogic: '<S176>/Logic' incorporates:
       *  Constant: '<S174>/Constant'
       */
      e_out = windEmulatorStep4_B.RelationalOperator_e;
      rowIdx = e_out;
      e_out = windEmulatorStep4_cal->Constant_Value_c;
      rowIdx = static_cast<int32_T>((static_cast<uint32_T>(rowIdx) << 1) + e_out);
      e_out = windEmulatorStep4_B.Memory_o;
      rowIdx = static_cast<int32_T>((static_cast<uint32_T>(rowIdx) << 1) + e_out);
      windEmulatorStep4_B.Logic_c[0U] = windEmulatorStep4_cal->Logic_table_h[
        static_cast<uint32_T>(rowIdx)];
      windEmulatorStep4_B.Logic_c[1U] = windEmulatorStep4_cal->Logic_table_h[
        static_cast<uint32_T>(rowIdx) + 8U];
      if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
        /* RateLimiter: '<S120>/Rate Limiter1' incorporates:
         *  Constant: '<S120>/shaftSpeedRefMin'
         */
        rateLimiterRate = windEmulatorStep4_cal->shaftSpeedRefMin_Value_k -
          windEmulatorStep4_DW.PrevY_fq;
        if (rateLimiterRate > windEmulatorStep4_cal->RateLimiter1_RisingLim_o *
            windEmulatorStep4_period) {
          /* RateLimiter: '<S120>/Rate Limiter1' */
          windEmulatorStep4_B.RateLimiter1_m =
            windEmulatorStep4_cal->RateLimiter1_RisingLim_o *
            windEmulatorStep4_period + windEmulatorStep4_DW.PrevY_fq;
        } else if (rateLimiterRate <
                   windEmulatorStep4_cal->RateLimiter1_FallingLim_e *
                   windEmulatorStep4_period) {
          /* RateLimiter: '<S120>/Rate Limiter1' */
          windEmulatorStep4_B.RateLimiter1_m =
            windEmulatorStep4_cal->RateLimiter1_FallingLim_e *
            windEmulatorStep4_period + windEmulatorStep4_DW.PrevY_fq;
        } else {
          /* RateLimiter: '<S120>/Rate Limiter1' */
          windEmulatorStep4_B.RateLimiter1_m =
            windEmulatorStep4_cal->shaftSpeedRefMin_Value_k;
        }

        windEmulatorStep4_DW.PrevY_fq = windEmulatorStep4_B.RateLimiter1_m;

        /* End of RateLimiter: '<S120>/Rate Limiter1' */

        /* Memory: '<S177>/Memory' */
        windEmulatorStep4_B.Memory_a =
          windEmulatorStep4_DW.Memory_PreviousInput_n;
      }

      /* Sum: '<S123>/Add1' */
      windEmulatorStep4_B.Add1_f =
        windEmulatorStep4_B.BusAssignment_c.genSpeedActual -
        windEmulatorStep4_B.RateLimiter1_m;

      /* Product: '<S175>/Product' incorporates:
       *  Constant: '<S175>/Constant1'
       */
      riseValLimit = -tmp_8;

      /* Product: '<S175>/Product' */
      windEmulatorStep4_B.ControlSignal31_m = riseValLimit *
        windEmulatorStep4_B.Add1_f;

      /* RelationalOperator: '<S175>/Relational Operator' */
      windEmulatorStep4_B.RelationalOperator_g =
        (windEmulatorStep4_B.RateLimiter1_m >=
         windEmulatorStep4_B.BusAssignment_c.genSpeedActual);

      /* CombinatorialLogic: '<S177>/Logic' incorporates:
       *  Constant: '<S175>/Constant'
       */
      e_out = windEmulatorStep4_B.RelationalOperator_g;
      rowIdx = e_out;
      e_out = windEmulatorStep4_cal->Constant_Value_ks;
      rowIdx = static_cast<int32_T>((static_cast<uint32_T>(rowIdx) << 1) + e_out);
      e_out = windEmulatorStep4_B.Memory_a;
      rowIdx = static_cast<int32_T>((static_cast<uint32_T>(rowIdx) << 1) + e_out);
      windEmulatorStep4_B.Logic_p[0U] = windEmulatorStep4_cal->Logic_table_n[
        static_cast<uint32_T>(rowIdx)];
      windEmulatorStep4_B.Logic_p[1U] = windEmulatorStep4_cal->Logic_table_n[
        static_cast<uint32_T>(rowIdx) + 8U];

      /* Switch: '<S123>/Switch' incorporates:
       *  Switch: '<S123>/Switch1'
       *  Switch: '<S175>/Switch'
       */
      if (windEmulatorStep4_B.Add_m > windEmulatorStep4_cal->Switch_Threshold_k)
      {
        /* Switch: '<S174>/Switch' */
        if (windEmulatorStep4_B.Logic_c[0]) {
          /* Saturate: '<S174>/Saturation' */
          riseValLimit = windEmulatorStep4_B.ControlSignal31_d;
          u1 = windEmulatorStep4_cal->Saturation_LowerSat_f;
          rateLimiterRate = windEmulatorStep4_cal->Saturation_UpperSat_f;
          if (riseValLimit > rateLimiterRate) {
            /* Saturate: '<S174>/Saturation' */
            windEmulatorStep4_B.Saturation_af = rateLimiterRate;
          } else if (riseValLimit < u1) {
            /* Saturate: '<S174>/Saturation' */
            windEmulatorStep4_B.Saturation_af = u1;
          } else {
            /* Saturate: '<S174>/Saturation' */
            windEmulatorStep4_B.Saturation_af = riseValLimit;
          }

          /* End of Saturate: '<S174>/Saturation' */

          /* Switch: '<S174>/Switch' */
          windEmulatorStep4_B.ControlSignal3_h =
            windEmulatorStep4_B.Saturation_af;
        } else {
          /* Switch: '<S174>/Switch' */
          windEmulatorStep4_B.ControlSignal3_h =
            windEmulatorStep4_B.ControlSignal31_d;
        }

        /* End of Switch: '<S174>/Switch' */

        /* Switch: '<S123>/Switch' */
        windEmulatorStep4_B.Switch_n = windEmulatorStep4_B.ControlSignal3_h;
      } else {
        if (windEmulatorStep4_B.Add1_f >
            windEmulatorStep4_cal->Switch1_Threshold_k) {
          /* Switch: '<S123>/Switch1' incorporates:
           *  Constant: '<S123>/Constant1'
           */
          windEmulatorStep4_B.Switch1_g =
            windEmulatorStep4_cal->Constant1_Value_j;
        } else {
          if (windEmulatorStep4_B.Logic_p[0]) {
            /* Saturate: '<S175>/Saturation' incorporates:
             *  Switch: '<S123>/Switch1'
             *  Switch: '<S175>/Switch'
             */
            riseValLimit = windEmulatorStep4_B.ControlSignal31_m;
            u1 = windEmulatorStep4_cal->Saturation_LowerSat_h;
            rateLimiterRate = windEmulatorStep4_cal->Saturation_UpperSat_m;
            if (riseValLimit > rateLimiterRate) {
              /* Saturate: '<S175>/Saturation' */
              windEmulatorStep4_B.Saturation_p = rateLimiterRate;
            } else if (riseValLimit < u1) {
              /* Saturate: '<S175>/Saturation' */
              windEmulatorStep4_B.Saturation_p = u1;
            } else {
              /* Saturate: '<S175>/Saturation' */
              windEmulatorStep4_B.Saturation_p = riseValLimit;
            }

            /* End of Saturate: '<S175>/Saturation' */

            /* Switch: '<S175>/Switch' incorporates:
             *  Switch: '<S123>/Switch1'
             */
            windEmulatorStep4_B.ControlSignal3 =
              windEmulatorStep4_B.Saturation_p;
          } else {
            /* Switch: '<S175>/Switch' incorporates:
             *  Switch: '<S123>/Switch1'
             */
            windEmulatorStep4_B.ControlSignal3 =
              windEmulatorStep4_B.ControlSignal31_m;
          }

          /* Switch: '<S123>/Switch1' */
          windEmulatorStep4_B.Switch1_g = windEmulatorStep4_B.ControlSignal3;
        }

        /* Switch: '<S123>/Switch' incorporates:
         *  Switch: '<S123>/Switch1'
         *  Switch: '<S175>/Switch'
         */
        windEmulatorStep4_B.Switch_n = windEmulatorStep4_B.Switch1_g;
      }

      /* End of Switch: '<S123>/Switch' */

      /* Gain: '<S123>/Gain' */
      windEmulatorStep4_B.Gain_f = tmp_m * windEmulatorStep4_B.Switch_n;

      /* RateLimiter: '<S120>/Rate Limiter' */
      if (windEmulatorStep4_DW.LastMajorTime_n == (rtInf)) {
        /* RateLimiter: '<S120>/Rate Limiter' */
        windEmulatorStep4_B.RateLimiter_aw = windEmulatorStep4_B.Gain_f;
      } else {
        u1 = windEmulatorStep4_M->Timing.t[0] -
          windEmulatorStep4_DW.LastMajorTime_n;
        riseValLimit = u1 * tmp_7;
        rateLimiterRate = windEmulatorStep4_B.Gain_f -
          windEmulatorStep4_DW.PrevY_e;
        if (rateLimiterRate > riseValLimit) {
          /* RateLimiter: '<S120>/Rate Limiter' */
          windEmulatorStep4_B.RateLimiter_aw = windEmulatorStep4_DW.PrevY_e +
            riseValLimit;
        } else {
          riseValLimit = -tmp_7;
          u1 *= riseValLimit;
          if (rateLimiterRate < u1) {
            /* RateLimiter: '<S120>/Rate Limiter' */
            windEmulatorStep4_B.RateLimiter_aw = windEmulatorStep4_DW.PrevY_e +
              u1;
          } else {
            /* RateLimiter: '<S120>/Rate Limiter' */
            windEmulatorStep4_B.RateLimiter_aw = windEmulatorStep4_B.Gain_f;
          }
        }
      }

      /* Sum: '<S122>/Sum' */
      windEmulatorStep4_B.wError_c =
        windEmulatorStep4_B.BusAssignment_c.genSpeedActual -
        windEmulatorStep4_B.BusAssignment_c.speedRef_rpm;
      if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
        /* Gain: '<S161>/Proportional Gain' */
        windEmulatorStep4_B.ProportionalGain_h = tmp_9->PG *
          windEmulatorStep4_B.wError_c;

        /* DiscreteIntegrator: '<S156>/Integrator' */
        if (windEmulatorStep4_B.BusAssignment_c.speedCtrlReset ||
            (windEmulatorStep4_DW.Integrator_PrevResetState_g != 0)) {
          windEmulatorStep4_DW.Integrator_DSTATE_e =
            windEmulatorStep4_cal->PIDController_InitialConditio_c;
        }

        /* DiscreteIntegrator: '<S156>/Integrator' */
        windEmulatorStep4_B.Integrator_l =
          windEmulatorStep4_DW.Integrator_DSTATE_e;

        /* Gain: '<S150>/Derivative Gain' */
        windEmulatorStep4_B.DerivativeGain_g =
          windEmulatorStep4_cal->PIDController_D_d *
          windEmulatorStep4_B.wError_c;

        /* DiscreteIntegrator: '<S151>/Filter' */
        if (windEmulatorStep4_B.BusAssignment_c.speedCtrlReset ||
            (windEmulatorStep4_DW.Filter_PrevResetState_g != 0)) {
          windEmulatorStep4_DW.Filter_DSTATE_b =
            windEmulatorStep4_cal->PIDController_InitialConditio_k;
        }

        /* DiscreteIntegrator: '<S151>/Filter' */
        windEmulatorStep4_B.Filter_j = windEmulatorStep4_DW.Filter_DSTATE_b;

        /* Sum: '<S151>/SumD' */
        windEmulatorStep4_B.SumD_c = windEmulatorStep4_B.DerivativeGain_g -
          windEmulatorStep4_B.Filter_j;

        /* Gain: '<S159>/Filter Coefficient' */
        windEmulatorStep4_B.FilterCoefficient_g =
          windEmulatorStep4_cal->PIDController_N_p * windEmulatorStep4_B.SumD_c;

        /* Sum: '<S166>/Sum' */
        windEmulatorStep4_B.Sum_n = (windEmulatorStep4_B.ProportionalGain_h +
          windEmulatorStep4_B.Integrator_l) +
          windEmulatorStep4_B.FilterCoefficient_g;

        /* RelationalOperator: '<S164>/LowerRelop1' incorporates:
         *  Constant: '<S122>/Constant'
         */
        windEmulatorStep4_B.LowerRelop1_h = (windEmulatorStep4_B.Sum_n > tmp_6);

        /* RelationalOperator: '<S164>/UpperRelop' incorporates:
         *  Constant: '<S122>/Constant1'
         */
        riseValLimit = -tmp_6;

        /* RelationalOperator: '<S164>/UpperRelop' */
        windEmulatorStep4_B.UpperRelop_m = (windEmulatorStep4_B.Sum_n <
          riseValLimit);

        /* Switch: '<S164>/Switch' */
        if (windEmulatorStep4_B.UpperRelop_m) {
          /* Switch: '<S164>/Switch' incorporates:
           *  Constant: '<S122>/Constant1'
           */
          windEmulatorStep4_B.Switch_c = -tmp_6;
        } else {
          /* Switch: '<S164>/Switch' */
          windEmulatorStep4_B.Switch_c = windEmulatorStep4_B.Sum_n;
        }

        /* Switch: '<S164>/Switch2' */
        if (windEmulatorStep4_B.LowerRelop1_h) {
          /* Switch: '<S164>/Switch2' incorporates:
           *  Constant: '<S122>/Constant'
           */
          windEmulatorStep4_B.Switch2_m = tmp_6;
        } else {
          /* Switch: '<S164>/Switch2' */
          windEmulatorStep4_B.Switch2_m = windEmulatorStep4_B.Switch_c;
        }

        /* End of Switch: '<S164>/Switch2' */

        /* Gain: '<S122>/Gain2' */
        windEmulatorStep4_B.ContolTorque_f = windEmulatorStep4_cal->Gain2_Gain_e
          * windEmulatorStep4_B.Switch2_m;
      }

      /* Switch generated from: '<S53>/Switch' incorporates:
       *  Constant: '<S120>/DeadBandController'
       *  Switch: '<S120>/Switch'
       */
      if (windEmulatorStep4_B.Abs >= tmp_h) {
        /* Switch: '<S61>/Switch' incorporates:
         *  Constant: '<S61>/DeadBandController'
         */
        if (tmp_g) {
          /* Switch: '<S61>/Switch' */
          windEmulatorStep4_B.Switch_cn = windEmulatorStep4_B.RateLimiter_b;
        } else {
          /* Switch: '<S61>/Switch' */
          windEmulatorStep4_B.Switch_cn = windEmulatorStep4_B.ContolTorque;
        }

        /* Switch generated from: '<S53>/Switch' */
        windEmulatorStep4_B.ControlTorqueLoad = windEmulatorStep4_B.Switch_cn;
      } else {
        if (tmp_g) {
          /* Switch: '<S120>/Switch' */
          windEmulatorStep4_B.Switch_f = windEmulatorStep4_B.RateLimiter_aw;
        } else {
          /* Switch: '<S120>/Switch' */
          windEmulatorStep4_B.Switch_f = windEmulatorStep4_B.ContolTorque_f;
        }

        /* Switch generated from: '<S53>/Switch' */
        windEmulatorStep4_B.ControlTorqueLoad = windEmulatorStep4_B.Switch_f;
      }

      if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
        /* Gain: '<S184>/Gain' */
        windEmulatorStep4_B.Pressure = windEmulatorStep4_cal->Gain_Gain_cc *
          windEmulatorStep4_B.OUTPUT_1_0[8];

        /* Gain: '<S6>/psi -> bar' */
        windEmulatorStep4_B.psibar = *get_psi2bar() *
          windEmulatorStep4_B.Pressure;

        /* Gain: '<S191>/Gain' */
        windEmulatorStep4_B.ShaftSpeedPump = tmp_k *
          windEmulatorStep4_B.OUTPUT_1_0[9];
      }

      /* Sum: '<S121>/Sum' */
      windEmulatorStep4_B.Sum_f = windEmulatorStep4_B.OUTPUT_1_0[10] -
        windEmulatorStep4_B.TorqueInputRef;
      if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
        /* DiscreteIntegrator: '<S121>/Discrete-Time Integrator' */
        windEmulatorStep4_B.DiscreteTimeIntegrator =
          windEmulatorStep4_DW.DiscreteTimeIntegrator_DSTATE;

        /* Gain: '<S121>/Gain1' */
        windEmulatorStep4_B.Gain1_b = tmp_f->IG *
          windEmulatorStep4_B.DiscreteTimeIntegrator;
      }

      /* Switch generated from: '<S53>/Switch' */
      if (windEmulatorStep4_B.Abs >= tmp_h) {
        /* Switch: '<S62>/Switch' */
        if (windEmulatorStep4_B.TorqueInputRef >=
            windEmulatorStep4_cal->Switch_Threshold) {
          /* Switch: '<S62>/Switch' incorporates:
           *  Constant: '<S62>/Constant'
           */
          windEmulatorStep4_B.Switch_ne = windEmulatorStep4_cal->Constant_Value;
        } else {
          /* Switch: '<S62>/Switch' incorporates:
           *  Constant: '<S62>/Constant1'
           */
          windEmulatorStep4_B.Switch_ne =
            windEmulatorStep4_cal->Constant1_Value_b;
        }

        /* End of Switch: '<S62>/Switch' */

        /* Switch generated from: '<S53>/Switch' */
        windEmulatorStep4_B.ControlSignal1 = windEmulatorStep4_B.Switch_ne;
      } else {
        /* Gain: '<S121>/Gain2' */
        windEmulatorStep4_B.Gain2_n = windEmulatorStep4_cal->Gain2_Gain *
          windEmulatorStep4_B.TorqueInputRef;

        /* Product: '<S121>/Divide' */
        windEmulatorStep4_B.Divide = windEmulatorStep4_B.Gain2_n /
          windEmulatorStep4_B.Pressure;

        /* Product: '<S121>/Product' incorporates:
         *  Constant: '<S121>/UnitsConversion'
         */
        riseValLimit = 6.283185307179586E+6 / (tmp_5 * 6894.75);

        /* Product: '<S121>/Product' */
        windEmulatorStep4_B.Product_aq = windEmulatorStep4_B.Divide *
          riseValLimit;

        /* Gain: '<S121>/Gain' */
        windEmulatorStep4_B.Gain_a = tmp_f->PG * windEmulatorStep4_B.Sum_f;

        /* Sum: '<S121>/Add' */
        windEmulatorStep4_B.Add_iu = windEmulatorStep4_B.Gain_a +
          windEmulatorStep4_B.Gain1_b;

        /* Sum: '<S121>/Add1' */
        windEmulatorStep4_B.Add1_l = windEmulatorStep4_B.Add_iu +
          windEmulatorStep4_B.Product_aq;

        /* Saturate: '<S121>/Saturation' */
        riseValLimit = windEmulatorStep4_B.Add1_l;
        u1 = windEmulatorStep4_cal->Saturation_LowerSat_b;
        rateLimiterRate = windEmulatorStep4_cal->Saturation_UpperSat_mj;
        if (riseValLimit > rateLimiterRate) {
          /* Saturate: '<S121>/Saturation' */
          windEmulatorStep4_B.Saturation_j = rateLimiterRate;
        } else if (riseValLimit < u1) {
          /* Saturate: '<S121>/Saturation' */
          windEmulatorStep4_B.Saturation_j = u1;
        } else {
          /* Saturate: '<S121>/Saturation' */
          windEmulatorStep4_B.Saturation_j = riseValLimit;
        }

        /* End of Saturate: '<S121>/Saturation' */

        /* Switch generated from: '<S53>/Switch' */
        windEmulatorStep4_B.ControlSignal1 = windEmulatorStep4_B.Saturation_j;
      }

      /* Abs: '<S53>/Abs2' */
      windEmulatorStep4_B.Abs2 = std::abs(windEmulatorStep4_B.TorqueInputRef);

      /* Switch: '<S53>/Switch2' */
      if (windEmulatorStep4_B.Abs2 >= tmp_h) {
        /* Gain: '<S53>/Gain' */
        riseValLimit = 1.0 / (tmp_5 * 1.0E-6 * 0.15915494309189535) / 6894.75;

        /* Gain: '<S53>/Gain' */
        windEmulatorStep4_B.Gain_k = riseValLimit *
          windEmulatorStep4_B.TorqueInputRef;

        /* Abs: '<S53>/Abs1' */
        windEmulatorStep4_B.Abs1 = std::abs(windEmulatorStep4_B.Gain_k);

        /* Switch: '<S53>/Switch2' */
        windEmulatorStep4_B.PressureRef = windEmulatorStep4_B.Abs1;
      } else {
        /* Switch: '<S53>/Switch2' incorporates:
         *  Constant: '<S52>/Constant1'
         */
        windEmulatorStep4_B.PressureRef = *get_minPressureRef_psi();
      }

      /* Sum: '<S60>/Sum' */
      windEmulatorStep4_B.Sum_me = windEmulatorStep4_B.Pressure -
        windEmulatorStep4_B.PressureRef;
      if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
        /* DiscreteIntegrator: '<S60>/Discrete-Time Integrator' */
        windEmulatorStep4_B.DiscreteTimeIntegrator_i =
          windEmulatorStep4_DW.DiscreteTimeIntegrator_DSTATE_n;

        /* Gain: '<S60>/Gain1' */
        windEmulatorStep4_B.Gain1_j = tmp_e->IG *
          windEmulatorStep4_B.DiscreteTimeIntegrator_i;

        /* DiscreteIntegrator: '<S119>/Discrete-Time Integrator' */
        windEmulatorStep4_B.DiscreteTimeIntegrator_e =
          windEmulatorStep4_DW.DiscreteTimeIntegrator_DSTATE_l;

        /* Gain: '<S119>/Gain1' */
        windEmulatorStep4_B.Gain1_h = tmp_e->IG *
          windEmulatorStep4_B.DiscreteTimeIntegrator_e;
      }

      /* Sum: '<S119>/Sum' */
      windEmulatorStep4_B.Sum_d = windEmulatorStep4_B.Pressure -
        windEmulatorStep4_B.PressureRef;

      /* Switch generated from: '<S53>/Switch' */
      if (windEmulatorStep4_B.Abs >= tmp_h) {
        /* Gain: '<S60>/Gain' */
        windEmulatorStep4_B.Gain_ge = tmp_e->PG * windEmulatorStep4_B.Sum_me;

        /* Sum: '<S60>/Add' */
        windEmulatorStep4_B.Add_i = windEmulatorStep4_B.Gain_ge +
          windEmulatorStep4_B.Gain1_j;

        /* Saturate: '<S60>/Saturation' */
        riseValLimit = windEmulatorStep4_B.Add_i;
        u1 = windEmulatorStep4_cal->Saturation_LowerSat_g;
        rateLimiterRate = windEmulatorStep4_cal->Saturation_UpperSat_d;
        if (riseValLimit > rateLimiterRate) {
          /* Saturate: '<S60>/Saturation' */
          windEmulatorStep4_B.Saturation_a = rateLimiterRate;
        } else if (riseValLimit < u1) {
          /* Saturate: '<S60>/Saturation' */
          windEmulatorStep4_B.Saturation_a = u1;
        } else {
          /* Saturate: '<S60>/Saturation' */
          windEmulatorStep4_B.Saturation_a = riseValLimit;
        }

        /* End of Saturate: '<S60>/Saturation' */

        /* Switch generated from: '<S53>/Switch' */
        windEmulatorStep4_B.ControlSignal2 = windEmulatorStep4_B.Saturation_a;
      } else {
        /* Gain: '<S119>/Gain' */
        windEmulatorStep4_B.Gain_j = tmp_e->PG * windEmulatorStep4_B.Sum_d;

        /* Sum: '<S119>/Add' */
        windEmulatorStep4_B.Add_k = windEmulatorStep4_B.Gain_j +
          windEmulatorStep4_B.Gain1_h;

        /* Saturate: '<S119>/Saturation' */
        riseValLimit = windEmulatorStep4_B.Add_k;
        u1 = windEmulatorStep4_cal->Saturation_LowerSat_e;
        rateLimiterRate = windEmulatorStep4_cal->Saturation_UpperSat_p;
        if (riseValLimit > rateLimiterRate) {
          /* Saturate: '<S119>/Saturation' */
          windEmulatorStep4_B.Saturation_f = rateLimiterRate;
        } else if (riseValLimit < u1) {
          /* Saturate: '<S119>/Saturation' */
          windEmulatorStep4_B.Saturation_f = u1;
        } else {
          /* Saturate: '<S119>/Saturation' */
          windEmulatorStep4_B.Saturation_f = riseValLimit;
        }

        /* End of Saturate: '<S119>/Saturation' */

        /* Switch generated from: '<S53>/Switch' */
        windEmulatorStep4_B.ControlSignal2 = windEmulatorStep4_B.Saturation_f;
      }

      if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
        /* Gain: '<S181>/Gain' */
        windEmulatorStep4_B.FlowMotor1 = windEmulatorStep4_cal->Gain_Gain_f *
          windEmulatorStep4_B.OUTPUT_1_0[2];
      }

      /* BusAssignment: '<S6>/Bus Assignment' incorporates:
       *  Constant: '<S6>/Constant'
       */
      windEmulatorStep4_B.BusAssignment_n = *get_hptoSignalStruct();

      /* BusAssignment: '<S6>/Bus Assignment' */
      windEmulatorStep4_B.BusAssignment_n.genTorqueCmd_Nm =
        windEmulatorStep4_B.ControlTorqueLoad;
      windEmulatorStep4_B.BusAssignment_n.pressure_bar =
        windEmulatorStep4_B.psibar;
      windEmulatorStep4_B.BusAssignment_n.hmOutputShafTorque_Nm =
        windEmulatorStep4_B.OUTPUT_1_0[3];
      windEmulatorStep4_B.BusAssignment_n.genShaftSpeed_rpm =
        windEmulatorStep4_B.BusAssignment_c.genSpeedActual;
      windEmulatorStep4_B.BusAssignment_n.excShaftSpeed_rpm =
        windEmulatorStep4_B.ShaftSpeedPump;
      windEmulatorStep4_B.BusAssignment_n.excShaftTorque_Nm =
        windEmulatorStep4_B.OUTPUT_1_0[10];
      windEmulatorStep4_B.BusAssignment_n.ctrlSignal1 =
        windEmulatorStep4_B.ControlSignal1;
      windEmulatorStep4_B.BusAssignment_n.ctrlSignal2 =
        windEmulatorStep4_B.ControlSignal2;
      windEmulatorStep4_B.BusAssignment_n.genPumpFlow_lpm =
        windEmulatorStep4_B.FlowMotor1;

      /* MultiPortSwitch: '<S3>/Multiport Switch' incorporates:
       *  DataTypeConversion: '<S3>/toExpTypeEnum'
       */
      switch (windEmulatorStep4_B.toExpTypeEnum) {
       case expTypeEnum_off:
        /* MultiPortSwitch: '<S3>/Multiport Switch' incorporates:
         *  Constant: '<S3>/zeroTorque'
         */
        windEmulatorStep4_B.MultiportSwitch_p =
          windEmulatorStep4_cal->zeroTorque_Value;
        break;

       case expTypeEnum_sid:
        /* MultiPortSwitch: '<S3>/Multiport Switch' */
        windEmulatorStep4_B.MultiportSwitch_p =
          windEmulatorStep4_B.BusAssignment_f.acs880Torque_Nm;
        break;

       case expTypeEnum_hil:
        /* MultiPortSwitch: '<S3>/Multiport Switch' */
        windEmulatorStep4_B.MultiportSwitch_p =
          windEmulatorStep4_B.BusAssignment_n.hmOutputShafTorque_Nm;
        break;

       default:
        /* MultiPortSwitch: '<S3>/Multiport Switch' incorporates:
         *  Constant: '<S3>/zeroTorque'
         */
        windEmulatorStep4_B.MultiportSwitch_p =
          windEmulatorStep4_cal->zeroTorque_Value;
        break;
      }

      /* End of MultiPortSwitch: '<S3>/Multiport Switch' */

      /* RateLimiter: '<S2>/acs880RateLim' */
      if (windEmulatorStep4_DW.LastMajorTime_d == (rtInf)) {
        /* RateLimiter: '<S2>/acs880RateLim' */
        windEmulatorStep4_B.acs880RateLim =
          windEmulatorStep4_B.MultiportSwitch_p;
      } else {
        u1 = windEmulatorStep4_M->Timing.t[0] -
          windEmulatorStep4_DW.LastMajorTime_d;
        riseValLimit = u1 * *get_acs880RateLimRising();
        rateLimiterRate = windEmulatorStep4_B.MultiportSwitch_p -
          windEmulatorStep4_DW.PrevY_b;
        if (rateLimiterRate > riseValLimit) {
          /* RateLimiter: '<S2>/acs880RateLim' */
          windEmulatorStep4_B.acs880RateLim = windEmulatorStep4_DW.PrevY_b +
            riseValLimit;
        } else {
          u1 *= *get_acs880RateLimFalling();
          if (rateLimiterRate < u1) {
            /* RateLimiter: '<S2>/acs880RateLim' */
            windEmulatorStep4_B.acs880RateLim = windEmulatorStep4_DW.PrevY_b +
              u1;
          } else {
            /* RateLimiter: '<S2>/acs880RateLim' */
            windEmulatorStep4_B.acs880RateLim =
              windEmulatorStep4_B.MultiportSwitch_p;
          }
        }
      }

      /* End of RateLimiter: '<S2>/acs880RateLim' */

      /* Saturate: '<S2>/Saturation' */
      riseValLimit = windEmulatorStep4_B.acs880RateLim;
      u1 = *get_acs880SetpointLimLower();
      rateLimiterRate = *get_acs880SetpointLimUpper();
      if (riseValLimit > rateLimiterRate) {
        /* Saturate: '<S2>/Saturation' */
        windEmulatorStep4_B.ACS880Setpoint = rateLimiterRate;
      } else if (riseValLimit < u1) {
        /* Saturate: '<S2>/Saturation' */
        windEmulatorStep4_B.ACS880Setpoint = u1;
      } else {
        /* Saturate: '<S2>/Saturation' */
        windEmulatorStep4_B.ACS880Setpoint = riseValLimit;
      }

      /* End of Saturate: '<S2>/Saturation' */

      /* Gain: '<S2>/Nm -> %' */
      riseValLimit = 1.0 / tmp_4 * 100.0;

      /* Gain: '<S2>/Nm -> %' */
      windEmulatorStep4_B.Nm = riseValLimit * windEmulatorStep4_B.ACS880Setpoint;

      /* BusAssignment: '<S2>/Bus Assignment' incorporates:
       *  Constant: '<S2>/Constant'
       */
      windEmulatorStep4_B.BusAssignment_k = *get_acs880CtrlStruct();

      /* BusAssignment: '<S2>/Bus Assignment' */
      windEmulatorStep4_B.BusAssignment_k.ctrlWord =
        windEmulatorStep4_B.ControlWord;
      windEmulatorStep4_B.BusAssignment_k.state =
        windEmulatorStep4_B.CastToDouble1_a;
      windEmulatorStep4_B.BusAssignment_k.torqueSetpoint_Nm =
        windEmulatorStep4_B.ACS880Setpoint;
      windEmulatorStep4_B.BusAssignment_k.torqueSetpoint_percent =
        windEmulatorStep4_B.Nm;
      if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
        /* ToAsyncQueueBlock generated from: '<S25>/acs880CtrlSignals' */
        slrtLogSignal
          (windEmulatorStep4_DW.TAQSigLogging_InsertedFor_acs_a.SLRTSigHandles,
           (((windEmulatorStep4_M->Timing.clockTick1+
              windEmulatorStep4_M->Timing.clockTickH1* 4294967296.0)) * 0.004));
      }

      /* MATLAB Function: '<S46>/parseCtrlWord' */
      windEmulatorStep4_parseCtrlWord
        (windEmulatorStep4_B.BusAssignment_k.ctrlWord,
         &windEmulatorStep4_B.sf_parseCtrlWord_h);

      /* Logic: '<S46>/enableOperation' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_B.enableOperation =
        ((windEmulatorStep4_B.sf_parseCtrlWord_h.enableOperation != 0) &&
         windEmulatorStep4_cal->Constant_Value_e);

      /* Logic: '<S46>/extCtrlLoc' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_B.extCtrlLoc_f =
        ((windEmulatorStep4_B.sf_parseCtrlWord_h.extCtrlLoc != 0) &&
         windEmulatorStep4_cal->Constant_Value_e);

      /* Logic: '<S46>/inching1' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_B.inching1 =
        ((windEmulatorStep4_B.sf_parseCtrlWord_h.inching1 != 0) &&
         windEmulatorStep4_cal->Constant_Value_e);

      /* Logic: '<S46>/inching2' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_B.inching2 =
        ((windEmulatorStep4_B.sf_parseCtrlWord_h.inching2 != 0) &&
         windEmulatorStep4_cal->Constant_Value_e);

      /* Logic: '<S46>/off1Ctrl' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_B.off1Ctrl =
        ((windEmulatorStep4_B.sf_parseCtrlWord_h.off1Ctrl != 0) &&
         windEmulatorStep4_cal->Constant_Value_e);

      /* Logic: '<S46>/off2Ctrl' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_B.off2Ctrl =
        ((windEmulatorStep4_B.sf_parseCtrlWord_h.off2Ctrl != 0) &&
         windEmulatorStep4_cal->Constant_Value_e);

      /* Logic: '<S46>/off3Ctrl' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_B.off3Ctrl =
        ((windEmulatorStep4_B.sf_parseCtrlWord_h.off3Ctrl != 0) &&
         windEmulatorStep4_cal->Constant_Value_e);

      /* Logic: '<S46>/rampHold' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_B.rampHold =
        ((windEmulatorStep4_B.sf_parseCtrlWord_h.rampHold != 0) &&
         windEmulatorStep4_cal->Constant_Value_e);

      /* Logic: '<S46>/rampInZero' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_B.rampInZero =
        ((windEmulatorStep4_B.sf_parseCtrlWord_h.rampInZero != 0) &&
         windEmulatorStep4_cal->Constant_Value_e);

      /* Logic: '<S46>/rampOutZero' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_B.rampOutZero =
        ((windEmulatorStep4_B.sf_parseCtrlWord_h.rampOutZero != 0) &&
         windEmulatorStep4_cal->Constant_Value_e);

      /* Logic: '<S46>/remoteCmd' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_B.remoteCmd =
        ((windEmulatorStep4_B.sf_parseCtrlWord_h.remoteCmd != 0) &&
         windEmulatorStep4_cal->Constant_Value_e);

      /* Logic: '<S46>/reset' incorporates:
       *  Constant: '<S46>/Constant'
       */
      windEmulatorStep4_B.reset = ((windEmulatorStep4_B.sf_parseCtrlWord_h.reset
        != 0) && windEmulatorStep4_cal->Constant_Value_e);

      /* Bias: '<S26>/ctrlWord' */
      windEmulatorStep4_B.ctrlWord = static_cast<uint16_T>(static_cast<uint32_T>
        (windEmulatorStep4_B.BusAssignment_k.ctrlWord) +
        windEmulatorStep4_cal->ctrlWord_Bias);

      /* Bias: '<S26>/state' */
      windEmulatorStep4_B.state = windEmulatorStep4_B.BusAssignment_k.state +
        windEmulatorStep4_cal->state_Bias;

      /* Gain: '<S26>/torqueSetpoint_Nm' */
      windEmulatorStep4_B.torqueSetpoint_Nm =
        windEmulatorStep4_cal->torqueSetpoint_Nm_Gain *
        windEmulatorStep4_B.BusAssignment_k.torqueSetpoint_Nm;

      /* Gain: '<S26>/torqueSetpoint_percent' */
      windEmulatorStep4_B.torqueSetpoint_percent =
        windEmulatorStep4_cal->torqueSetpoint_percent_Gain *
        windEmulatorStep4_B.BusAssignment_k.torqueSetpoint_percent;
      if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
        /* ToAsyncQueueBlock generated from: '<S23>/acs800Signals' */
        slrtLogSignal
          (windEmulatorStep4_DW.TAQSigLogging_InsertedFor_acs80.SLRTSigHandles,
           (((windEmulatorStep4_M->Timing.clockTick1+
              windEmulatorStep4_M->Timing.clockTickH1* 4294967296.0)) * 0.004));

        /* Gain: '<S43>/rpm -> rad//s' */
        windEmulatorStep4_B.rpmrads_o = tmp_a *
          windEmulatorStep4_B.BusAssignment_h.motorSpeed_rpm;

        /* Product: '<S43>/Product' */
        windEmulatorStep4_B.Product_o = windEmulatorStep4_B.rpmrads_o *
          windEmulatorStep4_B.BusAssignment_h.motorTorque_Nm;
        windEmulatorStep4_MovingAverage(windEmulatorStep4_B.Product_o,
          &windEmulatorStep4_B.MovingAverage,
          &windEmulatorStep4_DW.MovingAverage);

        /* Gain: '<S43>/shaftPowerAverage_W' */
        windEmulatorStep4_B.shaftPowerAverage_W_a =
          windEmulatorStep4_cal->shaftPowerAverage_W_Gain_l *
          windEmulatorStep4_B.MovingAverage.MovingAverage;

        /* Gain: '<S43>/shaftPower_W' */
        windEmulatorStep4_B.shaftPower_W_e =
          windEmulatorStep4_cal->shaftPower_W_Gain_b *
          windEmulatorStep4_B.Product_o;

        /* MATLAB Function: '<S44>/Parse Status Word' */
        windEmulatorSte_ParseStatusWord
          (windEmulatorStep4_B.BusAssignment_h.statusWord,
           &windEmulatorStep4_B.sf_ParseStatusWord);

        /* Logic: '<S44>/aboveLimit' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_B.aboveLimit_c =
          (windEmulatorStep4_cal->Constant_Value_j &&
           (windEmulatorStep4_B.sf_ParseStatusWord.above_limit != 0));

        /* Logic: '<S44>/alarm' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_B.alarm_c = (windEmulatorStep4_cal->Constant_Value_j &&
          (windEmulatorStep4_B.sf_ParseStatusWord.alarm != 0));

        /* Logic: '<S44>/atSetpoint' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_B.atSetpoint_a =
          (windEmulatorStep4_cal->Constant_Value_j &&
           (windEmulatorStep4_B.sf_ParseStatusWord.at_setpoint != 0));

        /* Logic: '<S44>/commErr' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_B.commErr_i = (windEmulatorStep4_cal->Constant_Value_j
          && (windEmulatorStep4_B.sf_ParseStatusWord.comm_err != 0));

        /* Logic: '<S44>/extCtrlLoc' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_B.extCtrlLoc_e =
          (windEmulatorStep4_cal->Constant_Value_j &&
           (windEmulatorStep4_B.sf_ParseStatusWord.ext_ctrl_loc != 0));

        /* Logic: '<S44>/extRunEnable' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_B.extRunEnable_k =
          (windEmulatorStep4_cal->Constant_Value_j &&
           (windEmulatorStep4_B.sf_ParseStatusWord.ext_run_enable != 0));

        /* Logic: '<S44>/mswB13' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_B.mswB13_a = (windEmulatorStep4_cal->Constant_Value_j &&
          (windEmulatorStep4_B.sf_ParseStatusWord.msw_b13 != 0));

        /* Logic: '<S44>/mswB14' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_B.mswB14_m = (windEmulatorStep4_cal->Constant_Value_j &&
          (windEmulatorStep4_B.sf_ParseStatusWord.msw_b14 != 0));

        /* Logic: '<S44>/off2' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_B.off2_n = (windEmulatorStep4_cal->Constant_Value_j &&
          (windEmulatorStep4_B.sf_ParseStatusWord.off2 != 0));

        /* Logic: '<S44>/off3' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_B.off3_p = (windEmulatorStep4_cal->Constant_Value_j &&
          (windEmulatorStep4_B.sf_ParseStatusWord.off3 != 0));

        /* Logic: '<S44>/rdyOn' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_B.rdyOn_b = (windEmulatorStep4_cal->Constant_Value_j &&
          (windEmulatorStep4_B.sf_ParseStatusWord.rdy_on != 0));

        /* Logic: '<S44>/rdyRef' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_B.rdyRef_f = (windEmulatorStep4_cal->Constant_Value_j &&
          (windEmulatorStep4_B.sf_ParseStatusWord.rdy_ref != 0));

        /* Logic: '<S44>/rdyRun' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_B.rdyRun_n = (windEmulatorStep4_cal->Constant_Value_j &&
          (windEmulatorStep4_B.sf_ParseStatusWord.rdy_run != 0));

        /* Logic: '<S44>/remote' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_B.remote_d = (windEmulatorStep4_cal->Constant_Value_j &&
          (windEmulatorStep4_B.sf_ParseStatusWord.remote != 0));

        /* Logic: '<S44>/switchOnInhibit' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_B.switchOnInhibit_p =
          (windEmulatorStep4_cal->Constant_Value_j &&
           (windEmulatorStep4_B.sf_ParseStatusWord.swc_on_inhib != 0));

        /* Logic: '<S44>/tripped' incorporates:
         *  Constant: '<S44>/Constant'
         */
        windEmulatorStep4_B.tripped_p = (windEmulatorStep4_cal->Constant_Value_j
          && (windEmulatorStep4_B.sf_ParseStatusWord.tripped != 0));

        /* Gain: '<S24>/dcBusVoltage_V' */
        windEmulatorStep4_B.dcBusVoltage_V =
          windEmulatorStep4_cal->dcBusVoltage_V_Gain *
          windEmulatorStep4_B.BusAssignment_h.dcBusVoltage_V;

        /* Gain: '<S24>/frequency_Hz' */
        windEmulatorStep4_B.frequency_Hz_j =
          windEmulatorStep4_cal->frequency_Hz_Gain_i *
          windEmulatorStep4_B.BusAssignment_h.frequency_Hz;

        /* Gain: '<S24>/motorSpeed_rpm' */
        windEmulatorStep4_B.motorSpeed_rpm_j =
          windEmulatorStep4_cal->motorSpeed_rpm_Gain_j *
          windEmulatorStep4_B.BusAssignment_h.motorSpeed_rpm;

        /* Gain: '<S24>/motorTorque_Nm' */
        windEmulatorStep4_B.motorTorque_Nm_p =
          windEmulatorStep4_cal->motorTorque_Nm_Gain_l *
          windEmulatorStep4_B.BusAssignment_h.motorTorque_Nm;

        /* Gain: '<S24>/shaftPower_W' */
        windEmulatorStep4_B.shaftPower_W_p =
          windEmulatorStep4_cal->shaftPower_W_Gain_l *
          windEmulatorStep4_B.BusAssignment_h.shaftPower_W;

        /* Bias: '<S24>/statusWord' */
        windEmulatorStep4_B.statusWord_m = static_cast<uint16_T>
          (static_cast<uint32_T>(windEmulatorStep4_B.BusAssignment_h.statusWord)
           + windEmulatorStep4_cal->statusWord_Bias_l);

        /* Gain: '<S24>/temperature' */
        windEmulatorStep4_B.temperature =
          windEmulatorStep4_cal->temperature_Gain *
          windEmulatorStep4_B.BusAssignment_h.temperature;

        /* Memory: '<S1>/Memory' */
        windEmulatorStep4_B.Memory_j =
          windEmulatorStep4_DW.Memory_PreviousInput_d;

        /* RelationalOperator: '<S1>/NotEqual' incorporates:
         *  Constant: '<S1>/powerUpButton'
         */
        windEmulatorStep4_B.NotEqual_f =
          (windEmulatorStep4_cal->powerUpButton_Value_b !=
           windEmulatorStep4_B.Memory_j);

        /* Memory: '<S1>/Memory1' */
        windEmulatorStep4_B.Memory1_f =
          windEmulatorStep4_DW.Memory1_PreviousInput_p;

        /* RelationalOperator: '<S1>/NotEqual1' incorporates:
         *  Constant: '<S1>/powerDownButton'
         */
        windEmulatorStep4_B.NotEqual1_m =
          (windEmulatorStep4_cal->powerDownButton_Value_k !=
           windEmulatorStep4_B.Memory1_f);

        /* Memory: '<S1>/Memory2' */
        windEmulatorStep4_B.Memory2_h =
          windEmulatorStep4_DW.Memory2_PreviousInput_h;

        /* RelationalOperator: '<S1>/NotEqual2' incorporates:
         *  Constant: '<S1>/resetFaultButton'
         */
        windEmulatorStep4_B.NotEqual2_h =
          (windEmulatorStep4_cal->resetFaultButton_Value_n !=
           windEmulatorStep4_B.Memory2_h);

        /* DataTypeConversion: '<S1>/Cast To Double' */
        windEmulatorStep4_B.CastToDouble_j = windEmulatorStep4_B.NotEqual2_h;

        /* Constant: '<S1>/ACS800CtrlMode' */
        windEmulatorStep4_B.ACS800CtrlMode = tmp_l;

        /* Chart: '<S16>/ABB Fieldbus Control' */
        if (windEmulatorStep4_DW.temporalCounter_i1_ly < 31U) {
          windEmulatorStep4_DW.temporalCounter_i1_ly = static_cast<uint8_T>
            (windEmulatorStep4_DW.temporalCounter_i1_ly + 1U);
        }

        windEmulatorStep4_DW.sfEvent_f = windEmulatorStep4_CALL_EVENT;
        if (windEmulatorStep4_DW.is_active_c9_windEmulatorStep4 == 0U) {
          windEmulatorStep4_DW.is_active_c9_windEmulatorStep4 = 1U;
          windEmulatorStep4_DW.is_active_UpdateStateMachine_o = 1U;
          windEmulatorStep4_DW.is_UpdateStateMachine_l =
            windEmulatorStep4_IN_initialize;
          windEmulatorStep4_cwInitialize();
          windEmulatorStep4_B.state_m = abbStateEnum_init;
          windEmulatorStep4_DW.is_active_UpdateControlWord_c = 1U;
        } else {
          windEmulatorS_swParseStatusWord();
          windEmulatorStep4_DW.cwRESET_b = windEmulatorStep4_B.CastToDouble_j;
          switch (windEmulatorStep4_DW.is_UpdateStateMachine_l) {
           case windEmulatorStep4_IN_DelayOFF1:
            if (windEmulatorStep4_DW.temporalCounter_i1_ly >= 25U) {
              windEmulatorStep4_DW.is_UpdateStateMachine_l =
                windEmula_IN_notReadyToSwitchOn;
              windEmulatorStep4_DW.cwOFF1_CONTROL_a = 0.0;
            } else {
              windEmulatorStep4_B.state_m = abbStateEnum_delayOff1;
            }
            break;

           case windEmulatorStep4_IN_initialize:
            windEmulatorStep4_DW.is_UpdateStateMachine_l =
              windEmula_IN_notReadyToSwitchOn;
            windEmulatorStep4_DW.cwOFF1_CONTROL_a = 0.0;
            break;

           case windEmula_IN_notReadyToSwitchOn:
            e_out = ((windEmulatorStep4_DW.swRDY_ON_k == 1.0) &&
                     (windEmulatorStep4_DW.swWARNING_j == 0.0));
            if (e_out) {
              windEmulatorStep4_DW.is_UpdateStateMachine_l =
                windEmulator_IN_readyToSwitchOn;
              windEmulatorStep4_DW.cwOFF1_CONTROL_a = 0.0;
            } else {
              windEmulatorStep4_B.state_m = abbStateEnum_notReadyToSwitchOn;
            }
            break;

           case windEmulat_IN_operationDisabled:
            e_out = (windEmulatorStep4_B.NotEqual_f &&
                     (windEmulatorStep4_DW.swRDY_RUN_i == 1.0));
            if (e_out) {
              windEmulatorStep4_DW.is_UpdateStateMachine_l =
                windEmulato_IN_operationEnabled;
              windEmulatorStep4_DW.cwOFF1_CONTROL_a = 1.0;
              windEmulatorStep4_DW.cwENABLE_OPERATION_k = 1.0;
            } else {
              e_out = (windEmulatorStep4_B.NotEqual1_m ||
                       windEmulatorStep4_B.NotEqual_i ||
                       (windEmulatorStep4_DW.swTRIPPED_h == 1.0) ||
                       (windEmulatorStep4_DW.swSWC_ON_INHIB_a == 1.0) ||
                       (windEmulatorStep4_DW.swWARNING_j == 1.0));
              if (e_out) {
                windEmulatorStep4_DW.is_UpdateStateMachine_l =
                  windEmulatorStep4_IN_DelayOFF1;
                windEmulatorStep4_DW.temporalCounter_i1_ly = 0U;
              } else {
                windEmulatorStep4_B.state_m = abbStateEnum_operationDisabled;
              }
            }
            break;

           case windEmulato_IN_operationEnabled:
            if (windEmulatorStep4_B.NotEqual1_m) {
              windEmulatorStep4_DW.cwENABLE_OPERATION_k = 0.0;
              windEmulatorStep4_DW.is_UpdateStateMachine_l =
                windEmulat_IN_operationDisabled;
              windEmulatorStep4_DW.cwOFF1_CONTROL_a = 1.0;
            } else {
              e_out = (windEmulatorStep4_B.NotEqual_i ||
                       (windEmulatorStep4_DW.swTRIPPED_h == 1.0) ||
                       (windEmulatorStep4_DW.swSWC_ON_INHIB_a == 1.0));
              if (e_out) {
                windEmulatorStep4_DW.cwENABLE_OPERATION_k = 0.0;
                windEmulatorStep4_DW.is_UpdateStateMachine_l =
                  windEmulatorStep4_IN_DelayOFF1;
                windEmulatorStep4_DW.temporalCounter_i1_ly = 0U;
              } else {
                windEmulatorStep4_B.state_m = abbStateEnum_operationEnabled;
              }
            }
            break;

           default:
            /* case IN_readyToSwitchOn: */
            e_out = (windEmulatorStep4_B.NotEqual_f &&
                     (windEmulatorStep4_DW.swREMOTE_l == 1.0));
            if (e_out) {
              windEmulatorStep4_DW.is_UpdateStateMachine_l =
                windEmulat_IN_operationDisabled;
              windEmulatorStep4_DW.cwOFF1_CONTROL_a = 1.0;
            } else {
              windEmulatorStep4_B.state_m = abbStateEnum_readyToSwitchOn;
            }
            break;
          }

          windEmulator_cwBuildControlWord();
        }

        /* End of Chart: '<S16>/ABB Fieldbus Control' */

        /* DataTypeConversion: '<S1>/Cast To Double1' */
        windEmulatorStep4_B.CastToDouble1_ia = windEmulatorStep4_B.state_m;
      }

      /* MultiPortSwitch: '<S3>/Multiport Switch1' incorporates:
       *  DataTypeConversion: '<S3>/toExpTypeEnum'
       */
      switch (windEmulatorStep4_B.toExpTypeEnum) {
       case expTypeEnum_off:
        /* MultiPortSwitch: '<S3>/Multiport Switch1' incorporates:
         *  Constant: '<S3>/zeroTorque'
         */
        windEmulatorStep4_B.MultiportSwitch1_m =
          windEmulatorStep4_cal->zeroTorque_Value;
        break;

       case expTypeEnum_sid:
        /* MultiPortSwitch: '<S3>/Multiport Switch1' */
        windEmulatorStep4_B.MultiportSwitch1_m =
          windEmulatorStep4_B.BusAssignment_f.acs800Torque_Nm;
        break;

       case expTypeEnum_hil:
        /* MultiPortSwitch: '<S3>/Multiport Switch1' */
        windEmulatorStep4_B.MultiportSwitch1_m =
          windEmulatorStep4_B.BusAssignment_n.genTorqueCmd_Nm;
        break;

       default:
        /* MultiPortSwitch: '<S3>/Multiport Switch1' incorporates:
         *  Constant: '<S3>/zeroTorque'
         */
        windEmulatorStep4_B.MultiportSwitch1_m =
          windEmulatorStep4_cal->zeroTorque_Value;
        break;
      }

      /* End of MultiPortSwitch: '<S3>/Multiport Switch1' */

      /* Gain: '<S1>/Nm -> %' */
      riseValLimit = 1.0 / tmp_4 * 100.0;

      /* Gain: '<S1>/Nm -> %' */
      windEmulatorStep4_B.Nm_j = riseValLimit *
        windEmulatorStep4_B.MultiportSwitch1_m;

      /* BusAssignment: '<S1>/Bus Assignment' incorporates:
       *  Constant: '<S1>/Constant'
       */
      windEmulatorStep4_B.BusAssignment_kc = *get_acs800CtrlStruct();

      /* BusAssignment: '<S1>/Bus Assignment' */
      windEmulatorStep4_B.BusAssignment_kc.ctrlWord =
        windEmulatorStep4_B.ControlWord_c;
      windEmulatorStep4_B.BusAssignment_kc.state =
        windEmulatorStep4_B.CastToDouble1_ia;
      windEmulatorStep4_B.BusAssignment_kc.torqueSetpoint_Nm =
        windEmulatorStep4_B.MultiportSwitch1_m;
      windEmulatorStep4_B.BusAssignment_kc.torqueSetpoint_percent =
        windEmulatorStep4_B.Nm_j;
      if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
        /* ToAsyncQueueBlock generated from: '<S21>/acs800CtrlSignals' */
        slrtLogSignal
          (windEmulatorStep4_DW.TAQSigLogging_InsertedFor_acs_l.SLRTSigHandles,
           (((windEmulatorStep4_M->Timing.clockTick1+
              windEmulatorStep4_M->Timing.clockTickH1* 4294967296.0)) * 0.004));
      }

      /* MATLAB Function: '<S41>/parseCtrlWord' */
      windEmulatorStep4_parseCtrlWord
        (windEmulatorStep4_B.BusAssignment_kc.ctrlWord,
         &windEmulatorStep4_B.sf_parseCtrlWord);

      /* Logic: '<S41>/enableOperation' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_B.enableOperation_c =
        ((windEmulatorStep4_B.sf_parseCtrlWord.enableOperation != 0) &&
         windEmulatorStep4_cal->Constant_Value_l);

      /* Logic: '<S41>/extCtrlLoc' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_B.extCtrlLoc_ec =
        ((windEmulatorStep4_B.sf_parseCtrlWord.extCtrlLoc != 0) &&
         windEmulatorStep4_cal->Constant_Value_l);

      /* Logic: '<S41>/inching1' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_B.inching1_o =
        ((windEmulatorStep4_B.sf_parseCtrlWord.inching1 != 0) &&
         windEmulatorStep4_cal->Constant_Value_l);

      /* Logic: '<S41>/inching2' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_B.inching2_j =
        ((windEmulatorStep4_B.sf_parseCtrlWord.inching2 != 0) &&
         windEmulatorStep4_cal->Constant_Value_l);

      /* Logic: '<S41>/off1Ctrl' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_B.off1Ctrl_a =
        ((windEmulatorStep4_B.sf_parseCtrlWord.off1Ctrl != 0) &&
         windEmulatorStep4_cal->Constant_Value_l);

      /* Logic: '<S41>/off2Ctrl' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_B.off2Ctrl_h =
        ((windEmulatorStep4_B.sf_parseCtrlWord.off2Ctrl != 0) &&
         windEmulatorStep4_cal->Constant_Value_l);

      /* Logic: '<S41>/off3Ctrl' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_B.off3Ctrl_a =
        ((windEmulatorStep4_B.sf_parseCtrlWord.off3Ctrl != 0) &&
         windEmulatorStep4_cal->Constant_Value_l);

      /* Logic: '<S41>/rampHold' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_B.rampHold_p =
        ((windEmulatorStep4_B.sf_parseCtrlWord.rampHold != 0) &&
         windEmulatorStep4_cal->Constant_Value_l);

      /* Logic: '<S41>/rampInZero' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_B.rampInZero_e =
        ((windEmulatorStep4_B.sf_parseCtrlWord.rampInZero != 0) &&
         windEmulatorStep4_cal->Constant_Value_l);

      /* Logic: '<S41>/rampOutZero' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_B.rampOutZero_p =
        ((windEmulatorStep4_B.sf_parseCtrlWord.rampOutZero != 0) &&
         windEmulatorStep4_cal->Constant_Value_l);

      /* Logic: '<S41>/remoteCmd' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_B.remoteCmd_p =
        ((windEmulatorStep4_B.sf_parseCtrlWord.remoteCmd != 0) &&
         windEmulatorStep4_cal->Constant_Value_l);

      /* Logic: '<S41>/reset' incorporates:
       *  Constant: '<S41>/Constant'
       */
      windEmulatorStep4_B.reset_b = ((windEmulatorStep4_B.sf_parseCtrlWord.reset
        != 0) && windEmulatorStep4_cal->Constant_Value_l);

      /* Bias: '<S22>/ctrlWord' */
      windEmulatorStep4_B.ctrlWord_g = static_cast<uint16_T>
        (static_cast<uint32_T>(windEmulatorStep4_B.BusAssignment_kc.ctrlWord) +
         windEmulatorStep4_cal->ctrlWord_Bias_a);

      /* Bias: '<S22>/state' */
      windEmulatorStep4_B.state_i = windEmulatorStep4_B.BusAssignment_kc.state +
        windEmulatorStep4_cal->state_Bias_l;

      /* Gain: '<S22>/torqueSetpoint_Nm' */
      windEmulatorStep4_B.torqueSetpoint_Nm_j =
        windEmulatorStep4_cal->torqueSetpoint_Nm_Gain_m *
        windEmulatorStep4_B.BusAssignment_kc.torqueSetpoint_Nm;

      /* Gain: '<S22>/torqueSetpoint_percent' */
      windEmulatorStep4_B.torqueSetpoint_percent_c =
        windEmulatorStep4_cal->torqueSetpoint_percent_Gain_k *
        windEmulatorStep4_B.BusAssignment_kc.torqueSetpoint_percent;
      if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
        /* ToAsyncQueueBlock generated from: '<S34>/hptoSignals' */
        slrtLogSignal
          (windEmulatorStep4_DW.TAQSigLogging_InsertedFor_hptoS.SLRTSigHandles,
           (((windEmulatorStep4_M->Timing.clockTick1+
              windEmulatorStep4_M->Timing.clockTickH1* 4294967296.0)) * 0.004));
      }

      /* Gain: '<S51>/bar->Pa' */
      windEmulatorStep4_B.barPa = *get_bar2pa() *
        windEmulatorStep4_B.BusAssignment_n.pressure_bar;

      /* Gain: '<S51>/l//m -> m3//s' */
      windEmulatorStep4_B.lmm3s = *get_lpm2m3ps() *
        windEmulatorStep4_B.BusAssignment_n.genPumpFlow_lpm;

      /* Product: '<S51>/Product' */
      windEmulatorStep4_B.Product_gm = windEmulatorStep4_B.barPa *
        windEmulatorStep4_B.lmm3s;
      windEmulatorStep4_MovingAverage(windEmulatorStep4_B.Product_gm,
        &windEmulatorStep4_B.MovingAverage_pn,
        &windEmulatorStep4_DW.MovingAverage_pn);

      /* Gain: '<S51>/rpm -> rad//s' */
      windEmulatorStep4_B.rpmrads_oh = tmp_a *
        windEmulatorStep4_B.BusAssignment_n.excShaftSpeed_rpm;

      /* Product: '<S51>/Product1' */
      windEmulatorStep4_B.Product1_i = windEmulatorStep4_B.rpmrads_oh *
        windEmulatorStep4_B.BusAssignment_n.excShaftTorque_Nm;
      windEmulatorStep4_MovingAverage(windEmulatorStep4_B.Product1_i,
        &windEmulatorStep4_B.MovingAverage1,
        &windEmulatorStep4_DW.MovingAverage1);

      /* Gain: '<S51>/excShaftPowerAverage_W' */
      windEmulatorStep4_B.excShaftPowerAverage_W =
        windEmulatorStep4_cal->excShaftPowerAverage_W_Gain *
        windEmulatorStep4_B.MovingAverage1.MovingAverage;

      /* Gain: '<S51>/excShaftPower_W' */
      windEmulatorStep4_B.excShaftPower_W =
        windEmulatorStep4_cal->excShaftPower_W_Gain *
        windEmulatorStep4_B.Product1_i;

      /* Gain: '<S51>/hydrPowerAverage_W' */
      windEmulatorStep4_B.hydrPowerAverage_W =
        windEmulatorStep4_cal->hydrPowerAverage_W_Gain *
        windEmulatorStep4_B.MovingAverage_pn.MovingAverage;

      /* Gain: '<S51>/hydrPower_W' */
      windEmulatorStep4_B.hydrPower_W = windEmulatorStep4_cal->hydrPower_W_Gain *
        windEmulatorStep4_B.Product_gm;

      /* Gain: '<S35>/ctrlSignal1' */
      windEmulatorStep4_B.ctrlSignal1 = windEmulatorStep4_cal->ctrlSignal1_Gain *
        windEmulatorStep4_B.BusAssignment_n.ctrlSignal1;

      /* Gain: '<S35>/ctrlSignal2' */
      windEmulatorStep4_B.ctrlSignal2 = windEmulatorStep4_cal->ctrlSignal2_Gain *
        windEmulatorStep4_B.BusAssignment_n.ctrlSignal2;

      /* Gain: '<S35>/excShaftSpeed_rpm' */
      windEmulatorStep4_B.excShaftSpeed_rpm =
        windEmulatorStep4_cal->excShaftSpeed_rpm_Gain *
        windEmulatorStep4_B.BusAssignment_n.excShaftSpeed_rpm;

      /* Gain: '<S35>/excShaftTorque_Nm' */
      windEmulatorStep4_B.excShaftTorque_Nm =
        windEmulatorStep4_cal->excShaftTorque_Nm_Gain *
        windEmulatorStep4_B.BusAssignment_n.excShaftTorque_Nm;

      /* Gain: '<S35>/genPumpFlow_lpm' */
      windEmulatorStep4_B.genPumpFlow_lpm =
        windEmulatorStep4_cal->genPumpFlow_lpm_Gain *
        windEmulatorStep4_B.BusAssignment_n.genPumpFlow_lpm;

      /* Gain: '<S35>/genShaftSpeed_rpm' */
      windEmulatorStep4_B.genShaftSpeed_rpm =
        windEmulatorStep4_cal->genShaftSpeed_rpm_Gain *
        windEmulatorStep4_B.BusAssignment_n.genShaftSpeed_rpm;

      /* Gain: '<S35>/genTorqueCmd_Nm' */
      windEmulatorStep4_B.genTorqueCmd_Nm =
        windEmulatorStep4_cal->genTorqueCmd_Nm_Gain *
        windEmulatorStep4_B.BusAssignment_n.genTorqueCmd_Nm;

      /* Gain: '<S35>/hmOutputShaftTorque_Nm' */
      windEmulatorStep4_B.hmOutputShaftTorque_Nm =
        windEmulatorStep4_cal->hmOutputShaftTorque_Nm_Gain *
        windEmulatorStep4_B.BusAssignment_n.hmOutputShafTorque_Nm;

      /* Gain: '<S35>/pressure_bar' */
      windEmulatorStep4_B.pressure_bar =
        windEmulatorStep4_cal->pressure_bar_Gain *
        windEmulatorStep4_B.BusAssignment_n.pressure_bar;
      if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
        /* ToAsyncQueueBlock generated from: '<S32>/hptoCtrl' */
        slrtLogSignal
          (windEmulatorStep4_DW.TAQSigLogging_InsertedFor_hptoC.SLRTSigHandles,
           (((windEmulatorStep4_M->Timing.clockTick1+
              windEmulatorStep4_M->Timing.clockTickH1* 4294967296.0)) * 0.004));

        /* Gain: '<S33>/excForce_N' */
        windEmulatorStep4_B.excForce_N = windEmulatorStep4_cal->excForce_N_Gain *
          windEmulatorStep4_B.BusAssignment_c.excForce_N;

        /* Gain: '<S33>/genSpeedActual' */
        windEmulatorStep4_B.genSpeedActual =
          windEmulatorStep4_cal->genSpeedActual_Gain *
          windEmulatorStep4_B.BusAssignment_c.genSpeedActual;

        /* Gain: '<S33>/speedCtrlReset' */
        windEmulatorStep4_B.speedCtrlReset = static_cast<uint8_T>
          (windEmulatorStep4_B.BusAssignment_c.speedCtrlReset ?
           static_cast<int32_T>(windEmulatorStep4_cal->speedCtrlReset_Gain) : 0);

        /* Gain: '<S33>/speedRef_rpm' */
        windEmulatorStep4_B.speedRef_rpm =
          windEmulatorStep4_cal->speedRef_rpm_Gain *
          windEmulatorStep4_B.BusAssignment_c.speedRef_rpm;

        /* ToAsyncQueueBlock generated from: '<S30>/expCtrlSignals' */
        slrtLogSignal
          (windEmulatorStep4_DW.TAQSigLogging_InsertedFor_expCt.SLRTSigHandles,
           (((windEmulatorStep4_M->Timing.clockTick1+
              windEmulatorStep4_M->Timing.clockTickH1* 4294967296.0)) * 0.004));

        /* Bias: '<S31>/expType' */
        windEmulatorStep4_B.expType = static_cast<uint16_T>(static_cast<uint32_T>
          (windEmulatorStep4_B.BusAssignment_b.expType) +
          windEmulatorStep4_cal->expType_Bias);

        /* Gain: '<S31>/ramp' */
        windEmulatorStep4_B.ramp = windEmulatorStep4_cal->ramp_Gain *
          windEmulatorStep4_B.BusAssignment_b.ramp;

        /* Logic: '<S31>/resetHilIntegrator' incorporates:
         *  Constant: '<S31>/Constant'
         */
        windEmulatorStep4_B.resetHilIntegrator =
          (windEmulatorStep4_cal->Constant_Value_kp &&
           windEmulatorStep4_B.BusAssignment_b.resetHilIntegrator);

        /* Logic: '<S31>/resetSidIntegrator' incorporates:
         *  Constant: '<S31>/Constant'
         */
        windEmulatorStep4_B.resetSidIntegrator =
          (windEmulatorStep4_cal->Constant_Value_kp &&
           windEmulatorStep4_B.BusAssignment_b.resetSidIntegrator);

        /* Bias: '<S31>/runCounter' */
        windEmulatorStep4_B.runCounter =
          windEmulatorStep4_B.BusAssignment_b.runCounter +
          windEmulatorStep4_cal->runCounter_Bias;

        /* Logic: '<S31>/runHil' incorporates:
         *  Constant: '<S31>/Constant'
         */
        windEmulatorStep4_B.runHil = (windEmulatorStep4_cal->Constant_Value_kp &&
          windEmulatorStep4_B.BusAssignment_b.runHil);

        /* Logic: '<S31>/runSid' incorporates:
         *  Constant: '<S31>/Constant'
         */
        windEmulatorStep4_B.runSid = (windEmulatorStep4_cal->Constant_Value_kp &&
          windEmulatorStep4_B.BusAssignment_b.runSid);

        /* Bias: '<S31>/stepCounter' */
        windEmulatorStep4_B.stepCounter =
          windEmulatorStep4_B.BusAssignment_b.stepCounter +
          windEmulatorStep4_cal->stepCounter_Bias;

        /* Gain: '<S31>/time' */
        windEmulatorStep4_B.time = windEmulatorStep4_cal->time_Gain *
          windEmulatorStep4_B.BusAssignment_b.time;

        /* S-Function (slecatpdorx): '<S12>/readTorqueInput' */
        {
          /*------------ S-Function Block: <S12>/readTorqueInput PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.readTorqueInput;
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
        windEmulatorStep4_B.CastToDouble_a = windEmulatorStep4_B.readTorqueInput;

        /* Gain: '<S12>/Gain' */
        windEmulatorStep4_B.Gain_nm = *get_futekTorqueScale() *
          windEmulatorStep4_B.CastToDouble_a;

        /* S-Function (slecatpdorx): '<S12>/readEncoderCounter' */
        {
          /*------------ S-Function Block: <S12>/readEncoderCounter PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.readEncoderCounter;
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

        /* Memory: '<S235>/lastRawCounts' */
        windEmulatorStep4_B.lastRawCounts =
          windEmulatorStep4_DW.lastRawCounts_PreviousInput;

        /* DataTypeConversion: '<S235>/Cast To Double' */
        windEmulatorStep4_B.CastToDouble_p = static_cast<int32_T>
          (windEmulatorStep4_B.readEncoderCounter);

        /* Sum: '<S235>/Add' */
        windEmulatorStep4_B.Add_d = windEmulatorStep4_B.lastRawCounts -
          windEmulatorStep4_B.CastToDouble_p;

        /* Abs: '<S235>/Abs' */
        rowIdx = windEmulatorStep4_B.Add_d;
        if (rowIdx < 0) {
          /* Abs: '<S235>/Abs' */
          windEmulatorStep4_B.Abs_a = -rowIdx;
        } else {
          /* Abs: '<S235>/Abs' */
          windEmulatorStep4_B.Abs_a = rowIdx;
        }

        /* End of Abs: '<S235>/Abs' */

        /* RelationalOperator: '<S235>/Relational Operator' incorporates:
         *  Constant: '<S235>/Constant'
         */
        windEmulatorStep4_B.RelationalOperator_l = (windEmulatorStep4_B.Abs_a >=
          windEmulatorStep4_cal->Constant_Value_p);

        /* Switch: '<S235>/Switch' */
        if (windEmulatorStep4_B.RelationalOperator_l) {
          /* Signum: '<S235>/Sign' */
          rowIdx = windEmulatorStep4_B.Add_d;
          if (rowIdx < 0) {
            /* Signum: '<S235>/Sign' */
            windEmulatorStep4_B.Sign = -1;
          } else {
            /* Signum: '<S235>/Sign' */
            windEmulatorStep4_B.Sign = (rowIdx > 0);
          }

          /* End of Signum: '<S235>/Sign' */

          /* Switch: '<S235>/Switch' */
          windEmulatorStep4_B.Switch_g0 = windEmulatorStep4_B.Sign;
        } else {
          /* Switch: '<S235>/Switch' incorporates:
           *  Constant: '<S235>/Constant1'
           */
          windEmulatorStep4_B.Switch_g0 =
            windEmulatorStep4_cal->Constant1_Value_jt;
        }

        /* End of Switch: '<S235>/Switch' */

        /* Memory: '<S235>/lastTurn' */
        windEmulatorStep4_B.lastTurn =
          windEmulatorStep4_DW.lastTurn_PreviousInput;

        /* Sum: '<S235>/Add1' */
        windEmulatorStep4_B.Add1_p = windEmulatorStep4_B.Switch_g0 +
          windEmulatorStep4_B.lastTurn;

        /* DataTypeConversion: '<S235>/Cast To Double3' */
        windEmulatorStep4_B.CastToDouble3_b = windEmulatorStep4_B.CastToDouble_p;

        /* Gain: '<S235>/encoderCountsToRad' */
        windEmulatorStep4_B.encoderCountsToRad = *get_absEncoderCountsToRad() *
          windEmulatorStep4_B.CastToDouble3_b;

        /* DataTypeConversion: '<S235>/Cast To Double1' */
        windEmulatorStep4_B.CastToDouble1_m = windEmulatorStep4_B.Add1_p;

        /* Gain: '<S235>/Gain' */
        windEmulatorStep4_B.Gain_b = windEmulatorStep4_cal->Gain_Gain_cf *
          windEmulatorStep4_B.CastToDouble1_m;

        /* Sum: '<S235>/Add2' */
        windEmulatorStep4_B.Add2 = windEmulatorStep4_B.encoderCountsToRad +
          windEmulatorStep4_B.Gain_b;

        /* SampleTimeMath: '<S236>/TSamp'
         *
         * About '<S236>/TSamp':
         *  y = u * K where K = 1 / ( w * Ts )
         */
        windEmulatorStep4_B.TSamp = windEmulatorStep4_B.Add2 *
          windEmulatorStep4_cal->TSamp_WtEt;

        /* UnitDelay: '<S236>/UD' */
        windEmulatorStep4_B.Uk1 = windEmulatorStep4_DW.UD_DSTATE;

        /* Sum: '<S236>/Diff' */
        windEmulatorStep4_B.Diff = windEmulatorStep4_B.TSamp -
          windEmulatorStep4_B.Uk1;

        /* Gain: '<S12>/rad//s->rpm' */
        windEmulatorStep4_B.radsrpm = tmp_k * windEmulatorStep4_B.Diff;

        /* S-Function (slecatpdorx): '<S12>/EtherCAT PDO Receive7' */
        {
          /*------------ S-Function Block: <S12>/EtherCAT PDO Receive7 PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)&windEmulatorStep4_B.encoderStatus
            [0];
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
        windEmulatorStep4_B.CastToDouble1_f = windEmulatorStep4_B.encoderStatus
          [0];

        /* DataTypeConversion: '<S12>/Cast To Double2' */
        windEmulatorStep4_B.CastToDouble2_d = windEmulatorStep4_B.encoderStatus
          [1];

        /* BusAssignment: '<S12>/Bus Assignment' incorporates:
         *  Constant: '<S12>/Constant'
         */
        windEmulatorStep4_B.BusAssignment_g = *get_shaftSignalStruct();

        /* BusAssignment: '<S12>/Bus Assignment' */
        windEmulatorStep4_B.BusAssignment_g.torqueActual_Nm =
          windEmulatorStep4_B.Gain_nm;
        windEmulatorStep4_B.BusAssignment_g.absEncoderCounts =
          windEmulatorStep4_B.readEncoderCounter;
        windEmulatorStep4_B.BusAssignment_g.absEncoderTurns =
          windEmulatorStep4_B.Add1_p;
        windEmulatorStep4_B.BusAssignment_g.absEncoderPosition_rad =
          windEmulatorStep4_B.Add2;
        windEmulatorStep4_B.BusAssignment_g.absEncoderSpeed_rpm =
          windEmulatorStep4_B.radsrpm;
        windEmulatorStep4_B.BusAssignment_g.absEncoderStatus1 =
          windEmulatorStep4_B.CastToDouble1_f;
        windEmulatorStep4_B.BusAssignment_g.absEncoderStatus2 =
          windEmulatorStep4_B.CastToDouble2_d;

        /* ToAsyncQueueBlock generated from: '<S38>/shaftSignals' */
        slrtLogSignal
          (windEmulatorStep4_DW.TAQSigLogging_InsertedFor_shaft.SLRTSigHandles,
           (((windEmulatorStep4_M->Timing.clockTick1+
              windEmulatorStep4_M->Timing.clockTickH1* 4294967296.0)) * 0.004));

        /* Bias: '<S39>/absEncoderCounts' */
        windEmulatorStep4_B.absEncoderCounts =
          windEmulatorStep4_B.BusAssignment_g.absEncoderCounts +
          windEmulatorStep4_cal->absEncoderCounts_Bias;

        /* Gain: '<S39>/absEncoderPosition_rad' */
        windEmulatorStep4_B.absEncoderPosition_rad =
          windEmulatorStep4_cal->absEncoderPosition_rad_Gain *
          windEmulatorStep4_B.BusAssignment_g.absEncoderPosition_rad;

        /* Gain: '<S39>/absEncoderSpeed_rpm' */
        windEmulatorStep4_B.absEncoderSpeed_rpm =
          windEmulatorStep4_cal->absEncoderSpeed_rpm_Gain *
          windEmulatorStep4_B.BusAssignment_g.absEncoderSpeed_rpm;

        /* Bias: '<S39>/absEncoderStatus1' */
        windEmulatorStep4_B.absEncoderStatus1 = static_cast<uint8_T>(
          static_cast<uint32_T>
          (windEmulatorStep4_B.BusAssignment_g.absEncoderStatus1) +
          windEmulatorStep4_cal->absEncoderStatus1_Bias);

        /* Bias: '<S39>/absEncoderStatus2' */
        windEmulatorStep4_B.absEncoderStatus2 = static_cast<uint8_T>(
          static_cast<uint32_T>
          (windEmulatorStep4_B.BusAssignment_g.absEncoderStatus2) +
          windEmulatorStep4_cal->absEncoderStatus2_Bias);

        /* Bias: '<S39>/absEncoderTurns' */
        windEmulatorStep4_B.absEncoderTurns =
          windEmulatorStep4_B.BusAssignment_g.absEncoderTurns +
          windEmulatorStep4_cal->absEncoderTurns_Bias;

        /* Gain: '<S39>/torqueActual_Nm' */
        windEmulatorStep4_B.torqueActual_Nm =
          windEmulatorStep4_cal->torqueActual_Nm_Gain *
          windEmulatorStep4_B.BusAssignment_g.torqueActual_Nm;

        /* S-Function (slecatpdorx): '<S9>/L1InaccurateURead' */
        {
          /*------------ S-Function Block: <S9>/L1InaccurateURead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L1InaccurateURead;
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
        windEmulatorStep4_B.L1InaccurateU =
          (windEmulatorStep4_cal->Constant1_Value_k &&
           windEmulatorStep4_B.L1InaccurateURead);

        /* S-Function (slecatpdorx): '<S9>/L1InaccurateIRead' */
        {
          /*------------ S-Function Block: <S9>/L1InaccurateIRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L1InaccurateIRead;
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
        windEmulatorStep4_B.L1InaccurateI =
          (windEmulatorStep4_cal->Constant1_Value_k &&
           windEmulatorStep4_B.L1InaccurateIRead);

        /* S-Function (slecatpdorx): '<S9>/L1VoltageRead' */
        {
          /*------------ S-Function Block: <S9>/L1VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)&windEmulatorStep4_B.L1VoltageRead;
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
        windEmulatorStep4_B.L1Voltage = static_cast<real_T>
          (windEmulatorStep4_cal->L1Voltage_Gain) *
          windEmulatorStep4_B.L1VoltageRead;

        /* S-Function (slecatpdorx): '<S9>/L1CurrentRead' */
        {
          /*------------ S-Function Block: <S9>/L1CurrentRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)&windEmulatorStep4_B.L1CurrentRead;
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
        windEmulatorStep4_B.L1Current = static_cast<real_T>
          (windEmulatorStep4_cal->L1Current_Gain) *
          windEmulatorStep4_B.L1CurrentRead;

        /* S-Function (slecatpdorx): '<S9>/L1PowFactorRead' */
        {
          /*------------ S-Function Block: <S9>/L1PowFactorRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L1PowFactorRead;
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
        windEmulatorStep4_B.L1PowFactor = static_cast<real_T>
          (windEmulatorStep4_cal->L1PowFactor_Gain) *
          windEmulatorStep4_B.L1PowFactorRead;

        /* S-Function (slecatpdorx): '<S9>/L1ActivePowRead' */
        {
          /*------------ S-Function Block: <S9>/L1ActivePowRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L1ActivePowRead;
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
        windEmulatorStep4_B.L1ActivePow = static_cast<real_T>
          (windEmulatorStep4_cal->L1ActivePow_Gain) *
          windEmulatorStep4_B.L1ActivePowRead;

        /* S-Function (slecatpdorx): '<S9>/L1THDuRead' */
        {
          /*------------ S-Function Block: <S9>/L1THDuRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)&windEmulatorStep4_B.L1THDuRead;
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
        windEmulatorStep4_B.L1THDu = static_cast<real_T>
          (windEmulatorStep4_cal->L1THDu_Gain) * windEmulatorStep4_B.L1THDuRead;

        /* S-Function (slecatpdorx): '<S9>/L1THDiRead' */
        {
          /*------------ S-Function Block: <S9>/L1THDiRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)&windEmulatorStep4_B.L1THDiRead;
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
        windEmulatorStep4_B.L1THDi = static_cast<real_T>
          (windEmulatorStep4_cal->L1THDi_Gain) * windEmulatorStep4_B.L1THDiRead;

        /* S-Function (slecatpdorx): '<S9>/L2InaccurateURead' */
        {
          /*------------ S-Function Block: <S9>/L2InaccurateURead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L2InaccurateURead;
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
        windEmulatorStep4_B.L2InaccurateU =
          (windEmulatorStep4_cal->Constant2_Value_e &&
           windEmulatorStep4_B.L2InaccurateURead);

        /* S-Function (slecatpdorx): '<S9>/L2InaccurateIRead' */
        {
          /*------------ S-Function Block: <S9>/L2InaccurateIRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L2InaccurateIRead;
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
        windEmulatorStep4_B.L2InaccurateI =
          (windEmulatorStep4_cal->Constant2_Value_e &&
           windEmulatorStep4_B.L2InaccurateIRead);

        /* S-Function (slecatpdorx): '<S9>/L2VoltageRead' */
        {
          /*------------ S-Function Block: <S9>/L2VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)&windEmulatorStep4_B.L2VoltageRead;
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
        windEmulatorStep4_B.L2Voltage = static_cast<real_T>
          (windEmulatorStep4_cal->L2Voltage_Gain) *
          windEmulatorStep4_B.L2VoltageRead;

        /* S-Function (slecatpdorx): '<S9>/L2CurrentRead' */
        {
          /*------------ S-Function Block: <S9>/L2CurrentRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)&windEmulatorStep4_B.L2CurrentRead;
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
        windEmulatorStep4_B.L2Current = static_cast<real_T>
          (windEmulatorStep4_cal->L2Current_Gain) *
          windEmulatorStep4_B.L2CurrentRead;

        /* S-Function (slecatpdorx): '<S9>/L2PowFactorRead' */
        {
          /*------------ S-Function Block: <S9>/L2PowFactorRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L2PowFactorRead;
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
        windEmulatorStep4_B.L2PowFactor = static_cast<real_T>
          (windEmulatorStep4_cal->L2PowFactor_Gain) *
          windEmulatorStep4_B.L2PowFactorRead;

        /* S-Function (slecatpdorx): '<S9>/L2ActivePowRead' */
        {
          /*------------ S-Function Block: <S9>/L2ActivePowRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L2ActivePowRead;
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
        windEmulatorStep4_B.L2ActivePow = static_cast<real_T>
          (windEmulatorStep4_cal->L2ActivePow_Gain) *
          windEmulatorStep4_B.L2ActivePowRead;

        /* S-Function (slecatpdorx): '<S9>/L2THDuRead' */
        {
          /*------------ S-Function Block: <S9>/L2THDuRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)&windEmulatorStep4_B.L2THDuRead;
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
        windEmulatorStep4_B.L2THDu = static_cast<real_T>
          (windEmulatorStep4_cal->L2THDu_Gain) * windEmulatorStep4_B.L2THDuRead;

        /* S-Function (slecatpdorx): '<S9>/L2THDiRead' */
        {
          /*------------ S-Function Block: <S9>/L2THDiRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)&windEmulatorStep4_B.L2THDiRead;
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
        windEmulatorStep4_B.L2THDi = static_cast<real_T>
          (windEmulatorStep4_cal->L2THDi_Gain) * windEmulatorStep4_B.L2THDiRead;

        /* S-Function (slecatpdorx): '<S9>/L3InaccurateURead' */
        {
          /*------------ S-Function Block: <S9>/L3InaccurateURead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L3InaccurateURead;
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
        windEmulatorStep4_B.L3InaccurateU =
          (windEmulatorStep4_cal->Constant3_Value_d &&
           windEmulatorStep4_B.L3InaccurateURead);

        /* S-Function (slecatpdorx): '<S9>/L3InaccurateIRead' */
        {
          /*------------ S-Function Block: <S9>/L3InaccurateIRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L3InaccurateIRead;
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
        windEmulatorStep4_B.L3InaccurateI =
          (windEmulatorStep4_cal->Constant3_Value_d &&
           windEmulatorStep4_B.L3InaccurateIRead);

        /* S-Function (slecatpdorx): '<S9>/L3VoltageRead' */
        {
          /*------------ S-Function Block: <S9>/L3VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)&windEmulatorStep4_B.L3VoltageRead;
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
        windEmulatorStep4_B.L3Voltage = static_cast<real_T>
          (windEmulatorStep4_cal->L3Voltage_Gain) *
          windEmulatorStep4_B.L3VoltageRead;

        /* S-Function (slecatpdorx): '<S9>/L3CurrentRead' */
        {
          /*------------ S-Function Block: <S9>/L3CurrentRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)&windEmulatorStep4_B.L3CurrentRead;
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
        windEmulatorStep4_B.L3Current = static_cast<real_T>
          (windEmulatorStep4_cal->L3Current_Gain) *
          windEmulatorStep4_B.L3CurrentRead;

        /* S-Function (slecatpdorx): '<S9>/L3PowFactorRead' */
        {
          /*------------ S-Function Block: <S9>/L3PowFactorRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L3PowFactorRead;
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
        windEmulatorStep4_B.L3PowFactor = static_cast<real_T>
          (windEmulatorStep4_cal->L3PowFactor_Gain) *
          windEmulatorStep4_B.L3PowFactorRead;

        /* S-Function (slecatpdorx): '<S9>/L3ActivePowRead' */
        {
          /*------------ S-Function Block: <S9>/L3ActivePowRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L3ActivePowRead;
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
        windEmulatorStep4_B.L3ActivePow = static_cast<real_T>
          (windEmulatorStep4_cal->L3ActivePow_Gain) *
          windEmulatorStep4_B.L3ActivePowRead;

        /* S-Function (slecatpdorx): '<S9>/L3THDuRead' */
        {
          /*------------ S-Function Block: <S9>/L3THDuRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)&windEmulatorStep4_B.L3THDuRead;
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
        windEmulatorStep4_B.L3THDu = static_cast<real_T>
          (windEmulatorStep4_cal->L3THDu_Gain) * windEmulatorStep4_B.L3THDuRead;

        /* S-Function (slecatpdorx): '<S9>/L3THDiRead' */
        {
          /*------------ S-Function Block: <S9>/L3THDiRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)&windEmulatorStep4_B.L3THDiRead;
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
        windEmulatorStep4_B.L3THDi = static_cast<real_T>
          (windEmulatorStep4_cal->L3THDi_Gain) * windEmulatorStep4_B.L3THDiRead;

        /* S-Function (slecatpdorx): '<S9>/FrequencyRead' */
        {
          /*------------ S-Function Block: <S9>/FrequencyRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)&windEmulatorStep4_B.FrequencyRead;
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
        windEmulatorStep4_B.totalFrequency = static_cast<real_T>
          (windEmulatorStep4_cal->totalFrequency_Gain) *
          windEmulatorStep4_B.FrequencyRead;

        /* S-Function (slecatpdorx): '<S9>/totalPowFactorRead' */
        {
          /*------------ S-Function Block: <S9>/totalPowFactorRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.totalPowFactorRead;
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
        windEmulatorStep4_B.totalPowFactor = static_cast<real_T>
          (windEmulatorStep4_cal->totalPowFactor_Gain) *
          windEmulatorStep4_B.totalPowFactorRead;

        /* S-Function (slecatpdorx): '<S9>/totalActivePowRead' */
        {
          /*------------ S-Function Block: <S9>/totalActivePowRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.totalActivePowRead;
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
        windEmulatorStep4_B.totalActivePow = static_cast<real_T>
          (windEmulatorStep4_cal->totalActivePow_Gain) *
          windEmulatorStep4_B.totalActivePowRead;

        /* S-Function (slecatpdorx): '<S9>/L1L2VoltageRead' */
        {
          /*------------ S-Function Block: <S9>/L1L2VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L1L2VoltageRead;
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
        windEmulatorStep4_B.L1L2Voltage = static_cast<real_T>
          (windEmulatorStep4_cal->L1L2Voltage_Gain) *
          windEmulatorStep4_B.L1L2VoltageRead;

        /* S-Function (slecatpdorx): '<S9>/L2L3VoltageRead' */
        {
          /*------------ S-Function Block: <S9>/L2L3VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L2L3VoltageRead;
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
        windEmulatorStep4_B.L2L3Voltage = static_cast<real_T>
          (windEmulatorStep4_cal->L2L3Voltage_Gain) *
          windEmulatorStep4_B.L2L3VoltageRead;

        /* S-Function (slecatpdorx): '<S9>/L3L1VoltageRead' */
        {
          /*------------ S-Function Block: <S9>/L3L1VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L3L1VoltageRead;
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
        windEmulatorStep4_B.L3L1Voltage = static_cast<real_T>
          (windEmulatorStep4_cal->L3L1Voltage_Gain) *
          windEmulatorStep4_B.L3L1VoltageRead;

        /* BusAssignment: '<S9>/Bus Assignment' */
        windEmulatorStep4_B.BusAssignment = *tmp_d;

        /* BusAssignment: '<S9>/Bus Assignment' */
        windEmulatorStep4_B.BusAssignment.L1InaccurateU =
          windEmulatorStep4_B.L1InaccurateU;
        windEmulatorStep4_B.BusAssignment.L1InaccurateI =
          windEmulatorStep4_B.L1InaccurateI;
        windEmulatorStep4_B.BusAssignment.L1Voltage =
          windEmulatorStep4_B.L1Voltage;
        windEmulatorStep4_B.BusAssignment.L1Current =
          windEmulatorStep4_B.L1Current;
        windEmulatorStep4_B.BusAssignment.L1PowFactor =
          windEmulatorStep4_B.L1PowFactor;
        windEmulatorStep4_B.BusAssignment.L1ActivePow =
          windEmulatorStep4_B.L1ActivePow;
        windEmulatorStep4_B.BusAssignment.L1THDu = windEmulatorStep4_B.L1THDu;
        windEmulatorStep4_B.BusAssignment.L1THDi = windEmulatorStep4_B.L1THDi;
        windEmulatorStep4_B.BusAssignment.L2InaccurateU =
          windEmulatorStep4_B.L2InaccurateU;
        windEmulatorStep4_B.BusAssignment.L2InaccurateI =
          windEmulatorStep4_B.L2InaccurateI;
        windEmulatorStep4_B.BusAssignment.L2Voltage =
          windEmulatorStep4_B.L2Voltage;
        windEmulatorStep4_B.BusAssignment.L2Current =
          windEmulatorStep4_B.L2Current;
        windEmulatorStep4_B.BusAssignment.L2PowFactor =
          windEmulatorStep4_B.L2PowFactor;
        windEmulatorStep4_B.BusAssignment.L2ActivePow =
          windEmulatorStep4_B.L2ActivePow;
        windEmulatorStep4_B.BusAssignment.L2THDu = windEmulatorStep4_B.L2THDu;
        windEmulatorStep4_B.BusAssignment.L2THDi = windEmulatorStep4_B.L2THDi;
        windEmulatorStep4_B.BusAssignment.L3InaccurateU =
          windEmulatorStep4_B.L3InaccurateU;
        windEmulatorStep4_B.BusAssignment.L3InaccurateI =
          windEmulatorStep4_B.L3InaccurateI;
        windEmulatorStep4_B.BusAssignment.L3Voltage =
          windEmulatorStep4_B.L3Voltage;
        windEmulatorStep4_B.BusAssignment.L3Current =
          windEmulatorStep4_B.L3Current;
        windEmulatorStep4_B.BusAssignment.L3PowFactor =
          windEmulatorStep4_B.L3PowFactor;
        windEmulatorStep4_B.BusAssignment.L3ActivePow =
          windEmulatorStep4_B.L3ActivePow;
        windEmulatorStep4_B.BusAssignment.L3THDu = windEmulatorStep4_B.L3THDu;
        windEmulatorStep4_B.BusAssignment.L3THDi = windEmulatorStep4_B.L3THDi;
        windEmulatorStep4_B.BusAssignment.Frequency =
          windEmulatorStep4_B.totalFrequency;
        windEmulatorStep4_B.BusAssignment.TotalPowFactor =
          windEmulatorStep4_B.totalPowFactor;
        windEmulatorStep4_B.BusAssignment.TotalActivePow =
          windEmulatorStep4_B.totalActivePow;
        windEmulatorStep4_B.BusAssignment.L1L2Voltage =
          windEmulatorStep4_B.L1L2Voltage;
        windEmulatorStep4_B.BusAssignment.L2L3Voltage =
          windEmulatorStep4_B.L2L3Voltage;
        windEmulatorStep4_B.BusAssignment.L3L1Volage =
          windEmulatorStep4_B.L3L1Voltage;

        /* ToAsyncQueueBlock generated from: '<S36>/invPowerAcs800' */
        slrtLogSignal
          (windEmulatorStep4_DW.TAQSigLogging_InsertedFor_invPo.SLRTSigHandles,
           (((windEmulatorStep4_M->Timing.clockTick1+
              windEmulatorStep4_M->Timing.clockTickH1* 4294967296.0)) * 0.004));

        /* S-Function (slecatpdorx): '<S11>/L1InaccurateURead' */
        {
          /*------------ S-Function Block: <S11>/L1InaccurateURead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L1InaccurateURead_l;
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
        windEmulatorStep4_B.L1InaccurateU_f =
          (windEmulatorStep4_cal->Constant1_Value_cc &&
           windEmulatorStep4_B.L1InaccurateURead_l);

        /* S-Function (slecatpdorx): '<S11>/L1InaccurateIRead' */
        {
          /*------------ S-Function Block: <S11>/L1InaccurateIRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L1InaccurateIRead_o;
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
        windEmulatorStep4_B.L1InaccurateI_f =
          (windEmulatorStep4_cal->Constant1_Value_cc &&
           windEmulatorStep4_B.L1InaccurateIRead_o);

        /* S-Function (slecatpdorx): '<S11>/L1VoltageRead' */
        {
          /*------------ S-Function Block: <S11>/L1VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L1VoltageRead_i;
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
        windEmulatorStep4_B.L1Voltage_f = static_cast<real_T>
          (windEmulatorStep4_cal->L1Voltage_Gain_h) *
          windEmulatorStep4_B.L1VoltageRead_i;

        /* S-Function (slecatpdorx): '<S11>/L1CurrentRead' */
        {
          /*------------ S-Function Block: <S11>/L1CurrentRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L1CurrentRead_c;
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
        windEmulatorStep4_B.L1Current_a = static_cast<real_T>
          (windEmulatorStep4_cal->L1Current_Gain_f) *
          windEmulatorStep4_B.L1CurrentRead_c;

        /* S-Function (slecatpdorx): '<S11>/L1PowFactorRead' */
        {
          /*------------ S-Function Block: <S11>/L1PowFactorRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L1PowFactorRead_k;
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
        windEmulatorStep4_B.L1PowFactor_l = static_cast<real_T>
          (windEmulatorStep4_cal->L1PowFactor_Gain_a) *
          windEmulatorStep4_B.L1PowFactorRead_k;

        /* S-Function (slecatpdorx): '<S11>/L1ActivePowRead' */
        {
          /*------------ S-Function Block: <S11>/L1ActivePowRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L1ActivePowRead_h;
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
        windEmulatorStep4_B.L1ActivePow_m = static_cast<real_T>
          (windEmulatorStep4_cal->L1ActivePow_Gain_o) *
          windEmulatorStep4_B.L1ActivePowRead_h;

        /* S-Function (slecatpdorx): '<S11>/L1THDuRead' */
        {
          /*------------ S-Function Block: <S11>/L1THDuRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)&windEmulatorStep4_B.L1THDuRead_a;
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
        windEmulatorStep4_B.L1THDu_l = static_cast<real_T>
          (windEmulatorStep4_cal->L1THDu_Gain_a) *
          windEmulatorStep4_B.L1THDuRead_a;

        /* S-Function (slecatpdorx): '<S11>/L1THDiRead' */
        {
          /*------------ S-Function Block: <S11>/L1THDiRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)&windEmulatorStep4_B.L1THDiRead_a;
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
        windEmulatorStep4_B.L1THDi_e = static_cast<real_T>
          (windEmulatorStep4_cal->L1THDi_Gain_f) *
          windEmulatorStep4_B.L1THDiRead_a;

        /* S-Function (slecatpdorx): '<S11>/L2InaccurateURead' */
        {
          /*------------ S-Function Block: <S11>/L2InaccurateURead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L2InaccurateURead_h;
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
        windEmulatorStep4_B.L2InaccurateU_p =
          (windEmulatorStep4_cal->Constant2_Value_eh &&
           windEmulatorStep4_B.L2InaccurateURead_h);

        /* S-Function (slecatpdorx): '<S11>/L2InaccurateIRead' */
        {
          /*------------ S-Function Block: <S11>/L2InaccurateIRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L2InaccurateIRead_f;
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
        windEmulatorStep4_B.L2InaccurateI_p =
          (windEmulatorStep4_cal->Constant2_Value_eh &&
           windEmulatorStep4_B.L2InaccurateIRead_f);

        /* S-Function (slecatpdorx): '<S11>/L2VoltageRead' */
        {
          /*------------ S-Function Block: <S11>/L2VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L2VoltageRead_p;
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
        windEmulatorStep4_B.L2Voltage_j = static_cast<real_T>
          (windEmulatorStep4_cal->L2Voltage_Gain_a) *
          windEmulatorStep4_B.L2VoltageRead_p;

        /* S-Function (slecatpdorx): '<S11>/L2CurrentRead' */
        {
          /*------------ S-Function Block: <S11>/L2CurrentRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L2CurrentRead_a;
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
        windEmulatorStep4_B.L2Current_h = static_cast<real_T>
          (windEmulatorStep4_cal->L2Current_Gain_j) *
          windEmulatorStep4_B.L2CurrentRead_a;

        /* S-Function (slecatpdorx): '<S11>/L2PowFactorRead' */
        {
          /*------------ S-Function Block: <S11>/L2PowFactorRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L2PowFactorRead_o;
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
        windEmulatorStep4_B.L2PowFactor_p = static_cast<real_T>
          (windEmulatorStep4_cal->L2PowFactor_Gain_p) *
          windEmulatorStep4_B.L2PowFactorRead_o;

        /* S-Function (slecatpdorx): '<S11>/L2ActivePowRead' */
        {
          /*------------ S-Function Block: <S11>/L2ActivePowRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L2ActivePowRead_e;
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
        windEmulatorStep4_B.L2ActivePow_b = static_cast<real_T>
          (windEmulatorStep4_cal->L2ActivePow_Gain_b) *
          windEmulatorStep4_B.L2ActivePowRead_e;

        /* S-Function (slecatpdorx): '<S11>/L2THDuRead' */
        {
          /*------------ S-Function Block: <S11>/L2THDuRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)&windEmulatorStep4_B.L2THDuRead_l;
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
        windEmulatorStep4_B.L2THDu_i = static_cast<real_T>
          (windEmulatorStep4_cal->L2THDu_Gain_i) *
          windEmulatorStep4_B.L2THDuRead_l;

        /* S-Function (slecatpdorx): '<S11>/L2THDiRead' */
        {
          /*------------ S-Function Block: <S11>/L2THDiRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)&windEmulatorStep4_B.L2THDiRead_j;
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
        windEmulatorStep4_B.L2THDi_m = static_cast<real_T>
          (windEmulatorStep4_cal->L2THDi_Gain_j) *
          windEmulatorStep4_B.L2THDiRead_j;

        /* S-Function (slecatpdorx): '<S11>/L3InaccurateURead' */
        {
          /*------------ S-Function Block: <S11>/L3InaccurateURead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L3InaccurateURead_o;
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
        windEmulatorStep4_B.L3InaccurateU_n =
          (windEmulatorStep4_cal->Constant3_Value_g &&
           windEmulatorStep4_B.L3InaccurateURead_o);

        /* S-Function (slecatpdorx): '<S11>/L3InaccurateIRead' */
        {
          /*------------ S-Function Block: <S11>/L3InaccurateIRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L3InaccurateIRead_p;
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
        windEmulatorStep4_B.L3InaccurateI_m =
          (windEmulatorStep4_cal->Constant3_Value_g &&
           windEmulatorStep4_B.L3InaccurateIRead_p);

        /* S-Function (slecatpdorx): '<S11>/L3VoltageRead' */
        {
          /*------------ S-Function Block: <S11>/L3VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L3VoltageRead_h;
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
        windEmulatorStep4_B.L3Voltage_n = static_cast<real_T>
          (windEmulatorStep4_cal->L3Voltage_Gain_h) *
          windEmulatorStep4_B.L3VoltageRead_h;

        /* S-Function (slecatpdorx): '<S11>/L3CurrentRead' */
        {
          /*------------ S-Function Block: <S11>/L3CurrentRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L3CurrentRead_h;
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
        windEmulatorStep4_B.L3Current_h = static_cast<real_T>
          (windEmulatorStep4_cal->L3Current_Gain_i) *
          windEmulatorStep4_B.L3CurrentRead_h;

        /* S-Function (slecatpdorx): '<S11>/L3PowFactorRead' */
        {
          /*------------ S-Function Block: <S11>/L3PowFactorRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L3PowFactorRead_o;
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
        windEmulatorStep4_B.L3PowFactor_o = static_cast<real_T>
          (windEmulatorStep4_cal->L3PowFactor_Gain_k) *
          windEmulatorStep4_B.L3PowFactorRead_o;

        /* S-Function (slecatpdorx): '<S11>/L3ActivePowRead' */
        {
          /*------------ S-Function Block: <S11>/L3ActivePowRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L3ActivePowRead_m;
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
        windEmulatorStep4_B.L3ActivePow_l = static_cast<real_T>
          (windEmulatorStep4_cal->L3ActivePow_Gain_f) *
          windEmulatorStep4_B.L3ActivePowRead_m;

        /* S-Function (slecatpdorx): '<S11>/L3THDuRead' */
        {
          /*------------ S-Function Block: <S11>/L3THDuRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)&windEmulatorStep4_B.L3THDuRead_b;
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
        windEmulatorStep4_B.L3THDu_i = static_cast<real_T>
          (windEmulatorStep4_cal->L3THDu_Gain_b) *
          windEmulatorStep4_B.L3THDuRead_b;

        /* S-Function (slecatpdorx): '<S11>/L3THDiRead' */
        {
          /*------------ S-Function Block: <S11>/L3THDiRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)&windEmulatorStep4_B.L3THDiRead_m;
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
        windEmulatorStep4_B.L3THDi_j = static_cast<real_T>
          (windEmulatorStep4_cal->L3THDi_Gain_p) *
          windEmulatorStep4_B.L3THDiRead_m;

        /* S-Function (slecatpdorx): '<S11>/FrequencyRead' */
        {
          /*------------ S-Function Block: <S11>/FrequencyRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.FrequencyRead_g;
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
        windEmulatorStep4_B.totalFrequency_g = static_cast<real_T>
          (windEmulatorStep4_cal->totalFrequency_Gain_b) *
          windEmulatorStep4_B.FrequencyRead_g;

        /* S-Function (slecatpdorx): '<S11>/totalPowFactorRead' */
        {
          /*------------ S-Function Block: <S11>/totalPowFactorRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.totalPowFactorRead_n;
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
        windEmulatorStep4_B.totalPowFactor_e = static_cast<real_T>
          (windEmulatorStep4_cal->totalPowFactor_Gain_f) *
          windEmulatorStep4_B.totalPowFactorRead_n;

        /* S-Function (slecatpdorx): '<S11>/totalActivePowRead' */
        {
          /*------------ S-Function Block: <S11>/totalActivePowRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.totalActivePowRead_o;
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
        windEmulatorStep4_B.totalActivePow_g = static_cast<real_T>
          (windEmulatorStep4_cal->totalActivePow_Gain_m) *
          windEmulatorStep4_B.totalActivePowRead_o;

        /* S-Function (slecatpdorx): '<S11>/L1L2VoltageRead' */
        {
          /*------------ S-Function Block: <S11>/L1L2VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L1L2VoltageRead_m;
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
        windEmulatorStep4_B.L1L2Voltage_n = static_cast<real_T>
          (windEmulatorStep4_cal->L1L2Voltage_Gain_j) *
          windEmulatorStep4_B.L1L2VoltageRead_m;

        /* S-Function (slecatpdorx): '<S11>/L2L3VoltageRead' */
        {
          /*------------ S-Function Block: <S11>/L2L3VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L2L3VoltageRead_e;
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
        windEmulatorStep4_B.L2L3Voltage_h = static_cast<real_T>
          (windEmulatorStep4_cal->L2L3Voltage_Gain_m) *
          windEmulatorStep4_B.L2L3VoltageRead_e;

        /* S-Function (slecatpdorx): '<S11>/L3L1VoltageRead' */
        {
          /*------------ S-Function Block: <S11>/L3L1VoltageRead PDO receive block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          uint8_T *sigOutputPtr = (uint8_T *)
            &windEmulatorStep4_B.L3L1VoltageRead_n;
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
        windEmulatorStep4_B.L3L1Voltage_m = static_cast<real_T>
          (windEmulatorStep4_cal->L3L1Voltage_Gain_k) *
          windEmulatorStep4_B.L3L1VoltageRead_n;

        /* BusAssignment: '<S11>/Bus Assignment' */
        windEmulatorStep4_B.BusAssignment_l = *tmp_d;

        /* BusAssignment: '<S11>/Bus Assignment' */
        windEmulatorStep4_B.BusAssignment_l.L1InaccurateU =
          windEmulatorStep4_B.L1InaccurateU_f;
        windEmulatorStep4_B.BusAssignment_l.L1InaccurateI =
          windEmulatorStep4_B.L1InaccurateI_f;
        windEmulatorStep4_B.BusAssignment_l.L1Voltage =
          windEmulatorStep4_B.L1Voltage_f;
        windEmulatorStep4_B.BusAssignment_l.L1Current =
          windEmulatorStep4_B.L1Current_a;
        windEmulatorStep4_B.BusAssignment_l.L1PowFactor =
          windEmulatorStep4_B.L1PowFactor_l;
        windEmulatorStep4_B.BusAssignment_l.L1ActivePow =
          windEmulatorStep4_B.L1ActivePow_m;
        windEmulatorStep4_B.BusAssignment_l.L1THDu =
          windEmulatorStep4_B.L1THDu_l;
        windEmulatorStep4_B.BusAssignment_l.L1THDi =
          windEmulatorStep4_B.L1THDi_e;
        windEmulatorStep4_B.BusAssignment_l.L2InaccurateU =
          windEmulatorStep4_B.L2InaccurateU_p;
        windEmulatorStep4_B.BusAssignment_l.L2InaccurateI =
          windEmulatorStep4_B.L2InaccurateI_p;
        windEmulatorStep4_B.BusAssignment_l.L2Voltage =
          windEmulatorStep4_B.L2Voltage_j;
        windEmulatorStep4_B.BusAssignment_l.L2Current =
          windEmulatorStep4_B.L2Current_h;
        windEmulatorStep4_B.BusAssignment_l.L2PowFactor =
          windEmulatorStep4_B.L2PowFactor_p;
        windEmulatorStep4_B.BusAssignment_l.L2ActivePow =
          windEmulatorStep4_B.L2ActivePow_b;
        windEmulatorStep4_B.BusAssignment_l.L2THDu =
          windEmulatorStep4_B.L2THDu_i;
        windEmulatorStep4_B.BusAssignment_l.L2THDi =
          windEmulatorStep4_B.L2THDi_m;
        windEmulatorStep4_B.BusAssignment_l.L3InaccurateU =
          windEmulatorStep4_B.L3InaccurateU_n;
        windEmulatorStep4_B.BusAssignment_l.L3InaccurateI =
          windEmulatorStep4_B.L3InaccurateI_m;
        windEmulatorStep4_B.BusAssignment_l.L3Voltage =
          windEmulatorStep4_B.L3Voltage_n;
        windEmulatorStep4_B.BusAssignment_l.L3Current =
          windEmulatorStep4_B.L3Current_h;
        windEmulatorStep4_B.BusAssignment_l.L3PowFactor =
          windEmulatorStep4_B.L3PowFactor_o;
        windEmulatorStep4_B.BusAssignment_l.L3ActivePow =
          windEmulatorStep4_B.L3ActivePow_l;
        windEmulatorStep4_B.BusAssignment_l.L3THDu =
          windEmulatorStep4_B.L3THDu_i;
        windEmulatorStep4_B.BusAssignment_l.L3THDi =
          windEmulatorStep4_B.L3THDi_j;
        windEmulatorStep4_B.BusAssignment_l.Frequency =
          windEmulatorStep4_B.totalFrequency_g;
        windEmulatorStep4_B.BusAssignment_l.TotalPowFactor =
          windEmulatorStep4_B.totalPowFactor_e;
        windEmulatorStep4_B.BusAssignment_l.TotalActivePow =
          windEmulatorStep4_B.totalActivePow_g;
        windEmulatorStep4_B.BusAssignment_l.L1L2Voltage =
          windEmulatorStep4_B.L1L2Voltage_n;
        windEmulatorStep4_B.BusAssignment_l.L2L3Voltage =
          windEmulatorStep4_B.L2L3Voltage_h;
        windEmulatorStep4_B.BusAssignment_l.L3L1Volage =
          windEmulatorStep4_B.L3L1Voltage_m;

        /* ToAsyncQueueBlock generated from: '<S37>/invPowerAcs880' */
        slrtLogSignal
          (windEmulatorStep4_DW.TAQSigLogging_InsertedFor_inv_p.SLRTSigHandles,
           (((windEmulatorStep4_M->Timing.clockTick1+
              windEmulatorStep4_M->Timing.clockTickH1* 4294967296.0)) * 0.004));

        /* DataTypeConversion: '<S238>/Cast To uint32' incorporates:
         *  Constant: '<S238>/sidType'
         */
        windEmulatorStep4_B.CastTouint32 = windEmulatorStep4_cal->sidType_Value;
      }

      /* Switch: '<S288>/Switch3' */
      if (windEmulatorStep4_B.Outofbounds) {
        /* Switch: '<S288>/Switch3' incorporates:
         *  Constant: '<S288>/Set bound'
         */
        windEmulatorStep4_B.Switch3 = windEmulatorStep4_cal->Setbound_Value;
      } else {
        /* Switch: '<S288>/Switch3' incorporates:
         *  Inport: '<Root>/inportCaseCounter'
         */
        windEmulatorStep4_B.Switch3 = windEmulatorStep4_U.inportCaseCounter;
      }

      /* End of Switch: '<S288>/Switch3' */

      /* Gain: '<S288>/caseCounterSignalsNow' */
      windEmulatorStep4_B.caseCounterSignalsNow =
        windEmulatorStep4_cal->caseCounterSignalsNow_Gain *
        windEmulatorStep4_B.Switch3;

      /* BusAssignment: '<S238>/Bus Assignment' incorporates:
       *  Constant: '<S238>/Constant'
       */
      windEmulatorStep4_B.BusAssignment_j = *get_sidInfoStruct();

      /* BusAssignment: '<S238>/Bus Assignment' */
      windEmulatorStep4_B.BusAssignment_j.acs800TorqueSetpoint_Nm =
        windEmulatorStep4_B.Product_g;
      windEmulatorStep4_B.BusAssignment_j.acs880SpeedSetpoint_rpm =
        windEmulatorStep4_B.Product1;
      windEmulatorStep4_B.BusAssignment_j.sidType =
        windEmulatorStep4_B.CastTouint32;
      windEmulatorStep4_B.BusAssignment_j.fromFileCaseCounter =
        windEmulatorStep4_B.caseCounterSignalsNow;
      if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
        /* ToAsyncQueueBlock generated from: '<S40>/sidInfoSignals' */
        slrtLogSignal
          (windEmulatorStep4_DW.TAQSigLogging_InsertedFor_sidIn.SLRTSigHandles,
           (((windEmulatorStep4_M->Timing.clockTick1+
              windEmulatorStep4_M->Timing.clockTickH1* 4294967296.0)) * 0.004));

        /* Constant: '<S5>/Constant' */
        windEmulatorStep4_B.Constant = windEmulatorStep4_cal->Constant_Value_m3;

        /* S-Function (slrealtimeenablelogging): '<S5>/Enable File Log' */

        /* Level2 S-Function Block: '<S5>/Enable File Log' (slrealtimeenablelogging) */
        {
          SimStruct *rts = windEmulatorStep4_M->childSfunctions[0];
          sfcnOutputs(rts,0);
        }

        /* Memory: '<S29>/Memory' */
        windEmulatorStep4_B.Memory = windEmulatorStep4_DW.Memory_PreviousInput;

        /* Sum: '<S29>/Sum' incorporates:
         *  Constant: '<S29>/loopAdd'
         */
        windEmulatorStep4_B.Sum_e = windEmulatorStep4_cal->loopAdd_Value +
          windEmulatorStep4_B.Memory;

        /* Gain: '<S29>/loopCounter' */
        windEmulatorStep4_B.loopCounter =
          windEmulatorStep4_cal->loopCounter_Gain * windEmulatorStep4_B.Sum_e;

        /* DataTypeConversion: '<S29>/Data Type Conversion' */
        windEmulatorStep4_B.DataTypeConversion_f =
          windEmulatorStep4_B.loopCounter;

        /* Gain: '<S29>/time_s' */
        windEmulatorStep4_B.time_s = *get_Ts() *
          windEmulatorStep4_B.DataTypeConversion_f;

        /* S-Function (slecatpdotx): '<S14>/ACS800CtrlWord' */
        {
          /*------------ S-Function Block: <S14>/ACS800CtrlWord PDO transmit block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          int_T i;
          uint8_T *sigInputPtr = (uint8_T *)
            &windEmulatorStep4_B.BusAssignment_kc.ctrlWord;
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

      /* Product: '<S290>/Product' incorporates:
       *  Constant: '<S290>/Constant EngineeringValue'
       *  Constant: '<S290>/Constant FieldbusValue'
       */
      windEmulatorStep4_B.Product_k =
        windEmulatorStep4_B.BusAssignment_kc.torqueSetpoint_percent *
        *get_acs800TorqueNomFb() / *get_acs800TorqueNomEng();

      /* DataTypeConversion: '<S14>/ACS800Ref2Int16' */
      u1 = std::floor(windEmulatorStep4_B.Product_k);
      if (rtIsNaN(u1) || rtIsInf(u1)) {
        u1 = 0.0;
      } else {
        u1 = std::fmod(u1, 65536.0);
      }

      /* DataTypeConversion: '<S14>/ACS800Ref2Int16' */
      windEmulatorStep4_B.ACS800Ref2Int16 = static_cast<int16_T>(u1 < 0.0 ?
        static_cast<int32_T>(static_cast<int16_T>(-static_cast<int16_T>(
        static_cast<uint16_T>(-u1)))) : static_cast<int32_T>(static_cast<int16_T>
        (static_cast<uint16_T>(u1))));
      if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
        /* S-Function (slecatpdotx): '<S14>/ACS800TorqueSetpoint' */
        {
          /*------------ S-Function Block: <S14>/ACS800TorqueSetpoint PDO transmit block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          int_T i;
          uint8_T *sigInputPtr = (uint8_T *)&windEmulatorStep4_B.ACS800Ref2Int16;
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

        /* Product: '<S289>/Product' incorporates:
         *  Constant: '<S14>/Constant'
         *  Constant: '<S289>/Constant EngineeringValue'
         *  Constant: '<S289>/Constant FieldbusValue'
         */
        windEmulatorStep4_B.Product_d = windEmulatorStep4_cal->Constant_Value_m *
          *get_acs800SpeedNomFb() / *get_acs800SpeedNomEng();

        /* DataTypeConversion: '<S14>/ACS800Ref1Int16' */
        u1 = std::floor(windEmulatorStep4_B.Product_d);
        if (rtIsNaN(u1) || rtIsInf(u1)) {
          u1 = 0.0;
        } else {
          u1 = std::fmod(u1, 65536.0);
        }

        /* DataTypeConversion: '<S14>/ACS800Ref1Int16' */
        windEmulatorStep4_B.ACS800Ref1Int16 = static_cast<int16_T>(u1 < 0.0 ?
          static_cast<int32_T>(static_cast<int16_T>(-static_cast<int16_T>(
          static_cast<uint16_T>(-u1)))) : static_cast<int32_T>(static_cast<
          int16_T>(static_cast<uint16_T>(u1))));

        /* S-Function (slecatpdotx): '<S14>/ACS800SpeedSetpoint' */
        {
          /*------------ S-Function Block: <S14>/ACS800SpeedSetpoint PDO transmit block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          int_T i;
          uint8_T *sigInputPtr = (uint8_T *)&windEmulatorStep4_B.ACS800Ref1Int16;
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
            &windEmulatorStep4_B.BusAssignment_k.ctrlWord;
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

      /* Product: '<S291>/Product' incorporates:
       *  Constant: '<S291>/Constant EngineeringValue'
       *  Constant: '<S291>/Constant FieldbusValue'
       */
      windEmulatorStep4_B.Product_l =
        windEmulatorStep4_B.BusAssignment_k.torqueSetpoint_percent *
        *get_acs880TorqueFieldbusScale() / *get_acs880TorqueSetpointScaling();

      /* DataTypeConversion: '<S15>/ACS880Ref2Int16' */
      u1 = std::floor(windEmulatorStep4_B.Product_l);
      if (rtIsNaN(u1) || rtIsInf(u1)) {
        u1 = 0.0;
      } else {
        u1 = std::fmod(u1, 65536.0);
      }

      /* DataTypeConversion: '<S15>/ACS880Ref2Int16' */
      windEmulatorStep4_B.ACS880Ref2Int16 = static_cast<int16_T>(u1 < 0.0 ?
        static_cast<int32_T>(static_cast<int16_T>(-static_cast<int16_T>(
        static_cast<uint16_T>(-u1)))) : static_cast<int32_T>(static_cast<int16_T>
        (static_cast<uint16_T>(u1))));
      if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
        /* S-Function (slecatpdotx): '<S15>/ACS880TorqueSetpoint' */
        {
          /*------------ S-Function Block: <S15>/ACS880TorqueSetpoint PDO transmit block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          int_T i;
          uint8_T *sigInputPtr = (uint8_T *)&windEmulatorStep4_B.ACS880Ref2Int16;
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
        u1 = std::floor(windEmulatorStep4_cal->Constant2_Value);
        if (rtIsNaN(u1) || rtIsInf(u1)) {
          u1 = 0.0;
        } else {
          u1 = std::fmod(u1, 65536.0);
        }

        /* DataTypeConversion: '<S15>/ACS880Ref1Int16' */
        windEmulatorStep4_B.ACS880Ref1Int16 = static_cast<int16_T>(u1 < 0.0 ?
          static_cast<int32_T>(static_cast<int16_T>(-static_cast<int16_T>(
          static_cast<uint16_T>(-u1)))) : static_cast<int32_T>(static_cast<
          int16_T>(static_cast<uint16_T>(u1))));

        /* S-Function (slecatpdotx): '<S15>/ACS880SpeedSetpoint' */
        {
          /*------------ S-Function Block: <S15>/ACS880SpeedSetpoint PDO transmit block  ------------*/
          static int counter= 0;
          int_T sigIdx;
          int_T bitOffset;
          int_T i;
          uint8_T *sigInputPtr = (uint8_T *)&windEmulatorStep4_B.ACS880Ref1Int16;
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

        /* Gain: '<S94>/Integral Gain' */
        windEmulatorStep4_B.IntegralGain = tmp_9->IG *
          windEmulatorStep4_B.wError;

        /* Gain: '<S153>/Integral Gain' */
        windEmulatorStep4_B.IntegralGain_h = tmp_9->IG *
          windEmulatorStep4_B.wError_c;
      }

      /* Switch: '<S53>/Switch1' */
      if (windEmulatorStep4_B.Abs >= tmp_h) {
        /* Switch: '<S53>/Switch1' incorporates:
         *  Constant: '<S53>/Constant'
         */
        windEmulatorStep4_B.SwitchLogic =
          windEmulatorStep4_cal->Constant_Value_i;
      } else {
        /* Switch: '<S53>/Switch1' incorporates:
         *  Constant: '<S53>/Constant1'
         */
        windEmulatorStep4_B.SwitchLogic =
          windEmulatorStep4_cal->Constant1_Value_p;
      }

      if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
        /* Gain: '<S179>/m3toL' */
        windEmulatorStep4_B.FlowPump1 = windEmulatorStep4_cal->m3toL_Gain *
          windEmulatorStep4_B.OUTPUT_1_0[0];

        /* Sum: '<S202>/Sum' */
        windEmulatorStep4_B.Sum_a =
          windEmulatorStep4_B.BusAssignment_c.excForce_N +
          windEmulatorStep4_B.OUTPUT_1_0[4];

        /* Gain: '<S180>/Gain' */
        windEmulatorStep4_B.FlowAccumulator = windEmulatorStep4_cal->Gain_Gain_p
          * windEmulatorStep4_B.OUTPUT_1_0[1];

        /* DataTypeConversion: '<S288>/vecIndex' */
        windEmulatorStep4_B.vecIndex = windEmulatorStep4_B.Mod;

        /* DataTypeConversion: '<S288>/Cast To Double' */
        windEmulatorStep4_B.CastToDouble_g = windEmulatorStep4_B.vecIndex;

        /* DataTypeConversion: '<S288>/Cast To Double1' incorporates:
         *  Constant: '<S288>/Length of input'
         */
        windEmulatorStep4_B.CastToDouble1_i =
          windEmulatorStep4_cal->Lengthofinput_Value;

        /* Product: '<S288>/Divide1' */
        windEmulatorStep4_B.Divide1 = windEmulatorStep4_B.CastToDouble_g /
          windEmulatorStep4_B.CastToDouble1_i;

        /* Bias: '<S288>/fileSamples' incorporates:
         *  Constant: '<S288>/Length of input'
         */
        windEmulatorStep4_B.fileSamples =
          windEmulatorStep4_cal->Lengthofinput_Value +
          windEmulatorStep4_cal->fileSamples_Bias;

        /* Gain: '<S288>/vecPercent' */
        windEmulatorStep4_B.vecPercent = windEmulatorStep4_cal->vecPercent_Gain *
          windEmulatorStep4_B.Divide1;
      }

      /* Product: '<S267>/IProd Out' incorporates:
       *  Constant: '<S13>/acs880SpeedIGain'
       */
      windEmulatorStep4_B.IProdOut = windEmulatorStep4_B.Sum *
        windEmulatorStep4_cal->acs880SpeedIGain_Value;

      /* user code (Output function Trailer) */
      {
        /*------------ S-Function Block: <Root>/EtherCAT Init Write Process Data ,Run Admin Tasks and then Write Acyclic Data------------*/
        xpcEtherCATWriteProcessData(0,NULL);
        xpcEtherCATExecAdminJobs(0);
        xpcEtherCATWriteAcyclicData(0);
      }
    }
  }

  if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
    NeslSimulationData *simulationData;
    NeslSimulator *simulator;
    NeuDiagnosticManager *diagnosticManager;
    real_T tmp_0[20];
    real_T time;
    int_T tmp_1[6];
    if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
      /* Update for Memory: '<S2>/Memory' incorporates:
       *  Constant: '<S2>/powerUpButton'
       */
      windEmulatorStep4_DW.Memory_PreviousInput_i =
        windEmulatorStep4_cal->powerUpButton_Value;

      /* Update for Memory: '<S2>/Memory1' incorporates:
       *  Constant: '<S2>/powerDownButton'
       */
      windEmulatorStep4_DW.Memory1_PreviousInput =
        windEmulatorStep4_cal->powerDownButton_Value;

      /* Update for Memory: '<S2>/Memory2' incorporates:
       *  Constant: '<S2>/resetFaultButton'
       */
      windEmulatorStep4_DW.Memory2_PreviousInput =
        windEmulatorStep4_cal->resetFaultButton_Value;

      /* Update for Memory: '<S4>/Memory' incorporates:
       *  Constant: '<S4>/eStopButton'
       */
      windEmulatorStep4_DW.Memory_PreviousInput_k =
        windEmulatorStep4_cal->eStopButton_Value;

      /* Update for Memory: '<S4>/Memory1' incorporates:
       *  Constant: '<S4>/startButton'
       */
      windEmulatorStep4_DW.Memory1_PreviousInput_d =
        windEmulatorStep4_cal->startButton_Value;

      /* Update for Memory: '<S4>/Memory2' incorporates:
       *  Constant: '<S4>/stopButton'
       */
      windEmulatorStep4_DW.Memory2_PreviousInput_l =
        windEmulatorStep4_cal->stopButton_Value;

      /* Update for Memory: '<S117>/Memory' */
      windEmulatorStep4_DW.Memory_PreviousInput_kk = windEmulatorStep4_B.Logic[0];

      /* Update for Memory: '<S118>/Memory' */
      windEmulatorStep4_DW.Memory_PreviousInput_h = windEmulatorStep4_B.Logic_g
        [0];
    }

    /* Update for RateLimiter: '<S288>/torqueSlewRate' */
    windEmulatorStep4_DW.PrevY = windEmulatorStep4_B.torqueSlewRate;
    windEmulatorStep4_DW.LastMajorTime = windEmulatorStep4_M->Timing.t[0];

    /* Update for RateLimiter: '<S288>/speedSlewRate' */
    windEmulatorStep4_DW.PrevY_a = windEmulatorStep4_B.speedSlewRate;
    windEmulatorStep4_DW.LastMajorTime_p = windEmulatorStep4_M->Timing.t[0];

    /* Update for RateLimiter: '<S61>/Rate Limiter' */
    windEmulatorStep4_DW.PrevY_l = windEmulatorStep4_B.RateLimiter_b;
    windEmulatorStep4_DW.LastMajorTime_m = windEmulatorStep4_M->Timing.t[0];
    if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
      NeuDiagnosticTree *diagnosticTree;
      int32_T tmp_2;
      boolean_T tmp;

      /* Update for DiscreteIntegrator: '<S97>/Integrator' */
      windEmulatorStep4_DW.Integrator_DSTATE +=
        windEmulatorStep4_cal->Integrator_gainval *
        windEmulatorStep4_B.IntegralGain;
      windEmulatorStep4_DW.Integrator_PrevResetState = static_cast<int8_T>
        (windEmulatorStep4_B.BusAssignment_c.speedCtrlReset);

      /* Update for DiscreteIntegrator: '<S92>/Filter' */
      windEmulatorStep4_DW.Filter_DSTATE +=
        windEmulatorStep4_cal->Filter_gainval *
        windEmulatorStep4_B.FilterCoefficient;
      windEmulatorStep4_DW.Filter_PrevResetState = static_cast<int8_T>
        (windEmulatorStep4_B.BusAssignment_c.speedCtrlReset);

      /* Update for SimscapeExecutionBlock: '<S222>/STATE_1' */
      simulationData = static_cast<NeslSimulationData *>
        (windEmulatorStep4_DW.STATE_1_SimData);
      time = windEmulatorStep4_M->Timing.t[0];
      simulationData->mData->mTime.mN = 1;
      simulationData->mData->mTime.mX = &time;
      simulationData->mData->mContStates.mN = 0;
      simulationData->mData->mContStates.mX = NULL;
      simulationData->mData->mDiscStates.mN = 22;
      simulationData->mData->mDiscStates.mX =
        &windEmulatorStep4_DW.STATE_1_Discrete[0];
      simulationData->mData->mModeVector.mN = 17;
      simulationData->mData->mModeVector.mX =
        &windEmulatorStep4_DW.STATE_1_Modes[0];
      tmp = false;
      simulationData->mData->mFoundZcEvents = tmp;
      simulationData->mData->mIsMajorTimeStep = true;
      tmp = false;
      simulationData->mData->mIsSolverAssertCheck = tmp;
      simulationData->mData->mIsSolverCheckingCIC = false;
      simulationData->mData->mIsComputingJacobian = false;
      simulationData->mData->mIsEvaluatingF0 = false;
      simulationData->mData->mIsSolverRequestingReset = false;
      simulationData->mData->mIsModeUpdateTimeStep = true;
      tmp_1[0] = 0;
      tmp_0[0] = windEmulatorStep4_B.INPUT_1_1_1[0];
      tmp_0[1] = windEmulatorStep4_B.INPUT_1_1_1[1];
      tmp_0[2] = windEmulatorStep4_B.INPUT_1_1_1[2];
      tmp_0[3] = windEmulatorStep4_B.INPUT_1_1_1[3];
      tmp_1[1] = 4;
      tmp_0[4] = windEmulatorStep4_B.INPUT_5_1_1[0];
      tmp_0[5] = windEmulatorStep4_B.INPUT_5_1_1[1];
      tmp_0[6] = windEmulatorStep4_B.INPUT_5_1_1[2];
      tmp_0[7] = windEmulatorStep4_B.INPUT_5_1_1[3];
      tmp_1[2] = 8;
      tmp_0[8] = windEmulatorStep4_B.INPUT_2_1_1[0];
      tmp_0[9] = windEmulatorStep4_B.INPUT_2_1_1[1];
      tmp_0[10] = windEmulatorStep4_B.INPUT_2_1_1[2];
      tmp_0[11] = windEmulatorStep4_B.INPUT_2_1_1[3];
      tmp_1[3] = 12;
      tmp_0[12] = windEmulatorStep4_B.INPUT_4_1_1[0];
      tmp_0[13] = windEmulatorStep4_B.INPUT_4_1_1[1];
      tmp_0[14] = windEmulatorStep4_B.INPUT_4_1_1[2];
      tmp_0[15] = windEmulatorStep4_B.INPUT_4_1_1[3];
      tmp_1[4] = 16;
      tmp_0[16] = windEmulatorStep4_B.INPUT_3_1_1[0];
      tmp_0[17] = windEmulatorStep4_B.INPUT_3_1_1[1];
      tmp_0[18] = windEmulatorStep4_B.INPUT_3_1_1[2];
      tmp_0[19] = windEmulatorStep4_B.INPUT_3_1_1[3];
      tmp_1[5] = 20;
      simulationData->mData->mInputValues.mN = 20;
      simulationData->mData->mInputValues.mX = &tmp_0[0];
      simulationData->mData->mInputOffsets.mN = 6;
      simulationData->mData->mInputOffsets.mX = &tmp_1[0];
      simulator = static_cast<NeslSimulator *>
        (windEmulatorStep4_DW.STATE_1_Simulator);
      diagnosticManager = static_cast<NeuDiagnosticManager *>
        (windEmulatorStep4_DW.STATE_1_DiagMgr);
      diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
      tmp_2 = ne_simulator_method(simulator, NESL_SIM_UPDATE, simulationData,
        diagnosticManager);
      if (tmp_2 != 0) {
        tmp = error_buffer_is_empty(rtmGetErrorStatus(windEmulatorStep4_M));
        if (tmp) {
          char *msg;
          msg = rtw_diagnostics_msg(diagnosticTree);
          rtmSetErrorStatus(windEmulatorStep4_M, msg);
        }
      }

      /* End of Update for SimscapeExecutionBlock: '<S222>/STATE_1' */

      /* Update for Memory: '<S176>/Memory' */
      windEmulatorStep4_DW.Memory_PreviousInput_g = windEmulatorStep4_B.Logic_c
        [0];

      /* Update for Memory: '<S177>/Memory' */
      windEmulatorStep4_DW.Memory_PreviousInput_n = windEmulatorStep4_B.Logic_p
        [0];

      /* Update for DiscreteIntegrator: '<S156>/Integrator' */
      windEmulatorStep4_DW.Integrator_DSTATE_e +=
        windEmulatorStep4_cal->Integrator_gainval_b *
        windEmulatorStep4_B.IntegralGain_h;
      windEmulatorStep4_DW.Integrator_PrevResetState_g = static_cast<int8_T>
        (windEmulatorStep4_B.BusAssignment_c.speedCtrlReset);

      /* Update for DiscreteIntegrator: '<S151>/Filter' */
      windEmulatorStep4_DW.Filter_DSTATE_b +=
        windEmulatorStep4_cal->Filter_gainval_d *
        windEmulatorStep4_B.FilterCoefficient_g;
      windEmulatorStep4_DW.Filter_PrevResetState_g = static_cast<int8_T>
        (windEmulatorStep4_B.BusAssignment_c.speedCtrlReset);

      /* Update for DiscreteIntegrator: '<S121>/Discrete-Time Integrator' */
      windEmulatorStep4_DW.DiscreteTimeIntegrator_DSTATE +=
        windEmulatorStep4_cal->DiscreteTimeIntegrator_gainval *
        windEmulatorStep4_B.Sum_f;

      /* Update for DiscreteIntegrator: '<S60>/Discrete-Time Integrator' */
      windEmulatorStep4_DW.DiscreteTimeIntegrator_DSTATE_n +=
        windEmulatorStep4_cal->DiscreteTimeIntegrator_gainva_l *
        windEmulatorStep4_B.Sum_me;

      /* Update for DiscreteIntegrator: '<S119>/Discrete-Time Integrator' */
      windEmulatorStep4_DW.DiscreteTimeIntegrator_DSTATE_l +=
        windEmulatorStep4_cal->DiscreteTimeIntegrator_gainva_b *
        windEmulatorStep4_B.Sum_d;
    }

    /* Update for RateLimiter: '<S120>/Rate Limiter' */
    windEmulatorStep4_DW.PrevY_e = windEmulatorStep4_B.RateLimiter_aw;
    windEmulatorStep4_DW.LastMajorTime_n = windEmulatorStep4_M->Timing.t[0];

    /* Update for RateLimiter: '<S2>/acs880RateLim' */
    windEmulatorStep4_DW.PrevY_b = windEmulatorStep4_B.acs880RateLim;
    windEmulatorStep4_DW.LastMajorTime_d = windEmulatorStep4_M->Timing.t[0];
    if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
      /* Update for Memory: '<S1>/Memory' incorporates:
       *  Constant: '<S1>/powerUpButton'
       */
      windEmulatorStep4_DW.Memory_PreviousInput_d =
        windEmulatorStep4_cal->powerUpButton_Value_b;

      /* Update for Memory: '<S1>/Memory1' incorporates:
       *  Constant: '<S1>/powerDownButton'
       */
      windEmulatorStep4_DW.Memory1_PreviousInput_p =
        windEmulatorStep4_cal->powerDownButton_Value_k;

      /* Update for Memory: '<S1>/Memory2' incorporates:
       *  Constant: '<S1>/resetFaultButton'
       */
      windEmulatorStep4_DW.Memory2_PreviousInput_h =
        windEmulatorStep4_cal->resetFaultButton_Value_n;

      /* Update for Memory: '<S235>/lastRawCounts' */
      windEmulatorStep4_DW.lastRawCounts_PreviousInput =
        windEmulatorStep4_B.CastToDouble_p;

      /* Update for Memory: '<S235>/lastTurn' */
      windEmulatorStep4_DW.lastTurn_PreviousInput = windEmulatorStep4_B.Add1_p;

      /* Update for UnitDelay: '<S236>/UD' */
      windEmulatorStep4_DW.UD_DSTATE = windEmulatorStep4_B.TSamp;

      /* Update for Memory: '<S29>/Memory' */
      windEmulatorStep4_DW.Memory_PreviousInput = windEmulatorStep4_B.Sum_e;
    }
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(windEmulatorStep4_M)) {
    rt_ertODEUpdateContinuousStates(&windEmulatorStep4_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++windEmulatorStep4_M->Timing.clockTick0)) {
      ++windEmulatorStep4_M->Timing.clockTickH0;
    }

    windEmulatorStep4_M->Timing.t[0] = rtsiGetSolverStopTime
      (&windEmulatorStep4_M->solverInfo);

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
      if (!(++windEmulatorStep4_M->Timing.clockTick1)) {
        ++windEmulatorStep4_M->Timing.clockTickH1;
      }

      windEmulatorStep4_M->Timing.t[1] = windEmulatorStep4_M->Timing.clockTick1 *
        windEmulatorStep4_M->Timing.stepSize1 +
        windEmulatorStep4_M->Timing.clockTickH1 *
        windEmulatorStep4_M->Timing.stepSize1 * 4294967296.0;
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void windEmulatorStep4_derivatives(void)
{
  XDot_windEmulatorStep4_T *_rtXdot;
  _rtXdot = ((XDot_windEmulatorStep4_T *) windEmulatorStep4_M->derivs);

  /* Derivatives for Integrator: '<S270>/Integrator' */
  if (!windEmulatorStep4_B.BusAssignment_b.resetSidIntegrator) {
    _rtXdot->Integrator_CSTATE = windEmulatorStep4_B.IProdOut;
  } else {
    /* level reset is active */
    _rtXdot->Integrator_CSTATE = 0.0;
  }

  /* End of Derivatives for Integrator: '<S270>/Integrator' */

  /* Derivatives for StateSpace: '<S211>/Internal' */
  _rtXdot->Internal_CSTATE[0] = 0.0;
  _rtXdot->Internal_CSTATE[1] = 0.0;
  _rtXdot->Internal_CSTATE[2] = 0.0;
  for (uint32_T ri = windEmulatorStep4_cal->Internal_A_jc[0U]; ri <
       windEmulatorStep4_cal->Internal_A_jc[1U]; ri++) {
    _rtXdot->Internal_CSTATE[windEmulatorStep4_cal->Internal_A_ir[ri]] +=
      windEmulatorStep4_cal->Internal_A_pr[ri] *
      windEmulatorStep4_X.Internal_CSTATE[0U];
  }

  for (uint32_T ri = windEmulatorStep4_cal->Internal_A_jc[1U]; ri <
       windEmulatorStep4_cal->Internal_A_jc[2U]; ri++) {
    _rtXdot->Internal_CSTATE[windEmulatorStep4_cal->Internal_A_ir[ri]] +=
      windEmulatorStep4_cal->Internal_A_pr[ri] *
      windEmulatorStep4_X.Internal_CSTATE[1U];
  }

  for (uint32_T ri = windEmulatorStep4_cal->Internal_A_jc[2U]; ri <
       windEmulatorStep4_cal->Internal_A_jc[3U]; ri++) {
    _rtXdot->Internal_CSTATE[windEmulatorStep4_cal->Internal_A_ir[ri]] +=
      windEmulatorStep4_cal->Internal_A_pr[ri] *
      windEmulatorStep4_X.Internal_CSTATE[2U];
  }

  for (uint32_T ri = windEmulatorStep4_cal->Internal_B_jc[0U]; ri <
       windEmulatorStep4_cal->Internal_B_jc[1U]; ri++) {
    _rtXdot->Internal_CSTATE[windEmulatorStep4_cal->Internal_B_ir] +=
      windEmulatorStep4_cal->Internal_B_pr * windEmulatorStep4_B.Sum_a;
  }

  /* End of Derivatives for StateSpace: '<S211>/Internal' */

  /* Derivatives for StateSpace: '<S227>/Internal' */
  _rtXdot->Internal_CSTATE_j = 0.0;
  for (uint32_T ri = windEmulatorStep4_cal->Internal_A_jc_i[0U]; ri <
       windEmulatorStep4_cal->Internal_A_jc_i[1U]; ri++) {
    _rtXdot->Internal_CSTATE_j += windEmulatorStep4_cal->Internal_A_pr_j *
      windEmulatorStep4_X.Internal_CSTATE_j;
  }

  for (uint32_T ri = windEmulatorStep4_cal->Internal_B_jc_i[0U]; ri <
       windEmulatorStep4_cal->Internal_B_jc_i[1U]; ri++) {
    _rtXdot->Internal_CSTATE_j += windEmulatorStep4_cal->Internal_B_pr_g *
      windEmulatorStep4_B.ControlSignal2;
  }

  /* End of Derivatives for StateSpace: '<S227>/Internal' */

  /* Derivatives for StateSpace: '<S223>/Internal' */
  _rtXdot->Internal_CSTATE_a = 0.0;
  for (uint32_T ri = windEmulatorStep4_cal->Internal_A_jc_e[0U]; ri <
       windEmulatorStep4_cal->Internal_A_jc_e[1U]; ri++) {
    _rtXdot->Internal_CSTATE_a += windEmulatorStep4_cal->Internal_A_pr_e *
      windEmulatorStep4_X.Internal_CSTATE_a;
  }

  for (uint32_T ri = windEmulatorStep4_cal->Internal_B_jc_h[0U]; ri <
       windEmulatorStep4_cal->Internal_B_jc_h[1U]; ri++) {
    _rtXdot->Internal_CSTATE_a += windEmulatorStep4_cal->Internal_B_pr_k *
      windEmulatorStep4_B.ControlSignal1;
  }

  /* End of Derivatives for StateSpace: '<S223>/Internal' */
}

/* Model initialize function */
void windEmulatorStep4_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&windEmulatorStep4_M->solverInfo,
                          &windEmulatorStep4_M->Timing.simTimeStep);
    rtsiSetTPtr(&windEmulatorStep4_M->solverInfo, &rtmGetTPtr
                (windEmulatorStep4_M));
    rtsiSetStepSizePtr(&windEmulatorStep4_M->solverInfo,
                       &windEmulatorStep4_M->Timing.stepSize0);
    rtsiSetdXPtr(&windEmulatorStep4_M->solverInfo, &windEmulatorStep4_M->derivs);
    rtsiSetContStatesPtr(&windEmulatorStep4_M->solverInfo, (real_T **)
                         &windEmulatorStep4_M->contStates);
    rtsiSetNumContStatesPtr(&windEmulatorStep4_M->solverInfo,
      &windEmulatorStep4_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&windEmulatorStep4_M->solverInfo,
      &windEmulatorStep4_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&windEmulatorStep4_M->solverInfo,
      &windEmulatorStep4_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&windEmulatorStep4_M->solverInfo,
      &windEmulatorStep4_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&windEmulatorStep4_M->solverInfo, (&rtmGetErrorStatus
      (windEmulatorStep4_M)));
    rtsiSetRTModelPtr(&windEmulatorStep4_M->solverInfo, windEmulatorStep4_M);
  }

  rtsiSetSimTimeStep(&windEmulatorStep4_M->solverInfo, MAJOR_TIME_STEP);
  rtsiSetIsMinorTimeStepWithModeChange(&windEmulatorStep4_M->solverInfo, false);
  windEmulatorStep4_M->intgData.y = windEmulatorStep4_M->odeY;
  windEmulatorStep4_M->intgData.f[0] = windEmulatorStep4_M->odeF[0];
  windEmulatorStep4_M->intgData.f[1] = windEmulatorStep4_M->odeF[1];
  windEmulatorStep4_M->intgData.f[2] = windEmulatorStep4_M->odeF[2];
  windEmulatorStep4_M->intgData.f[3] = windEmulatorStep4_M->odeF[3];
  windEmulatorStep4_M->contStates = ((X_windEmulatorStep4_T *)
    &windEmulatorStep4_X);
  rtsiSetSolverData(&windEmulatorStep4_M->solverInfo, static_cast<void *>
                    (&windEmulatorStep4_M->intgData));
  rtsiSetSolverName(&windEmulatorStep4_M->solverInfo,"ode4");
  windEmulatorStep4_M->solverInfoPtr = (&windEmulatorStep4_M->solverInfo);

  /* Initialize timing info */
  {
    int_T *mdlTsMap = windEmulatorStep4_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;

    /* polyspace +2 MISRA2012:D4.1 [Justified:Low] "windEmulatorStep4_M points to
       static memory which is guaranteed to be non-NULL" */
    windEmulatorStep4_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    windEmulatorStep4_M->Timing.sampleTimes =
      (&windEmulatorStep4_M->Timing.sampleTimesArray[0]);
    windEmulatorStep4_M->Timing.offsetTimes =
      (&windEmulatorStep4_M->Timing.offsetTimesArray[0]);

    /* task periods */
    windEmulatorStep4_M->Timing.sampleTimes[0] = (0.0);
    windEmulatorStep4_M->Timing.sampleTimes[1] = (0.004);

    /* task offsets */
    windEmulatorStep4_M->Timing.offsetTimes[0] = (0.0);
    windEmulatorStep4_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(windEmulatorStep4_M, &windEmulatorStep4_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = windEmulatorStep4_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    windEmulatorStep4_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(windEmulatorStep4_M, -1);
  windEmulatorStep4_M->Timing.stepSize0 = 0.004;
  windEmulatorStep4_M->Timing.stepSize1 = 0.004;
  windEmulatorStep4_M->solverInfoPtr = (&windEmulatorStep4_M->solverInfo);
  windEmulatorStep4_M->Timing.stepSize = (0.004);
  rtsiSetFixedStepSize(&windEmulatorStep4_M->solverInfo, 0.004);
  rtsiSetSolverMode(&windEmulatorStep4_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  (void) std::memset((static_cast<void *>(&windEmulatorStep4_B)), 0,
                     sizeof(B_windEmulatorStep4_T));

  {
    windEmulatorStep4_B.state_j = abbStateEnum_undefined;
    windEmulatorStep4_B.state_m = abbStateEnum_undefined;
    windEmulatorStep4_B.expType_a = expTypeEnum_off;
    windEmulatorStep4_B.toExpTypeEnum = expTypeEnum_off;
  }

  /* states (continuous) */
  {
    (void) std::memset(static_cast<void *>(&windEmulatorStep4_X), 0,
                       sizeof(X_windEmulatorStep4_T));
  }

  /* states (dwork) */
  (void) std::memset(static_cast<void *>(&windEmulatorStep4_DW), 0,
                     sizeof(DW_windEmulatorStep4_T));

  /* external inputs */
  (void)std::memset(&windEmulatorStep4_U, 0, sizeof(ExtU_windEmulatorStep4_T));

  /* child S-Function registration */
  {
    RTWSfcnInfo *sfcnInfo = &windEmulatorStep4_M->NonInlinedSFcns.sfcnInfo;
    windEmulatorStep4_M->sfcnInfo = (sfcnInfo);
    rtssSetErrorStatusPtr(sfcnInfo, (&rtmGetErrorStatus(windEmulatorStep4_M)));
    windEmulatorStep4_M->Sizes.numSampTimes = (2);
    rtssSetNumRootSampTimesPtr(sfcnInfo,
      &windEmulatorStep4_M->Sizes.numSampTimes);
    windEmulatorStep4_M->NonInlinedSFcns.taskTimePtrs[0] = &(rtmGetTPtr
      (windEmulatorStep4_M)[0]);
    windEmulatorStep4_M->NonInlinedSFcns.taskTimePtrs[1] = &(rtmGetTPtr
      (windEmulatorStep4_M)[1]);
    rtssSetTPtrPtr(sfcnInfo,windEmulatorStep4_M->NonInlinedSFcns.taskTimePtrs);
    rtssSetTStartPtr(sfcnInfo, &rtmGetTStart(windEmulatorStep4_M));
    rtssSetTFinalPtr(sfcnInfo, &rtmGetTFinal(windEmulatorStep4_M));
    rtssSetTimeOfLastOutputPtr(sfcnInfo, &rtmGetTimeOfLastOutput
      (windEmulatorStep4_M));
    rtssSetStepSizePtr(sfcnInfo, &windEmulatorStep4_M->Timing.stepSize);
    rtssSetStopRequestedPtr(sfcnInfo, &rtmGetStopRequested(windEmulatorStep4_M));
    rtssSetDerivCacheNeedsResetPtr(sfcnInfo,
      &windEmulatorStep4_M->derivCacheNeedsReset);
    rtssSetZCCacheNeedsResetPtr(sfcnInfo,
      &windEmulatorStep4_M->zCCacheNeedsReset);
    rtssSetContTimeOutputInconsistentWithStateAtMajorStepPtr(sfcnInfo,
      &windEmulatorStep4_M->CTOutputIncnstWithState);
    rtssSetSampleHitsPtr(sfcnInfo, &windEmulatorStep4_M->Timing.sampleHits);
    rtssSetPerTaskSampleHitsPtr(sfcnInfo,
      &windEmulatorStep4_M->Timing.perTaskSampleHits);
    rtssSetSimModePtr(sfcnInfo, &windEmulatorStep4_M->simMode);
    rtssSetSolverInfoPtr(sfcnInfo, &windEmulatorStep4_M->solverInfoPtr);
  }

  windEmulatorStep4_M->Sizes.numSFcns = (1);

  /* register each child */
  {
    (void) std::memset(static_cast<void *>
                       (&windEmulatorStep4_M->NonInlinedSFcns.childSFunctions[0]),
                       0,
                       1*sizeof(SimStruct));
    windEmulatorStep4_M->childSfunctions =
      (&windEmulatorStep4_M->NonInlinedSFcns.childSFunctionPtrs[0]);
    windEmulatorStep4_M->childSfunctions[0] =
      (&windEmulatorStep4_M->NonInlinedSFcns.childSFunctions[0]);

    /* Level2 S-Function Block: windEmulatorStep4/<S5>/Enable File Log (slrealtimeenablelogging) */
    {
      SimStruct *rts = windEmulatorStep4_M->childSfunctions[0];

      /* timing info */
      time_T *sfcnPeriod = windEmulatorStep4_M->NonInlinedSFcns.Sfcn0.sfcnPeriod;
      time_T *sfcnOffset = windEmulatorStep4_M->NonInlinedSFcns.Sfcn0.sfcnOffset;
      int_T *sfcnTsMap = windEmulatorStep4_M->NonInlinedSFcns.Sfcn0.sfcnTsMap;
      (void) std::memset(static_cast<void*>(sfcnPeriod), 0,
                         sizeof(time_T)*1);
      (void) std::memset(static_cast<void*>(sfcnOffset), 0,
                         sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      {
        ssSetBlkInfo2Ptr(rts, &windEmulatorStep4_M->NonInlinedSFcns.blkInfo2[0]);
      }

      _ssSetBlkInfo2PortInfo2Ptr(rts,
        &windEmulatorStep4_M->NonInlinedSFcns.inputOutputPortInfo2[0]);

      /* Set up the mdlInfo pointer */
      ssSetRTWSfcnInfo(rts, windEmulatorStep4_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts, &windEmulatorStep4_M->NonInlinedSFcns.methods2[0]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts, &windEmulatorStep4_M->NonInlinedSFcns.methods3[0]);
      }

      /* Allocate memory of model methods 4 */
      {
        ssSetModelMethods4(rts, &windEmulatorStep4_M->NonInlinedSFcns.methods4[0]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts, &windEmulatorStep4_M->NonInlinedSFcns.statesInfo2
                         [0]);
        ssSetPeriodicStatesInfo(rts,
          &windEmulatorStep4_M->NonInlinedSFcns.periodicStatesInfo[0]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &windEmulatorStep4_M->NonInlinedSFcns.Sfcn0.inputPortInfo[0]);
        ssSetPortInfoForInputs(rts,
          &windEmulatorStep4_M->NonInlinedSFcns.Sfcn0.inputPortInfo[0]);
        _ssSetPortInfo2ForInputUnits(rts,
          &windEmulatorStep4_M->NonInlinedSFcns.Sfcn0.inputPortUnits[0]);
        ssSetInputPortUnit(rts, 0, 0);
        _ssSetPortInfo2ForInputCoSimAttribute(rts,
          &windEmulatorStep4_M->NonInlinedSFcns.Sfcn0.inputPortCoSimAttribute[0]);
        ssSetInputPortIsContinuousQuantity(rts, 0, 0);

        /* port 0 */
        {
          ssSetInputPortRequiredContiguous(rts, 0, 1);
          ssSetInputPortSignal(rts, 0, &windEmulatorStep4_B.Constant);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidthAsInt(rts, 0, 1);
        }
      }

      /* path info */
      ssSetModelName(rts, "Enable File Log");
      ssSetPath(rts, "windEmulatorStep4/fileAndUI/Enable File Log");
      ssSetRTModel(rts,windEmulatorStep4_M);
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
    NeslRtpManager *manager;
    NeslSimulationData *tmp;
    NeslSimulator *simulator;
    NeuDiagnosticManager *diagnosticManager;
    NeuDiagnosticTree *diagnosticTree;
    char *msg;
    real_T tmp_0;
    int32_T tmp_1;
    boolean_T tmp_2;
    boolean_T zcDisabled;

    /* Start for Constant: '<S2>/ACS880CtrlMode' */
    tmp_2 = *get_ctrlModeTorque();

    /* Start for S-Function (slecatinit): '<Root>/EtherCAT Init' */
    slrealtime::StartCallbackService::registerCB( std::bind
      ( Root_EtherCATInit_callback, nullptr ), 10 );

    /* Start for ToAsyncQueueBlock generated from: '<S27>/acs880Signals' */
    windEmulatorStep4_DW.TAQSigLogging_InsertedFor_acs88.SLRTSigHandles =
      slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
      (&windEmulatorStep4_B.BusAssignment_a));

    /* Start for Constant: '<S2>/ACS880CtrlMode' */
    windEmulatorStep4_B.ACS880CtrlMode = tmp_2;

    /* Start for Constant: '<S4>/expType' */
    windEmulatorStep4_B.expType_a = windEmulatorStep4_cal->expType_Value;

    /* Start for Constant: '<S4>/expRunTime' */
    windEmulatorStep4_B.expRunTime = windEmulatorStep4_cal->expRunTime_Value;

    /* Start for SimscapeRtp: '<S188>/RTP_1' */
    manager = nesl_lease_rtp_manager(
      "windEmulatorStep4/hptoSim/hptoModel/HPTO/Solver Configuration_1", 0);
    zcDisabled = pointer_is_null(manager);
    if (zcDisabled) {
      windEmulatorStep4_1e9c788f_1_gateway();
      manager = nesl_lease_rtp_manager(
        "windEmulatorStep4/hptoSim/hptoModel/HPTO/Solver Configuration_1", 0);
    }

    windEmulatorStep4_DW.RTP_1_RtpManager = (void *)manager;
    windEmulatorStep4_DW.RTP_1_SetParametersNeeded = true;

    /* End of Start for SimscapeRtp: '<S188>/RTP_1' */

    /* Start for SimscapeExecutionBlock: '<S222>/STATE_1' */
    simulator = nesl_lease_simulator(
      "windEmulatorStep4/hptoSim/hptoModel/HPTO/Solver Configuration_1", 0, 0);
    windEmulatorStep4_DW.STATE_1_Simulator = (void *)simulator;
    zcDisabled = pointer_is_null(windEmulatorStep4_DW.STATE_1_Simulator);
    if (zcDisabled) {
      windEmulatorStep4_1e9c788f_1_gateway();
      simulator = nesl_lease_simulator(
        "windEmulatorStep4/hptoSim/hptoModel/HPTO/Solver Configuration_1", 0, 0);
      windEmulatorStep4_DW.STATE_1_Simulator = (void *)simulator;
    }

    tmp = nesl_create_simulation_data();
    windEmulatorStep4_DW.STATE_1_SimData = (void *)tmp;
    diagnosticManager = rtw_create_diagnostics();
    windEmulatorStep4_DW.STATE_1_DiagMgr = (void *)diagnosticManager;
    modelParameters.mSolverType = NE_SOLVER_TYPE_ODE;
    modelParameters.mSolverAbsTol = 0.001;
    modelParameters.mSolverRelTol = 0.001;
    modelParameters.mSolverModifyAbsTol = NE_MODIFY_ABS_TOL_NO;
    modelParameters.mStartTime = 0.0;
    modelParameters.mLoadInitialState = false;
    modelParameters.mUseSimState = false;
    modelParameters.mLinTrimCompile = false;
    modelParameters.mLoggingMode = SSC_LOGGING_NONE;
    modelParameters.mRTWModifiedTimeStamp = 7.04312912E+8;
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
      (windEmulatorStep4_DW.STATE_1_Simulator);
    diagnosticManager = static_cast<NeuDiagnosticManager *>
      (windEmulatorStep4_DW.STATE_1_DiagMgr);
    diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
    tmp_1 = nesl_initialize_simulator(simulator, &modelParameters,
      diagnosticManager);
    if (tmp_1 != 0) {
      zcDisabled = error_buffer_is_empty(rtmGetErrorStatus(windEmulatorStep4_M));
      if (zcDisabled) {
        msg = rtw_diagnostics_msg(diagnosticTree);
        rtmSetErrorStatus(windEmulatorStep4_M, msg);
      }
    }

    /* End of Start for SimscapeExecutionBlock: '<S222>/STATE_1' */

    /* Start for SimscapeExecutionBlock: '<S222>/OUTPUT_1_0' */
    simulator = nesl_lease_simulator(
      "windEmulatorStep4/hptoSim/hptoModel/HPTO/Solver Configuration_1", 1, 0);
    windEmulatorStep4_DW.OUTPUT_1_0_Simulator = (void *)simulator;
    zcDisabled = pointer_is_null(windEmulatorStep4_DW.OUTPUT_1_0_Simulator);
    if (zcDisabled) {
      windEmulatorStep4_1e9c788f_1_gateway();
      simulator = nesl_lease_simulator(
        "windEmulatorStep4/hptoSim/hptoModel/HPTO/Solver Configuration_1", 1, 0);
      windEmulatorStep4_DW.OUTPUT_1_0_Simulator = (void *)simulator;
    }

    tmp = nesl_create_simulation_data();
    windEmulatorStep4_DW.OUTPUT_1_0_SimData = (void *)tmp;
    diagnosticManager = rtw_create_diagnostics();
    windEmulatorStep4_DW.OUTPUT_1_0_DiagMgr = (void *)diagnosticManager;
    modelParameters_0.mSolverType = NE_SOLVER_TYPE_ODE;
    modelParameters_0.mSolverAbsTol = 0.001;
    modelParameters_0.mSolverRelTol = 0.001;
    modelParameters_0.mSolverModifyAbsTol = NE_MODIFY_ABS_TOL_NO;
    modelParameters_0.mStartTime = 0.0;
    modelParameters_0.mLoadInitialState = false;
    modelParameters_0.mUseSimState = false;
    modelParameters_0.mLinTrimCompile = false;
    modelParameters_0.mLoggingMode = SSC_LOGGING_NONE;
    modelParameters_0.mRTWModifiedTimeStamp = 7.04312912E+8;
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
      (windEmulatorStep4_DW.OUTPUT_1_0_Simulator);
    diagnosticManager = static_cast<NeuDiagnosticManager *>
      (windEmulatorStep4_DW.OUTPUT_1_0_DiagMgr);
    diagnosticTree = neu_diagnostic_manager_get_initial_tree(diagnosticManager);
    tmp_1 = nesl_initialize_simulator(simulator, &modelParameters_0,
      diagnosticManager);
    if (tmp_1 != 0) {
      zcDisabled = error_buffer_is_empty(rtmGetErrorStatus(windEmulatorStep4_M));
      if (zcDisabled) {
        msg = rtw_diagnostics_msg(diagnosticTree);
        rtmSetErrorStatus(windEmulatorStep4_M, msg);
      }
    }

    /* End of Start for SimscapeExecutionBlock: '<S222>/OUTPUT_1_0' */

    /* Start for ToAsyncQueueBlock generated from: '<S25>/acs880CtrlSignals' */
    windEmulatorStep4_DW.TAQSigLogging_InsertedFor_acs_a.SLRTSigHandles =
      slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
      (&windEmulatorStep4_B.BusAssignment_k));

    /* Start for ToAsyncQueueBlock generated from: '<S23>/acs800Signals' */
    windEmulatorStep4_DW.TAQSigLogging_InsertedFor_acs80.SLRTSigHandles =
      slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
      (&windEmulatorStep4_B.BusAssignment_h));

    /* Start for Constant: '<S1>/ACS800CtrlMode' */
    windEmulatorStep4_B.ACS800CtrlMode = tmp_2;

    /* Start for ToAsyncQueueBlock generated from: '<S21>/acs800CtrlSignals' */
    windEmulatorStep4_DW.TAQSigLogging_InsertedFor_acs_l.SLRTSigHandles =
      slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
      (&windEmulatorStep4_B.BusAssignment_kc));

    /* Start for ToAsyncQueueBlock generated from: '<S34>/hptoSignals' */
    windEmulatorStep4_DW.TAQSigLogging_InsertedFor_hptoS.SLRTSigHandles =
      slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
      (&windEmulatorStep4_B.BusAssignment_n));

    /* Start for ToAsyncQueueBlock generated from: '<S32>/hptoCtrl' */
    windEmulatorStep4_DW.TAQSigLogging_InsertedFor_hptoC.SLRTSigHandles =
      slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
      (&windEmulatorStep4_B.BusAssignment_c));

    /* Start for ToAsyncQueueBlock generated from: '<S30>/expCtrlSignals' */
    windEmulatorStep4_DW.TAQSigLogging_InsertedFor_expCt.SLRTSigHandles =
      slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
      (&windEmulatorStep4_B.BusAssignment_b));

    /* Start for ToAsyncQueueBlock generated from: '<S38>/shaftSignals' */
    windEmulatorStep4_DW.TAQSigLogging_InsertedFor_shaft.SLRTSigHandles =
      slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
      (&windEmulatorStep4_B.BusAssignment_g));

    /* Start for ToAsyncQueueBlock generated from: '<S36>/invPowerAcs800' */
    windEmulatorStep4_DW.TAQSigLogging_InsertedFor_invPo.SLRTSigHandles =
      slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
      (&windEmulatorStep4_B.BusAssignment));

    /* Start for ToAsyncQueueBlock generated from: '<S37>/invPowerAcs880' */
    windEmulatorStep4_DW.TAQSigLogging_InsertedFor_inv_p.SLRTSigHandles =
      slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
      (&windEmulatorStep4_B.BusAssignment_l));

    /* Start for ToAsyncQueueBlock generated from: '<S40>/sidInfoSignals' */
    windEmulatorStep4_DW.TAQSigLogging_InsertedFor_sidIn.SLRTSigHandles =
      slrtRegisterSignalToLoggingService(reinterpret_cast<uintptr_t>
      (&windEmulatorStep4_B.BusAssignment_j));

    /* Start for Constant: '<S5>/Constant' */
    windEmulatorStep4_B.Constant = windEmulatorStep4_cal->Constant_Value_m3;

    /* Start for S-Function (slrealtimeenablelogging): '<S5>/Enable File Log' */
    /* Level2 S-Function Block: '<S5>/Enable File Log' (slrealtimeenablelogging) */
    {
      SimStruct *rts = windEmulatorStep4_M->childSfunctions[0];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }
  }

  windEmulatorStep4_PrevZCX.Integrator_Reset_ZCE = UNINITIALIZED_ZCSIG;

  /* InitializeConditions for Memory: '<S2>/Memory' */
  windEmulatorStep4_DW.Memory_PreviousInput_i =
    windEmulatorStep4_cal->Memory_InitialCondition_a;

  /* InitializeConditions for Memory: '<S2>/Memory1' */
  windEmulatorStep4_DW.Memory1_PreviousInput =
    windEmulatorStep4_cal->Memory1_InitialCondition;

  /* InitializeConditions for Memory: '<S2>/Memory2' */
  windEmulatorStep4_DW.Memory2_PreviousInput =
    windEmulatorStep4_cal->Memory2_InitialCondition;

  /* InitializeConditions for Memory: '<S4>/Memory' */
  windEmulatorStep4_DW.Memory_PreviousInput_k =
    windEmulatorStep4_cal->Memory_InitialCondition_d;

  /* InitializeConditions for Memory: '<S4>/Memory1' */
  windEmulatorStep4_DW.Memory1_PreviousInput_d =
    windEmulatorStep4_cal->Memory1_InitialCondition_h;

  /* InitializeConditions for Memory: '<S4>/Memory2' */
  windEmulatorStep4_DW.Memory2_PreviousInput_l =
    windEmulatorStep4_cal->Memory2_InitialCondition_i;

  /* InitializeConditions for RateLimiter: '<S288>/torqueSlewRate' */
  windEmulatorStep4_DW.LastMajorTime = (rtInf);

  /* InitializeConditions for RateLimiter: '<S288>/speedSlewRate' */
  windEmulatorStep4_DW.LastMajorTime_p = (rtInf);

  /* InitializeConditions for Integrator: '<S270>/Integrator' */
  windEmulatorStep4_X.Integrator_CSTATE =
    windEmulatorStep4_cal->PIDController_InitialConditio_o;

  /* InitializeConditions for RateLimiter: '<S7>/Rate Limiter' */
  windEmulatorStep4_DW.PrevY_f = windEmulatorStep4_cal->RateLimiter_IC;

  /* InitializeConditions for Memory: '<S117>/Memory' */
  windEmulatorStep4_DW.Memory_PreviousInput_kk =
    windEmulatorStep4_cal->SRFlipFlop_initial_condition;

  /* InitializeConditions for RateLimiter: '<S61>/Rate Limiter1' */
  windEmulatorStep4_DW.PrevY_m = windEmulatorStep4_cal->RateLimiter1_IC;

  /* InitializeConditions for Memory: '<S118>/Memory' */
  windEmulatorStep4_DW.Memory_PreviousInput_h =
    windEmulatorStep4_cal->SRFlipFlop_initial_condition_d;

  /* InitializeConditions for RateLimiter: '<S61>/Rate Limiter' */
  windEmulatorStep4_DW.LastMajorTime_m = (rtInf);

  /* InitializeConditions for DiscreteIntegrator: '<S97>/Integrator' */
  windEmulatorStep4_DW.Integrator_DSTATE =
    windEmulatorStep4_cal->PIDController_InitialConditio_a;
  windEmulatorStep4_DW.Integrator_PrevResetState = 0;

  /* InitializeConditions for DiscreteIntegrator: '<S92>/Filter' */
  windEmulatorStep4_DW.Filter_DSTATE =
    windEmulatorStep4_cal->PIDController_InitialConditionF;
  windEmulatorStep4_DW.Filter_PrevResetState = 0;

  /* InitializeConditions for RateLimiter: '<S178>/Rate Limiter' */
  windEmulatorStep4_DW.PrevY_k = windEmulatorStep4_cal->RateLimiter_IC_g;

  /* InitializeConditions for StateSpace: '<S211>/Internal' */
  windEmulatorStep4_X.Internal_CSTATE[0] =
    windEmulatorStep4_cal->Internal_InitialCondition;
  windEmulatorStep4_X.Internal_CSTATE[1] =
    windEmulatorStep4_cal->Internal_InitialCondition;
  windEmulatorStep4_X.Internal_CSTATE[2] =
    windEmulatorStep4_cal->Internal_InitialCondition;

  /* InitializeConditions for StateSpace: '<S227>/Internal' */
  windEmulatorStep4_X.Internal_CSTATE_j =
    windEmulatorStep4_cal->Internal_InitialCondition_p;

  /* InitializeConditions for StateSpace: '<S223>/Internal' */
  windEmulatorStep4_X.Internal_CSTATE_a =
    windEmulatorStep4_cal->Internal_InitialCondition_j;

  /* InitializeConditions for Memory: '<S176>/Memory' */
  windEmulatorStep4_DW.Memory_PreviousInput_g =
    windEmulatorStep4_cal->SRFlipFlop_initial_condition_k;

  /* InitializeConditions for RateLimiter: '<S120>/Rate Limiter1' */
  windEmulatorStep4_DW.PrevY_fq = windEmulatorStep4_cal->RateLimiter1_IC_g;

  /* InitializeConditions for Memory: '<S177>/Memory' */
  windEmulatorStep4_DW.Memory_PreviousInput_n =
    windEmulatorStep4_cal->SRFlipFlop_initial_condition_j;

  /* InitializeConditions for RateLimiter: '<S120>/Rate Limiter' */
  windEmulatorStep4_DW.LastMajorTime_n = (rtInf);

  /* InitializeConditions for DiscreteIntegrator: '<S156>/Integrator' */
  windEmulatorStep4_DW.Integrator_DSTATE_e =
    windEmulatorStep4_cal->PIDController_InitialConditio_c;
  windEmulatorStep4_DW.Integrator_PrevResetState_g = 0;

  /* InitializeConditions for DiscreteIntegrator: '<S151>/Filter' */
  windEmulatorStep4_DW.Filter_DSTATE_b =
    windEmulatorStep4_cal->PIDController_InitialConditio_k;
  windEmulatorStep4_DW.Filter_PrevResetState_g = 0;

  /* InitializeConditions for DiscreteIntegrator: '<S121>/Discrete-Time Integrator' */
  windEmulatorStep4_DW.DiscreteTimeIntegrator_DSTATE =
    windEmulatorStep4_cal->DiscreteTimeIntegrator_IC;

  /* InitializeConditions for DiscreteIntegrator: '<S60>/Discrete-Time Integrator' */
  windEmulatorStep4_DW.DiscreteTimeIntegrator_DSTATE_n =
    windEmulatorStep4_cal->DiscreteTimeIntegrator_IC_e;

  /* InitializeConditions for DiscreteIntegrator: '<S119>/Discrete-Time Integrator' */
  windEmulatorStep4_DW.DiscreteTimeIntegrator_DSTATE_l =
    windEmulatorStep4_cal->DiscreteTimeIntegrator_IC_l;

  /* InitializeConditions for RateLimiter: '<S2>/acs880RateLim' */
  windEmulatorStep4_DW.LastMajorTime_d = (rtInf);

  /* InitializeConditions for Memory: '<S1>/Memory' */
  windEmulatorStep4_DW.Memory_PreviousInput_d =
    windEmulatorStep4_cal->Memory_InitialCondition_l;

  /* InitializeConditions for Memory: '<S1>/Memory1' */
  windEmulatorStep4_DW.Memory1_PreviousInput_p =
    windEmulatorStep4_cal->Memory1_InitialCondition_a;

  /* InitializeConditions for Memory: '<S1>/Memory2' */
  windEmulatorStep4_DW.Memory2_PreviousInput_h =
    windEmulatorStep4_cal->Memory2_InitialCondition_p;

  /* InitializeConditions for Memory: '<S235>/lastRawCounts' */
  windEmulatorStep4_DW.lastRawCounts_PreviousInput =
    windEmulatorStep4_cal->lastRawCounts_InitialCondition;

  /* InitializeConditions for Memory: '<S235>/lastTurn' */
  windEmulatorStep4_DW.lastTurn_PreviousInput =
    windEmulatorStep4_cal->lastTurn_InitialCondition;

  /* InitializeConditions for UnitDelay: '<S236>/UD' */
  windEmulatorStep4_DW.UD_DSTATE =
    windEmulatorStep4_cal->posToVel_ICPrevScaledInput;

  /* InitializeConditions for Memory: '<S29>/Memory' */
  windEmulatorStep4_DW.Memory_PreviousInput =
    windEmulatorStep4_cal->Memory_InitialCondition;

  /* SystemInitialize for Chart: '<S18>/ABB Fieldbus Control' */
  windEmulatorStep4_DW.sfEvent_g = windEmulatorStep4_CALL_EVENT;
  windEmulatorStep4_DW.is_active_UpdateControlWord = 0U;
  windEmulatorStep4_DW.is_active_UpdateStateMachine = 0U;
  windEmulatorStep4_DW.is_UpdateStateMachine = windEmulator_IN_NO_ACTIVE_CHILD;
  windEmulatorStep4_DW.temporalCounter_i1_l = 0U;
  windEmulatorStep4_DW.is_active_c7_windEmulatorStep4 = 0U;
  windEmulatorStep4_DW.swRDY_ON = 0.0;
  windEmulatorStep4_DW.swRDY_RUN = 0.0;
  windEmulatorStep4_DW.swRDY_REF = 0.0;
  windEmulatorStep4_DW.swTRIPPED = 0.0;
  windEmulatorStep4_DW.swOFF_2_STA = 0.0;
  windEmulatorStep4_DW.swOFF_3_STA = 0.0;
  windEmulatorStep4_DW.swSWC_ON_INHIB = 0.0;
  windEmulatorStep4_DW.swAT_SETPOINT = 0.0;
  windEmulatorStep4_DW.swEXT_RUN_ENABLE = 0.0;
  windEmulatorStep4_DW.cwOFF2_CONTROL = 0.0;
  windEmulatorStep4_DW.cwOFF3_CONTROL = 0.0;
  windEmulatorStep4_DW.cwENABLE_OPERATION = 0.0;
  windEmulatorStep4_DW.cwRAMP_OUT_ZERO = 0.0;
  windEmulatorStep4_DW.cwRAMP_HOLD = 0.0;
  windEmulatorStep4_DW.cwRAMP_IN_ZERO = 0.0;
  windEmulatorStep4_DW.cwRESET = 0.0;
  windEmulatorStep4_DW.cwREMOTE_CMD = 0.0;
  windEmulatorStep4_DW.cwOFF1_CONTROL = 0.0;
  windEmulatorStep4_DW.swEXT_CTRL_LOC = 0.0;
  windEmulatorStep4_DW.swREMOTE = 0.0;
  windEmulatorStep4_DW.swWARNING = 0.0;
  windEmulatorStep4_DW.swABOVE_LIMIT = 0.0;
  windEmulatorStep4_DW.swMSW_B13 = 0.0;
  windEmulatorStep4_DW.swMSW_B14 = 0.0;
  windEmulatorStep4_DW.swCOMM_ERR = 0.0;
  windEmulatorStep4_B.ControlWord = 0U;
  windEmulatorStep4_B.state_j = abbStateEnum_undefined;

  /* SystemInitialize for Chart: '<S4>/FexcRamp' */
  windEmulatorStep4_DW.sfEvent = windEmulatorStep4_CALL_EVENT;
  windEmulatorStep4_DW.temporalCounter_i1 = 0U;
  windEmulatorStep4_DW.is_active_c3_windEmulatorStep4 = 0U;
  windEmulatorStep4_DW.is_c3_windEmulatorStep4 = windEmulator_IN_NO_ACTIVE_CHILD;
  windEmulatorStep4_DW.rampLast = 0.0;
  windEmulatorStep4_DW.rampUpTime = 0.0;
  windEmulatorStep4_DW.runTime = 0.0;
  windEmulatorStep4_B.time_d = 0.0;
  windEmulatorStep4_B.ramp_i = 0.0;
  windEmulatorStep4_B.runCounter_m = 0U;
  windEmulatorStep4_B.stepCounter_a = 0U;
  windEmulatorStep4_B.resetHilIntegrator_i = false;
  windEmulatorStep4_B.resetSidIntegrator_i = false;

  /* SystemInitialize for Chart: '<S16>/ABB Fieldbus Control' */
  windEmulatorStep4_DW.sfEvent_f = windEmulatorStep4_CALL_EVENT;
  windEmulatorStep4_DW.is_active_UpdateControlWord_c = 0U;
  windEmulatorStep4_DW.is_active_UpdateStateMachine_o = 0U;
  windEmulatorStep4_DW.is_UpdateStateMachine_l = windEmulator_IN_NO_ACTIVE_CHILD;
  windEmulatorStep4_DW.temporalCounter_i1_ly = 0U;
  windEmulatorStep4_DW.is_active_c9_windEmulatorStep4 = 0U;
  windEmulatorStep4_DW.swRDY_ON_k = 0.0;
  windEmulatorStep4_DW.swRDY_RUN_i = 0.0;
  windEmulatorStep4_DW.swRDY_REF_o = 0.0;
  windEmulatorStep4_DW.swTRIPPED_h = 0.0;
  windEmulatorStep4_DW.swOFF_2_STA_e = 0.0;
  windEmulatorStep4_DW.swOFF_3_STA_a = 0.0;
  windEmulatorStep4_DW.swSWC_ON_INHIB_a = 0.0;
  windEmulatorStep4_DW.swAT_SETPOINT_b = 0.0;
  windEmulatorStep4_DW.swEXT_RUN_ENABLE_g = 0.0;
  windEmulatorStep4_DW.cwOFF2_CONTROL_f = 0.0;
  windEmulatorStep4_DW.cwOFF3_CONTROL_c = 0.0;
  windEmulatorStep4_DW.cwENABLE_OPERATION_k = 0.0;
  windEmulatorStep4_DW.cwRAMP_OUT_ZERO_a = 0.0;
  windEmulatorStep4_DW.cwRAMP_HOLD_n = 0.0;
  windEmulatorStep4_DW.cwRAMP_IN_ZERO_k = 0.0;
  windEmulatorStep4_DW.cwRESET_b = 0.0;
  windEmulatorStep4_DW.cwREMOTE_CMD_a = 0.0;
  windEmulatorStep4_DW.cwOFF1_CONTROL_a = 0.0;
  windEmulatorStep4_DW.swEXT_CTRL_LOC_m = 0.0;
  windEmulatorStep4_DW.swREMOTE_l = 0.0;
  windEmulatorStep4_DW.swWARNING_j = 0.0;
  windEmulatorStep4_DW.swABOVE_LIMIT_c = 0.0;
  windEmulatorStep4_DW.swMSW_B13_o = 0.0;
  windEmulatorStep4_DW.swMSW_B14_d = 0.0;
  windEmulatorStep4_DW.swCOMM_ERR_g = 0.0;
  windEmulatorStep4_B.ControlWord_c = 0U;
  windEmulatorStep4_B.state_m = abbStateEnum_undefined;
  windEmulator_MovingAverage_Init(&windEmulatorStep4_DW.MovingAverage_p);
  windEmulator_MovingAverage_Init(&windEmulatorStep4_DW.MovingAverage);
  windEmulator_MovingAverage_Init(&windEmulatorStep4_DW.MovingAverage_pn);
  windEmulator_MovingAverage_Init(&windEmulatorStep4_DW.MovingAverage1);
}

/* Model terminate function */
void windEmulatorStep4_terminate(void)
{
  NeslSimulationData *simulationData;
  NeuDiagnosticManager *diagnosticManager;
  windEmulator_MovingAverage_Term(&windEmulatorStep4_DW.MovingAverage_p);

  /* Terminate for SimscapeExecutionBlock: '<S222>/STATE_1' */
  diagnosticManager = static_cast<NeuDiagnosticManager *>
    (windEmulatorStep4_DW.STATE_1_DiagMgr);
  neu_destroy_diagnostic_manager(diagnosticManager);
  simulationData = static_cast<NeslSimulationData *>
    (windEmulatorStep4_DW.STATE_1_SimData);
  nesl_destroy_simulation_data(simulationData);
  nesl_erase_simulator("windEmulatorStep4/hptoSim/hptoModel/HPTO/Solver Configuration_1");
  nesl_destroy_registry();

  /* Terminate for SimscapeExecutionBlock: '<S222>/OUTPUT_1_0' */
  diagnosticManager = static_cast<NeuDiagnosticManager *>
    (windEmulatorStep4_DW.OUTPUT_1_0_DiagMgr);
  neu_destroy_diagnostic_manager(diagnosticManager);
  simulationData = static_cast<NeslSimulationData *>
    (windEmulatorStep4_DW.OUTPUT_1_0_SimData);
  nesl_destroy_simulation_data(simulationData);
  nesl_erase_simulator("windEmulatorStep4/hptoSim/hptoModel/HPTO/Solver Configuration_1");
  nesl_destroy_registry();
  windEmulator_MovingAverage_Term(&windEmulatorStep4_DW.MovingAverage);
  windEmulator_MovingAverage_Term(&windEmulatorStep4_DW.MovingAverage_pn);
  windEmulator_MovingAverage_Term(&windEmulatorStep4_DW.MovingAverage1);

  /* Terminate for S-Function (slrealtimeenablelogging): '<S5>/Enable File Log' */
  /* Level2 S-Function Block: '<S5>/Enable File Log' (slrealtimeenablelogging) */
  {
    SimStruct *rts = windEmulatorStep4_M->childSfunctions[0];
    sfcnTerminate(rts);
  }

  /* user code (Terminate function Trailer) */

  /*------------ S-Function Block: <Root>/EtherCAT Init Process Shutdown Network ------------*/
  {
    int_T status;
    status = xpcEtherCATstop(0, 1000 );/* 1 second timeout */
  }
}
