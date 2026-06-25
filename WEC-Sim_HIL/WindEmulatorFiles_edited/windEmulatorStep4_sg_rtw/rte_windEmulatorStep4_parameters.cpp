#include "rte_windEmulatorStep4_parameters.h"
#include "windEmulatorStep4.h"
#include "windEmulatorStep4_cal.h"

RTE_Param_Service_T RTE_Param_Service = {
  {
    0,
    0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0,
    0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0,
    0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0
  },

  {
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0
  },

  {
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0U
  },

  {
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0U
  },

  {
    0.0,
    0U,
    0,
    0.0,
    0.0,
    0U,
    0U
  },

  {
    0.0,
    0.0,
    0.0,
    0
  },

  {
    0U,
    0,
    0,
    0,
    0,
    0.0,
    0.0,
    0U,
    0U
  },

  {
    0.0,
    0.0,
    0U,
    0.0
  },

  {
    0U,
    0,
    0.0,
    0.0
  },

  {
    0U,
    0,
    0.0,
    0.0
  },

  {
    0.2,
    0.0
  },

  {
    0.7,
    0.3
  },

  {
    0.05,
    0.0
  },

  {
    0.0,
    0.0
  },
  1.0,
  40.0,
  4.0,
  75.0,
  0.004,
  5.9921124526782858E-6,
  1.0,
  0.01,
  14.9109,
  1500.0,
  20000.0,
  0.075,
  0.1,
  100.0,
  10000.0,
  0.008035933,
  0.003,
  0.01,
  1.0,
  7.457,
  -10000.0,
  10000.0,
  80.494,
  -78.48165,
  78.48165,
  -78.48165,
  78.48165,
  0.1,
  10000.0,
  0.0080494,
  100.0,
  100000.0,
  0.005,
  300.0,
  250.0,
  25.0,
  -9.3132257461547852E-8,
  78.48165,
  1.6666666666666667E-5,
  1000.0,
  43.893388185128963,
  0.068947572931783066,
  9.5492965855137211,
  5.0,
  0.10471975511965977,
  209.43951023931956,
  true,
  false
};

RTE_Param_Service_T *RTE_Param_Service_ptr = &RTE_Param_Service;
invPowerBus* get_invPowerStruct(void)
{
  return &RTE_Param_Service_ptr->invPowerStruct;
}

hptoSignalBus* get_hptoSignalStruct(void)
{
  return &RTE_Param_Service_ptr->hptoSignalStruct;
}

acs800SignalBus* get_acs800SignalStruct(void)
{
  return &RTE_Param_Service_ptr->acs800SignalStruct;
}

acs880SignalBus* get_acs880SignalStruct(void)
{
  return &RTE_Param_Service_ptr->acs880SignalStruct;
}

shaftSignalBus* get_shaftSignalStruct(void)
{
  return &RTE_Param_Service_ptr->shaftSignalStruct;
}

hptoCtrlBus* get_hptoCtrlStruct(void)
{
  return &RTE_Param_Service_ptr->hptoCtrlStruct;
}

expCtrlBus* get_expCtrlStruct(void)
{
  return &RTE_Param_Service_ptr->expCtrlStruct;
}

sidInfoBus* get_sidInfoStruct(void)
{
  return &RTE_Param_Service_ptr->sidInfoStruct;
}

acs800CtrlBus* get_acs800CtrlStruct(void)
{
  return &RTE_Param_Service_ptr->acs800CtrlStruct;
}

acs880CtrlBus* get_acs880CtrlStruct(void)
{
  return &RTE_Param_Service_ptr->acs880CtrlStruct;
}

struct_KY1U3Kyrwv5e6VnUIBWG5G* get_PressureControl(void)
{
  return &RTE_Param_Service_ptr->PressureControl;
}

struct_KY1U3Kyrwv5e6VnUIBWG5G* get_SpeedControl(void)
{
  return &RTE_Param_Service_ptr->SpeedControl;
}

struct_KY1U3Kyrwv5e6VnUIBWG5G* get_TorqueInputControl(void)
{
  return &RTE_Param_Service_ptr->TorqueInputControl;
}

sidCtrlBus* get_sidCtrlStruct(void)
{
  return &RTE_Param_Service_ptr->sidCtrlStruct;
}

real_T* get_A(void)
{
  return &RTE_Param_Service_ptr->A;
}

real_T* get_Dm_max(void)
{
  return &RTE_Param_Service_ptr->Dm_max;
}

real_T* get_T(void)
{
  return &RTE_Param_Service_ptr->T;
}

real_T* get_TorqueLoadMax(void)
{
  return &RTE_Param_Service_ptr->TorqueLoadMax;
}

real_T* get_Ts(void)
{
  return &RTE_Param_Service_ptr->Ts;
}

real_T* get_absEncoderCountsToRad(void)
{
  return &RTE_Param_Service_ptr->absEncoderCountsToRad;
}

real_T* get_acs800DcBusVoltsScaling(void)
{
  return &RTE_Param_Service_ptr->acs800DcBusVoltsScaling;
}

real_T* get_acs800FreqScaling(void)
{
  return &RTE_Param_Service_ptr->acs800FreqScaling;
}

real_T* get_acs800PowerScaling(void)
{
  return &RTE_Param_Service_ptr->acs800PowerScaling;
}

real_T* get_acs800SpeedNomEng(void)
{
  return &RTE_Param_Service_ptr->acs800SpeedNomEng;
}

real_T* get_acs800SpeedNomFb(void)
{
  return &RTE_Param_Service_ptr->acs800SpeedNomFb;
}

real_T* get_acs800SpeedScaling(void)
{
  return &RTE_Param_Service_ptr->acs800SpeedScaling;
}

real_T* get_acs800TempScaling(void)
{
  return &RTE_Param_Service_ptr->acs800TempScaling;
}

real_T* get_acs800TorqueNomEng(void)
{
  return &RTE_Param_Service_ptr->acs800TorqueNomEng;
}

real_T* get_acs800TorqueNomFb(void)
{
  return &RTE_Param_Service_ptr->acs800TorqueNomFb;
}

real_T* get_acs800TorqueScaling(void)
{
  return &RTE_Param_Service_ptr->acs800TorqueScaling;
}

real_T* get_acs880FreqScaling(void)
{
  return &RTE_Param_Service_ptr->acs880FreqScaling;
}

real_T* get_acs880MotorCurrentScaling(void)
{
  return &RTE_Param_Service_ptr->acs880MotorCurrentScaling;
}

real_T* get_acs880MotorVoltsScaling(void)
{
  return &RTE_Param_Service_ptr->acs880MotorVoltsScaling;
}

real_T* get_acs880PowerScaling(void)
{
  return &RTE_Param_Service_ptr->acs880PowerScaling;
}

real_T* get_acs880RateLimFalling(void)
{
  return &RTE_Param_Service_ptr->acs880RateLimFalling;
}

real_T* get_acs880RateLimRising(void)
{
  return &RTE_Param_Service_ptr->acs880RateLimRising;
}

real_T* get_acs880RatedTorque(void)
{
  return &RTE_Param_Service_ptr->acs880RatedTorque;
}

real_T* get_acs880SetpointLimLower(void)
{
  return &RTE_Param_Service_ptr->acs880SetpointLimLower;
}

real_T* get_acs880SetpointLimUpper(void)
{
  return &RTE_Param_Service_ptr->acs880SetpointLimUpper;
}

real_T* get_acs880SpeedPILimLo(void)
{
  return &RTE_Param_Service_ptr->acs880SpeedPILimLo;
}

real_T* get_acs880SpeedPILimUp(void)
{
  return &RTE_Param_Service_ptr->acs880SpeedPILimUp;
}

real_T* get_acs880SpeedScaling(void)
{
  return &RTE_Param_Service_ptr->acs880SpeedScaling;
}

real_T* get_acs880TorqueFieldbusScale(void)
{
  return &RTE_Param_Service_ptr->acs880TorqueFieldbusScale;
}

real_T* get_acs880TorqueScaling(void)
{
  return &RTE_Param_Service_ptr->acs880TorqueScaling;
}

real_T* get_acs880TorqueSetpointScaling(void)
{
  return &RTE_Param_Service_ptr->acs880TorqueSetpointScaling;
}

real_T* get_bar2pa(void)
{
  return &RTE_Param_Service_ptr->bar2pa;
}

real_T* get_belowMinPGain(void)
{
  return &RTE_Param_Service_ptr->belowMinPGain;
}

real_T* get_deadbandTorqueSlewRate(void)
{
  return &RTE_Param_Service_ptr->deadbandTorqueSlewRate;
}

real_T* get_fromFileSpeedSlewRate(void)
{
  return &RTE_Param_Service_ptr->fromFileSpeedSlewRate;
}

real_T* get_fromFileTorqueSlewRate(void)
{
  return &RTE_Param_Service_ptr->fromFileTorqueSlewRate;
}

real_T* get_futekTorqueScale(void)
{
  return &RTE_Param_Service_ptr->futekTorqueScale;
}

real_T* get_genMaxTorque(void)
{
  return &RTE_Param_Service_ptr->genMaxTorque;
}

real_T* get_lpm2m3ps(void)
{
  return &RTE_Param_Service_ptr->lpm2m3ps;
}

real_T* get_minPressureRef_psi(void)
{
  return &RTE_Param_Service_ptr->minPressureRef_psi;
}

real_T* get_minTorqueRef_Nm(void)
{
  return &RTE_Param_Service_ptr->minTorqueRef_Nm;
}

real_T* get_psi2bar(void)
{
  return &RTE_Param_Service_ptr->psi2bar;
}

real_T* get_radps2rpm(void)
{
  return &RTE_Param_Service_ptr->radps2rpm;
}

real_T* get_rampTime(void)
{
  return &RTE_Param_Service_ptr->rampTime;
}

real_T* get_rpm2radps(void)
{
  return &RTE_Param_Service_ptr->rpm2radps;
}

real_T* get_targetOmega(void)
{
  return &RTE_Param_Service_ptr->targetOmega;
}

boolean_T* get_ctrlModeTorque(void)
{
  return &RTE_Param_Service_ptr->ctrlModeTorque;
}

boolean_T* get_deadBandController(void)
{
  return &RTE_Param_Service_ptr->deadBandController;
}

extern windEmulatorStep4_cal_type windEmulatorStep4_cal_impl;
extern RTE_Param_Service_T RTE_Param_Service;
namespace slrealtime
{
  /* Description of SEGMENTS */
  SegmentVector segmentInfo {
    { (void*)&RTE_Param_Service, (void**)&RTE_Param_Service_ptr, sizeof
      (RTE_Param_Service_T), 2 },

    { (void*)&windEmulatorStep4_cal_impl, (void**)&windEmulatorStep4_cal, sizeof
      (windEmulatorStep4_cal_type), 2 }
  };

  SegmentVector &getSegmentVector(void)
  {
    return segmentInfo;
  }
}                                      // slrealtime
