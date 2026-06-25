#include "windEmulatorStep4_cal.h"
#include "windEmulatorStep4.h"

/* Storage class 'PageSwitching' */
windEmulatorStep4_cal_type windEmulatorStep4_cal_impl = {
  /* Computed Parameter: Internal_A_pr
   * Referenced by: '<S211>/Internal'
   */
  { -118.13271821199052, 8.0, -11.623239835183748, 8.0, -15.144377718883565 },

  /* Mask Parameter: PIDController_D
   * Referenced by: '<S91>/Derivative Gain'
   */
  0.0,

  /* Mask Parameter: PIDController_D_d
   * Referenced by: '<S150>/Derivative Gain'
   */
  0.0,

  /* Mask Parameter: posToVel_ICPrevScaledInput
   * Referenced by: '<S236>/UD'
   */
  0.0,

  /* Mask Parameter: PIDController_InitialConditionF
   * Referenced by: '<S92>/Filter'
   */
  0.0,

  /* Mask Parameter: PIDController_InitialConditio_k
   * Referenced by: '<S151>/Filter'
   */
  0.0,

  /* Mask Parameter: PIDController_InitialConditio_o
   * Referenced by: '<S270>/Integrator'
   */
  0.0,

  /* Mask Parameter: PIDController_InitialConditio_a
   * Referenced by: '<S97>/Integrator'
   */
  0.0,

  /* Mask Parameter: PIDController_InitialConditio_c
   * Referenced by: '<S156>/Integrator'
   */
  0.0,

  /* Mask Parameter: Ramp_InitialOutput
   * Referenced by: '<S58>/Constant1'
   */
  0.0,

  /* Mask Parameter: PIDController_N
   * Referenced by: '<S100>/Filter Coefficient'
   */
  100.0,

  /* Mask Parameter: PIDController_N_p
   * Referenced by: '<S159>/Filter Coefficient'
   */
  100.0,

  /* Mask Parameter: Ramp_slope
   * Referenced by: '<S58>/Step'
   */
  0.1,

  /* Mask Parameter: Ramp_start
   * Referenced by:
   *   '<S58>/Constant'
   *   '<S58>/Step'
   */
  10.0,

  /* Expression: 0
   * Referenced by: '<S115>/Saturation'
   */
  0.0,

  /* Expression: -1
   * Referenced by: '<S115>/Saturation'
   */
  -1.0,

  /* Expression: 1
   * Referenced by: '<S116>/Saturation'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S116>/Saturation'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S64>/Switch1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S174>/Saturation'
   */
  0.0,

  /* Expression: -1
   * Referenced by: '<S174>/Saturation'
   */
  -1.0,

  /* Expression: 1
   * Referenced by: '<S175>/Saturation'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S175>/Saturation'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S123>/Switch1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S62>/Switch'
   */
  0.0,

  /* Expression: -1
   * Referenced by: '<S121>/Gain2'
   */
  -1.0,

  /* Expression: 1
   * Referenced by: '<S121>/Saturation'
   */
  1.0,

  /* Expression: -1
   * Referenced by: '<S121>/Saturation'
   */
  -1.0,

  /* Expression: 1
   * Referenced by: '<S60>/Saturation'
   */
  1.0,

  /* Expression: -1
   * Referenced by: '<S60>/Saturation'
   */
  -1.0,

  /* Expression: 1
   * Referenced by: '<S119>/Saturation'
   */
  1.0,

  /* Expression: -1
   * Referenced by: '<S119>/Saturation'
   */
  -1.0,

  /* Expression: 0
   * Referenced by: '<S7>/Constant1'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S288>/fromFileSpeedNow_rpm'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S288>/fromFileTorqueNow_Nm'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S238>/Constant1'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S48>/shaftPowerAverage_W'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S48>/shaftPower_W'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S28>/frequency_Hz'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S28>/motorCurrent_A'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S28>/motorSpeed_rpm'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S28>/motorTorque_Nm'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S28>/motorVoltage_V'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S28>/shaftPower_W'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S4>/expRunTime'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S3>/zeroTorque'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S238>/Constant3'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S238>/manualTorqueSetpoint_Nm'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S288>/Set bound'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S238>/manualSpeedSetpoint_rpm'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S13>/acs880SpeedPGain'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S7>/speedReference'
   */
  0.0,

  /* Expression: 50
   * Referenced by: '<S7>/Rate Limiter'
   */
  50.0,

  /* Expression: -50
   * Referenced by: '<S7>/Rate Limiter'
   */
  -50.0,

  /* Expression: 0
   * Referenced by: '<S7>/Rate Limiter'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S7>/excForceFreq_Hz'
   */
  0.0,

  /* Expression: 2*pi
   * Referenced by: '<S7>/f->w'
   */
  6.2831853071795862,

  /* Expression: 0
   * Referenced by: '<S7>/excForceAmp_N'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S7>/excForceAmpNow_N'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S64>/Constant1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S61>/shaftSpeedRefMin'
   */
  0.0,

  /* Expression: 50
   * Referenced by: '<S61>/Rate Limiter1'
   */
  50.0,

  /* Expression: -50
   * Referenced by: '<S61>/Rate Limiter1'
   */
  -50.0,

  /* Expression: 0
   * Referenced by: '<S61>/Rate Limiter1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S64>/Switch'
   */
  0.0,

  /* Computed Parameter: Integrator_gainval
   * Referenced by: '<S97>/Integrator'
   */
  0.004,

  /* Computed Parameter: Filter_gainval
   * Referenced by: '<S92>/Filter'
   */
  0.004,

  /* Expression: -1
   * Referenced by: '<S63>/Gain2'
   */
  -1.0,

  /* Expression: 0
   * Referenced by: '<S58>/Step'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S53>/Saturation'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S53>/Saturation'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S59>/kDamping'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S59>/kDampingNow'
   */
  1.0,

  /* Expression: 200
   * Referenced by: '<S178>/Rate Limiter'
   */
  200.0,

  /* Expression: -200
   * Referenced by: '<S178>/Rate Limiter'
   */
  -200.0,

  /* Expression: 0
   * Referenced by: '<S178>/Rate Limiter'
   */
  0.0,

  /* Computed Parameter: Internal_B_pr
   * Referenced by: '<S211>/Internal'
   */
  0.03125,

  /* Computed Parameter: Internal_C_pr
   * Referenced by: '<S211>/Internal'
   */
  0.042279268804058107,

  /* Expression: 0.0
   * Referenced by: '<S211>/Internal'
   */
  0.0,

  /* Computed Parameter: Internal_A_pr_j
   * Referenced by: '<S227>/Internal'
   */
  -628.31853071795865,

  /* Computed Parameter: Internal_B_pr_g
   * Referenced by: '<S227>/Internal'
   */
  32.0,

  /* Computed Parameter: Internal_C_pr_a
   * Referenced by: '<S227>/Internal'
   */
  19.634954084936208,

  /* Expression: 0.0
   * Referenced by: '<S227>/Internal'
   */
  0.0,

  /* Expression: 0.0075
   * Referenced by: '<S190>/Gain'
   */
  0.0075,

  /* Computed Parameter: Internal_A_pr_e
   * Referenced by: '<S223>/Internal'
   */
  -628.31853071795865,

  /* Computed Parameter: Internal_B_pr_k
   * Referenced by: '<S223>/Internal'
   */
  32.0,

  /* Computed Parameter: Internal_C_pr_m
   * Referenced by: '<S223>/Internal'
   */
  19.634954084936208,

  /* Expression: 0.0
   * Referenced by: '<S223>/Internal'
   */
  0.0,

  /* Expression: 0.0075
   * Referenced by: '<S189>/Gain'
   */
  0.0075,

  /* Expression: 5e-6
   * Referenced by: '<S54>/Subsystem_around_RTP_D290B913_fluid_volume'
   */
  5.0E-6,

  /* Expression: 400
   * Referenced by: '<S178>/Subsystem_around_RTP_D2E1D090_liquid_pressure'
   */
  400.0,

  /* Expression: 4.9e-3
   * Referenced by: '<S178>/Subsystem_around_RTP_D2E1D090_liquid_volume'
   */
  0.0049,

  /* Expression: 0
   * Referenced by: '<S59>/kSpring'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S59>/kSpringNow'
   */
  1.0,

  /* Expression: -1
   * Referenced by: '<S55>/Gain'
   */
  -1.0,

  /* Expression: 0
   * Referenced by: '<S123>/Constant1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S120>/shaftSpeedRefMin'
   */
  0.0,

  /* Expression: 50
   * Referenced by: '<S120>/Rate Limiter1'
   */
  50.0,

  /* Expression: -50
   * Referenced by: '<S120>/Rate Limiter1'
   */
  -50.0,

  /* Expression: 0
   * Referenced by: '<S120>/Rate Limiter1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S123>/Switch'
   */
  0.0,

  /* Computed Parameter: Integrator_gainval_b
   * Referenced by: '<S156>/Integrator'
   */
  0.004,

  /* Computed Parameter: Filter_gainval_d
   * Referenced by: '<S151>/Filter'
   */
  0.004,

  /* Expression: -1
   * Referenced by: '<S122>/Gain2'
   */
  -1.0,

  /* Expression: 1/6894.75
   * Referenced by: '<S184>/Gain'
   */
  0.00014503789114906271,

  /* Expression: -1
   * Referenced by: '<S62>/Constant'
   */
  -1.0,

  /* Expression: 1
   * Referenced by: '<S62>/Constant1'
   */
  1.0,

  /* Computed Parameter: DiscreteTimeIntegrator_gainval
   * Referenced by: '<S121>/Discrete-Time Integrator'
   */
  0.004,

  /* Expression: 0
   * Referenced by: '<S121>/Discrete-Time Integrator'
   */
  0.0,

  /* Computed Parameter: DiscreteTimeIntegrator_gainva_l
   * Referenced by: '<S60>/Discrete-Time Integrator'
   */
  0.004,

  /* Expression: 0
   * Referenced by: '<S60>/Discrete-Time Integrator'
   */
  0.0,

  /* Computed Parameter: DiscreteTimeIntegrator_gainva_b
   * Referenced by: '<S119>/Discrete-Time Integrator'
   */
  0.004,

  /* Expression: 0
   * Referenced by: '<S119>/Discrete-Time Integrator'
   */
  0.0,

  /* Expression: 1000*60
   * Referenced by: '<S181>/Gain'
   */
  60000.0,

  /* Expression: 1
   * Referenced by: '<S26>/torqueSetpoint_Nm'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S26>/torqueSetpoint_percent'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S43>/shaftPowerAverage_W'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S43>/shaftPower_W'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S24>/dcBusVoltage_V'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S24>/frequency_Hz'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S24>/motorSpeed_rpm'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S24>/motorTorque_Nm'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S24>/shaftPower_W'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S24>/temperature'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S22>/torqueSetpoint_Nm'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S22>/torqueSetpoint_percent'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S51>/excShaftPowerAverage_W'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S51>/excShaftPower_W'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S51>/hydrPowerAverage_W'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S51>/hydrPower_W'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S35>/ctrlSignal1'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S35>/ctrlSignal2'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S35>/excShaftSpeed_rpm'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S35>/excShaftTorque_Nm'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S35>/genPumpFlow_lpm'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S35>/genShaftSpeed_rpm'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S35>/genTorqueCmd_Nm'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S35>/hmOutputShaftTorque_Nm'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S35>/pressure_bar'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S33>/excForce_N'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S33>/genSpeedActual'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S33>/speedRef_rpm'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S31>/ramp'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S31>/time'
   */
  1.0,

  /* Expression: 2*pi
   * Referenced by: '<S235>/Gain'
   */
  6.2831853071795862,

  /* Computed Parameter: TSamp_WtEt
   * Referenced by: '<S236>/TSamp'
   */
  250.0,

  /* Expression: 1
   * Referenced by: '<S39>/absEncoderPosition_rad'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S39>/absEncoderSpeed_rpm'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S39>/torqueActual_Nm'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S288>/caseCounterSignalsNow'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S14>/Constant'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S15>/Constant2'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S53>/Constant'
   */
  1.0,

  /* Expression: -1
   * Referenced by: '<S53>/Constant1'
   */
  -1.0,

  /* Expression: 1000*60
   * Referenced by: '<S179>/m3toL'
   */
  60000.0,

  /* Expression: 1000*60
   * Referenced by: '<S180>/Gain'
   */
  60000.0,

  /* Expression: 0
   * Referenced by: '<S13>/acs880SpeedIGain'
   */
  0.0,

  /* Expression: 100
   * Referenced by: '<S288>/vecPercent'
   */
  100.0,

  /* Computed Parameter: L1Voltage_Gain
   * Referenced by: '<S9>/L1Voltage'
   */
  1.0F,

  /* Computed Parameter: L1Current_Gain
   * Referenced by: '<S9>/L1Current'
   */
  1.0F,

  /* Computed Parameter: L1PowFactor_Gain
   * Referenced by: '<S9>/L1PowFactor'
   */
  1.0F,

  /* Computed Parameter: L1ActivePow_Gain
   * Referenced by: '<S9>/L1ActivePow'
   */
  1.0F,

  /* Computed Parameter: L1THDu_Gain
   * Referenced by: '<S9>/L1THDu'
   */
  1.0F,

  /* Computed Parameter: L1THDi_Gain
   * Referenced by: '<S9>/L1THDi'
   */
  1.0F,

  /* Computed Parameter: L2Voltage_Gain
   * Referenced by: '<S9>/L2Voltage'
   */
  1.0F,

  /* Computed Parameter: L2Current_Gain
   * Referenced by: '<S9>/L2Current'
   */
  1.0F,

  /* Computed Parameter: L2PowFactor_Gain
   * Referenced by: '<S9>/L2PowFactor'
   */
  1.0F,

  /* Computed Parameter: L2ActivePow_Gain
   * Referenced by: '<S9>/L2ActivePow'
   */
  1.0F,

  /* Computed Parameter: L2THDu_Gain
   * Referenced by: '<S9>/L2THDu'
   */
  1.0F,

  /* Computed Parameter: L2THDi_Gain
   * Referenced by: '<S9>/L2THDi'
   */
  1.0F,

  /* Computed Parameter: L3Voltage_Gain
   * Referenced by: '<S9>/L3Voltage'
   */
  1.0F,

  /* Computed Parameter: L3Current_Gain
   * Referenced by: '<S9>/L3Current'
   */
  1.0F,

  /* Computed Parameter: L3PowFactor_Gain
   * Referenced by: '<S9>/L3PowFactor'
   */
  1.0F,

  /* Computed Parameter: L3ActivePow_Gain
   * Referenced by: '<S9>/L3ActivePow'
   */
  1.0F,

  /* Computed Parameter: L3THDu_Gain
   * Referenced by: '<S9>/L3THDu'
   */
  1.0F,

  /* Computed Parameter: L3THDi_Gain
   * Referenced by: '<S9>/L3THDi'
   */
  1.0F,

  /* Computed Parameter: totalFrequency_Gain
   * Referenced by: '<S9>/totalFrequency'
   */
  1.0F,

  /* Computed Parameter: totalPowFactor_Gain
   * Referenced by: '<S9>/totalPowFactor'
   */
  1.0F,

  /* Computed Parameter: totalActivePow_Gain
   * Referenced by: '<S9>/totalActivePow'
   */
  1.0F,

  /* Computed Parameter: L1L2Voltage_Gain
   * Referenced by: '<S9>/L1L2Voltage'
   */
  1.0F,

  /* Computed Parameter: L2L3Voltage_Gain
   * Referenced by: '<S9>/L2L3Voltage'
   */
  1.0F,

  /* Computed Parameter: L3L1Voltage_Gain
   * Referenced by: '<S9>/L3L1Voltage'
   */
  1.0F,

  /* Computed Parameter: L1Voltage_Gain_h
   * Referenced by: '<S11>/L1Voltage'
   */
  1.0F,

  /* Computed Parameter: L1Current_Gain_f
   * Referenced by: '<S11>/L1Current'
   */
  1.0F,

  /* Computed Parameter: L1PowFactor_Gain_a
   * Referenced by: '<S11>/L1PowFactor'
   */
  1.0F,

  /* Computed Parameter: L1ActivePow_Gain_o
   * Referenced by: '<S11>/L1ActivePow'
   */
  1.0F,

  /* Computed Parameter: L1THDu_Gain_a
   * Referenced by: '<S11>/L1THDu'
   */
  1.0F,

  /* Computed Parameter: L1THDi_Gain_f
   * Referenced by: '<S11>/L1THDi'
   */
  1.0F,

  /* Computed Parameter: L2Voltage_Gain_a
   * Referenced by: '<S11>/L2Voltage'
   */
  1.0F,

  /* Computed Parameter: L2Current_Gain_j
   * Referenced by: '<S11>/L2Current'
   */
  1.0F,

  /* Computed Parameter: L2PowFactor_Gain_p
   * Referenced by: '<S11>/L2PowFactor'
   */
  1.0F,

  /* Computed Parameter: L2ActivePow_Gain_b
   * Referenced by: '<S11>/L2ActivePow'
   */
  1.0F,

  /* Computed Parameter: L2THDu_Gain_i
   * Referenced by: '<S11>/L2THDu'
   */
  1.0F,

  /* Computed Parameter: L2THDi_Gain_j
   * Referenced by: '<S11>/L2THDi'
   */
  1.0F,

  /* Computed Parameter: L3Voltage_Gain_h
   * Referenced by: '<S11>/L3Voltage'
   */
  1.0F,

  /* Computed Parameter: L3Current_Gain_i
   * Referenced by: '<S11>/L3Current'
   */
  1.0F,

  /* Computed Parameter: L3PowFactor_Gain_k
   * Referenced by: '<S11>/L3PowFactor'
   */
  1.0F,

  /* Computed Parameter: L3ActivePow_Gain_f
   * Referenced by: '<S11>/L3ActivePow'
   */
  1.0F,

  /* Computed Parameter: L3THDu_Gain_b
   * Referenced by: '<S11>/L3THDu'
   */
  1.0F,

  /* Computed Parameter: L3THDi_Gain_p
   * Referenced by: '<S11>/L3THDi'
   */
  1.0F,

  /* Computed Parameter: totalFrequency_Gain_b
   * Referenced by: '<S11>/totalFrequency'
   */
  1.0F,

  /* Computed Parameter: totalPowFactor_Gain_f
   * Referenced by: '<S11>/totalPowFactor'
   */
  1.0F,

  /* Computed Parameter: totalActivePow_Gain_m
   * Referenced by: '<S11>/totalActivePow'
   */
  1.0F,

  /* Computed Parameter: L1L2Voltage_Gain_j
   * Referenced by: '<S11>/L1L2Voltage'
   */
  1.0F,

  /* Computed Parameter: L2L3Voltage_Gain_m
   * Referenced by: '<S11>/L2L3Voltage'
   */
  1.0F,

  /* Computed Parameter: L3L1Voltage_Gain_k
   * Referenced by: '<S11>/L3L1Voltage'
   */
  1.0F,

  /* Computed Parameter: Constant1_Value_jt
   * Referenced by: '<S235>/Constant1'
   */
  0,

  /* Computed Parameter: state_Bias
   * Referenced by: '<S26>/state'
   */
  0,

  /* Computed Parameter: state_Bias_l
   * Referenced by: '<S22>/state'
   */
  0,

  /* Computed Parameter: lastRawCounts_InitialCondition
   * Referenced by: '<S235>/lastRawCounts'
   */
  0,

  /* Computed Parameter: Constant_Value_p
   * Referenced by: '<S235>/Constant'
   */
  524288,

  /* Computed Parameter: lastTurn_InitialCondition
   * Referenced by: '<S235>/lastTurn'
   */
  0,

  /* Computed Parameter: absEncoderTurns_Bias
   * Referenced by: '<S39>/absEncoderTurns'
   */
  0,

  /* Computed Parameter: Internal_A_ir
   * Referenced by: '<S211>/Internal'
   */
  { 0U, 1U, 0U, 2U, 0U },

  /* Computed Parameter: Internal_A_jc
   * Referenced by: '<S211>/Internal'
   */
  { 0U, 2U, 4U, 5U },

  /* Computed Parameter: Internal_B_jc
   * Referenced by: '<S211>/Internal'
   */
  { 0U, 1U },

  /* Computed Parameter: Internal_C_jc
   * Referenced by: '<S211>/Internal'
   */
  { 0U, 0U, 1U, 1U },

  /* Computed Parameter: Internal_A_jc_i
   * Referenced by: '<S227>/Internal'
   */
  { 0U, 1U },

  /* Computed Parameter: Internal_B_jc_i
   * Referenced by: '<S227>/Internal'
   */
  { 0U, 1U },

  /* Computed Parameter: Internal_C_jc_k
   * Referenced by: '<S227>/Internal'
   */
  { 0U, 1U },

  /* Computed Parameter: Internal_A_jc_e
   * Referenced by: '<S223>/Internal'
   */
  { 0U, 1U },

  /* Computed Parameter: Internal_B_jc_h
   * Referenced by: '<S223>/Internal'
   */
  { 0U, 1U },

  /* Computed Parameter: Internal_C_jc_h
   * Referenced by: '<S223>/Internal'
   */
  { 0U, 1U },

  /* Computed Parameter: Constant_Value_im
   * Referenced by: '<S288>/Constant'
   */
  1U,

  /* Computed Parameter: Lengthofinput_Value
   * Referenced by: '<S288>/Length of input'
   */
  560000U,

  /* Computed Parameter: Internal_B_ir
   * Referenced by: '<S211>/Internal'
   */
  0U,

  /* Computed Parameter: Internal_C_ir
   * Referenced by: '<S211>/Internal'
   */
  0U,

  /* Computed Parameter: Internal_A_ir_f
   * Referenced by: '<S227>/Internal'
   */
  0U,

  /* Computed Parameter: Internal_B_ir_l
   * Referenced by: '<S227>/Internal'
   */
  0U,

  /* Computed Parameter: Internal_C_ir_n
   * Referenced by: '<S227>/Internal'
   */
  0U,

  /* Computed Parameter: Internal_A_ir_j
   * Referenced by: '<S223>/Internal'
   */
  0U,

  /* Computed Parameter: Internal_B_ir_a
   * Referenced by: '<S223>/Internal'
   */
  0U,

  /* Computed Parameter: Internal_C_ir_e
   * Referenced by: '<S223>/Internal'
   */
  0U,

  /* Computed Parameter: runCounter_Bias
   * Referenced by: '<S31>/runCounter'
   */
  0U,

  /* Computed Parameter: stepCounter_Bias
   * Referenced by: '<S31>/stepCounter'
   */
  0U,

  /* Computed Parameter: absEncoderCounts_Bias
   * Referenced by: '<S39>/absEncoderCounts'
   */
  0U,

  /* Computed Parameter: loopAdd_Value
   * Referenced by: '<S29>/loopAdd'
   */
  1U,

  /* Computed Parameter: Memory_InitialCondition
   * Referenced by: '<S29>/Memory'
   */
  0U,

  /* Computed Parameter: loopCounter_Gain
   * Referenced by: '<S29>/loopCounter'
   */
  1U,

  /* Computed Parameter: fileSamples_Bias
   * Referenced by: '<S288>/fileSamples'
   */
  0U,

  /* Expression: expTypeEnum.off
   * Referenced by: '<S4>/expType'
   */
  expTypeEnum_off,

  /* Expression: expTypeEnum.hil
   * Referenced by: '<S4>/expModeHil'
   */
  expTypeEnum_hil,

  /* Expression: expTypeEnum.sid
   * Referenced by: '<S4>/expModeSid'
   */
  expTypeEnum_sid,

  /* Expression: sidTypeEnum.off
   * Referenced by: '<S238>/sidType'
   */
  sidTypeEnum_off,

  /* Computed Parameter: statusWord_Bias
   * Referenced by: '<S28>/statusWord'
   */
  0U,

  /* Computed Parameter: ctrlWord_Bias
   * Referenced by: '<S26>/ctrlWord'
   */
  0U,

  /* Computed Parameter: statusWord_Bias_l
   * Referenced by: '<S24>/statusWord'
   */
  0U,

  /* Computed Parameter: ctrlWord_Bias_a
   * Referenced by: '<S22>/ctrlWord'
   */
  0U,

  /* Computed Parameter: expType_Bias
   * Referenced by: '<S31>/expType'
   */
  0U,

  /* Computed Parameter: speedCtrlReset_Gain
   * Referenced by: '<S33>/speedCtrlReset'
   */
  128U,

  /* Computed Parameter: absEncoderStatus1_Bias
   * Referenced by: '<S39>/absEncoderStatus1'
   */
  0U,

  /* Computed Parameter: absEncoderStatus2_Bias
   * Referenced by: '<S39>/absEncoderStatus2'
   */
  0U,

  /* Computed Parameter: Logic_table
   * Referenced by: '<S117>/Logic'
   */
  { false, true, false, false, true, true, false, false, true, false, true, true,
    false, false, false, false },

  /* Computed Parameter: Logic_table_o
   * Referenced by: '<S118>/Logic'
   */
  { false, true, false, false, true, true, false, false, true, false, true, true,
    false, false, false, false },

  /* Computed Parameter: Logic_table_h
   * Referenced by: '<S176>/Logic'
   */
  { false, true, false, false, true, true, false, false, true, false, true, true,
    false, false, false, false },

  /* Computed Parameter: Logic_table_n
   * Referenced by: '<S177>/Logic'
   */
  { false, true, false, false, true, true, false, false, true, false, true, true,
    false, false, false, false },

  /* Mask Parameter: SRFlipFlop_initial_condition
   * Referenced by: '<S117>/Memory'
   */
  false,

  /* Mask Parameter: SRFlipFlop_initial_condition_d
   * Referenced by: '<S118>/Memory'
   */
  false,

  /* Mask Parameter: SRFlipFlop_initial_condition_k
   * Referenced by: '<S176>/Memory'
   */
  false,

  /* Mask Parameter: SRFlipFlop_initial_condition_j
   * Referenced by: '<S177>/Memory'
   */
  false,

  /* Expression: true
   * Referenced by: '<S49>/Constant'
   */
  true,

  /* Expression: false
   * Referenced by: '<S2>/powerUpButton'
   */
  false,

  /* Computed Parameter: Memory_InitialCondition_a
   * Referenced by: '<S2>/Memory'
   */
  false,

  /* Expression: false
   * Referenced by: '<S2>/powerDownButton'
   */
  false,

  /* Computed Parameter: Memory1_InitialCondition
   * Referenced by: '<S2>/Memory1'
   */
  false,

  /* Expression: false
   * Referenced by: '<S2>/resetFaultButton'
   */
  false,

  /* Computed Parameter: Memory2_InitialCondition
   * Referenced by: '<S2>/Memory2'
   */
  false,

  /* Expression: false
   * Referenced by: '<S4>/eStopButton'
   */
  false,

  /* Computed Parameter: Memory_InitialCondition_d
   * Referenced by: '<S4>/Memory'
   */
  false,

  /* Expression: false
   * Referenced by: '<S4>/startButton'
   */
  false,

  /* Computed Parameter: Memory1_InitialCondition_h
   * Referenced by: '<S4>/Memory1'
   */
  false,

  /* Expression: false
   * Referenced by: '<S4>/stopButton'
   */
  false,

  /* Computed Parameter: Memory2_InitialCondition_i
   * Referenced by: '<S4>/Memory2'
   */
  false,

  /* Computed Parameter: Constant_Value_b
   * Referenced by: '<S115>/Constant'
   */
  false,

  /* Computed Parameter: Constant_Value_k
   * Referenced by: '<S116>/Constant'
   */
  false,

  /* Computed Parameter: Constant_Value_c
   * Referenced by: '<S174>/Constant'
   */
  false,

  /* Computed Parameter: Constant_Value_ks
   * Referenced by: '<S175>/Constant'
   */
  false,

  /* Expression: true
   * Referenced by: '<S46>/Constant'
   */
  true,

  /* Expression: true
   * Referenced by: '<S44>/Constant'
   */
  true,

  /* Expression: false
   * Referenced by: '<S1>/powerUpButton'
   */
  false,

  /* Computed Parameter: Memory_InitialCondition_l
   * Referenced by: '<S1>/Memory'
   */
  false,

  /* Expression: false
   * Referenced by: '<S1>/powerDownButton'
   */
  false,

  /* Computed Parameter: Memory1_InitialCondition_a
   * Referenced by: '<S1>/Memory1'
   */
  false,

  /* Expression: false
   * Referenced by: '<S1>/resetFaultButton'
   */
  false,

  /* Computed Parameter: Memory2_InitialCondition_p
   * Referenced by: '<S1>/Memory2'
   */
  false,

  /* Expression: true
   * Referenced by: '<S41>/Constant'
   */
  true,

  /* Expression: true
   * Referenced by: '<S31>/Constant'
   */
  true,

  /* Expression: true
   * Referenced by: '<S9>/Constant1'
   */
  true,

  /* Expression: true
   * Referenced by: '<S9>/Constant2'
   */
  true,

  /* Expression: true
   * Referenced by: '<S9>/Constant3'
   */
  true,

  /* Expression: true
   * Referenced by: '<S11>/Constant1'
   */
  true,

  /* Expression: true
   * Referenced by: '<S11>/Constant2'
   */
  true,

  /* Expression: true
   * Referenced by: '<S11>/Constant3'
   */
  true,

  /* Expression: true
   * Referenced by: '<S5>/Constant'
   */
  true
};

windEmulatorStep4_cal_type *windEmulatorStep4_cal = &windEmulatorStep4_cal_impl;
