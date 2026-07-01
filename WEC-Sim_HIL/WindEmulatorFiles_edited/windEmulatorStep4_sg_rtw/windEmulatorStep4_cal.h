#ifndef windEmulatorStep4_cal_h_
#define windEmulatorStep4_cal_h_
#include "rtwtypes.h"
#include "expType.h"
#include "sidType.h"

/* Storage class 'PageSwitching', for system '<Root>' */
struct windEmulatorStep4_cal_type {
  real_T Internal_A_pr[5];             /* Computed Parameter: Internal_A_pr
                                        * Referenced by: '<S219>/Internal'
                                        */
  real_T PIDController_D;              /* Mask Parameter: PIDController_D
                                        * Referenced by: '<S93>/Derivative Gain'
                                        */
  real_T PIDController_D_d;            /* Mask Parameter: PIDController_D_d
                                        * Referenced by: '<S156>/Derivative Gain'
                                        */
  real_T posToVel_ICPrevScaledInput;
                                   /* Mask Parameter: posToVel_ICPrevScaledInput
                                    * Referenced by: '<S241>/UD'
                                    */
  real_T PIDController_InitialConditionF;
                              /* Mask Parameter: PIDController_InitialConditionF
                               * Referenced by: '<S95>/Filter'
                               */
  real_T PIDController_InitialConditio_k;
                              /* Mask Parameter: PIDController_InitialConditio_k
                               * Referenced by: '<S158>/Filter'
                               */
  real_T PIDController_InitialConditio_o;
                              /* Mask Parameter: PIDController_InitialConditio_o
                               * Referenced by: '<S278>/Integrator'
                               */
  real_T PIDController_InitialConditio_a;
                              /* Mask Parameter: PIDController_InitialConditio_a
                               * Referenced by: '<S100>/Integrator'
                               */
  real_T PIDController_InitialConditio_c;
                              /* Mask Parameter: PIDController_InitialConditio_c
                               * Referenced by: '<S163>/Integrator'
                               */
  real_T Ramp_InitialOutput;           /* Mask Parameter: Ramp_InitialOutput
                                        * Referenced by: '<S58>/Constant1'
                                        */
  real_T PIDController_N;              /* Mask Parameter: PIDController_N
                                        * Referenced by: '<S103>/Filter Coefficient'
                                        */
  real_T PIDController_N_p;            /* Mask Parameter: PIDController_N_p
                                        * Referenced by: '<S166>/Filter Coefficient'
                                        */
  real_T Ramp_slope;                   /* Mask Parameter: Ramp_slope
                                        * Referenced by: '<S58>/Step'
                                        */
  real_T Ramp_start;                   /* Mask Parameter: Ramp_start
                                        * Referenced by:
                                        *   '<S58>/Constant'
                                        *   '<S58>/Step'
                                        */
  real_T Saturation_UpperSat;          /* Expression: 0
                                        * Referenced by: '<S119>/Saturation'
                                        */
  real_T Saturation_LowerSat;          /* Expression: -1
                                        * Referenced by: '<S119>/Saturation'
                                        */
  real_T Saturation_UpperSat_g;        /* Expression: 1
                                        * Referenced by: '<S120>/Saturation'
                                        */
  real_T Saturation_LowerSat_n;        /* Expression: 0
                                        * Referenced by: '<S120>/Saturation'
                                        */
  real_T Switch1_Threshold;            /* Expression: 0
                                        * Referenced by: '<S64>/Switch1'
                                        */
  real_T Saturation_UpperSat_f;        /* Expression: 0
                                        * Referenced by: '<S182>/Saturation'
                                        */
  real_T Saturation_LowerSat_f;        /* Expression: -1
                                        * Referenced by: '<S182>/Saturation'
                                        */
  real_T Saturation_UpperSat_m;        /* Expression: 1
                                        * Referenced by: '<S183>/Saturation'
                                        */
  real_T Saturation_LowerSat_h;        /* Expression: 0
                                        * Referenced by: '<S183>/Saturation'
                                        */
  real_T Switch1_Threshold_k;          /* Expression: 0
                                        * Referenced by: '<S127>/Switch1'
                                        */
  real_T Switch_Threshold;             /* Expression: 0
                                        * Referenced by: '<S62>/Switch'
                                        */
  real_T Gain2_Gain;                   /* Expression: -1
                                        * Referenced by: '<S125>/Gain2'
                                        */
  real_T Saturation_UpperSat_mj;       /* Expression: 1
                                        * Referenced by: '<S125>/Saturation'
                                        */
  real_T Saturation_LowerSat_b;        /* Expression: -1
                                        * Referenced by: '<S125>/Saturation'
                                        */
  real_T Saturation_UpperSat_d;        /* Expression: 1
                                        * Referenced by: '<S60>/Saturation'
                                        */
  real_T Saturation_LowerSat_g;        /* Expression: -1
                                        * Referenced by: '<S60>/Saturation'
                                        */
  real_T Saturation_UpperSat_p;        /* Expression: 1
                                        * Referenced by: '<S123>/Saturation'
                                        */
  real_T Saturation_LowerSat_e;        /* Expression: -1
                                        * Referenced by: '<S123>/Saturation'
                                        */
  real_T Constant1_Value;              /* Expression: 0
                                        * Referenced by: '<S7>/Constant1'
                                        */
  real_T fromFileSpeedNow_rpm_Gain;    /* Expression: 1
                                        * Referenced by: '<S297>/fromFileSpeedNow_rpm'
                                        */
  real_T fromFileTorqueNow_Nm_Gain;    /* Expression: 1
                                        * Referenced by: '<S297>/fromFileTorqueNow_Nm'
                                        */
  real_T Constant1_Value_e;            /* Expression: 0
                                        * Referenced by: '<S243>/Constant1'
                                        */
  real_T shaftPowerAverage_W_Gain;     /* Expression: 1
                                        * Referenced by: '<S48>/shaftPowerAverage_W'
                                        */
  real_T shaftPower_W_Gain;            /* Expression: 1
                                        * Referenced by: '<S48>/shaftPower_W'
                                        */
  real_T frequency_Hz_Gain;            /* Expression: 1
                                        * Referenced by: '<S28>/frequency_Hz'
                                        */
  real_T motorCurrent_A_Gain;          /* Expression: 1
                                        * Referenced by: '<S28>/motorCurrent_A'
                                        */
  real_T motorSpeed_rpm_Gain;          /* Expression: 1
                                        * Referenced by: '<S28>/motorSpeed_rpm'
                                        */
  real_T motorTorque_Nm_Gain;          /* Expression: 1
                                        * Referenced by: '<S28>/motorTorque_Nm'
                                        */
  real_T motorVoltage_V_Gain;          /* Expression: 1
                                        * Referenced by: '<S28>/motorVoltage_V'
                                        */
  real_T shaftPower_W_Gain_e;          /* Expression: 1
                                        * Referenced by: '<S28>/shaftPower_W'
                                        */
  real_T expRunTime_Value;             /* Expression: 0
                                        * Referenced by: '<S4>/expRunTime'
                                        */
  real_T zeroTorque_Value;             /* Expression: 0
                                        * Referenced by: '<S3>/zeroTorque'
                                        */
  real_T Constant3_Value;              /* Expression: 0
                                        * Referenced by: '<S243>/Constant3'
                                        */
  real_T manualTorqueSetpoint_Nm_Value;/* Expression: 0
                                        * Referenced by: '<S243>/manualTorqueSetpoint_Nm'
                                        */
  real_T Setbound_Value;               /* Expression: 0
                                        * Referenced by: '<S297>/Set bound'
                                        */
  real_T manualSpeedSetpoint_rpm_Value;/* Expression: 0
                                        * Referenced by: '<S243>/manualSpeedSetpoint_rpm'
                                        */
  real_T acs880SpeedPGain_Value;       /* Expression: 0
                                        * Referenced by: '<S13>/acs880SpeedPGain'
                                        */
  real_T speedReference_Value;         /* Expression: 0
                                        * Referenced by: '<S7>/speedReference'
                                        */
  real_T RateLimiter_RisingLim;        /* Expression: 50
                                        * Referenced by: '<S7>/Rate Limiter'
                                        */
  real_T RateLimiter_FallingLim;       /* Expression: -50
                                        * Referenced by: '<S7>/Rate Limiter'
                                        */
  real_T RateLimiter_IC;               /* Expression: 0
                                        * Referenced by: '<S7>/Rate Limiter'
                                        */
  real_T excForceFreq_Hz_Value;        /* Expression: 0
                                        * Referenced by: '<S7>/excForceFreq_Hz'
                                        */
  real_T fw_Gain;                      /* Expression: 2*pi
                                        * Referenced by: '<S7>/f->w'
                                        */
  real_T excForceAmp_N_Value;          /* Expression: 0
                                        * Referenced by: '<S7>/excForceAmp_N'
                                        */
  real_T excForceAmpNow_N_Gain;        /* Expression: 1
                                        * Referenced by: '<S7>/excForceAmpNow_N'
                                        */
  real_T Constant1_Value_c;            /* Expression: 0
                                        * Referenced by: '<S64>/Constant1'
                                        */
  real_T shaftSpeedRefMin_Value;       /* Expression: 0
                                        * Referenced by: '<S61>/shaftSpeedRefMin'
                                        */
  real_T RateLimiter1_RisingLim;       /* Expression: 50
                                        * Referenced by: '<S61>/Rate Limiter1'
                                        */
  real_T RateLimiter1_FallingLim;      /* Expression: -50
                                        * Referenced by: '<S61>/Rate Limiter1'
                                        */
  real_T RateLimiter1_IC;              /* Expression: 0
                                        * Referenced by: '<S61>/Rate Limiter1'
                                        */
  real_T Switch_Threshold_j;           /* Expression: 0
                                        * Referenced by: '<S64>/Switch'
                                        */
  real_T Integrator_gainval;           /* Computed Parameter: Integrator_gainval
                                        * Referenced by: '<S100>/Integrator'
                                        */
  real_T Filter_gainval;               /* Computed Parameter: Filter_gainval
                                        * Referenced by: '<S95>/Filter'
                                        */
  real_T Gain2_Gain_o;                 /* Expression: -1
                                        * Referenced by: '<S63>/Gain2'
                                        */
  real_T Step_Y0;                      /* Expression: 0
                                        * Referenced by: '<S58>/Step'
                                        */
  real_T Saturation_UpperSat_da;       /* Expression: 1
                                        * Referenced by: '<S53>/Saturation'
                                        */
  real_T Saturation_LowerSat_em;       /* Expression: 0
                                        * Referenced by: '<S53>/Saturation'
                                        */
  real_T kDamping_Value;               /* Expression: 0
                                        * Referenced by: '<S59>/kDamping'
                                        */
  real_T kDampingNow_Gain;             /* Expression: 1
                                        * Referenced by: '<S59>/kDampingNow'
                                        */
  real_T RateLimiter_RisingLim_j;      /* Expression: 200
                                        * Referenced by: '<S186>/Rate Limiter'
                                        */
  real_T RateLimiter_FallingLim_e;     /* Expression: -200
                                        * Referenced by: '<S186>/Rate Limiter'
                                        */
  real_T RateLimiter_IC_g;             /* Expression: 0
                                        * Referenced by: '<S186>/Rate Limiter'
                                        */
  real_T Internal_B_pr;                /* Computed Parameter: Internal_B_pr
                                        * Referenced by: '<S219>/Internal'
                                        */
  real_T Internal_C_pr;                /* Computed Parameter: Internal_C_pr
                                        * Referenced by: '<S219>/Internal'
                                        */
  real_T Internal_InitialCondition;    /* Expression: xinit
                                        * Referenced by: '<S219>/Internal'
                                        */
  real_T Internal_A_pr_j;              /* Computed Parameter: Internal_A_pr_j
                                        * Referenced by: '<S233>/Internal'
                                        */
  real_T Internal_B_pr_g;              /* Computed Parameter: Internal_B_pr_g
                                        * Referenced by: '<S233>/Internal'
                                        */
  real_T Internal_C_pr_a;              /* Computed Parameter: Internal_C_pr_a
                                        * Referenced by: '<S233>/Internal'
                                        */
  real_T Internal_InitialCondition_p;  /* Expression: xinit
                                        * Referenced by: '<S233>/Internal'
                                        */
  real_T Gain_Gain;                    /* Expression: 0.0075
                                        * Referenced by: '<S198>/Gain'
                                        */
  real_T Internal_A_pr_e;              /* Computed Parameter: Internal_A_pr_e
                                        * Referenced by: '<S230>/Internal'
                                        */
  real_T Internal_B_pr_k;              /* Computed Parameter: Internal_B_pr_k
                                        * Referenced by: '<S230>/Internal'
                                        */
  real_T Internal_C_pr_m;              /* Computed Parameter: Internal_C_pr_m
                                        * Referenced by: '<S230>/Internal'
                                        */
  real_T Internal_InitialCondition_j;  /* Expression: xinit
                                        * Referenced by: '<S230>/Internal'
                                        */
  real_T Gain_Gain_j;                  /* Expression: 0.0075
                                        * Referenced by: '<S197>/Gain'
                                        */
  real_T RTP_D290B913_fluid_volume_Value;/* Expression: 5e-6
                                          * Referenced by: '<S54>/Subsystem_around_RTP_D290B913_fluid_volume'
                                          */
  real_T RTP_D2E1D090_liquid_pressure_Va;/* Expression: 400
                                          * Referenced by: '<S186>/Subsystem_around_RTP_D2E1D090_liquid_pressure'
                                          */
  real_T RTP_D2E1D090_liquid_volume_Valu;/* Expression: 4.9e-3
                                          * Referenced by: '<S186>/Subsystem_around_RTP_D2E1D090_liquid_volume'
                                          */
  real_T kSpring_Value;                /* Expression: 0
                                        * Referenced by: '<S59>/kSpring'
                                        */
  real_T kSpringNow_Gain;              /* Expression: 1
                                        * Referenced by: '<S59>/kSpringNow'
                                        */
  real_T Gain_Gain_c;                  /* Expression: -1
                                        * Referenced by: '<S55>/Gain'
                                        */
  real_T Constant1_Value_j;            /* Expression: 0
                                        * Referenced by: '<S127>/Constant1'
                                        */
  real_T shaftSpeedRefMin_Value_k;     /* Expression: 0
                                        * Referenced by: '<S124>/shaftSpeedRefMin'
                                        */
  real_T RateLimiter1_RisingLim_o;     /* Expression: 50
                                        * Referenced by: '<S124>/Rate Limiter1'
                                        */
  real_T RateLimiter1_FallingLim_e;    /* Expression: -50
                                        * Referenced by: '<S124>/Rate Limiter1'
                                        */
  real_T RateLimiter1_IC_g;            /* Expression: 0
                                        * Referenced by: '<S124>/Rate Limiter1'
                                        */
  real_T Switch_Threshold_k;           /* Expression: 0
                                        * Referenced by: '<S127>/Switch'
                                        */
  real_T Integrator_gainval_b;       /* Computed Parameter: Integrator_gainval_b
                                      * Referenced by: '<S163>/Integrator'
                                      */
  real_T Filter_gainval_d;             /* Computed Parameter: Filter_gainval_d
                                        * Referenced by: '<S158>/Filter'
                                        */
  real_T Gain2_Gain_e;                 /* Expression: -1
                                        * Referenced by: '<S126>/Gain2'
                                        */
  real_T Gain_Gain_cc;                 /* Expression: 1/6894.75
                                        * Referenced by: '<S192>/Gain'
                                        */
  real_T Constant_Value;               /* Expression: -1
                                        * Referenced by: '<S62>/Constant'
                                        */
  real_T Constant1_Value_b;            /* Expression: 1
                                        * Referenced by: '<S62>/Constant1'
                                        */
  real_T DiscreteTimeIntegrator_gainval;
                           /* Computed Parameter: DiscreteTimeIntegrator_gainval
                            * Referenced by: '<S125>/Discrete-Time Integrator'
                            */
  real_T DiscreteTimeIntegrator_IC;    /* Expression: 0
                                        * Referenced by: '<S125>/Discrete-Time Integrator'
                                        */
  real_T DiscreteTimeIntegrator_gainva_l;
                          /* Computed Parameter: DiscreteTimeIntegrator_gainva_l
                           * Referenced by: '<S60>/Discrete-Time Integrator'
                           */
  real_T DiscreteTimeIntegrator_IC_e;  /* Expression: 0
                                        * Referenced by: '<S60>/Discrete-Time Integrator'
                                        */
  real_T DiscreteTimeIntegrator_gainva_b;
                          /* Computed Parameter: DiscreteTimeIntegrator_gainva_b
                           * Referenced by: '<S123>/Discrete-Time Integrator'
                           */
  real_T DiscreteTimeIntegrator_IC_l;  /* Expression: 0
                                        * Referenced by: '<S123>/Discrete-Time Integrator'
                                        */
  real_T Gain_Gain_f;                  /* Expression: 1000*60
                                        * Referenced by: '<S189>/Gain'
                                        */
  real_T torqueSetpoint_Nm_Gain;       /* Expression: 1
                                        * Referenced by: '<S26>/torqueSetpoint_Nm'
                                        */
  real_T torqueSetpoint_percent_Gain;  /* Expression: 1
                                        * Referenced by: '<S26>/torqueSetpoint_percent'
                                        */
  real_T shaftPowerAverage_W_Gain_l;   /* Expression: 1
                                        * Referenced by: '<S43>/shaftPowerAverage_W'
                                        */
  real_T shaftPower_W_Gain_b;          /* Expression: 1
                                        * Referenced by: '<S43>/shaftPower_W'
                                        */
  real_T dcBusVoltage_V_Gain;          /* Expression: 1
                                        * Referenced by: '<S24>/dcBusVoltage_V'
                                        */
  real_T frequency_Hz_Gain_i;          /* Expression: 1
                                        * Referenced by: '<S24>/frequency_Hz'
                                        */
  real_T motorSpeed_rpm_Gain_j;        /* Expression: 1
                                        * Referenced by: '<S24>/motorSpeed_rpm'
                                        */
  real_T motorTorque_Nm_Gain_l;        /* Expression: 1
                                        * Referenced by: '<S24>/motorTorque_Nm'
                                        */
  real_T shaftPower_W_Gain_l;          /* Expression: 1
                                        * Referenced by: '<S24>/shaftPower_W'
                                        */
  real_T temperature_Gain;             /* Expression: 1
                                        * Referenced by: '<S24>/temperature'
                                        */
  real_T torqueSetpoint_Nm_Gain_m;     /* Expression: 1
                                        * Referenced by: '<S22>/torqueSetpoint_Nm'
                                        */
  real_T torqueSetpoint_percent_Gain_k;/* Expression: 1
                                        * Referenced by: '<S22>/torqueSetpoint_percent'
                                        */
  real_T excShaftPowerAverage_W_Gain;  /* Expression: 1
                                        * Referenced by: '<S51>/excShaftPowerAverage_W'
                                        */
  real_T excShaftPower_W_Gain;         /* Expression: 1
                                        * Referenced by: '<S51>/excShaftPower_W'
                                        */
  real_T hydrPowerAverage_W_Gain;      /* Expression: 1
                                        * Referenced by: '<S51>/hydrPowerAverage_W'
                                        */
  real_T hydrPower_W_Gain;             /* Expression: 1
                                        * Referenced by: '<S51>/hydrPower_W'
                                        */
  real_T ctrlSignal1_Gain;             /* Expression: 1
                                        * Referenced by: '<S35>/ctrlSignal1'
                                        */
  real_T ctrlSignal2_Gain;             /* Expression: 1
                                        * Referenced by: '<S35>/ctrlSignal2'
                                        */
  real_T excShaftSpeed_rpm_Gain;       /* Expression: 1
                                        * Referenced by: '<S35>/excShaftSpeed_rpm'
                                        */
  real_T excShaftTorque_Nm_Gain;       /* Expression: 1
                                        * Referenced by: '<S35>/excShaftTorque_Nm'
                                        */
  real_T genPumpFlow_lpm_Gain;         /* Expression: 1
                                        * Referenced by: '<S35>/genPumpFlow_lpm'
                                        */
  real_T genShaftSpeed_rpm_Gain;       /* Expression: 1
                                        * Referenced by: '<S35>/genShaftSpeed_rpm'
                                        */
  real_T genTorqueCmd_Nm_Gain;         /* Expression: 1
                                        * Referenced by: '<S35>/genTorqueCmd_Nm'
                                        */
  real_T hmOutputShaftTorque_Nm_Gain;  /* Expression: 1
                                        * Referenced by: '<S35>/hmOutputShaftTorque_Nm'
                                        */
  real_T pressure_bar_Gain;            /* Expression: 1
                                        * Referenced by: '<S35>/pressure_bar'
                                        */
  real_T excForce_N_Gain;              /* Expression: 1
                                        * Referenced by: '<S33>/excForce_N'
                                        */
  real_T genSpeedActual_Gain;          /* Expression: 1
                                        * Referenced by: '<S33>/genSpeedActual'
                                        */
  real_T speedRef_rpm_Gain;            /* Expression: 1
                                        * Referenced by: '<S33>/speedRef_rpm'
                                        */
  real_T ramp_Gain;                    /* Expression: 1
                                        * Referenced by: '<S31>/ramp'
                                        */
  real_T time_Gain;                    /* Expression: 1
                                        * Referenced by: '<S31>/time'
                                        */
  real_T Gain_Gain_cf;                 /* Expression: 2*pi
                                        * Referenced by: '<S240>/Gain'
                                        */
  real_T TSamp_WtEt;                   /* Computed Parameter: TSamp_WtEt
                                        * Referenced by: '<S241>/TSamp'
                                        */
  real_T absEncoderPosition_rad_Gain;  /* Expression: 1
                                        * Referenced by: '<S39>/absEncoderPosition_rad'
                                        */
  real_T absEncoderSpeed_rpm_Gain;     /* Expression: 1
                                        * Referenced by: '<S39>/absEncoderSpeed_rpm'
                                        */
  real_T torqueActual_Nm_Gain;         /* Expression: 1
                                        * Referenced by: '<S39>/torqueActual_Nm'
                                        */
  real_T caseCounterSignalsNow_Gain;   /* Expression: 1
                                        * Referenced by: '<S297>/caseCounterSignalsNow'
                                        */
  real_T Constant_Value_m;             /* Expression: 0
                                        * Referenced by: '<S14>/Constant'
                                        */
  real_T Constant2_Value;              /* Expression: 0
                                        * Referenced by: '<S15>/Constant2'
                                        */
  real_T Constant_Value_i;             /* Expression: 1
                                        * Referenced by: '<S53>/Constant'
                                        */
  real_T Constant1_Value_p;            /* Expression: -1
                                        * Referenced by: '<S53>/Constant1'
                                        */
  real_T m3toL_Gain;                   /* Expression: 1000*60
                                        * Referenced by: '<S187>/m3toL'
                                        */
  real_T Gain_Gain_p;                  /* Expression: 1000*60
                                        * Referenced by: '<S188>/Gain'
                                        */
  real_T acs880SpeedIGain_Value;       /* Expression: 0
                                        * Referenced by: '<S13>/acs880SpeedIGain'
                                        */
  real_T vecPercent_Gain;              /* Expression: 100
                                        * Referenced by: '<S297>/vecPercent'
                                        */
  real32_T L1Voltage_Gain;             /* Computed Parameter: L1Voltage_Gain
                                        * Referenced by: '<S9>/L1Voltage'
                                        */
  real32_T L1Current_Gain;             /* Computed Parameter: L1Current_Gain
                                        * Referenced by: '<S9>/L1Current'
                                        */
  real32_T L1PowFactor_Gain;           /* Computed Parameter: L1PowFactor_Gain
                                        * Referenced by: '<S9>/L1PowFactor'
                                        */
  real32_T L1ActivePow_Gain;           /* Computed Parameter: L1ActivePow_Gain
                                        * Referenced by: '<S9>/L1ActivePow'
                                        */
  real32_T L1THDu_Gain;                /* Computed Parameter: L1THDu_Gain
                                        * Referenced by: '<S9>/L1THDu'
                                        */
  real32_T L1THDi_Gain;                /* Computed Parameter: L1THDi_Gain
                                        * Referenced by: '<S9>/L1THDi'
                                        */
  real32_T L2Voltage_Gain;             /* Computed Parameter: L2Voltage_Gain
                                        * Referenced by: '<S9>/L2Voltage'
                                        */
  real32_T L2Current_Gain;             /* Computed Parameter: L2Current_Gain
                                        * Referenced by: '<S9>/L2Current'
                                        */
  real32_T L2PowFactor_Gain;           /* Computed Parameter: L2PowFactor_Gain
                                        * Referenced by: '<S9>/L2PowFactor'
                                        */
  real32_T L2ActivePow_Gain;           /* Computed Parameter: L2ActivePow_Gain
                                        * Referenced by: '<S9>/L2ActivePow'
                                        */
  real32_T L2THDu_Gain;                /* Computed Parameter: L2THDu_Gain
                                        * Referenced by: '<S9>/L2THDu'
                                        */
  real32_T L2THDi_Gain;                /* Computed Parameter: L2THDi_Gain
                                        * Referenced by: '<S9>/L2THDi'
                                        */
  real32_T L3Voltage_Gain;             /* Computed Parameter: L3Voltage_Gain
                                        * Referenced by: '<S9>/L3Voltage'
                                        */
  real32_T L3Current_Gain;             /* Computed Parameter: L3Current_Gain
                                        * Referenced by: '<S9>/L3Current'
                                        */
  real32_T L3PowFactor_Gain;           /* Computed Parameter: L3PowFactor_Gain
                                        * Referenced by: '<S9>/L3PowFactor'
                                        */
  real32_T L3ActivePow_Gain;           /* Computed Parameter: L3ActivePow_Gain
                                        * Referenced by: '<S9>/L3ActivePow'
                                        */
  real32_T L3THDu_Gain;                /* Computed Parameter: L3THDu_Gain
                                        * Referenced by: '<S9>/L3THDu'
                                        */
  real32_T L3THDi_Gain;                /* Computed Parameter: L3THDi_Gain
                                        * Referenced by: '<S9>/L3THDi'
                                        */
  real32_T totalFrequency_Gain;       /* Computed Parameter: totalFrequency_Gain
                                       * Referenced by: '<S9>/totalFrequency'
                                       */
  real32_T totalPowFactor_Gain;       /* Computed Parameter: totalPowFactor_Gain
                                       * Referenced by: '<S9>/totalPowFactor'
                                       */
  real32_T totalActivePow_Gain;       /* Computed Parameter: totalActivePow_Gain
                                       * Referenced by: '<S9>/totalActivePow'
                                       */
  real32_T L1L2Voltage_Gain;           /* Computed Parameter: L1L2Voltage_Gain
                                        * Referenced by: '<S9>/L1L2Voltage'
                                        */
  real32_T L2L3Voltage_Gain;           /* Computed Parameter: L2L3Voltage_Gain
                                        * Referenced by: '<S9>/L2L3Voltage'
                                        */
  real32_T L3L1Voltage_Gain;           /* Computed Parameter: L3L1Voltage_Gain
                                        * Referenced by: '<S9>/L3L1Voltage'
                                        */
  real32_T L1Voltage_Gain_h;           /* Computed Parameter: L1Voltage_Gain_h
                                        * Referenced by: '<S11>/L1Voltage'
                                        */
  real32_T L1Current_Gain_f;           /* Computed Parameter: L1Current_Gain_f
                                        * Referenced by: '<S11>/L1Current'
                                        */
  real32_T L1PowFactor_Gain_a;         /* Computed Parameter: L1PowFactor_Gain_a
                                        * Referenced by: '<S11>/L1PowFactor'
                                        */
  real32_T L1ActivePow_Gain_o;         /* Computed Parameter: L1ActivePow_Gain_o
                                        * Referenced by: '<S11>/L1ActivePow'
                                        */
  real32_T L1THDu_Gain_a;              /* Computed Parameter: L1THDu_Gain_a
                                        * Referenced by: '<S11>/L1THDu'
                                        */
  real32_T L1THDi_Gain_f;              /* Computed Parameter: L1THDi_Gain_f
                                        * Referenced by: '<S11>/L1THDi'
                                        */
  real32_T L2Voltage_Gain_a;           /* Computed Parameter: L2Voltage_Gain_a
                                        * Referenced by: '<S11>/L2Voltage'
                                        */
  real32_T L2Current_Gain_j;           /* Computed Parameter: L2Current_Gain_j
                                        * Referenced by: '<S11>/L2Current'
                                        */
  real32_T L2PowFactor_Gain_p;         /* Computed Parameter: L2PowFactor_Gain_p
                                        * Referenced by: '<S11>/L2PowFactor'
                                        */
  real32_T L2ActivePow_Gain_b;         /* Computed Parameter: L2ActivePow_Gain_b
                                        * Referenced by: '<S11>/L2ActivePow'
                                        */
  real32_T L2THDu_Gain_i;              /* Computed Parameter: L2THDu_Gain_i
                                        * Referenced by: '<S11>/L2THDu'
                                        */
  real32_T L2THDi_Gain_j;              /* Computed Parameter: L2THDi_Gain_j
                                        * Referenced by: '<S11>/L2THDi'
                                        */
  real32_T L3Voltage_Gain_h;           /* Computed Parameter: L3Voltage_Gain_h
                                        * Referenced by: '<S11>/L3Voltage'
                                        */
  real32_T L3Current_Gain_i;           /* Computed Parameter: L3Current_Gain_i
                                        * Referenced by: '<S11>/L3Current'
                                        */
  real32_T L3PowFactor_Gain_k;         /* Computed Parameter: L3PowFactor_Gain_k
                                        * Referenced by: '<S11>/L3PowFactor'
                                        */
  real32_T L3ActivePow_Gain_f;         /* Computed Parameter: L3ActivePow_Gain_f
                                        * Referenced by: '<S11>/L3ActivePow'
                                        */
  real32_T L3THDu_Gain_b;              /* Computed Parameter: L3THDu_Gain_b
                                        * Referenced by: '<S11>/L3THDu'
                                        */
  real32_T L3THDi_Gain_p;              /* Computed Parameter: L3THDi_Gain_p
                                        * Referenced by: '<S11>/L3THDi'
                                        */
  real32_T totalFrequency_Gain_b;   /* Computed Parameter: totalFrequency_Gain_b
                                     * Referenced by: '<S11>/totalFrequency'
                                     */
  real32_T totalPowFactor_Gain_f;   /* Computed Parameter: totalPowFactor_Gain_f
                                     * Referenced by: '<S11>/totalPowFactor'
                                     */
  real32_T totalActivePow_Gain_m;   /* Computed Parameter: totalActivePow_Gain_m
                                     * Referenced by: '<S11>/totalActivePow'
                                     */
  real32_T L1L2Voltage_Gain_j;         /* Computed Parameter: L1L2Voltage_Gain_j
                                        * Referenced by: '<S11>/L1L2Voltage'
                                        */
  real32_T L2L3Voltage_Gain_m;         /* Computed Parameter: L2L3Voltage_Gain_m
                                        * Referenced by: '<S11>/L2L3Voltage'
                                        */
  real32_T L3L1Voltage_Gain_k;         /* Computed Parameter: L3L1Voltage_Gain_k
                                        * Referenced by: '<S11>/L3L1Voltage'
                                        */
  int32_T Constant1_Value_jt;          /* Computed Parameter: Constant1_Value_jt
                                        * Referenced by: '<S240>/Constant1'
                                        */
  int32_T state_Bias;                  /* Computed Parameter: state_Bias
                                        * Referenced by: '<S26>/state'
                                        */
  int32_T state_Bias_l;                /* Computed Parameter: state_Bias_l
                                        * Referenced by: '<S22>/state'
                                        */
  int32_T lastRawCounts_InitialCondition;
                           /* Computed Parameter: lastRawCounts_InitialCondition
                            * Referenced by: '<S240>/lastRawCounts'
                            */
  int32_T Constant_Value_p;            /* Computed Parameter: Constant_Value_p
                                        * Referenced by: '<S240>/Constant'
                                        */
  int32_T lastTurn_InitialCondition;
                                /* Computed Parameter: lastTurn_InitialCondition
                                 * Referenced by: '<S240>/lastTurn'
                                 */
  int32_T absEncoderTurns_Bias;      /* Computed Parameter: absEncoderTurns_Bias
                                      * Referenced by: '<S39>/absEncoderTurns'
                                      */
  uint32_T Internal_A_ir[5];           /* Computed Parameter: Internal_A_ir
                                        * Referenced by: '<S219>/Internal'
                                        */
  uint32_T Internal_A_jc[4];           /* Computed Parameter: Internal_A_jc
                                        * Referenced by: '<S219>/Internal'
                                        */
  uint32_T Internal_B_jc[2];           /* Computed Parameter: Internal_B_jc
                                        * Referenced by: '<S219>/Internal'
                                        */
  uint32_T Internal_C_jc[4];           /* Computed Parameter: Internal_C_jc
                                        * Referenced by: '<S219>/Internal'
                                        */
  uint32_T Internal_A_jc_i[2];         /* Computed Parameter: Internal_A_jc_i
                                        * Referenced by: '<S233>/Internal'
                                        */
  uint32_T Internal_B_jc_i[2];         /* Computed Parameter: Internal_B_jc_i
                                        * Referenced by: '<S233>/Internal'
                                        */
  uint32_T Internal_C_jc_k[2];         /* Computed Parameter: Internal_C_jc_k
                                        * Referenced by: '<S233>/Internal'
                                        */
  uint32_T Internal_A_jc_e[2];         /* Computed Parameter: Internal_A_jc_e
                                        * Referenced by: '<S230>/Internal'
                                        */
  uint32_T Internal_B_jc_h[2];         /* Computed Parameter: Internal_B_jc_h
                                        * Referenced by: '<S230>/Internal'
                                        */
  uint32_T Internal_C_jc_h[2];         /* Computed Parameter: Internal_C_jc_h
                                        * Referenced by: '<S230>/Internal'
                                        */
  uint32_T Constant_Value_im;          /* Computed Parameter: Constant_Value_im
                                        * Referenced by: '<S297>/Constant'
                                        */
  uint32_T Lengthofinput_Value;       /* Computed Parameter: Lengthofinput_Value
                                       * Referenced by: '<S297>/Length of input'
                                       */
  uint32_T Internal_B_ir;              /* Computed Parameter: Internal_B_ir
                                        * Referenced by: '<S219>/Internal'
                                        */
  uint32_T Internal_C_ir;              /* Computed Parameter: Internal_C_ir
                                        * Referenced by: '<S219>/Internal'
                                        */
  uint32_T Internal_A_ir_f;            /* Computed Parameter: Internal_A_ir_f
                                        * Referenced by: '<S233>/Internal'
                                        */
  uint32_T Internal_B_ir_l;            /* Computed Parameter: Internal_B_ir_l
                                        * Referenced by: '<S233>/Internal'
                                        */
  uint32_T Internal_C_ir_n;            /* Computed Parameter: Internal_C_ir_n
                                        * Referenced by: '<S233>/Internal'
                                        */
  uint32_T Internal_A_ir_j;            /* Computed Parameter: Internal_A_ir_j
                                        * Referenced by: '<S230>/Internal'
                                        */
  uint32_T Internal_B_ir_a;            /* Computed Parameter: Internal_B_ir_a
                                        * Referenced by: '<S230>/Internal'
                                        */
  uint32_T Internal_C_ir_e;            /* Computed Parameter: Internal_C_ir_e
                                        * Referenced by: '<S230>/Internal'
                                        */
  uint32_T runCounter_Bias;            /* Computed Parameter: runCounter_Bias
                                        * Referenced by: '<S31>/runCounter'
                                        */
  uint32_T stepCounter_Bias;           /* Computed Parameter: stepCounter_Bias
                                        * Referenced by: '<S31>/stepCounter'
                                        */
  uint32_T absEncoderCounts_Bias;   /* Computed Parameter: absEncoderCounts_Bias
                                     * Referenced by: '<S39>/absEncoderCounts'
                                     */
  uint32_T loopAdd_Value;              /* Computed Parameter: loopAdd_Value
                                        * Referenced by: '<S29>/loopAdd'
                                        */
  uint32_T Memory_InitialCondition;
                                  /* Computed Parameter: Memory_InitialCondition
                                   * Referenced by: '<S29>/Memory'
                                   */
  uint32_T loopCounter_Gain;           /* Computed Parameter: loopCounter_Gain
                                        * Referenced by: '<S29>/loopCounter'
                                        */
  uint32_T fileSamples_Bias;           /* Computed Parameter: fileSamples_Bias
                                        * Referenced by: '<S297>/fileSamples'
                                        */
  expTypeEnum expType_Value;           /* Expression: expTypeEnum.off
                                        * Referenced by: '<S4>/expType'
                                        */
  expTypeEnum expModeHil_Value;        /* Expression: expTypeEnum.hil
                                        * Referenced by: '<S4>/expModeHil'
                                        */
  expTypeEnum expModeSid_Value;        /* Expression: expTypeEnum.sid
                                        * Referenced by: '<S4>/expModeSid'
                                        */
  sidTypeEnum sidType_Value;           /* Expression: sidTypeEnum.off
                                        * Referenced by: '<S243>/sidType'
                                        */
  uint16_T statusWord_Bias;            /* Computed Parameter: statusWord_Bias
                                        * Referenced by: '<S28>/statusWord'
                                        */
  uint16_T ctrlWord_Bias;              /* Computed Parameter: ctrlWord_Bias
                                        * Referenced by: '<S26>/ctrlWord'
                                        */
  uint16_T statusWord_Bias_l;          /* Computed Parameter: statusWord_Bias_l
                                        * Referenced by: '<S24>/statusWord'
                                        */
  uint16_T ctrlWord_Bias_a;            /* Computed Parameter: ctrlWord_Bias_a
                                        * Referenced by: '<S22>/ctrlWord'
                                        */
  uint16_T expType_Bias;               /* Computed Parameter: expType_Bias
                                        * Referenced by: '<S31>/expType'
                                        */
  uint8_T speedCtrlReset_Gain;        /* Computed Parameter: speedCtrlReset_Gain
                                       * Referenced by: '<S33>/speedCtrlReset'
                                       */
  uint8_T absEncoderStatus1_Bias;  /* Computed Parameter: absEncoderStatus1_Bias
                                    * Referenced by: '<S39>/absEncoderStatus1'
                                    */
  uint8_T absEncoderStatus2_Bias;  /* Computed Parameter: absEncoderStatus2_Bias
                                    * Referenced by: '<S39>/absEncoderStatus2'
                                    */
  boolean_T Logic_table[16];           /* Computed Parameter: Logic_table
                                        * Referenced by: '<S121>/Logic'
                                        */
  boolean_T Logic_table_o[16];         /* Computed Parameter: Logic_table_o
                                        * Referenced by: '<S122>/Logic'
                                        */
  boolean_T Logic_table_h[16];         /* Computed Parameter: Logic_table_h
                                        * Referenced by: '<S184>/Logic'
                                        */
  boolean_T Logic_table_n[16];         /* Computed Parameter: Logic_table_n
                                        * Referenced by: '<S185>/Logic'
                                        */
  boolean_T SRFlipFlop_initial_condition;
                                 /* Mask Parameter: SRFlipFlop_initial_condition
                                  * Referenced by: '<S121>/Memory'
                                  */
  boolean_T SRFlipFlop_initial_condition_d;
                               /* Mask Parameter: SRFlipFlop_initial_condition_d
                                * Referenced by: '<S122>/Memory'
                                */
  boolean_T SRFlipFlop_initial_condition_k;
                               /* Mask Parameter: SRFlipFlop_initial_condition_k
                                * Referenced by: '<S184>/Memory'
                                */
  boolean_T SRFlipFlop_initial_condition_j;
                               /* Mask Parameter: SRFlipFlop_initial_condition_j
                                * Referenced by: '<S185>/Memory'
                                */
  boolean_T Constant_Value_o;          /* Expression: true
                                        * Referenced by: '<S49>/Constant'
                                        */
  boolean_T powerUpButton_Value;       /* Expression: false
                                        * Referenced by: '<S2>/powerUpButton'
                                        */
  boolean_T Memory_InitialCondition_a;
                                /* Computed Parameter: Memory_InitialCondition_a
                                 * Referenced by: '<S2>/Memory'
                                 */
  boolean_T powerDownButton_Value;     /* Expression: false
                                        * Referenced by: '<S2>/powerDownButton'
                                        */
  boolean_T Memory1_InitialCondition;
                                 /* Computed Parameter: Memory1_InitialCondition
                                  * Referenced by: '<S2>/Memory1'
                                  */
  boolean_T resetFaultButton_Value;    /* Expression: false
                                        * Referenced by: '<S2>/resetFaultButton'
                                        */
  boolean_T Memory2_InitialCondition;
                                 /* Computed Parameter: Memory2_InitialCondition
                                  * Referenced by: '<S2>/Memory2'
                                  */
  boolean_T eStopButton_Value;         /* Expression: false
                                        * Referenced by: '<S4>/eStopButton'
                                        */
  boolean_T Memory_InitialCondition_d;
                                /* Computed Parameter: Memory_InitialCondition_d
                                 * Referenced by: '<S4>/Memory'
                                 */
  boolean_T startButton_Value;         /* Expression: false
                                        * Referenced by: '<S4>/startButton'
                                        */
  boolean_T Memory1_InitialCondition_h;
                               /* Computed Parameter: Memory1_InitialCondition_h
                                * Referenced by: '<S4>/Memory1'
                                */
  boolean_T stopButton_Value;          /* Expression: false
                                        * Referenced by: '<S4>/stopButton'
                                        */
  boolean_T Memory2_InitialCondition_i;
                               /* Computed Parameter: Memory2_InitialCondition_i
                                * Referenced by: '<S4>/Memory2'
                                */
  boolean_T Constant_Value_b;          /* Computed Parameter: Constant_Value_b
                                        * Referenced by: '<S119>/Constant'
                                        */
  boolean_T Constant_Value_k;          /* Computed Parameter: Constant_Value_k
                                        * Referenced by: '<S120>/Constant'
                                        */
  boolean_T Constant_Value_c;          /* Computed Parameter: Constant_Value_c
                                        * Referenced by: '<S182>/Constant'
                                        */
  boolean_T Constant_Value_ks;         /* Computed Parameter: Constant_Value_ks
                                        * Referenced by: '<S183>/Constant'
                                        */
  boolean_T Constant_Value_e;          /* Expression: true
                                        * Referenced by: '<S46>/Constant'
                                        */
  boolean_T Constant_Value_j;          /* Expression: true
                                        * Referenced by: '<S44>/Constant'
                                        */
  boolean_T powerUpButton_Value_b;     /* Expression: false
                                        * Referenced by: '<S1>/powerUpButton'
                                        */
  boolean_T Memory_InitialCondition_l;
                                /* Computed Parameter: Memory_InitialCondition_l
                                 * Referenced by: '<S1>/Memory'
                                 */
  boolean_T powerDownButton_Value_k;   /* Expression: false
                                        * Referenced by: '<S1>/powerDownButton'
                                        */
  boolean_T Memory1_InitialCondition_a;
                               /* Computed Parameter: Memory1_InitialCondition_a
                                * Referenced by: '<S1>/Memory1'
                                */
  boolean_T resetFaultButton_Value_n;  /* Expression: false
                                        * Referenced by: '<S1>/resetFaultButton'
                                        */
  boolean_T Memory2_InitialCondition_p;
                               /* Computed Parameter: Memory2_InitialCondition_p
                                * Referenced by: '<S1>/Memory2'
                                */
  boolean_T Constant_Value_l;          /* Expression: true
                                        * Referenced by: '<S41>/Constant'
                                        */
  boolean_T Constant_Value_kp;         /* Expression: true
                                        * Referenced by: '<S31>/Constant'
                                        */
  boolean_T Constant1_Value_k;         /* Expression: true
                                        * Referenced by: '<S9>/Constant1'
                                        */
  boolean_T Constant2_Value_e;         /* Expression: true
                                        * Referenced by: '<S9>/Constant2'
                                        */
  boolean_T Constant3_Value_d;         /* Expression: true
                                        * Referenced by: '<S9>/Constant3'
                                        */
  boolean_T Constant1_Value_cc;        /* Expression: true
                                        * Referenced by: '<S11>/Constant1'
                                        */
  boolean_T Constant2_Value_eh;        /* Expression: true
                                        * Referenced by: '<S11>/Constant2'
                                        */
  boolean_T Constant3_Value_g;         /* Expression: true
                                        * Referenced by: '<S11>/Constant3'
                                        */
  boolean_T Constant_Value_m3;         /* Expression: true
                                        * Referenced by: '<S5>/Constant'
                                        */
};

/* Storage class 'PageSwitching' */
extern windEmulatorStep4_cal_type windEmulatorStep4_cal_impl;
extern windEmulatorStep4_cal_type *windEmulatorStep4_cal;

#endif                                 /* windEmulatorStep4_cal_h_ */
