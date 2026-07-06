#include "windEmulatorStep4_WECSim_cal.h"
#include "windEmulatorStep4_WECSim.h"

/* Storage class 'PageSwitching' */
windEmulatorStep4_WECS_cal_type windEmulatorStep4_WECS_cal_impl = {
  /* Start of '<S139>/Nonlinear Wave Elevation' */
  {
    /* Expression: body.centerGravity
     * Referenced by: '<S144>/Center of Gravity'
     */
    { 0.0, 0.0, -3.9 },

    /* Expression: zeros(1,body.dof-3)
     * Referenced by: '<S144>/Constant'
     */
    { 0.0, 0.0, 0.0 },

    /* Expression: 0
     * Referenced by: '<S159>/zero'
     */
    0.0
  }
  ,

  /* End of '<S139>/Nonlinear Wave Elevation' */

  /* Start of '<S60>/Nonlinear Wave Elevation' */
  {
    /* Expression: body.centerGravity
     * Referenced by: '<S65>/Center of Gravity'
     */
    { 0.0, 0.0, -10.9 },

    /* Expression: zeros(1,body.dof-3)
     * Referenced by: '<S65>/Constant'
     */
    { 0.0, 0.0, 0.0 },

    /* Expression: 0
     * Referenced by: '<S80>/zero'
     */
    0.0
  }
  ,

  /* End of '<S60>/Nonlinear Wave Elevation' */

  /* Expression: body.hydroForce
   * Referenced by: '<S146>/Constant'
   */
  {
    {
      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        317844.0, 0.11142198, -0.062912511, 0.0, 0.0, 0.0, 0.11142198, 6562105.2,
        0.033379506, 0.0, 0.0, 0.0, -0.062912511, 0.033379506, -1933551.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
      0.0,

      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
      297.376,

      { 0.0, 0.0, -4.591477 },

      { 0.0, 0.0, -3.9 },

      {
        { 1420336.4556687877, 0.039768779733214833, 42201.8395382624,
          -0.60229974985914891, 1108703.9041504369, 4.0348415103940658 },

        { 434240.43249219697, -0.063881123329988923, 10522.100368383344,
          -0.041323664857100859, 360239.73574879719, 1.2598444843165717 },

        { 0.0, 0.0, 0.0, -0.0, 0.0, 0.0 }
      },

      { -563567.17026318051, -0.083004503715479558, 0.014694211351193975,
        -0.13317392630358385, -1121313.4113269693, 3.4140650849654759,
        -0.20162563218497784, -870800.01595735759, -0.0015155520333447747,
        -3312.86723444801, -0.36391108616245388, 1.100905427628196,
        -0.12859779859779924, 0.0060547084937924954, -825834.87162085623,
        0.24432671188708519, 1.165906852314706, -2.8072354008603857,
        0.79640235258951453, -3101.4724700953689, 0.14528840988108829, 0.0,
        -1.998105390002531, -1.6133709882810949, -1122010.987871229,
        0.00095553741703220782, 2.1841332415998984, 0.0, 0.0, 3.7440855687089059,
        4.7199544645628881, -0.85752125699387338, -0.28226193164378072, 0.0, 0.0,
        0.0 },

      { 1672778.5029394468, -0.09590099083845717, -0.502807128760516,
        -0.16987872068873214, 1313790.1286297347, -15.550725607622594,
        -0.032646965968858313, 20184.550554646645, 0.0024374308375090946,
        -42156.2648733619, 0.18824075359010295, -0.42356350152292688,
        -0.096848862215034068, -0.0001027202328730797, 1789.5268171706973,
        0.0030528342958933163, -0.18819653977614487, -1.3319556587579227,
        -0.065646338768427015, -42085.629605227, -0.0030994749931039572,
        90944.5186375705, -0.91760785863806171, 3.6547269467100705,
        1315587.995824272, -0.042930015557612965, -0.3804142565103511,
        -0.24318017051428911, 1033261.1616709508, -9.467900312769796,
        -21.096315731307119, 0.23785776209866979, 0.20603808240157168,
        5.4082899438748857, -9.3828352096765446, 28694056.040865388 },

      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

      {
        { 340513.6528733772, -0.083004503715479558, 0.014694211351193975,
          -0.13317392630358385, -1121313.4113269693, 3.4140650849654759,
          -0.20162563218497784, 33280.807179200165, -0.0015155520333447747,
          -3312.86723444801, -0.36391108616245388, 1.100905427628196,
          -0.12859779859779924, 0.0060547084937924954, 78245.951515701468,
          0.24432671188708519, 1.165906852314706, -2.8072354008603857,
          0.79640235258951453, -3101.4724700953689, 0.14528840988108829,
          1685242.5288239564, -0.068665257274994118, 3.4804947139078406,
          -1122010.987871229, 0.00095553741703220782, 2.1841332415998984,
          1.9294401327275368, 6991389.641449254, 14.411996204924794,
          4.7199544645628881, -0.85752125699387338, -0.28226193164378072,
          5.0938657021889355, 10.667910636215888, 32921632.662225153 }
      },
      127000.0,
      1031080.8231365577,

      { 3535242.5288239564, 8841389.641449254, 34771632.662225157 },

      { 1.9294401327275368, 5.0938657021889355, 10.667910636215888 }
    }
  },

  /* Expression: body.hydroForce
   * Referenced by: '<S67>/Constant'
   */
  {
    {
      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        -0.0053355609, -0.042099615000000007, -7.3090386000000009e-5, 0.0, 0.0,
        0.0, -0.042099615000000007, 358457.4, 0.12571515, 0.0, 0.0, 0.0,
        -7.3090386000000009e-5, 0.12571515, 358457.4, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0 },
      0.0,

      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
      48.708,

      { 0.0, 0.0, -10.149818 },

      { 0.0, 0.0, -10.9 },

      {
        { 73832.469720141467, 0.0065858536568465884, -42032.637316529632,
          0.092871726806195834, 34832.23798307174, -0.69489183096855778 },

        { 22572.851281314761, 0.00028261400698524955, -10479.913516362081,
          0.008688899402696619, 11317.680184571971, -0.317595151623588 },

        { 0.0, 0.0, -0.0, 0.0, 0.0, -0.0 }
      },

      { -292036.06959670153, 0.0087173819022238891, 0.20837976072843875,
        -0.019804418818860779, 58716.996737042427, -1.2127532666740526,
        -0.013293868476024415, -394871.90862085985, 6.3356416800664466e-5,
        10820.82059688866, -0.0078668120188822282, 0.19213952271597071,
        0.14742875129112853, 0.010193851034339993, -311479.537406777,
        0.47412429984806004, 0.10962536884749394, 0.065624711615548059,
        0.32727289895830997, -12928.26780038833, -0.011894814858454013, 0.0,
        0.12232713772663142, 0.86166496350276978, 79129.779836656133,
        0.0054197022435921519, 0.20614216826764326, 0.0, 0.0,
        -1.5751490994698774, 0.076719824341061171, -0.073710095088816985,
        0.33582461699243754, 0.0, 0.0, 0.0 },

      { 4509.7052122507739, -0.00014285363716791483, 0.013596886509649171,
        -0.00844175886269857, 2136.5447294023584, -0.056272913775421658,
        -0.0019945886697622221, 55.234748163962628, -0.0002182066688891862,
        1299.948720346915, -0.0010744525498388739, 0.00522613712468587,
        0.0052574896933206619, -0.00017510740733657677, 127.1395295938388,
        -0.0060069175071425108, 0.0023976048025386347, 0.00715530498306989,
        0.00082674487263708827, 58.452315215733087, 0.00044220000442507295,
        1592.6492626141051, 0.00045848338972609529, -0.053119136574052955,
        2948.2434908032506, -0.00010120434144842498, 0.008471421473174378,
        -0.0059204505939566424, 1396.7985286082092, -0.059954367918402224,
        -0.041179548140901724, -0.001987043044361154, 0.0078717913742945757,
        -0.17051471875016605, -0.014595916386968305, 75706.192759707585 },

      { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

      {
        { 107318.93665303382, 0.0087173819022238891, 0.20837976072843875,
          -0.019804418818860779, 58716.996737042427, -1.2127532666740526,
          -0.013293868476024415, 4483.09762887549, 6.3356416800664466e-5,
          10820.82059688866, -0.0078668120188822282, 0.19213952271597071,
          0.14742875129112853, 0.010193851034339993, 87875.468842958377,
          0.47412429984806004, 0.10962536884749394, 0.065624711615548059,
          0.32727289895830997, -12928.26780038833, -0.011894814858454013,
          2171872.2052828129, 0.20534544204913072, -3.0865398577015868,
          79129.779836656133, 0.0054197022435921519, 0.20614216826764326,
          0.083018304322499292, 56491.766378381777, -1.3147829338724752,
          0.076719824341061171, -0.073710095088816985, 0.33582461699243754,
          -3.9482048212043566, 0.26036616559740206, 2179049.8418775173 }
      },
      999.0,
      400354.00624973536,

      { 2172871.2052828129, 57490.766378381777, 2180048.8418775173 },

      { 0.083018304322499292, -3.9482048212043566, 0.26036616559740206 }
    }
  },

  /* Computed Parameter: Internal_A_pr
   * Referenced by: '<S531>/Internal'
   */
  { -118.13271821199052, 8.0, -11.623239835183748, 8.0, -15.144377718883565 },

  /* Expression: [0 0 0 0 0 0]
   * Referenced by: '<S58>/Constant'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  /* Expression: [0 0 0 0 0 0]
   * Referenced by: '<S58>/Constant1'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  /* Expression: [0 0 0 0 0 0]
   * Referenced by: '<S58>/Constant2'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  /* Expression: body.centerGravity
   * Referenced by: '<S126>/Center of Gravity'
   */
  { 0.0, 0.0, -10.9 },

  /* Expression: zeros(1,body.dof-3)
   * Referenced by: '<S126>/Constant1'
   */
  { 0.0, 0.0, 0.0 },

  /* Expression: zeros(1,body.dof)
   * Referenced by: '<S128>/Constant'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  /* Expression: zeros(1,body.dof)
   * Referenced by: '<S129>/Constant'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  /* Expression: zeros(1,body.dof)
   * Referenced by: '<S75>/Constant1'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  /* Expression: [0 0 0]
   * Referenced by: '<S75>/Constant2'
   */
  { 0.0, 0.0, 0.0 },

  /* Expression: body.centerGravity
   * Referenced by: '<S75>/Center of Gravity'
   */
  { 0.0, 0.0, -10.9 },

  /* Expression: zeros(1,body.dof-3)
   * Referenced by: '<S63>/Constant'
   */
  { 0.0, 0.0, 0.0 },

  /* Expression: zeros(1,body.dof)
   * Referenced by: '<S78>/Constant'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  /* Expression: body.centerGravity
   * Referenced by: '<S205>/Center of Gravity'
   */
  { 0.0, 0.0, -3.9 },

  /* Expression: zeros(1,body.dof-3)
   * Referenced by: '<S205>/Constant1'
   */
  { 0.0, 0.0, 0.0 },

  /* Expression: zeros(1,body.dof)
   * Referenced by: '<S207>/Constant'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  /* Expression: zeros(1,body.dof)
   * Referenced by: '<S208>/Constant'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  /* Expression: zeros(1,body.dof)
   * Referenced by: '<S154>/Constant1'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  /* Expression: [0 0 0]
   * Referenced by: '<S154>/Constant2'
   */
  { 0.0, 0.0, 0.0 },

  /* Expression: body.centerGravity
   * Referenced by: '<S154>/Center of Gravity'
   */
  { 0.0, 0.0, -3.9 },

  /* Expression: zeros(1,body.dof-3)
   * Referenced by: '<S142>/Constant'
   */
  { 0.0, 0.0, 0.0 },

  /* Expression: zeros(1,body.dof)
   * Referenced by: '<S157>/Constant'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  /* Expression: [0 0 0 0 0 0]
   * Referenced by: '<S58>/Constant3'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  /* Expression: [0 0 0 0 0 0]
   * Referenced by: '<S55>/Constant'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  /* Expression: [0 0 0 0 0 0]
   * Referenced by: '<S55>/Constant1'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  /* Expression: [0 0 0 0 0 0]
   * Referenced by: '<S55>/Constant2'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  /* Mask Parameter: PIDController_D
   * Referenced by: '<S405>/Derivative Gain'
   */
  0.0,

  /* Mask Parameter: PIDController_D_d
   * Referenced by: '<S468>/Derivative Gain'
   */
  0.0,

  /* Mask Parameter: DiscretePIDController_D
   * Referenced by: '<S291>/Derivative Gain'
   */
  0.0,

  /* Mask Parameter: DiscretePIDController_Different
   * Referenced by: '<S293>/UD'
   */
  0.0,

  /* Mask Parameter: posToVel_ICPrevScaledInput
   * Referenced by: '<S553>/UD'
   */
  0.0,

  /* Mask Parameter: PIDController_InitialConditionF
   * Referenced by: '<S407>/Filter'
   */
  0.0,

  /* Mask Parameter: PIDController_InitialConditio_k
   * Referenced by: '<S470>/Filter'
   */
  0.0,

  /* Mask Parameter: PIDController_InitialConditio_o
   * Referenced by: '<S590>/Integrator'
   */
  0.0,

  /* Mask Parameter: PIDController_InitialConditio_a
   * Referenced by: '<S412>/Integrator'
   */
  0.0,

  /* Mask Parameter: PIDController_InitialConditio_c
   * Referenced by: '<S475>/Integrator'
   */
  0.0,

  /* Mask Parameter: DiscretePIDController_InitialCo
   * Referenced by: '<S300>/Integrator'
   */
  0.0,

  /* Mask Parameter: Ramp_InitialOutput
   * Referenced by: '<S370>/Constant1'
   */
  0.0,

  /* Mask Parameter: PIDController_N
   * Referenced by: '<S415>/Filter Coefficient'
   */
  100.0,

  /* Mask Parameter: PIDController_N_p
   * Referenced by: '<S478>/Filter Coefficient'
   */
  100.0,

  /* Mask Parameter: Ramp_slope
   * Referenced by: '<S370>/Step'
   */
  0.1,

  /* Mask Parameter: Ramp_start
   * Referenced by:
   *   '<S370>/Constant'
   *   '<S370>/Step'
   */
  10.0,

  /* Expression: 1
   * Referenced by: '<S124>/Sine Wave Function'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S124>/Sine Wave Function'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S124>/Sine Wave Function'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S124>/Sine Wave Function'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S203>/Sine Wave Function'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S203>/Sine Wave Function'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S203>/Sine Wave Function'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S203>/Sine Wave Function'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S431>/Saturation'
   */
  0.0,

  /* Expression: -1
   * Referenced by: '<S431>/Saturation'
   */
  -1.0,

  /* Expression: 1
   * Referenced by: '<S432>/Saturation'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S432>/Saturation'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S376>/Switch1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S494>/Saturation'
   */
  0.0,

  /* Expression: -1
   * Referenced by: '<S494>/Saturation'
   */
  -1.0,

  /* Expression: 1
   * Referenced by: '<S495>/Saturation'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S495>/Saturation'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S439>/Switch1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S374>/Switch'
   */
  0.0,

  /* Expression: -1
   * Referenced by: '<S437>/Gain2'
   */
  -1.0,

  /* Expression: 1
   * Referenced by: '<S437>/Saturation'
   */
  1.0,

  /* Expression: -1
   * Referenced by: '<S437>/Saturation'
   */
  -1.0,

  /* Expression: 1
   * Referenced by: '<S372>/Saturation'
   */
  1.0,

  /* Expression: -1
   * Referenced by: '<S372>/Saturation'
   */
  -1.0,

  /* Expression: 1
   * Referenced by: '<S435>/Saturation'
   */
  1.0,

  /* Expression: -1
   * Referenced by: '<S435>/Saturation'
   */
  -1.0,

  /* Expression: 0
   * Referenced by: '<S7>/Constant1'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S609>/fromFileSpeedNow_rpm'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S609>/fromFileTorqueNow_Nm'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S555>/Constant1'
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
   * Referenced by: '<S555>/Constant3'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S555>/manualTorqueSetpoint_Nm'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S609>/Set bound'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S555>/manualSpeedSetpoint_rpm'
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
   * Referenced by: '<S376>/Constant1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S373>/shaftSpeedRefMin'
   */
  0.0,

  /* Expression: 50
   * Referenced by: '<S373>/Rate Limiter1'
   */
  50.0,

  /* Expression: -50
   * Referenced by: '<S373>/Rate Limiter1'
   */
  -50.0,

  /* Expression: 0
   * Referenced by: '<S373>/Rate Limiter1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S376>/Switch'
   */
  0.0,

  /* Computed Parameter: Integrator_gainval
   * Referenced by: '<S412>/Integrator'
   */
  0.004,

  /* Computed Parameter: Filter_gainval
   * Referenced by: '<S407>/Filter'
   */
  0.004,

  /* Expression: -1
   * Referenced by: '<S375>/Gain2'
   */
  -1.0,

  /* Expression: 0
   * Referenced by: '<S370>/Step'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S365>/Saturation'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S365>/Saturation'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S371>/kDamping'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S371>/kDampingNow'
   */
  1.0,

  /* Expression: 200
   * Referenced by: '<S498>/Rate Limiter'
   */
  200.0,

  /* Expression: -200
   * Referenced by: '<S498>/Rate Limiter'
   */
  -200.0,

  /* Expression: 0
   * Referenced by: '<S498>/Rate Limiter'
   */
  0.0,

  /* Computed Parameter: Internal_B_pr
   * Referenced by: '<S531>/Internal'
   */
  0.03125,

  /* Computed Parameter: Internal_C_pr
   * Referenced by: '<S531>/Internal'
   */
  0.042279268804058107,

  /* Expression: xinit
   * Referenced by: '<S531>/Internal'
   */
  0.0,

  /* Computed Parameter: Internal_A_pr_j
   * Referenced by: '<S545>/Internal'
   */
  -628.31853071795865,

  /* Computed Parameter: Internal_B_pr_g
   * Referenced by: '<S545>/Internal'
   */
  32.0,

  /* Computed Parameter: Internal_C_pr_a
   * Referenced by: '<S545>/Internal'
   */
  19.634954084936208,

  /* Expression: xinit
   * Referenced by: '<S545>/Internal'
   */
  0.0,

  /* Expression: 0.0075
   * Referenced by: '<S510>/Gain'
   */
  0.0075,

  /* Computed Parameter: Internal_A_pr_e
   * Referenced by: '<S542>/Internal'
   */
  -628.31853071795865,

  /* Computed Parameter: Internal_B_pr_k
   * Referenced by: '<S542>/Internal'
   */
  32.0,

  /* Computed Parameter: Internal_C_pr_m
   * Referenced by: '<S542>/Internal'
   */
  19.634954084936208,

  /* Expression: xinit
   * Referenced by: '<S542>/Internal'
   */
  0.0,

  /* Expression: 0.0075
   * Referenced by: '<S509>/Gain'
   */
  0.0075,

  /* Expression: 5e-6
   * Referenced by: '<S366>/Subsystem_around_RTP_D290B913_fluid_volume'
   */
  5.0E-6,

  /* Expression: 400
   * Referenced by: '<S498>/Subsystem_around_RTP_D2E1D090_liquid_pressure'
   */
  400.0,

  /* Expression: 4.9e-3
   * Referenced by: '<S498>/Subsystem_around_RTP_D2E1D090_liquid_volume'
   */
  0.0049,

  /* Expression: 0
   * Referenced by: '<S371>/kSpring'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S371>/kSpringNow'
   */
  1.0,

  /* Expression: -1
   * Referenced by: '<S367>/Gain'
   */
  -1.0,

  /* Expression: 0
   * Referenced by: '<S439>/Constant1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S436>/shaftSpeedRefMin'
   */
  0.0,

  /* Expression: 50
   * Referenced by: '<S436>/Rate Limiter1'
   */
  50.0,

  /* Expression: -50
   * Referenced by: '<S436>/Rate Limiter1'
   */
  -50.0,

  /* Expression: 0
   * Referenced by: '<S436>/Rate Limiter1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S439>/Switch'
   */
  0.0,

  /* Computed Parameter: Integrator_gainval_b
   * Referenced by: '<S475>/Integrator'
   */
  0.004,

  /* Computed Parameter: Filter_gainval_d
   * Referenced by: '<S470>/Filter'
   */
  0.004,

  /* Expression: -1
   * Referenced by: '<S438>/Gain2'
   */
  -1.0,

  /* Expression: 1/6894.75
   * Referenced by: '<S504>/Gain'
   */
  0.00014503789114906271,

  /* Expression: -1
   * Referenced by: '<S374>/Constant'
   */
  -1.0,

  /* Expression: 1
   * Referenced by: '<S374>/Constant1'
   */
  1.0,

  /* Computed Parameter: DiscreteTimeIntegrator_gainval
   * Referenced by: '<S437>/Discrete-Time Integrator'
   */
  0.004,

  /* Expression: 0
   * Referenced by: '<S437>/Discrete-Time Integrator'
   */
  0.0,

  /* Computed Parameter: DiscreteTimeIntegrator_gainva_l
   * Referenced by: '<S372>/Discrete-Time Integrator'
   */
  0.004,

  /* Expression: 0
   * Referenced by: '<S372>/Discrete-Time Integrator'
   */
  0.0,

  /* Computed Parameter: DiscreteTimeIntegrator_gainva_b
   * Referenced by: '<S435>/Discrete-Time Integrator'
   */
  0.004,

  /* Expression: 0
   * Referenced by: '<S435>/Discrete-Time Integrator'
   */
  0.0,

  /* Expression: 1000*60
   * Referenced by: '<S501>/Gain'
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
   * Referenced by: '<S552>/Gain'
   */
  6.2831853071795862,

  /* Computed Parameter: TSamp_WtEt
   * Referenced by: '<S553>/TSamp'
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
   * Referenced by: '<S609>/caseCounterSignalsNow'
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
   * Referenced by: '<S365>/Constant'
   */
  1.0,

  /* Expression: -1
   * Referenced by: '<S365>/Constant1'
   */
  -1.0,

  /* Expression: 1000*60
   * Referenced by: '<S499>/m3toL'
   */
  60000.0,

  /* Expression: 1000*60
   * Referenced by: '<S500>/Gain'
   */
  60000.0,

  /* Expression: 1
   * Referenced by: '<S6>/Constant1'
   */
  1.0,

  /* Expression: -1
   * Referenced by: '<S59>/Constant3'
   */
  -1.0,

  /* Expression: 1
   * Referenced by: '<S59>/Constant4'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S59>/Switch'
   */
  0.0,

  /* Expression: 100
   * Referenced by: '<S59>/Step'
   */
  100.0,

  /* Expression: 0.0075
   * Referenced by: '<S59>/Step'
   */
  0.0075,

  /* Expression: 0.0075
   * Referenced by: '<S59>/Step'
   */
  0.0075,

  /* Expression: 0.0
   * Referenced by: '<S59>/Delay One Step'
   */
  0.0,

  /* Computed Parameter: Integrator_gainval_n
   * Referenced by: '<S300>/Integrator'
   */
  0.004,

  /* Computed Parameter: Tsamp_WtEt
   * Referenced by: '<S295>/Tsamp'
   */
  250.0,

  /* Expression: 1
   * Referenced by: '<S59>/Gain'
   */
  1.0,

  /* Expression: body.yaw.option
   * Referenced by: '<S70>/Constant'
   */
  0.0,

  /* Expression: body.yaw.option
   * Referenced by: '<S133>/Constant'
   */
  0.0,

  /* Expression: 10e-8
   * Referenced by: '<S60>/Transport Delay'
   */
  1.0E-7,

  /* Expression: 0
   * Referenced by: '<S60>/Transport Delay'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S124>/Constant2'
   */
  1.0,

  /* Expression: pi
   * Referenced by: '<S124>/Constant4'
   */
  3.1415926535897931,

  /* Expression: simu.rampTime
   * Referenced by: '<S124>/Ramp Time'
   */
  100.0,

  /* Expression: 3*pi/2
   * Referenced by: '<S124>/Constant3'
   */
  4.71238898038469,

  /* Expression: 1/2
   * Referenced by: '<S124>/Constant1'
   */
  0.5,

  /* Expression: simu.rampTime
   * Referenced by: '<S68>/Ramp Function Time'
   */
  100.0,

  /* Expression: 1
   * Referenced by: '<S68>/Constant'
   */
  1.0,

  /* Expression: Nwave.amplitude(1,:).*Nwave.amplitude(1,:)
   * Referenced by: '<S126>/Wave Amplitude1'
   */
  4.0,

  /* Expression: Nwave.amplitude(1,:)
   * Referenced by: '<S126>/Wave Amplitude'
   */
  2.0,

  /* Expression: Nwave.omega(1,:)
   * Referenced by: '<S126>/Wave Frequency'
   */
  1.5707963267948966,

  /* Expression: pi/2
   * Referenced by: '<S126>/Constant'
   */
  1.5707963267948966,

  /* Expression: body.largeXYDisplacement.option
   * Referenced by: '<S126>/Displacement Phase Enable1'
   */
  0.0,

  /* Expression: Nwave.omega
   * Referenced by: '<S126>/Wave Frequency1'
   */
  1.5707963267948966,

  /* Expression: Nwave.wavenumber
   * Referenced by: '<S126>/Wave Number'
   */
  0.25352817244519071,

  /* Expression: Nwave.direction
   * Referenced by: '<S126>/Wave direction'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S126>/Sine Wave Function1'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S126>/Sine Wave Function1'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S126>/Sine Wave Function1'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S126>/Sine Wave Function1'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S126>/Sine Wave Function'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S126>/Sine Wave Function'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S126>/Sine Wave Function'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S126>/Sine Wave Function'
   */
  0.0,

  /* Expression: simu.gravity
   * Referenced by: '<S75>/Gravity'
   */
  9.81,

  /* Expression: simu.rho
   * Referenced by: '<S75>/Water Density'
   */
  1000.0,

  /* Expression: body.yaw.option
   * Referenced by: '<S149>/Constant'
   */
  0.0,

  /* Expression: body.yaw.option
   * Referenced by: '<S212>/Constant'
   */
  0.0,

  /* Expression: 10e-8
   * Referenced by: '<S139>/Transport Delay'
   */
  1.0E-7,

  /* Expression: 0
   * Referenced by: '<S139>/Transport Delay'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S203>/Constant2'
   */
  1.0,

  /* Expression: pi
   * Referenced by: '<S203>/Constant4'
   */
  3.1415926535897931,

  /* Expression: simu.rampTime
   * Referenced by: '<S203>/Ramp Time'
   */
  100.0,

  /* Expression: 3*pi/2
   * Referenced by: '<S203>/Constant3'
   */
  4.71238898038469,

  /* Expression: 1/2
   * Referenced by: '<S203>/Constant1'
   */
  0.5,

  /* Expression: simu.rampTime
   * Referenced by: '<S147>/Ramp Function Time'
   */
  100.0,

  /* Expression: 1
   * Referenced by: '<S147>/Constant'
   */
  1.0,

  /* Expression: Nwave.amplitude(1,:).*Nwave.amplitude(1,:)
   * Referenced by: '<S205>/Wave Amplitude1'
   */
  4.0,

  /* Expression: Nwave.amplitude(1,:)
   * Referenced by: '<S205>/Wave Amplitude'
   */
  2.0,

  /* Expression: Nwave.omega(1,:)
   * Referenced by: '<S205>/Wave Frequency'
   */
  1.5707963267948966,

  /* Expression: pi/2
   * Referenced by: '<S205>/Constant'
   */
  1.5707963267948966,

  /* Expression: body.largeXYDisplacement.option
   * Referenced by: '<S205>/Displacement Phase Enable1'
   */
  0.0,

  /* Expression: Nwave.omega
   * Referenced by: '<S205>/Wave Frequency1'
   */
  1.5707963267948966,

  /* Expression: Nwave.wavenumber
   * Referenced by: '<S205>/Wave Number'
   */
  0.25352817244519071,

  /* Expression: Nwave.direction
   * Referenced by: '<S205>/Wave direction'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S205>/Sine Wave Function1'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S205>/Sine Wave Function1'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S205>/Sine Wave Function1'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S205>/Sine Wave Function1'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S205>/Sine Wave Function'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S205>/Sine Wave Function'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S205>/Sine Wave Function'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S205>/Sine Wave Function'
   */
  0.0,

  /* Expression: simu.gravity
   * Referenced by: '<S154>/Gravity'
   */
  9.81,

  /* Expression: simu.rho
   * Referenced by: '<S154>/Water Density'
   */
  1000.0,

  /* Expression: -1
   * Referenced by: '<S58>/Gain6'
   */
  -1.0,

  /* Expression: -1
   * Referenced by: '<S58>/Gain7'
   */
  -1.0,

  /* Expression: -1
   * Referenced by: '<S58>/Gain4'
   */
  -1.0,

  /* Expression: -1
   * Referenced by: '<S58>/Gain5'
   */
  -1.0,

  /* Expression: -1
   * Referenced by: '<S59>/Constant2'
   */
  -1.0,

  /* Expression: 1
   * Referenced by: '<S59>/Constant5'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S59>/Switch1'
   */
  0.0,

  /* Expression: 10e-8
   * Referenced by: '<S64>/Transport Delay'
   */
  1.0E-7,

  /* Expression: 0
   * Referenced by: '<S64>/Transport Delay'
   */
  0.0,

  /* Expression: 10e-8
   * Referenced by: '<S143>/Transport Delay'
   */
  1.0E-7,

  /* Expression: 0
   * Referenced by: '<S143>/Transport Delay'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S13>/acs880SpeedIGain'
   */
  0.0,

  /* Expression: 100
   * Referenced by: '<S609>/vecPercent'
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
   * Referenced by: '<S552>/Constant1'
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
   * Referenced by: '<S552>/lastRawCounts'
   */
  0,

  /* Computed Parameter: Constant_Value_pf
   * Referenced by: '<S552>/Constant'
   */
  524288,

  /* Computed Parameter: lastTurn_InitialCondition
   * Referenced by: '<S552>/lastTurn'
   */
  0,

  /* Computed Parameter: absEncoderTurns_Bias
   * Referenced by: '<S39>/absEncoderTurns'
   */
  0,

  /* Computed Parameter: Internal_A_ir
   * Referenced by: '<S531>/Internal'
   */
  { 0U, 1U, 0U, 2U, 0U },

  /* Computed Parameter: Internal_A_jc
   * Referenced by: '<S531>/Internal'
   */
  { 0U, 2U, 4U, 5U },

  /* Computed Parameter: Internal_B_jc
   * Referenced by: '<S531>/Internal'
   */
  { 0U, 1U },

  /* Computed Parameter: Internal_C_jc
   * Referenced by: '<S531>/Internal'
   */
  { 0U, 0U, 1U, 1U },

  /* Computed Parameter: Internal_A_jc_i
   * Referenced by: '<S545>/Internal'
   */
  { 0U, 1U },

  /* Computed Parameter: Internal_B_jc_i
   * Referenced by: '<S545>/Internal'
   */
  { 0U, 1U },

  /* Computed Parameter: Internal_C_jc_k
   * Referenced by: '<S545>/Internal'
   */
  { 0U, 1U },

  /* Computed Parameter: Internal_A_jc_e
   * Referenced by: '<S542>/Internal'
   */
  { 0U, 1U },

  /* Computed Parameter: Internal_B_jc_h
   * Referenced by: '<S542>/Internal'
   */
  { 0U, 1U },

  /* Computed Parameter: Internal_C_jc_h
   * Referenced by: '<S542>/Internal'
   */
  { 0U, 1U },

  /* Computed Parameter: Constant_Value_im
   * Referenced by: '<S609>/Constant'
   */
  1U,

  /* Computed Parameter: Lengthofinput_Value
   * Referenced by: '<S609>/Length of input'
   */
  560000U,

  /* Computed Parameter: Internal_B_ir
   * Referenced by: '<S531>/Internal'
   */
  0U,

  /* Computed Parameter: Internal_C_ir
   * Referenced by: '<S531>/Internal'
   */
  0U,

  /* Computed Parameter: Internal_A_ir_f
   * Referenced by: '<S545>/Internal'
   */
  0U,

  /* Computed Parameter: Internal_B_ir_l
   * Referenced by: '<S545>/Internal'
   */
  0U,

  /* Computed Parameter: Internal_C_ir_n
   * Referenced by: '<S545>/Internal'
   */
  0U,

  /* Computed Parameter: Internal_A_ir_j
   * Referenced by: '<S542>/Internal'
   */
  0U,

  /* Computed Parameter: Internal_B_ir_a
   * Referenced by: '<S542>/Internal'
   */
  0U,

  /* Computed Parameter: Internal_C_ir_e
   * Referenced by: '<S542>/Internal'
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
   * Referenced by: '<S609>/fileSamples'
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
   * Referenced by: '<S555>/sidType'
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
   * Referenced by: '<S433>/Logic'
   */
  { false, true, false, false, true, true, false, false, true, false, true, true,
    false, false, false, false },

  /* Computed Parameter: Logic_table_o
   * Referenced by: '<S434>/Logic'
   */
  { false, true, false, false, true, true, false, false, true, false, true, true,
    false, false, false, false },

  /* Computed Parameter: Logic_table_h
   * Referenced by: '<S496>/Logic'
   */
  { false, true, false, false, true, true, false, false, true, false, true, true,
    false, false, false, false },

  /* Computed Parameter: Logic_table_n
   * Referenced by: '<S497>/Logic'
   */
  { false, true, false, false, true, true, false, false, true, false, true, true,
    false, false, false, false },

  /* Mask Parameter: SRFlipFlop_initial_condition
   * Referenced by: '<S433>/Memory'
   */
  false,

  /* Mask Parameter: SRFlipFlop_initial_condition_d
   * Referenced by: '<S434>/Memory'
   */
  false,

  /* Mask Parameter: SRFlipFlop_initial_condition_k
   * Referenced by: '<S496>/Memory'
   */
  false,

  /* Mask Parameter: SRFlipFlop_initial_condition_j
   * Referenced by: '<S497>/Memory'
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

  /* Computed Parameter: Constant_Value_bl
   * Referenced by: '<S431>/Constant'
   */
  false,

  /* Computed Parameter: Constant_Value_k
   * Referenced by: '<S432>/Constant'
   */
  false,

  /* Computed Parameter: Constant_Value_cl
   * Referenced by: '<S494>/Constant'
   */
  false,

  /* Computed Parameter: Constant_Value_ks
   * Referenced by: '<S495>/Constant'
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

windEmulatorStep4_WECS_cal_type *windEmulatorStep4_WECSim_cal =
  &windEmulatorStep4_WECS_cal_impl;
