#ifndef RTW_HEADER_abbState_h_
#define RTW_HEADER_abbState_h_
#include "rtwtypes.h"

typedef int32_T abbStateEnum;

/* enum abbStateEnum */
const abbStateEnum abbStateEnum_undefined = 0;/* Default value */
const abbStateEnum abbStateEnum_init = 1;
const abbStateEnum abbStateEnum_notReadyToSwitchOn = 2;
const abbStateEnum abbStateEnum_readyToSwitchOn = 3;
const abbStateEnum abbStateEnum_operationDisabled = 4;
const abbStateEnum abbStateEnum_operationEnabled = 5;
const abbStateEnum abbStateEnum_delayOff1 = 6;

#endif                                 /* RTW_HEADER_abbState_h_ */
