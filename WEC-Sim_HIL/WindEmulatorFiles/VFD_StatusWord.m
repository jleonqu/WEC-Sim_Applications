classdef VFD_StatusWord
    properties
        RDY_ON
        RDY_RUN
        RDY_REF
        TRIPPED
        OFF_2_STA
        OFF_3_STA
        SWC_ON_INHIB
        ALARM
        AT_SETPOINT
        REMOTE
        ABOVE_LIMIT
        EXT_CTRL_LOC
        EXT_RUN_ENABLE
        MSW_B13
        MSW_B14
        COMM_ERR
    end
    
    methods
        function this = VFD_StatusWord(StatusWord)
            if nargin==0
                this.RDY_ON = 0;
                this.RDY_RUN = 0;
                this.RDY_REF = 0;
                this.TRIPPED = 0;
                this.OFF_2_STA = 1;
                this.OFF_3_STA = 1;
                this.SWC_ON_INHIB = 0;
                this.ALARM = 0;
                this.AT_SETPOINT = 0;
                this.REMOTE = 1;
                this.ABOVE_LIMIT = 0;
                this.EXT_CTRL_LOC = 1;
                this.EXT_RUN_ENABLE = 0;
                this.MSW_B13 = 0;
                this.MSW_B14 = 0;
                this.COMM_ERR = 0;
            else
                this.RDY_ON = bitget(StatusWord,1);
                this.RDY_RUN = bitget(StatusWord,2);
                this.RDY_REF = bitget(StatusWord,3);
                this.TRIPPED = bitget(StatusWord,4);
                this.OFF_2_STA = bitget(StatusWord,5);
                this.OFF_3_STA = bitget(StatusWord,6);
                this.SWC_ON_INHIB = bitget(StatusWord,7);
                this.ALARM = bitget(StatusWord,8);
                this.AT_SETPOINT = bitget(StatusWord,9);
                this.REMOTE = bitget(StatusWord,10);
                this.ABOVE_LIMIT = bitget(StatusWord,11);
                this.EXT_CTRL_LOC = bitget(StatusWord,12);
                this.EXT_RUN_ENABLE = bitget(StatusWord,13);
                this.MSW_B13 = bitget(StatusWord,14);
                this.MSW_B14 = bitget(StatusWord,15);
                this.COMM_ERR = bitget(StatusWord,16);
            end
        end
        
        function StatusWord = GetStatusWord(this)
            StatusWord = uint16(0);
            if this.RDY_ON
                StatusWord = bitset(StatusWord,1);
            end
            if this.RDY_RUN
                StatusWord = bitset(StatusWord,2);
            end
            if this.RDY_REF
                StatusWord = bitset(StatusWord,3);
            end
            if this.TRIPPED
                StatusWord = bitset(StatusWord,4);
            end
            if this.OFF_2_STA
                StatusWord = bitset(StatusWord,5);
            end
            if this.OFF_3_STA
                StatusWord = bitset(StatusWord,6);
            end
            if this.SWC_ON_INHIB
                StatusWord = bitset(StatusWord,7);
            end
            if this.ALARM
                StatusWord = bitset(StatusWord,8);
            end
            if this.AT_SETPOINT
                StatusWord = bitset(StatusWord,9);
            end
            if this.REMOTE
                StatusWord = bitset(StatusWord,10);
            end
            if this.ABOVE_LIMIT
                StatusWord = bitset(StatusWord,11);
            end
            if this.EXT_CTRL_LOC
                StatusWord = bitset(StatusWord,12);
            end
            if this.EXT_RUN_ENABLE
                StatusWord = bitset(StatusWord,13);
            end
            if this.MSW_B13
                StatusWord = bitset(StatusWord,14);
            end
            if this.MSW_B14
                StatusWord = bitset(StatusWord,15);
            end
            if this.COMM_ERR
                StatusWord = bitset(StatusWord,16);
            end
        end
    end
    
end

