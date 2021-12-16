function [ STIM, RESP, OUTCOME, RULE, EVENT ] = MP_getPresentationCodes(presCodeSet)
% % bandit_getPresentationCodes %
%
%PURPOSE: To read and parse Presentation logfile for further analysis.
%AUTHORS: AC Kwan, 170518.
%
%INPUT VARIABLES
%   presCodeSEt:    Allow toggling between different types of event code
%                   definitions
%
%OUTPUT VARIABLES
%   STIM:     fields containing stimulus-related eventCode defined in Presentation
%   RESP:     fields containing response-related eventCode defined in Presentation
%   OUTCOME:  fields containing outcome-related eventCode defined in Presentation
%   RULE:     fields containing rule-related eventCode defined in Presentation
%   EVENT:    fields containing other event-related eventCode defined in Presentation

%% for matching pennies task (Hongli, 05/2008)
if presCodeSet == 1
    STIM.GO = 0;
    
    RESP.LEFT=2;
    RESP.RIGHT=3;
    
    OUTCOME.REWARDLEFT = 100;
    OUTCOME.REWARDRIGHT = 111;
    OUTCOME.NOREWARDLEFT = 101;
    OUTCOME.NOREWARDRIGHT = 110;
    OUTCOME.REWARDMANUAL = 10;
    OUTCOME.MISS = 77;      %miss
    
    RULE.STARTEXPTLEFT = 51;
    RULE.STARTEXPTRIGHT = 52;
    
    EVENT.TRIGGER = 99;
    % EVENT.INTERPULSE = 100;  %time between pulses of water reward
    % no interpulse interval anymore
    EVENT.ENDEXPT = 19;      %inter-trial period
    
elseif presCodeSet == 2
    %% for dynamic foraging task (Huriye, 05/2008)
    STIM.GO = 21;
    
    RESP.LEFT=2;
    RESP.RIGHT=3;
    
    OUTCOME.REWARDLEFT = 5;
    OUTCOME.REWARDRIGHT = 6;
    OUTCOME.NOREWARDLEFT = 75;
    OUTCOME.NOREWARDRIGHT = 76; 
    OUTCOME.MISS = 8;      %miss
    
    RULE.L70R10 = 41;
    RULE.L10R70 = 42;
    
    EVENT.ENDEXPT = 90;      %inter-trial period

elseif presCodeSet == 3
    STIM.GO = 5;

    RESP.LEFT=2;
    RESP.RIGHT=3;

%    OUTCOME.REWARD = 6;    %reward
%    OUTCOME.NOREWARD = 7;  %no reward
    OUTCOME.REWARDLEFT = 6;    
    OUTCOME.REWARDRIGHT = 66;    %place holder
    OUTCOME.NOREWARDLEFT = 7;  
    OUTCOME.NOREWARDRIGHT = 77;  %place holder
    OUTCOME.MISS = 8;      %miss

    RULE.L70R10 = 1;
    RULE.L10R70 = 2;

    EVENT.INTERPULSE = 100;  %time between pulses of water reward
    EVENT.ENDEXPT = 90;      %inter-trial period
    
else
    disp('Warning: Code set is invalid for MP_getPresentationCodes');
    
end