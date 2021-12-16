function [ STIM, RESP, OUTCOME, RULE, EVENT ] = value_getPresentationCodes(presCodeSet)
% % value_getPresentationCodes %
%
%PURPOSE: To read and parse Presentation logfile for further analysis.
%AUTHORS: H Atilgan and AC Kwan, 191210.
%
%INPUT VARIABLES
%   presCodeSEt:    Allow toggling between different types of event code
%                   definitions depending on the Scenario file used (i.e.
%                   task)
%
%OUTPUT VARIABLES
%   STIM:     fields containing stimulus-related eventCode defined in Presentation
%   RESP:     fields containing response-related eventCode defined in Presentation
%   OUTCOME:  fields containing outcome-related eventCode defined in Presentation
%   RULE:     fields containing rule-related eventCode defined in Presentation
%   EVENT:    fields containing other event-related eventCode defined in Presentation

if presCodeSet == 3 || presCodeSet == 8  
    %% Reversal Version 9/30/18 
        STIM.GO = 21;
    
        RESP.MANUAL=1; %when the experimenter manually press a key to deliver a reward to the animal.
        RESP.LEFT=2;
        RESP.RIGHT=3;
    
        OUTCOME.REWARDLEFT = 5;
        OUTCOME.REWARDRIGHT = 6;
        OUTCOME.REWARDMANUAL = 7; %when the experimenter manually press a key to deliver a reward to the animal, 
                %delivered either to left or right port (opposite to the direction of animal's last response)
                %for practical purposes, considering this as an 'Animal did not response', i.e. Miss trial
        OUTCOME.MISS = 8;       
        OUTCOME.NOREWARDLEFT = 75;
        OUTCOME.NOREWARDRIGHT = 76;
    
        RULE.L70R10 = 41;
        RULE.L10R70 = 42;
    
        EVENT.WAITCUE = 90;    %no-lick period, waiting for the next cue
        
elseif presCodeSet == 21 || presCodeSet == 22
    %% Reversal + OptoLaser
     
        STIM.GO = 21;
    
        RESP.MANUAL=1;
        RESP.LEFT=2;
        RESP.RIGHT=3;
    
        OUTCOME.REWARDLEFT = 5;
        OUTCOME.REWARDRIGHT = 6;
        OUTCOME.REWARDMANUAL = 7;
        OUTCOME.NOREWARDLEFT = 75;
        OUTCOME.NOREWARDRIGHT = 76;
        OUTCOME.MISS = 8;      %miss
    
        RULE.L70R10 = 41;
        RULE.L10R70 = 42;
    
        EVENT.WAITCUE = 90;       %inter-trial nolickperiod
        EVENT.SESSIONSTART = 61;
        EVENT.SESSIONEND   = 60;
        EVENT.LASERON      = 63;
        EVENT.LASEROFF     = 62;
        EVENT.REGION       = 1000;

 elseif presCodeSet == 6
    %% Dynamic foraging with 6 sets of reward probabilities
    STIM.GO = 21;
    
    RESP.MANUAL=1;
    RESP.LEFT=2;
    RESP.RIGHT=3;
    
    OUTCOME.REWARDLEFT = 5;
    OUTCOME.REWARDRIGHT = 6;
    OUTCOME.REWARDMANUAL = 7;
    OUTCOME.NOREWARDLEFT = 75;
    OUTCOME.NOREWARDRIGHT = 76;
    OUTCOME.MISS = 8;      %miss
    
    RULE.L70R30 = 41;
    RULE.L70R10 = 42;
    RULE.L30R10 = 43;
    RULE.L30R70 = 44;
    RULE.L10R70 = 45;
    RULE.L10R30 = 46;
    
    EVENT.WAITCUE = 90;  
 
elseif presCodeSet == 31 
    %% Reversal + NE RECORDING
    STIM.GO = 21;
    
    RESP.MANUAL=1;
    RESP.LEFT=2;
    RESP.RIGHT=3;
    
    OUTCOME.REWARDLEFT = 5;
    OUTCOME.REWARDRIGHT = 6;    %place holder
    OUTCOME.REWARDMANUAL = 7;
    OUTCOME.NOREWARDLEFT = 75;
    OUTCOME.NOREWARDRIGHT = 76;  %place holder
    OUTCOME.MISS = 8;      %miss
    
    RULE.L70R10 = 41;
    RULE.L10R70 = 42;
    
    EVENT.WAITCUE = 90;       %inter-trial period
    
else
    disp('Warning: Code set is invalid for value_getPresentationCodes');
    
end