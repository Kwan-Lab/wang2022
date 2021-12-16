function [ sessionData, trialData ] = MP_getSessionData( logData )
% % MP_getSessionData %
%PURPOSE: Retrieve session data for a two-choice value-based decision task.
%AUTHORS: AC Kwan 170518
%
%INPUT ARGUMENTS
%   logdata:    Structure obtained with a call to MP_parseLogfile().
%   
%OUTPUT VARIABLES
%   sessionData:    Structure containing these fields:
%                   {subject, dateTime, nTrials, *lickTimes, *nSwitch}. 
%                   * lickTimes([1 2]):=[left right] lick times.
%   trialData:      Fields: 
%                   {startTimes, cueTimes, outcomeTimes, *cue, *response, *outcome}
%               *cue, response, and outcome for each trial contains
%               correspoinding eventCode from NBS Presentation

%COPY FROM LOGDATA
sessionData.subject = logData.subject;
sessionData.dateTime = logData.dateTime;

%SESSION DATA <<logData.header: 'Subject' 'Trial' 'Event Type' 'Code' 'Time'>>
TYPE = logData.values{3}; %Intersectional approach necessary, because values 2,3 were reused; 
CODE = logData.values{4}; %change in future Presentation scripts: with unique codes, only CODE would be needed to parse the logfile...

% Get the event codes that are used in this log file, and then based on the set of codes used,
% make an educated guess on the correct event code set
tempidx=(strcmp(TYPE,'Nothing') | strcmp(TYPE,'Sound')); %do not consider RESPONSE or MANUAL
codeUsed = unique(CODE(tempidx));         %List the set of event codes found in this logfile

trialData.presCodeSet =1;   %which is the event code set used?
% for l = 1:3                    %test 3 potential event code sets
%     [ STIM, RESP, OUTCOME, RULE, EVENT ] = MP_getPresentationCodes(l);
%     cueCodes = cell2mat(struct2cell(STIM)); %get expected list of stimulus-associated codes as vector
%     outcomeCodes = cell2mat(struct2cell(OUTCOME)); %get expected outcome-associated codes as vector
%     ruleCodes = cell2mat(struct2cell(RULE));       %get expected rule-associated codes as vector
%     eventCodes = cell2mat(struct2cell(EVENT));     %get expected event-associated codes as vector
%     
%     if sum(ismember(codeUsed,[cueCodes; outcomeCodes; ruleCodes; eventCodes])==0)==0     %are all the used code a subset of the expected codes?
%         trialData.presCodeSet = l;
%     end
% end
if isnan(trialData.presCodeSet)     %did not recognize the event code set
    disp('Codes that appeared in this log file');
    codeUsed
    error('ERROR in MP_getSessionData: there are unrecognized event codes found in the log file.');
else                    %otherwise, recognize the event code set and will use it!
    [ STIM, RESP, OUTCOME, RULE, EVENT ] = MP_getPresentationCodes(trialData.presCodeSet);
    
    %interpulse is not used in later presentation files, 100 is actually an
    %outcome code
end

% Get rid of any CODE associated with the last, unfinished trial
lastCode=find(CODE==EVENT.ENDEXPT,1,'last'); 
TYPE = TYPE(1:lastCode);
CODE = CODE(1:lastCode);

% How many trials?
if trialData.presCodeSet == 1
    sessionData.nTrials = sum(CODE==RULE.STARTEXPTLEFT | CODE == RULE.STARTEXPTRIGHT);
elseif trialData.presCodeSet == 2
    sessionData.nTrials = sum(CODE==RULE.L70R10 | CODE == RULE.L10R70);
end

%sessionData.nRules = numel(fieldnames(RULE));
%sessionData.rule_labels = fieldnames(RULE);

% Set up the time axis, and identify lick times
eventCodes = cell2mat(struct2cell(RULE));  %looking for the first of any startexpt Events
time_0 = logData.values{5}(find(ismember(CODE,eventCodes),1,'first'));
time = logData.values{5}-time_0;   %time starts at first instance of startExpt
time = double(time)/10000;         %time as double in seconds
sessionData.lickTimes{1} = time(strcmp(TYPE,'Response') & CODE==RESP.LEFT);    %left licktimes
sessionData.lickTimes{2} = time(strcmp(TYPE,'Response') & CODE==RESP.RIGHT);   %right licktimes

% Stimlus trials - which occurred and when?
cueCodes = cell2mat(struct2cell(STIM)); %stimlus-associated codes as vector
trialData.cue =  CODE(strcmp(TYPE,'Sound') & ismember(CODE,cueCodes));
trialData.cueTimes = time(strcmp(TYPE,'Sound') & CODE==STIM.GO);

% Rule - which was specified and when?
ruleCodes = cell2mat(struct2cell(RULE)); %outcome-associated codes as vector
trialData.rule =  CODE((strcmp(TYPE,'Nothing') | strcmp(TYPE,'Manual')) & ismember(CODE,ruleCodes));
trialData.ruleTimes =  time((strcmp(TYPE,'Nothing') | strcmp(TYPE,'Manual')) & ismember(CODE,ruleCodes));

% Outcome trials - which occurred and when?
outcomeCodes = cell2mat(struct2cell(OUTCOME)); %outcome-associated codes as vector
trialData.outcome =  CODE((strcmp(TYPE,'Nothing') | strcmp(TYPE,'Sound')) & ismember(CODE,outcomeCodes));
trialData.outcomeTimes = time((strcmp(TYPE,'Nothing') | strcmp(TYPE,'Sound')) & ismember(CODE,outcomeCodes));

% number of ITIs - which occurred and when, % 19 stands for ITIs
trialData.iti =  CODE((strcmp(TYPE,'Nothing') ) & ismember(CODE,19));
trialData.itiTimes = time((strcmp(TYPE,'Nothing') ) & ismember(CODE,19));

% Check consistency
if numel(trialData.outcomeTimes) ~= numel(trialData.cueTimes)
    if numel(trialData.cueTimes) ~= 0  % use Rule for trial start      
        disp(['#outcomeTimes = ' int2str(numel(trialData.outcomeTimes)) '; #cueTimes= ' int2str(numel(trialData.cueTimes))]);
        error('ERROR in MP_getSessionData: check #1');
    else
        trialData.cueTimes = trialData.ruleTimes; % in newer verison of MP presentation file, waitlick and go cue are combined
        trialData.cue = trialData.rule;
    end
end
if numel(trialData.outcomeTimes) ~= numel(trialData.ruleTimes)
    disp(['#outcomeTimes = ' int2str(numel(trialData.outcomeTimes)) '; #ruleTimes= ' int2str(numel(trialData.ruleTimes))]);
    error('ERROR in MP_getSessionData: check #2');
end

%% if with pupil recording, get the ele pulse time
if ismember(99, CODE)
    trialData.trigger =  CODE((strcmp(TYPE,'Nothing') ) & ismember(CODE,99));
    trialData.triggerTimes = time((strcmp(TYPE,'Nothing') ) & ismember(CODE,99));
end

if ismember(20, CODE)
    trialData.save = CODE((strcmp(TYPE,'Nothing') ) & ismember(CODE,20));
    trialData.saveTimes = time((strcmp(TYPE,'Nothing') ) & ismember(CODE,20));
end

if ismember(21, CODE)
    trialData.IR = CODE((strcmp(TYPE,'Nothing') ) & ismember(CODE,21));
    trialData.IRTimes = time((strcmp(TYPE,'Nothing') ) & ismember(CODE,21));
end


%% Responses - when?
respIdx = find(strcmp(TYPE,'Response'));
respTimes = time(respIdx);

% What is the time of first lick after the cue
trialData.response = zeros(sessionData.nTrials,1,'uint32');   %Direction of the first lick
trialData.rt = nan(sessionData.nTrials,1);                    %Time of the first lick

idx = find(trialData.outcome~=OUTCOME.MISS); %Idx all non-miss trials
for i = 1:numel(idx)
     %Direction and time of the first lick
     temp = find(respTimes>trialData.cueTimes(idx(i)),1,'first');
     trialData.response(idx(i)) = CODE(respIdx(temp));
     trialData.rt(idx(i)) = respTimes(temp)-trialData.cueTimes(idx(i));
end

% What the lick times associated with each trial?
for i = 1:sessionData.nTrials
     %Times of all licks from prior trial to next trial
     if i>1
         time1=trialData.cueTimes(i-1);
     else
         time1=0;
     end
     
     if i<sessionData.nTrials
         time2=trialData.cueTimes(i+1);
     else
         time2=time(end);
     end
     
     temp=sessionData.lickTimes{1};
     temp=temp(temp>=time1 & temp<=time2); %save only those licks within range
     trialData.leftlickTimes{i}=temp'-trialData.cueTimes(i); %make lick time relative to go cue

     temp=sessionData.lickTimes{2};
     temp=temp(temp>=time1 & temp<=time2); %save only those licks within range
     trialData.rightlickTimes{i}=temp'-trialData.cueTimes(i); %make lick time relative to go cue
end

%if interpulse events are present
%that means to deliver the reward, we did TTL->interpulse event->TTL->reward
%event, i.e. deliver 2 consecutive short pulses (instead of 1 long pulse)

if isfield(EVENT,'INTERPULSE')
    trialData.nInterpulse = sum(strcmp(TYPE,'Nothing') & CODE==EVENT.INTERPULSE);
    if trialData.nInterpulse > 0
        rewardCodes = [OUTCOME.REWARDLEFT OUTCOME.REWARDRIGHT];
        rewardTrials = ismember(trialData.outcome,rewardCodes);
        trialData.outcomeTimes(rewardTrials) = trialData.outcomeTimes(rewardTrials) - 0.1;  %subtract 100 ms, considering the interpulse event
    end
end

end

