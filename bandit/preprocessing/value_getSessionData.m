function [ sessionData, trialData ] = value_getSessionData( logData, phase )
% % value_getSessionData %
%PURPOSE: Retrieve session data for a two-choice decision task.
%AUTHORS: H Atilgan and AC Kwan 190509
%
%INPUT ARGUMENTS
%   logdata:    Structure obtained with a call to parseLogfile().
%   phase:      The task for this logfile, so we know how to read out the
%               Code associated with the Event Type
%   
%OUTPUT VARIABLES
%   sessionData:    Structure containing these fields:
%                   {subject, dateTime, nTrials, *lickTimes, *nSwitch}. 
%                   * lickTimes([1 2]):=[left right] lick times.
%   trialData:      Fields: 
%                   {startTimes, cueTimes, outcomeTimes, *cue, *response, *outcome}
%                   *cue, response, and outcome for each trial contains
%                   correspoinding eventCode from NBS Presentation

%% Parse the data extracted from logfile
%COPY FROM LOGDATA
sessionData.subject = logData.subject;
sessionData.dateTime = logData.dateTime;
trialData.presCodeSet = phase;

%SESSION DATA <<logData.header: 'Subject' 'Trial' 'Event Type' 'Code' 'Time'>>
% to infer what has occurred, will use both Event Type and Code
TYPE = logData.values{3}; 
CODE = str2double(logData.values{4});
TIME = str2double(logData.values{5});

%Is there any Event Type = Port? (this occcurs if INPUT is unchecked
%accidentally in the Presentation software during experiment), remove
%those entries
tempidx=strcmp(TYPE,'Port');
if any(tempidx)
    error('Port code exists - check this part of the code and then remove error');
    TYPE=TYPE(~tempidx);
    CODE=CODE(~tempidx);
    TIME=TIME(~tempidx);
end
    
% Report all the unique event codes in this log file - this is useful for diagnostic
tempidx=(strcmp(TYPE,'Nothing') | strcmp(TYPE,'Sound')); %do not consider RESPONSE or MANUAL
codeUsed = unique(CODE(tempidx));         %List the set of event codes found in this logfile

[ STIM, RESP, OUTCOME, RULE, EVENT ] = value_getPresentationCodes(trialData.presCodeSet);

% Get rid of any CODE before the first trial
ruleCodes = cell2mat(struct2cell(RULE));  %looking for the first Event, which is a Rule Event
firstCode = find(ismember(CODE,99),1,'first');  % the trigger code is the first code
if isempty(firstCode)
    % if there is no trigger code
    firstCode = find(ismember(CODE,ruleCodes),1,'first');
end
TYPE = TYPE(firstCode:end);
CODE = CODE(firstCode:end);
TIME = TIME(firstCode:end);

% Get rid of any CODE beyond the last full trial
lastCode=find(CODE==EVENT.WAITCUE,1,'last'); 
TYPE = TYPE(1:lastCode);
CODE = CODE(1:lastCode);
TIME = TIME(1:lastCode);

% Set up the time axis, and identify lick times
t = TIME-TIME(1);           %time starts at zero
t = double(t)/10000;        %time as double in seconds
sessionData.lickTimes{1} = t(strcmp(TYPE,'Response') & CODE==RESP.LEFT);    %left licktimes
sessionData.lickTimes{2} = t(strcmp(TYPE,'Response') & CODE==RESP.RIGHT);   %right licktimes

% Stimlus trials - which occurred and when?
cueCodes = cell2mat(struct2cell(STIM)); %stimlus-associated codes as vector
trialData.cue =  CODE(strcmp(TYPE,'Sound') & ismember(CODE,cueCodes));
trialData.cueTimes = t(strcmp(TYPE,'Sound') & CODE==STIM.GO);

% nolick
nolickCodes = 90;
trialData.nolick = CODE(strcmp(TYPE,'Nothing') & ismember(CODE, nolickCodes));
trialData.nolickTimes = t(strcmp(TYPE, 'Nothing') & CODE == nolickCodes);

    
% triggerTimes
if ismember(99, CODE)
    trialData.trigger =  CODE((strcmp(TYPE,'Nothing') ) & ismember(CODE,99));
    trialData.triggerTimes = t((strcmp(TYPE,'Nothing') ) & ismember(CODE,99));
end

% Outcome trials - which occurred and when?
outcomeCodes = cell2mat(struct2cell(OUTCOME)); %outcome-associated codes as vector
trialData.outcome =  CODE((strcmp(TYPE,'Nothing') | strcmp(TYPE,'Sound')) & ismember(CODE,outcomeCodes));
trialData.outcomeTimes = t((strcmp(TYPE,'Nothing') | strcmp(TYPE,'Sound')) & ismember(CODE,outcomeCodes));

% Rule - which was specified and when?
ruleCodes = cell2mat(struct2cell(RULE)); %outcome-associated codes as vector
trialData.rule =  CODE(strcmp(TYPE,'Nothing') & ismember(CODE,ruleCodes));
trialData.ruleTimes =  t(strcmp(TYPE,'Nothing') & ismember(CODE,ruleCodes));

% How many trials?
sessionData.nTrials = numel(trialData.rule);
sessionData.nRules = numel(fieldnames(RULE));
sessionData.rule_labels = fieldnames(RULE);

%% Check consistency

if numel(trialData.outcomeTimes) ~= numel(trialData.cueTimes) || ...
        numel(trialData.outcomeTimes) ~= numel(trialData.ruleTimes) || ...
        numel(trialData.cueTimes) ~= numel(trialData.ruleTimes)
    disp(['    #ruleTimes = ' int2str(numel(trialData.ruleTimes)) '; #cueTimes= ' int2str(numel(trialData.cueTimes)) '; #outcomeTimes= ' int2str(numel(trialData.outcomeTimes))]);
    
    % 191202: Sometimes the occurrence of a cue is not reported by the NBS
    % Presentation software. That is, there is no code associated with the
    % cue that should have occurred. This appears to happen when animal
    % responds immediately after cue onset, so the cue duration is very
    % short (<5 ms?). In later implementation, we use stimulus_time_in = 10
    % to force the cue to have a duration of at least 10 ms, so that should 
    % always be logged
    
    if numel(trialData.ruleTimes) > numel(trialData.cueTimes)   
        disp(['    applying fix to the bug of missing cue']);
        
        % check timing from rule to cue for each trial
        for j = 1:numel(trialData.cueTimes)
            temp = trialData.cueTimes(j) - trialData.ruleTimes;
            diff_rulecue(j) = min(temp(temp>0)); %expect cue to follow rule with a small positive difference
        end
        
        % check for cue in each trial with a rule
        for j = 1:numel(trialData.ruleTimes)-1
            num_cue(j) = sum((trialData.cueTimes > trialData.ruleTimes(j)) & (trialData.cueTimes < trialData.ruleTimes(j+1)));
        end
        idxMissingCue = find(num_cue==0); %trial with no cue event

        trialData.cueTimes = sort([trialData.cueTimes; trialData.ruleTimes(idxMissingCue) + min(diff_rulecue)]); %estimate the cue times and add them
        trialData.cue = [trialData.cue; repmat(trialData.cue(1),numel(idxMissingCue),1)];   %go cue, all the cues are the same indicating start of trial
        
        if numel(trialData.outcomeTimes) ~= numel(trialData.cueTimes) || ...
                numel(trialData.outcomeTimes) ~= numel(trialData.ruleTimes) || ...
                numel(trialData.cueTimes) ~= numel(trialData.ruleTimes)
            error('ERROR in value_getSessionData: still inconsistent after the missing-cue fix');
        end
                
    else
        error('ERROR in value_getSessionData: check consistency');
        
    end
end

trialData.numNolicks = zeros(length(trialData.cue),1);
for tt = 1:length(trialData.numNolicks)-1
    trialData.numNolicks(tt) = sum(trialData.nolickTimes >= trialData.cueTimes(tt) & trialData.nolickTimes < trialData.cueTimes(tt+1));
end
%% Responses - when?
respIdx = find(strcmp(TYPE,'Response') & (CODE==RESP.LEFT | CODE==RESP.RIGHT));   %only want lick responses, not experimenter's manual responses
respTimes = t(respIdx);

% What is the time of first lick after the cue
trialData.response = zeros(sessionData.nTrials,1,'uint32');   %Direction of the first lick
trialData.rt = nan(sessionData.nTrials,1);                    %Time of the first lick

idx = find(all([trialData.outcome~=OUTCOME.MISS trialData.outcome~=OUTCOME.REWARDMANUAL],2)); %Idx all non-miss trials
for i = 1:numel(idx)
     %Direction and time of the first lick
     temp = find(respTimes>trialData.cueTimes(idx(i)),1,'first');
     if ~isempty(temp)                          %if there were a lick (and there should be, given the outcome)
         trialData.response(idx(i)) = CODE(respIdx(temp));
         trialData.rt(idx(i)) = respTimes(temp)-trialData.cueTimes(idx(i));
     end
end

% What the lick times associated with each trial?
trialData.leftlickTimes = cell(sessionData.nTrials,1);
trialData.rightlickTimes = cell(sessionData.nTrials,1);
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
         time2=t(end);
     end
     
     temp=sessionData.lickTimes{1};
     temp=temp(temp>=time1 & temp<=time2); %save only those licks within range
     trialData.leftlickTimes{i}=temp'-trialData.cueTimes(i); %make lick time relative to go cue

     temp=sessionData.lickTimes{2};
     temp=temp(temp>=time1 & temp<=time2); %save only those licks within range
     trialData.rightlickTimes{i}=temp'-trialData.cueTimes(i); %make lick time relative to go cue
end

%inter-trial interval (from outcome to next cue)
trialData.iti = [trialData.cueTimes(2:end)-trialData.outcomeTimes(1:end-1); NaN];

end

