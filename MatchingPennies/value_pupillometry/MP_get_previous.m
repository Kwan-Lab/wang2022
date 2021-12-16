function params1 = MP_get_previous(prcIndex, trials, trialData)

choice = NaN(size(trials.left));
choice(trialData.response == 2) = 0;
choice(trialData.response == 3) = 1;

reward = NaN(size(trials.left));
reward(trialData.response ~= 0 & trials.reward == 0) = 0;
reward(trialData.response ~= 0 & trials.reward == 1) = 1;

choice_01 = [choice(2:end);NaN];  % next trial masks
choice_1 = [NaN;choice(1:end-1)];  % previous trial masks
choice_2 = [NaN;NaN;choice(1:end-2)];

reward_01 = [reward(2:end);NaN];  % next trial masks
reward_1 = [NaN;reward(1:end-1)];  % previous trial masks
reward_2 = [NaN;NaN;reward(1:end-2)];

%params1.trigEvent=NaN(size(trials.left)-1);

    
    params1.trigEvent = choice_01(prcIndex);  % next trial
    params1.trigEvent_1 = choice(prcIndex);  % current trial
    params1.trigEvent_2 = choice_1(prcIndex); % previous 1 trial
    params1.trigEvent_3 = choice_2(prcIndex); % previous 2 trial
    
    %second predictor is outcome; dummy-code: reward=1, error=0, miss=NaN
    %params1.trigEvent2=NaN(size(trials.go)-1);
    params1.trigEvent2 = reward_01(prcIndex);    % next trial
    params1.trigEvent2_1 = reward(prcIndex);     % current trial
    params1.trigEvent2_2 = reward_1(prcIndex); % previous 1 trial
    params1.trigEvent2_3 = reward_2(prcIndex); % previous 2 trial
    
%interaction
params1.trigEvent3 = params1.trigEvent.*params1.trigEvent2;
params1.trigEvent3_1 = params1.trigEvent_1 .* params1.trigEvent2_1;
params1.trigEvent3_2 = params1.trigEvent_2 .* params1.trigEvent2_2;
params1.trigEvent3_3 = params1.trigEvent_3 .* params1.trigEvent2_3;

params1.xtitle = 'Time from cue (s)';
params1.window = [-1:0.1:5];
params1.nback = 0;       %how many trials back to regress against
params1.pvalThresh = 0.01;   %p-value for coefficient be considered significant
params1.interaction = false;
params1.trigTime = trialData.cueTimes(prcIndex);
end