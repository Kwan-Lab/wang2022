function MP_pupilMLR_StaySwitch(dataIndex)

% regression with only C(n+1), separated by stay or switch trials:
%   if the choice signal is movement related:
%      1) in stay trial, the model R = beta*C(n-1? is equivalent to R =
%      beta*C(n), we should see the results are the same (always the same,
%      doesn't matter whether it is actual choice or movement
%      2) in switch trial, the model R = beta*C(n-1) is different, if it is
%      movement related, we should see opposite effects since the movement
%      is the opposite right now, however, if it is choice signal, it
%      should still be the same.

% C(n-1) follows the same reasoning

%% the above reasoning is not going to work: in switch trials the choice preditor itself is opposite

% how about remove the effect of C(n) first, then regress against
% switch/stay trials?
nFiles = size(dataIndex,1);

for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},[dataIndex.LogFileName{ii}(1:end-4),'_beh_cut.mat']));
    load(fullfile(fn_beh.folder,fn_beh.name));
    

    % load pupil files
    date = num2str(dataIndex.DateNumber(ii));
    pup_name = fullfile(dataIndex.BehPath{ii}, ['*',date(1:6),'*_pup.mat']);
    fn_pup = dir(pup_name);
    if length(fn_pup) == 1
        
        savepupilfigpath = fullfile(dataIndex.BehPath{ii},[date(1:6),'_figs-pupil']);
        if ~exist(savepupilfigpath,'dir')
            mkdir(savepupilfigpath);
        end
        
        cd(savepupilfigpath);
        
    
        load(fullfile(fn_pup.folder,fn_pup.name));
        
        % make folders to save analysis and plots
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
        if ~exist(savematpath,'dir')
            mkdir(savematpath);
        end
        savepupilfigpath = fullfile(dataIndex.BehPath{ii},[date(1:6),'_figs-pupil']);
        if ~exist(savepupilfigpath,'dir')
            mkdir(savepupilfigpath);
        end
        
        cd(savepupilfigpath);
        
      
        saveRegName = fullfile(savematpath,[fn_beh.name(1:end-7),'regCR_changeSS.mat']);
%% control regression with C(n)
        params=[];
        
        if trialData.cutPoint ~= 0
            % cut the trial
            fn = fieldnames(trials);
            for tt = 1:length(fn)
                trials.(fn{tt}) = trials.(fn{tt})(1:trialData.cutPoint);
            end
            cueTime = trialData.cueTimes(1:trialData.cutPoint);
            response = trialData.response(1:trialData.cutPoint);
        else
            cueTime = trialData.cueTimes;
            response = trialData.response;
        end
        
        params.trigTime = cueTime;
        
        params.trigEvent=NaN(size(trials.left));
        params.trigEvent(trials.left) = 0;
        params.trigEvent(trials.right) = 1;

        fieldname={'go'};
        trialMask = getMask(trials,fieldname);
       
        
        params.xtitle = 'Time from cue (s)';
        params.window = [-3:0.1:5];
        params.nback = 0;       %how many trials back to regress against
        params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
        params.interaction = true;
            
        
        reg_cr_Cn=linear_regr( pupil.dia, pupil.t, [params.trigEvent], params.trigTime, trialMask, params );
        
            
%% regression with c(n+1)
        params.trigEvent=NaN(size(trials.left));
        params.trigEvent(trials.left) = 0;
        params.trigEvent(trials.right) = 1;
% separate trials based on stay (C(n+1) = C(n)) or switch

        stayTrials = [response(2:end)==response(1:end-1);0];
        switchTrials = [((response(2:end) ~= response(1:end-1)) & response(1:end-1)~=0 & response(2:end)~=0);0];
        
        trialIndex = 1:length(response);
        stayIndex = trialIndex(logical(stayTrials));
        switchIndex = trialIndex(logical(switchTrials));
        
        paramsStay.trigTime = cueTime(logical(stayTrials));
        paramsSwitch.trigTime = cueTime(logical(switchTrials));
        
        % regress for C(n) as control
        paramsStay.trigEvent = params.trigEvent(logical(stayTrials));
        
        paramsSwitch.trigEvent = params.trigEvent(logical(switchTrials));
        
        params.xtitle = 'Time from cue (s)';
        params.window = [-3:0.1:5];
        params.nback = 0;       %how many trials back to regress against
        params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
        params.interaction = true;
        
        reg_cr_CnStay=linear_regr( pupil.dia, pupil.t, [paramsStay.trigEvent], paramsStay.trigTime, trialMask, params );
        reg_cr_CnSwitch=linear_regr( pupil.dia, pupil.t, [paramsSwitch.trigEvent], paramsSwitch.trigTime, trialMask, params );
        
        % regress for C(n+1)
%         paramsStay.trigEvent1 = NaN(size(trials.left));
%         paramsStay.trigEvent1(stayIndex) = params.trigEvent(stayIndex+1);
%         paramsSwitch.trigEvent1 = NaN(size(trials.left));
%         paramsSwitch.trigEvent1(switchIndex)=params.trigEvent(switchIndex+1);
        paramsStay.trigEvent1 = params.trigEvent(stayIndex+1);
        paramsSwitch.trigEvent1=params.trigEvent(switchIndex+1);
        reg_cr_CfutureStay=linear_regr( pupil.dia, pupil.t, [paramsStay.trigEvent1], paramsStay.trigTime, trialMask, params );
        reg_cr_CfutureSwitch=linear_regr( pupil.dia, pupil.t, [paramsSwitch.trigEvent1], paramsSwitch.trigTime, trialMask, params );
        
        
%% regression with C(n-1)

% separate trials based on stay (C(n+1) = C(n)) or switch

        stayTrialsback = [0;response(2:end)==response(1:end-1)];
        switchTrialsback = [0;((response(2:end) ~= response(1:end-1)) & response(1:end-1)~=0 & response(2:end)~=0)];
        
        trialIndex = 1:length(response);
        staybackIndex = trialIndex(logical(stayTrialsback));
        switchbackIndex = trialIndex(logical(switchTrialsback));
        
        paramsStay.trigTime = cueTime(logical(stayTrials));
        paramsSwitch.trigTime = cueTime(logical(switchTrials));
        
        % regress for C(n) as control
        paramsStay.trigEvent = params.trigEvent(logical(stayTrials));
        
        paramsSwitch.trigEvent = params.trigEvent(logical(switchTrials));
        
        params.xtitle = 'Time from cue (s)';
        params.window = [-3:0.1:5];
        params.nback = 0;       %how many trials back to regress against
        params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
        params.interaction = true;
        
        reg_cr_CnStay=linear_regr( pupil.dia, pupil.t, [paramsStay.trigEvent], paramsStay.trigTime, trialMask, params );
        reg_cr_CnSwitch=linear_regr( pupil.dia, pupil.t, [paramsSwitch.trigEvent], paramsSwitch.trigTime, trialMask, params );
        
        % regress for C(n+1) 
        paramsStay.trigEvent1 = params.trigEvent(stayIndex+1);
        paramsSwitch.trigEvent1 = params.trigEvent(switchIndex+1);
        
        reg_cr_CfutureStay=linear_regr( pupil.dia, pupil.t, [paramsStay.trigEvent1], paramsStay.trigTime, trialMask, params );
        reg_cr_CfutureSwitch=linear_regr( pupil.dia, pupil.t, [paramsSwitch.trigEvent1], paramsSwitch.trigTime, trialMask, params );
        
        

    end
end