function MP_pupilReward(dataIndex)

%% reference to Sul et al.2011
% running multilinear regression 
% z-score = b0 + b1C(n) + b2C(n-1) + b3C(n-1) + b4R(n) + b5R(n-1) +
% b6R(n-2) + b7X(n) + b8X(n-1) + b9X(n-2);  
% C: choice, R: reward, X:interaction



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
        
      
        saveRegName = fullfile(savematpath,[fn_beh.name(1:end-7),'regReward.mat']);
        saveRegName_change = fullfile(savematpath,[fn_beh.name(1:end-7),'regReward_change.mat']);
        
        %% linear regression with cumulative number of rewards/average reward rate

% try only cumulative number of rewards/average reward rate first
% compare these different regresson with mean squared error
        params=[];
        
    
        params.trigTime = trialData.cueTimes;

        % cumulative number of rewards
        params.trigEvent=NaN(size(trials.left));
        for kk = 1:length(trials.left)
            params.trigEvent(kk) = sum(trials.reward(1:kk));
        end
        %params.trigEvent = (params.trigEvent-nanmean(params.trigEvent))/nanstd(params.trigEvent);
        % try normalization
        params.trigEvent = (params.trigEvent - min(params.trigEvent))/max(params.trigEvent);
        
        %cumulavtive average reward rate
        params.trigEvent2=NaN(size(trials.go));
        for kk = 1:length(trials.left)
            params.trigEvent2(kk) = sum(trials.reward(1:kk))/kk;
        end
        params.trigEvent2 = (params.trigEvent2-nanmean(params.trigEvent2))/nanstd(params.trigEvent2);
        
        % running average reward rate, try 10/20/40 trials
        params.trigEvent3 = NaN(size(trials.go));
        for kk = 1:length(trials.left)
            if kk <= 10
                params.trigEvent3(kk) = sum(trials.reward(1:kk))/kk;
            else
                params.trigEvent3(kk) = sum(trials.reward(kk-10:kk))/10;
            end
        end
        params.trigEvent3 = (params.trigEvent3-nanmean(params.trigEvent3))/nanstd(params.trigEvent3);
        
        params.trigEvent4 = NaN(size(trials.go));
        for kk = 1:length(trials.left)
            if kk <= 20
                params.trigEvent4(kk) = sum(trials.reward(1:kk))/kk;
            else
                params.trigEvent4(kk) = sum(trials.reward(kk-20:kk))/20;
            end
        end
        params.trigEvent4 = (params.trigEvent4-nanmean(params.trigEvent4))/nanstd(params.trigEvent4);
        
        params.trigEvent5 = NaN(size(trials.go));
        for kk = 1:length(trials.left)
            if kk <= 40
                params.trigEvent5(kk) = sum(trials.reward(1:kk))/kk;
            else
                params.trigEvent5(kk) = sum(trials.reward(kk-40:kk))/40;
            end
        end
        params.trigEvent5 = (params.trigEvent5-nanmean(params.trigEvent5))/nanstd(params.trigEvent5);
        
        fieldname={'go'};
        trialMask = getMask(trials,fieldname);
        
        
        params.xtitle = 'Time from cue (s)';
        params.window = [-3:0.1:5];
        params.nback = 0;       %how many trials back to regress against
        params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
        params.interaction = true;
        
        
        reg_crR1=linear_regr( pupil.dia, pupil.t, [params.trigEvent], params.trigTime, trialMask, params );
        reg_crR2=linear_regr( pupil.dia, pupil.t, [params.trigEvent2], params.trigTime, trialMask, params );
        reg_crR3=linear_regr( pupil.dia, pupil.t, [params.trigEvent3], params.trigTime, trialMask, params );
        reg_crR4=linear_regr( pupil.dia, pupil.t, [params.trigEvent4], params.trigTime, trialMask, params );
        reg_crR5=linear_regr( pupil.dia, pupil.t, [params.trigEvent5], params.trigTime, trialMask, params );
        
        
        %regression for pupil change
        params.window = [-3:0.1:5];
        reg_crR1_change=linear_regr( pupil.resp, pupil.respT, [params.trigEvent], params.trigTime, trialMask, params );
        reg_crR2_change=linear_regr( pupil.resp, pupil.respT, [params.trigEvent2], params.trigTime, trialMask, params );
        reg_crR3_change=linear_regr( pupil.resp, pupil.respT, [params.trigEvent3], params.trigTime, trialMask, params );
        reg_crR4_change=linear_regr( pupil.resp, pupil.respT, [params.trigEvent4], params.trigTime, trialMask, params );
        reg_crR5_change=linear_regr( pupil.resp, pupil.respT, [params.trigEvent5], params.trigTime, trialMask, params );
        
        % how to compare?
        % R/amount of variance explained by factors -- anova?
        % MP_plot_regrcoef_pupil(reg_cr_change,params.pvalThresh,tlabel,params.xtitle);
        
  
        %% save the results
        save(saveRegName, 'reg_crR1','reg_crR2','reg_crR3','reg_crR4','reg_crR5');
        save(saveRegName_change,'reg_crR1_change','reg_crR2_change','reg_crR3_change','reg_crR4_change','reg_crR5_change');

    end
end