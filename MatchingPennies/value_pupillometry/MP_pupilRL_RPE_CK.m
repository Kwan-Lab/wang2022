function MP_pupilRL_RPE_CK(dataIndex)

nFiles = size(dataIndex,1);
%% reference: Sul et al.
% pupil = b0 + b1c(n) + b2r(n) + b3x(n) + b4c(n-1) + b5r(n-1) + b6x(n-1) +
% b7deltaQ(n) + b8RPE(n) + b9deltaK(n) + b10CKE(n) + b11averageReward +
% b12cumulativeReward(n)

% no A(t) first
%% fit the model per session
for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},[dataIndex.LogFileName{ii}(1:end-4),'_beh_cut.mat']));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    fn_latent = fullfile(dataIndex.BehPath{ii}, [dataIndex.LogFileName{ii}(1:end-4),'_FQRPECKlatentV.mat']);
    load(fn_latent);
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
        %saveMLRmatpath = fullfile(savematpath,[fn_beh.name(1:end-7),'regRL_RPE_CK.mat']);
        saveMLRmatpath_change = fullfile(savematpath,[fn_beh.name(1:end-7),'regRL_RPE_change_CK.mat']);  % regression for pupil change
        savepupilfigpath = fullfile(dataIndex.BehPath{ii},[date(1:6),'_figs-pupil']);
        if ~exist(savepupilfigpath,'dir')
            mkdir(savepupilfigpath);
        end
        
        cd(savepupilfigpath);
        
        % simulate to get latent varibles
        
        
        params=[];
    
        %first predictor is action value left
        % find the positive and negative RPE (reward > 0)
        
        %% linear regression for whole session

     
        % dummycode left: 0, right 1
        choice = NaN(size(trials.left));
        choice(trials.left) = 0;
        choice(trials.right) = 1;
        % dummycode left: 0, right 1
        %C(n)
        params.trigEvent = choice;
        
        
        params.trigEvent2= [NaN; params.trigEvent(1:end-1,1)]; % c(n-1)
         
        %params.trigEvent4=[stats_sim.ql(2:end)-stats_sim.qr(2:end);NaN];
        reward = NaN(size(trials.go));
        reward(trials.reward) = 1;
        reward(trials.noreward) = 0;
        params.trigEvent3= [NaN;reward(1:end-1)];   % r(n-1)
        
        % dQ
        params.trigEvent4 = stats_new.ql-stats_new.qr; % delta q
        %RPE
        params.trigEvent5= stats_new.rpe;
        % dK
        params.trigEvent6 = stats_new.ckl-stats_new.ckr;
        % CKE
        % 1-choice choice
        params.trigEvent7 = NaN(size(trials.go));
        params.trigEvent7(stats_new.c==-1) = 1-stats_new.ckl(stats_new.c==-1); % choosing left
        params.trigEvent7(stats_new.c==1) = 1-stats_new.ckr(stats_new.c==1);
        
        params.trigEvent8 = NaN(size(trials.go));
        for kk = 1:length(trials.left)
            if kk <= 20
                params.trigEvent8(kk) = sum(trials.reward(1:kk))/kk;
            else
                params.trigEvent8(kk) = sum(trials.reward(kk-19:kk))/20;
            end
        end
        
        params.trigEvent9=NaN(size(trials.left));
        for kk = 1:length(trials.left)
            params.trigEvent9(kk) = sum(trials.reward(1:kk));
        end
        params.trigEvent9 = params.trigEvent9/sum(trials.reward);
        
         % make matrix
        RL_event = concat_event(params);
        
        %third predictor is rpe
        params.trigTime = trialData.cueTimes;
         params.xtitle = 'Time from cue (s)';
        params.window = [-1:0.1:5];
        params.nback = 0;       %how many trials back to regress against
        params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
        params.interaction = false;
        params.ifplot = 0;
        %only perform analysis on this subset of trials
        fieldname={'go'};
        trialMask = getMask(trials,fieldname);
        
        
        
        % pupil change
        reg_cr_RPE_change = linear_regr_PR(pupil.resp, pupil.respT, RL_event, params.trigTime, trialMask, params, trialData.cueTimes );

%% separate trials into pos/neg RPE trials - distinguish RPE and current reward
%% examine the coefficients - 
% positive RPE trials
        posIndex = stats.r > 0;
        negIndex = stats.r == 0;
        
        % C(n)
        paramsPos.trigEvent = params.trigEvent(posIndex);
        % C(n-1)
        paramsPos.trigEvent2= params.trigEvent2(posIndex);
        %R(n-1)
        paramsPos.trigEvent3= params.trigEvent3(posIndex);
        %dQ
        paramsPos.trigEvent4 = params.trigEvent4(posIndex); % delta q
        %RPE
        paramsPos.trigEvent5= params.trigEvent5(posIndex);
        %dK
        paramsPos.trigEvent6 = params.trigEvent6(posIndex);
        % CKE
        paramsPos.trigEvent7 = params.trigEvent7(posIndex);
        %average reward
        paramsPos.trigEvent8 = params.trigEvent8(posIndex);
        %cumulative reward
        paramsPos.trigEvent9=params.trigEvent9(posIndex);
          % make matrix
        RL_event_pos = concat_event(paramsPos);
        
        paramsPos.trigTime = trialData.cueTimes(posIndex);
        paramsPos.window = [-1:0.1:5];
        paramsPos.nback = 0;
        paramsPos.pvalThresh = 0.01;   %p-value for coefficient be considered significant
        paramsPos.interaction = false;
        paramsPos.ifplot = 0;
        % when align pupil signal to cue

        
        %only perform analysis on this subset of trials
        fieldname={'go'};
        trialMask = getMask(trials,fieldname);
        
        tlabel={'C(n)','C(n-1)','R(n-1)','QL-QR', 'posRPE', 'Average reward', 'Cumulative reward'};
       
        reg_cr_RPEpos_change = linear_regr_PR(pupil.resp, pupil.respT, RL_event_pos, paramsPos.trigTime, trialMask, paramsPos, trialData.cueTimes );
        
        %% negative RPE trials
       % C(n)
        paramsNeg.trigEvent = params.trigEvent(negIndex);
        % C(n-1)
        paramsNeg.trigEvent2= params.trigEvent2(negIndex);
        %R(n-1)
        paramsNeg.trigEvent3= params.trigEvent3(negIndex);
        %dQ
        paramsNeg.trigEvent4 = params.trigEvent4(negIndex); % delta q
        %RPE
        paramsNeg.trigEvent5= params.trigEvent5(negIndex);
        %dK
        paramsNeg.trigEvent6 = params.trigEvent6(negIndex);
        % CKE
        paramsNeg.trigEvent7 = params.trigEvent7(negIndex);
        %average reward
        paramsNeg.trigEvent8 = params.trigEvent8(negIndex);
        %cumulative reward
        paramsNeg.trigEvent9=params.trigEvent9(negIndex);
        
         % make matrix
        RL_event_neg = concat_event(paramsNeg);
        paramsNeg.trigTime = trialData.cueTimes(negIndex);
        paramsNeg.window = [-1:0.1:5];
        paramsNeg.nback = 0;
        paramsNeg.pvalThresh = 0.01;   %p-value for coefficient be considered significant
        paramsNeg.interaction = false;
        paramsNeg.ifplot = 0;
        % when align pupil signal to cue
        
        %only perform analysis on this subset of trials
        fieldname={'go'};
        trialMask = getMask(trials,fieldname);
        
        tlabel={'C(n)','C(n-1)','R(n-1)','QL-QR', 'negRPE', 'Average reward', 'Cumulative reward'};
        

        reg_cr_RPEneg_change = linear_regr_PR(pupil.resp, pupil.respT, RL_event_neg, paramsNeg.trigTime, trialMask, paramsNeg, trialData.cueTimes );

        
        


%% save the results 

        save(saveMLRmatpath_change, 'reg_cr_RPE_change', 'reg_cr_RPEpos_change','reg_cr_RPEneg_change');
        % clear params variable
        
        close all;
    end
end
end