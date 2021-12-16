function bandit_pupilRL_RPE_CK(dataIndex)

nFiles = size(dataIndex,1);
%% reference: Sul et al.
% s(t) = a0 + a1C(t) + a2R(t) + a3X(t) + a4deltaQ(t) + a5Qc(t) + A(t)
% no A(t) first
%% fit the model per session
for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},['bandit',dataIndex.LogFileName{ii}(end-30:end-4),'_beh.mat']));
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
        stats = stats_new;
       
        %% linear regression for whole session
        
        choice = NaN(size(trials.left));
        choice(trialData.response == 2) = 0;
        choice(trialData.response == 3) = 1;
        
        reward = NaN(size(trials.left));
        reward(trialData.response ~= 0 & trials.reward == 0) = 0;
        reward(trialData.response ~= 0 & trials.reward == 1) = 1;
        
        params.trigEvent = choice;
        % dummycode left: 0, right 1
        
        params.trigEvent2= [NaN; choice(1:end-1,1)]; % c(n-1)
         
        %params.trigEvent4=[stats_sim.ql(2:end)-stats_sim.qr(2:end);NaN];
        params.trigEvent3= [NaN;reward(1:end-1)];   % r(n-1)
        
        % dQ
        params.trigEvent4 = stats_new.ql-stats_new.qr; % delta q
        %params.trigEvent5=[stats_sim.qr(2:end) + stats_sim.ql(2:end);NaN];
        params.trigEvent5= stats_new.rpe;
        
        % delta K
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
        params.ifplot = 1;
        %only perform analysis on this subset of trials
        fieldname={'go'};
        trialMask = getMask(trials,fieldname);
        
        tlabel={'C(n)','QL-QR','C(n-1)','R(n-1)', 'RPE', 'Average reward','Cumulative reward'};
         %reg_cr_RPE=linear_regr_PR( pupil.dia, pupil.t, RL_event, params.trigTime, trialMask, params, trialData.cueTimes );
%         reg_cr_RPEctrl = linear_regr_ctrl( pupil.dia, pupil.t, RL_event, params.trigTime, trialMask, params, tlabel);
         params.xtitle = {'Time from cue (s)'};
    
        % pupil change
        reg_cr_RPE_change = linear_regr_PR(pupil.resp, pupil.respT, RL_event, params.trigTime, trialMask, params, trialData.cueTimes );
        
 
        %% pos/neg RPE trials
        % separate trials into positive/negative trials
        % outcome of current trial is trivial
        % RPE = rn - chosenQ
        % rn = 1, RPE > 0; rn = 0, RPE < 0
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
        
        params.trigTime = trialData.cueTimes(posIndex);
        params.window = [-1:0.1:5];
        % when align pupil signal to cue

        
        %only perform analysis on this subset of trials
        fieldname={'go'};
        trialMask = getMask(trials,fieldname);
        
        tlabel={'C(n)','QL-QR','C(n-1)','R(n-1)', 'posRPE','Average reward','Cumulative reward'};

        params.xtitle = {'Time from cue (s)'};

        reg_cr_RPEpos_change = linear_regr_PR(pupil.resp, pupil.respT, RL_event_pos, params.trigTime, trialMask, params, trialData.cueTimes );
        
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
        params.trigTime = trialData.cueTimes(negIndex);
        params.window = [-1:0.1:5];
        % when align pupil signal to cue
        
        %only perform analysis on this subset of trials
        fieldname={'go'};
        trialMask = getMask(trials,fieldname);
        
        tlabel={'C(n)','QL-QR','C(n-1)','R(n-1)', 'negRPE','Cumulative reward'};
 
        params.xtitle = {'Time from cue (s)'};

        reg_cr_RPEneg_change = linear_regr_PR(pupil.resp, pupil.respT, RL_event_neg, params.trigTime, trialMask, params, trialData.cueTimes );
     
        %% save the regression results
        save(saveMLRmatpath_change, 'reg_cr_RPE_change', 'reg_cr_RPEpos_change','reg_cr_RPEneg_change');

        
        close all;
    end
end
end