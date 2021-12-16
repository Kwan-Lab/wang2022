function bandit_pupilRL_MLR_CK(dataIndex)

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
    stats = stats_new;
    % load pupil files
    date = num2str(dataIndex.DateNumber(ii));
    pup_name = fullfile(dataIndex.BehPath{ii}, ['*',date(1:6),'*_pup.mat']);
    fn_pup = dir(pup_name);
    if length(fn_pup) == 1
        
        savepupilfigpath = fullfile(dataIndex.BehPath{ii},[date(1:6),'_figs-pupil']);
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
      
        saveMLRmatpath_change = fullfile(savematpath,[fn_beh.name(1:end-7),'regRL_change_CK.mat']);  % regression for pupil change
        if ~exist(savepupilfigpath,'dir')
            mkdir(savepupilfigpath);
        end
        
        cd(savepupilfigpath);
        load(fullfile(fn_pup.folder,fn_pup.name));
        
        % make folders to save analysis and plots
        
        if ~exist(savematpath,'dir')
            mkdir(savematpath);
        end
        
        params=[];
        choice = NaN(size(trials.left));
        choice(trialData.response == 2) = 0;
        choice(trialData.response == 3) = 1;
        
        reward = NaN(size(trials.left));
        reward(trialData.response ~= 0 & trials.reward == 0) = 0;
        reward(trialData.response ~= 0 & trials.reward == 1) = 1;
        
        % dummycode left: 0, right 1
        %C(n)
        params.trigEvent = choice;
        %R(n)
        params.trigEvent2 = reward;
        %C(n)*R(n)
        params.trigEvent3 = params.trigEvent .* params.trigEvent2; % interaction term
        % C(n-1)
        params.trigEvent4= [NaN; choice(1:end-1)];
        % R(n-1)
        params.trigEvent5 = [NaN;reward(1:end-1)];
        % C(n-1)*R(n-1)
        params.trigEvent6 = params.trigEvent4 .* params.trigEvent5;
        
        % delta Q
        params.trigEvent7=stats_new.ql-stats_new.qr;
       
        % chosen Q
        params.trigEvent8 = NaN(length(stats_new.ql),1);
        params.trigEvent8(stats.c(:,1)==-1) = stats_new.ql(stats.c(:,1)==-1);
        params.trigEvent8(stats.c(:,1) == 1) = stats_new.qr(stats.c(:,1) == 1);
        
        % delta K
        params.trigEvent9 = stats_new.ckl - stats_new.ckr;
        
        % chosen K
        params.trigEvent10 = NaN(length(stats_new.ql),1);
        params.trigEvent10(stats.c(:,1)==-1) = stats_new.ckl(stats.c(:,1)==-1);
        params.trigEvent10(stats.c(:,1) == 1) = stats_new.ckr(stats.c(:,1) == 1);
        
        % average reward rate on 20 trials window
        params.trigEvent11 = NaN(size(trials.go));
        for kk = 1:length(trials.left)
            if kk <= 20
                params.trigEvent11(kk) = sum(trials.reward(1:kk))/kk;
            else
                params.trigEvent11(kk) = sum(trials.reward(kk-19:kk))/20;
            end
        end
        
        % cumulative rewards
        params.trigEvent12=NaN(size(trials.left));
        for kk = 1:length(trials.left)
            params.trigEvent12(kk) = sum(trials.reward(1:kk));
        end
        params.trigEvent12 = params.trigEvent12/ sum(trials.reward);

        % make matrix
        RL_event = concat_event(params);
        
        % when align pupil signal to cue
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
        tlabel={'C(n)','R(n)','C(n)xR(n)','C(n-1)','R(n-1)', 'C(n-1)xR(n-1)','QL-QR', 'ChosenQ', 'Reward Rate', 'Cumulavtive reward'};
        
        params.xtitle = {'Time from cue (s)'};
         % regression for pupil change
        reg_cr_change = linear_regr_PR(pupil.resp, pupil.respT, RL_event, params.trigTime,   trialMask, params, trialData.cueTimes );
        
        
%%  save results
            save(saveMLRmatpath_change, 'reg_cr_change');
            close all;
    end
end
end