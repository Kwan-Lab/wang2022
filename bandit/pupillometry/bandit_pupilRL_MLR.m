function bandit_pupilRL_MLR(dataIndex)

nFiles = size(dataIndex,1);
%% reference: Sul et al.
% s(t) = a0 + a1C(t) + a2R(t) + a3X(t) + a4deltaQ(t) + a5Qc(t) + A(t)
% no A(t) first
%% fit the model per session
for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},['bandit',dataIndex.LogFileName{ii}(end-30:end-4),'_beh.mat']));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    fn_latent = fullfile(dataIndex.BehPath{ii}, [dataIndex.LogFileName{ii}(1:end-4),'_FQlatentV.mat']);
    load(fn_latent);
    stats = stats_new;
    % load pupil files
    date = num2str(dataIndex.DateNumber(ii));
    pup_name = fullfile(dataIndex.BehPath{ii}, ['*',date(1:6),'*_pup.mat']);
    fn_pup = dir(pup_name);
    if length(fn_pup) == 1
        
        savepupilfigpath = fullfile(dataIndex.BehPath{ii},[date(1:6),'_figs-pupil']);
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
        saveMLRmatpath = fullfile(savematpath,[fn_beh.name(1:end-7),'regRL.mat']);
        saveMLRmatpath_change = fullfile(savematpath,[fn_beh.name(1:end-7),'regRL_change.mat']);  % regression for pupil change
        %saveMLRmatpath_outcome = fullfile(dataIndex.BehPath{ii},[fn_beh.name(1:end-7),'regRL_lag0_outcome_cut_fitall.mat']);
        if ~exist(savepupilfigpath,'dir')
            mkdir(savepupilfigpath);
        end
        
        cd(savepupilfigpath);
        
    
        load(fullfile(fn_pup.folder,fn_pup.name));
        
        % make folders to save analysis and plots
        
        if ~exist(savematpath,'dir')
            mkdir(savematpath);
        end
        
        % simulate to get latent varibles
        
        
        params=[];
    
        %first predictor is action value left
        params.trigEvent = stats.c(:,1);
        % dummycode left: 0, right 1
        %C(n)
        params.trigEvent(params.trigEvent == -1) = 0;
        %R(n)
        params.trigEvent2 = stats.r(:,1);
        %C(n)*R(n)
        params.trigEvent3 = params.trigEvent .* params.trigEvent2; % interaction term
        % C(n-1)
        params.trigEvent4= [NaN; params.trigEvent(1:end-1,1)];
        % R(n-1)
        params.trigEvent5 = [NaN;stats.r(1:end-1)];
        % C(n-1)*R(n-1)
        params.trigEvent6 = params.trigEvent4 .* params.trigEvent5;
        
        params.trigEvent7=[stats_new.ql(2:end)-stats_new.qr(2:end);NaN];
        % delta Q
        %params.trigEvent4=(stats_new.ql-stats_new.qr);
        
        % chosen Q
        params.trigEvent8 = NaN(length(stats_new.ql),1);
        params.trigEvent8(stats.c(:,1)==-1) = stats_new.ql(stats.c(:,1)==-1);
        params.trigEvent8(stats.c(:,1) == 1) = stats_new.qr(stats.c(:,1) == 1);
        
        % average reward rate on 20 trials window
        params.trigEvent9 = NaN(size(trials.go));
        for kk = 1:length(trials.left)
            if kk <= 20
                params.trigEvent9(kk) = sum(trials.reward(1:kk))/kk;
            else
                params.trigEvent9(kk) = sum(trials.reward(kk-20:kk))/20;
            end
        end
        
        params.trigEvent10=NaN(size(trials.left));
        for kk = 1:length(trials.left)
            params.trigEvent10(kk) = sum(trials.reward(1:kk));
        end
        params.trigEvent10 = (params.trigEvent10-min(params.trigEvent10))/(max(params.trigEvent10)-min(params.trigEvent10));
        
            
        % make matrix
        RL_event = concat_event(params);
        
        % when align pupil signal to cue
        params.trigTime = trialData.cueTimes;
        params.xtitle = 'Time from cue (s)';
        params.window = [-3:0.1:5];
        params.nback = 0;       %how many trials back to regress against
        params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
        params.interaction = false;
        params.ifplot = 1;
        %only perform analysis on this subset of trials
        fieldname={'go'};
        trialMask = getMask(trials,fieldname);
        tlabel={'C(n)','R(n)','C(n)xR(n)','C(n-1)','R(n-1)', 'C(n-1)xR(n-1)','QL-QR', 'ChosenQ', 'Reward Rate', 'Cumulavtive reward'};
        
         reg_cr=linear_regr( pupil.dia, pupil.t, RL_event, params.trigTime, trialMask, params );
%         reg_cr_ctrl = linear_regr_ctrl( pupil.dia, pupil.t, RL_event, params.trigTime, trialMask, params, tlabel);
         
%         
        params.xtitle = {'Time from cue (s)'};

         %MP_plot_regrcoef_pupil(reg_cr,params.pvalThresh,tlabel,params.xtitle);
 
         %print(gcf,'-dpng','MLR-RL');    %png format
         %saveas(gcf, 'MLR-RL', 'fig');
        
        
         % regression for pupil change
        params.window = [-3:0.1:5];
        reg_cr_change = linear_regr(pupil.resp, pupil.respT, RL_event, params.trigTime,   trialMask, params );
        
         %% running control multilinear regression
        % shuffle every factor one by one, keeping other factors intact
        
        % check
       
            % construct every regression factor
            params_ctrl = params;
            
            % concatenate it into a matrix
            params_ctrlMat = [];
            fields = fieldnames(params_ctrl);
            for jj = 1:length(fields)
                if contains(fields{jj},'trigEvent')
                    params_ctrlMat = [params_ctrlMat params_ctrl.(fields{jj})];
                end
            end
            
            % shul
            fieldname={'go'};
            trialMask = getMask(trials,fieldname);
       
        
            params.xtitle = 'Time from cue (s)';
            params.window = [-3:0.1:5];
            params.nback = 0;       %how many trials back to regress against
            params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
            params.interaction = false;
            params.trigTime = trialData.cueTimes;
            params.ifplot = 0;
            %only perform analysis on this subset of trials
            
            % iterate through all 9 factors, shuffle name one by one to get
            % the control regression for every factor
            tlabel={'C(n)','R(n)','C(n)xR(n)','C(n-1)','R(n-1)', 'C(n-1)xR(n-1)','QL-QR', 'ChosenQ', 'Reward Rate', 'Cumulavtive reward'};
        
        
            reg_cr_ctrl = linear_regr_ctrl(pupil.dia, pupil.t, params_ctrlMat, params.trigTime, trialMask, params, tlabel);
            reg_cr_change_ctrl = linear_regr_ctrl(pupil.resp, pupil.respT, params_ctrlMat,  params.trigTime, trialMask, params, tlabel );
        
            save(saveMLRmatpath, 'reg_cr','reg_cr_ctrl');
            save(saveMLRmatpath_change, 'reg_cr_change', 'reg_cr_change_ctrl');
            close all;
    end
end
end