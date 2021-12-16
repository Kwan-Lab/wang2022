function MP_pupilRL_delta_MLR_withCK(dataIndex)

% separate the trials into two part:
% 1. choose the high value side
% 2. choose the low value side
% value = beta*dQ + beta_K*dK
nFiles = size(dataIndex,1);
%% reference: Sul et al.
% s(t) = a0 + a1C(t) + a2R(t) + a3X(t) + a4deltaQ(t) + a5Qc(t) + A(t)
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
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
        saveMLRmatpath = fullfile(savematpath,[fn_beh.name(1:end-7),'regRL_change_CK_value.mat']);
        % saveMLRmatpath_outcome = fullfile(dataIndex.BehPath{ii},[fn_beh.name(1:end-7),'regRL_lag0_outcome_cut_fitall.mat']);
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
        
        
        
        
          value = stats_new.beta*(stats_new.ql-stats_new.qr)+stats_new.betac*(stats_new.ckl-stats_new.ckr);
        
            params.trigEvent=NaN(size(trials.left));
            params.trigEvent(trials.left) = 0;
            params.trigEvent(trials.right) = 1;
            
            params.trigEvent2 = [NaN;params.trigEvent(1:end-1)];  %C(n-1)
            params.trigEvent3 = [NaN; NaN; params.trigEvent(1:end-2)];  %c(n-2);
            %second predictor is outcome; dummy-code: reward=1, error=0, miss=NaN
            params.trigEvent4=NaN(size(trials.go));
            params.trigEvent4(trials.reward) = 1;
            params.trigEvent4(trials.noreward) = 0;
           
            params.trigEvent5 = [NaN; params.trigEvent4(1:end-1)]; % R(n-1)
            params.trigEvent6 = [NaN; NaN; params.trigEvent4(1:end-2)];  %R(n-2)
            params.trigEvent7 = params.trigEvent .* params.trigEvent4;  %X(n)
            params.trigEvent8 = params.trigEvent2 .* params.trigEvent5;  %X(n-1);
            params.trigEvent9 = params.trigEvent3 .* params.trigEvent6;  %X(n-1);
            
            params_Mat = [];
            fields = fieldnames(params);
            for jj = 1:length(fields)
                params_Mat = [params_Mat params.(fields{jj})];
            end
            % params.trigEvent2(trials.omitreward) = 0;
            % params.trigEvent2(trials.doublereward) = 2;
            trials.highvalue = (trials.left & value >= 0) | (trials.right & value < 0);
            trials.lowvalue = (trials.left & value <= 0) | (trials.right & value > 0);
            
            
            
            params.trigTime = trialData.cueTimes;
            params.xtitle = 'Time from cue (s)';
            params.window = [-3:0.1:5];
            params.nback = 0;       %how many trials back to regress against
            params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
            params.interaction = false;

           tlabel={'C(n)','C(n-1)','C(n-2)','R(n)','R(n-1)','R(n-2)','C(n)xR(n)','C(n-1)xR(n-1)','C(n-2)xR(n-2)'};

           % MP_plot_regrcoef_pupil(reg_cr,params.pvalThresh,tlabel,params.xtitle);

            %print(gcf,'-dpng','MLR-choiceoutcome_cut_interaction');    %png format
            % saveas(gcf, 'MLR-choiceoutcome_cut_interaction', 'fig');
      
            %regression for pupil change
           params.window = [-3:0.1:5];
           
           % high value regression
           fieldname={'highvalue'};
           trialMask = getMask(trials,fieldname);
           reg_cr_change_high =linear_regr( pupil.resp, pupil.respT, [params.trigEvent params.trigEvent2], params.trigTime, trialMask, params );
          % MP_plot_regrcoef_pupil(reg_cr_change,params.pvalThresh,tlabel,params.xtitle);
          
            fieldname={'lowvalue'};
           trialMask = getMask(trials,fieldname);
           reg_cr_change_low =linear_regr( pupil.resp, pupil.respT, [params.trigEvent params.trigEvent2], params.trigTime, trialMask, params );
       
      
        %% running control multilinear regression
        % shuffle every factor one by one, keeping other factors intact
        
        % check
       
            % construct every regression factor
            params_ctrl.trigEvent1 = params.trigEvent;  % C(n)
            params_ctrl.trigEvent2 = [NaN;params_ctrl.trigEvent1(1:end-1)];  %C(n-1)
            params_ctrl.trigEvent3 = [NaN; NaN; params_ctrl.trigEvent1(1:end-2)];  %c(n-2);
            params_ctrl.trigEvent4 = params.trigEvent2;  % R(n);
            params_ctrl.trigEvent5 = [NaN; params_ctrl.trigEvent4(1:end-1)]; % R(n-1)
            params_ctrl.trigEvent6 = [NaN; NaN; params_ctrl.trigEvent4(1:end-2)];  %R(n-2)
            params_ctrl.trigEvent7 = params_ctrl.trigEvent1 .* params_ctrl.trigEvent4;  %X(n)
            params_ctrl.trigEvent8 = params_ctrl.trigEvent2 .* params_ctrl.trigEvent5;  %X(n-1);
            params_ctrl.trigEvent9 = params_ctrl.trigEvent3 .* params_ctrl.trigEvent6;  %X(n-1);
            
            % concatenate it into a matrix
            params_ctrlMat = [];
            fields = fieldnames(params_ctrl);
            for jj = 1:length(fields)
                params_ctrlMat = [params_ctrlMat params_ctrl.(fields{jj})];
            end
            
           
        
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
             % shul
           
       
            tlabel={'C(n)','C(n-1)','C(n-2)','R(n)','R(n-1)','R(n-2)','C(n)xR(n)','C(n-1)xR(n-1)','C(n-2)xR(n-2)'};
            
            fieldname={'highvalue'};
            trialMask = getMask(trials,fieldname);
            reg_cr_change_high_ctrl = linear_regr_ctrl(pupil.resp, pupil.respT, params_ctrlMat, params.trigTime, trialMask, params, tlabel);
            
            fieldname={'lowvalue'};
            trialMask = getMask(trials,fieldname);
            reg_cr_change_low_ctrl = linear_regr_ctrl(pupil.resp, pupil.respT, params_ctrlMat, params.trigTime, trialMask, params, tlabel);
            
            save(saveMLRmatpath, 'reg_cr_change_high','reg_cr_change_high_ctrl','reg_cr_change_low','reg_cr_change_low_ctrl');
            
            close all;
    end
end
end