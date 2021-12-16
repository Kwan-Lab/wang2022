function bandit_pupilRL_RPE(dataIndex)

nFiles = size(dataIndex,1);
%% reference: Sul et al.
% s(t) = a0 + a1C(t) + a2R(t) + a3X(t) + a4deltaQ(t) + a5Qc(t) + A(t)
% no A(t) first
%% fit the model per session
for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},['bandit',dataIndex.LogFileName{ii}(end-30:end-4),'_beh.mat']));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    fn_latent = fullfile(dataIndex.BehPath{ii}, [dataIndex.LogFileName{ii}(1:end-4),'_latentV.mat']);
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
        saveMLRmatpath = fullfile(savematpath,[fn_beh.name(1:end-7),'regRL_RPE.mat']);
        saveMLRmatpath_change = fullfile(savematpath,[fn_beh.name(1:end-7),'regRL_RPE_change.mat']);  % regression for pupil change
        savepupilfigpath = fullfile(dataIndex.BehPath{ii},[date(1:6),'_figs-pupil']);
        if ~exist(savepupilfigpath,'dir')
            mkdir(savepupilfigpath);
        end
        
        cd(savepupilfigpath);
        
        % simulate to get latent varibles
        
        
        params=[];
        stats = stats_new;
        %first predictor is action value left
        % find the positive and negative RPE (reward > 0)
        
        %% linear regression for whole session
        posIndex = stats.r > 0;
        negIndex = stats.r == 0;
        params.trigEvent = stats.c(:,1);
        % dummycode left: 0, right 1
        params.trigEvent(params.trigEvent == -1) = 0;
        
        params.trigEvent2 = stats_new.ql-stats_new.qr; % delta q
        
        params.trigEvent3= [NaN; params.trigEvent(1:end-1,1)]; % c(n-1)
         
        %params.trigEvent4=[stats_sim.ql(2:end)-stats_sim.qr(2:end);NaN];
        params.trigEvent4= [NaN;stats.r(1:end-1)];   % r(n-1)
        %second predictor is action value right
        %params.trigEvent5=[stats_sim.qr(2:end) + stats_sim.ql(2:end);NaN];
        params.trigEvent5= stats_new.rpe;
         params.trigEvent6 = NaN(size(trials.go));
        for kk = 1:length(trials.left)
            if kk <= 20
                params.trigEvent6(kk) = sum(trials.reward(1:kk))/kk;
            else
                params.trigEvent6(kk) = sum(trials.reward(kk-20:kk))/20;
            end
        end
        
        params.trigEvent7=NaN(size(trials.left));
        for kk = 1:length(trials.left)
            params.trigEvent7(kk) = sum(trials.reward(1:kk));
        end
        params.trigEvent7 = (params.trigEvent7-min(params.trigEvent7))/(max(params.trigEvent7)-min(params.trigEvent7));
        
            
         % make matrix
        RL_event = concat_event(params);
        
        %third predictor is rpe
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
        
        tlabel={'C(n)','QL-QR','C(n-1)','R(n-1)', 'RPE', 'Average reward','Cumulative reward'};
         reg_cr_RPE=linear_regr( pupil.dia, pupil.t, RL_event, params.trigTime, trialMask, params );
%         reg_cr_RPEctrl = linear_regr_ctrl( pupil.dia, pupil.t, RL_event, params.trigTime, trialMask, params, tlabel);
         params.xtitle = {'Time from cue (s)'};
         %MP_plot_regrRL_pupil(reg_cr_RPE,params.pvalThresh,tlabel,params.xtitle);
         %print(gcf,'-dpng','MLR-RL_RPE_fitall');    %png format
         %saveas(gcf, 'MLR-RL_RPE_fitall', 'fig');
        
        
        % pupil change
        params.window = [-3:0.1:5];
        reg_cr_RPE_change = linear_regr(pupil.resp, pupil.respT, RL_event, params.trigTime, trialMask, params );
        
        %% RPE controls    
         params.ifplot = 0;
         tlabel={'C(n)','QL-QR','C(n-1)','R(n-1)', 'RPE', 'Average reward', 'Cumulative reward'};
         reg_cr_RPE_ctrl = linear_regr_ctrl(pupil.dia, pupil.t, RL_event, params.trigTime, trialMask, params, tlabel);
         reg_cr_RPE_change_ctrl = linear_regr_ctrl(pupil.resp, pupil.respT, RL_event,  params.trigTime, trialMask, params, tlabel );

        %% positive RPE trials
        paramsPos.trigEvent = stats.c(posIndex,1);
        % dummycode left: 0, right 1
        paramsPos.trigEvent(paramsPos.trigEvent == -1) = 0;
        
        paramsPos.trigEvent2 = stats_new.ql(posIndex)-stats_new.qr(posIndex); % delta q
        
        preChoice = [NaN; params.trigEvent(1:end-1,1)];
        paramsPos.trigEvent3= preChoice(posIndex); % c(n-1)
         
        %params.trigEvent4=[stats_sim.ql(2:end)-stats_sim.qr(2:end);NaN];
        preReward = [NaN;stats.r(1:end-1)];
        paramsPos.trigEvent4= preReward(posIndex);   % r(n-1)
        %second predictor is action value right
        %params.trigEvent5=[stats_sim.qr(2:end) + stats_sim.ql(2:end);NaN];
        paramsPos.trigEvent5= stats_new.rpe(posIndex);
        paramsPos.trigEvent6 = NaN(size(trials.go));
        for kk = 1:length(trials.left)
            if kk <= 20
                paramsPos.trigEvent6(kk) = sum(trials.reward(1:kk))/kk;
            else
                paramsPos.trigEvent6(kk) = sum(trials.reward(kk-20:kk))/20;
            end
        end
        paramsPos.trigEvent6 =  paramsPos.trigEvent6(posIndex);
        
        paramsPos.trigEvent7=NaN(size(trials.left));
        for kk = 1:length(trials.left)
            paramsPos.trigEvent7(kk) = sum(trials.reward(1:kk));
        end
        paramsPos.trigEvent7 = (paramsPos.trigEvent7-min(paramsPos.trigEvent7))/(max(paramsPos.trigEvent7)-min(paramsPos.trigEvent7));
        paramsPos.trigEvent7 =  paramsPos.trigEvent7(posIndex);
        
          % make matrix
        RL_event_pos = concat_event(paramsPos);
        
        params.trigTime = trialData.cueTimes(posIndex);
        params.window = [-3:0.1:5];
        % when align pupil signal to cue

        
        %only perform analysis on this subset of trials
        fieldname={'go'};
        trialMask = getMask(trials,fieldname);
        
        tlabel={'C(n)','QL-QR','C(n-1)','R(n-1)', 'posRPE','Average reward','Cumulative reward'};
        reg_cr_pos=linear_regr( pupil.dia, pupil.t, RL_event_pos, params.trigTime, trialMask, params );
%         % reg_cr_posctrl=linear_regr_ctrl( pupil.dia, pupil.t, RL_event_pos, params.trigTime, trialMask, params, tlabel );
%         
%         
         params.xtitle = {'Time from cue (s)'};
         %MP_plot_regrRL_pupil(reg_cr_pos,params.pvalThresh,tlabel,params.xtitle);
%         
%        
        %print(gcf,'-dpng','MLR-RL_posRPE_fitall');    %png format
        %saveas(gcf, 'MLR-RL_posRPE_fitall', 'fig');
        
          params.window = [-3:0.1:5];
        reg_cr_RPEpos_change = linear_regr(pupil.resp, pupil.respT, RL_event_pos, params.trigTime, trialMask, params );
        
               %% pos RPE controls
         params.ifplot = 0;
         tlabel={'C(n)','QL-QR','C(n-1)','R(n-1)', 'Positive RPE', 'Average reward', 'Cumulative reward'};
         reg_cr_RPEpos_ctrl = linear_regr_ctrl(pupil.dia, pupil.t, RL_event_pos, params.trigTime, trialMask, params, tlabel);
         reg_cr_RPEpos_change_ctrl = linear_regr_ctrl(pupil.resp, pupil.respT, RL_event_pos,  params.trigTime, trialMask, params, tlabel );

        %% negative RPE trials
        paramsNeg.trigEvent = stats.c(negIndex,1);
        % dummycode left: 0, right 1
        paramsNeg.trigEvent(paramsNeg.trigEvent == -1) = 0;
        
        paramsNeg.trigEvent2 = stats_new.ql(negIndex)-stats_new.qr(negIndex); % delta q
        

        paramsNeg.trigEvent3= preChoice(negIndex); % c(n-1)
         
        %params.trigEvent4=[stats_sim.ql(2:end)-stats_sim.qr(2:end);NaN];

        paramsNeg.trigEvent4= preReward(negIndex);   % r(n-1)
        %second predictor is action value right
        %params.trigEvent5=[stats_sim.qr(2:end) + stats_sim.ql(2:end);NaN];
        paramsNeg.trigEvent5= stats_new.rpe(negIndex);
         paramsNeg.trigEvent6 = NaN(size(trials.go));
        for kk = 1:length(trials.left)
            if kk <= 20
                paramsNeg.trigEvent6(kk) = sum(trials.reward(1:kk))/kk;
            else
                paramsNeg.trigEvent6(kk) = sum(trials.reward(kk-20:kk))/20;
            end
        end
        paramsNeg.trigEvent6 =  paramsNeg.trigEvent6(negIndex);
        paramsNeg.trigEvent7=NaN(size(trials.left));
        for kk = 1:length(trials.left)
            paramsNeg.trigEvent7(kk) = sum(trials.reward(1:kk));
        end
        paramsNeg.trigEvent7 = (paramsNeg.trigEvent7-min(paramsNeg.trigEvent7))/(max(paramsNeg.trigEvent7)-min(paramsNeg.trigEvent7));
        paramsNeg.trigEvent7 =  paramsNeg.trigEvent7(negIndex);
        
          % make matrix
        RL_event_neg = concat_event(paramsNeg);
        params.trigTime = trialData.cueTimes(negIndex);
        params.window = [-3:0.1:5];
        % when align pupil signal to cue
        
        %only perform analysis on this subset of trials
        fieldname={'go'};
        trialMask = getMask(trials,fieldname);
        
        tlabel={'C(n)','QL-QR','C(n-1)','R(n-1)', 'negRPE','Cumulative reward'};
        reg_cr_neg=linear_regr( pupil.dia, pupil.t, RL_event_neg, params.trigTime, trialMask, params );
%         % reg_cr_negctrl=linear_regr_ctrl( pupil.dia, pupil.t, RL_event_neg, params.trigTime, trialMask, params, tlabel );
%         
%         
         params.xtitle = {'Time from cue (s)'};
         %MP_plot_regrRL_pupil(reg_cr_neg,params.pvalThresh,tlabel,params.xtitle);
%         
%        
         
%print(gcf,'-dpng','MLR-RL_negRPE_fitall');    %png format
 %        saveas(gcf, 'MLR-RL_RPE_negfitall', 'fig');
        
          params.window = [-3:0.1:5];
        reg_cr_RPEneg_change = linear_regr(pupil.resp, pupil.respT, RL_event_neg, params.trigTime, trialMask, params );
       %% running control for neg RPE regression
         params.ifplot = 0;
         tlabel={'C(n)','QL-QR','C(n-1)','R(n-1)', 'Negative RPE', 'Average reward', 'Cumulative reward'};
         reg_cr_RPEneg_ctrl = linear_regr_ctrl(pupil.dia, pupil.t, RL_event_neg, params.trigTime, trialMask, params, tlabel);
         reg_cr_RPEneg_change_ctrl = linear_regr_ctrl(pupil.resp, pupil.respT, RL_event_neg,  params.trigTime, trialMask, params, tlabel );
   
        %% updated choosen value

        params.trigEvent5= NaN(length(stats_new.rpe),1);
        params.trigEvent5(stats.c(:,1)==-1) = stats_new.ql(stats.c(:,1)==-1)+stats_new.alpha*stats_new.rpe(stats.c(:,1)==-1);
        params.trigEvent5(stats.c(:,1)==1) = stats_new.ql(stats.c(:,1)==1)+stats_new.alpha*stats_new.rpe(stats.c(:,1)==1);
        params.trigEvent6=NaN(size(trials.left));
        for kk = 1:length(trials.left)
            params.trigEvent6(kk) = sum(trials.reward(1:kk));
        end
        params.trigEvent6 = (params.trigEvent6-min(params.trigEvent6))/(max(params.trigEvent6)-min(params.trigEvent6));
        
        params.trigTime = trialData.cueTimes;
        RL_event_chosenQ = RL_event;
        RL_event_chosenQ(:,5) = params.trigEvent5; % chosen the last factor only
        
        params.window = [-3:0.1:5];
         
       
        %only perform analysis on this subset of trials
        fieldname={'go'};
        trialMask = getMask(trials,fieldname);
        
        tlabel={'C(n)','QL-QR','C(n-1)','R(n-1)', 'updatedQ','Average reward','Cumulative reward'};
         reg_cr_chosenQ=linear_regr( pupil.dia, pupil.t, RL_event_chosenQ, params.trigTime, trialMask, params );
%         reg_cr_updatedQctrl = linear_regr_ctrl( pupil.dia, pupil.t, RL_event_chosenQ, params.trigTime, trialMask, params, tlabel );
%         
%         
         params.xtitle = {'Time from cue (s)'};
         %MP_plot_regrRL_pupil(reg_cr_updatedQ,params.pvalThresh,tlabel,params.xtitle);
%         
%        
         %print(gcf,'-dpng','MLR-RL_updatedQ_fitall');    %png format
         %saveas(gcf, 'MLR-RL_updatdQ_fitall', 'fig');
%         
          params.window = [-3:0.1:5];
        reg_cr_chosenQ_change = linear_regr(pupil.resp, pupil.respT, RL_event_chosenQ, params.trigTime, trialMask, params );
        
        % plot RPE and updated Q together
%         figure;
%         
%         plot(reg_cr_RPE.regr_time,reg_cr_RPE.coeff(:,end),'r.-','MarkerSize',30);
%         hold on;
%         plot(reg_cr_updatedQ.regr_time,reg_cr_updatedQ.coeff(:,end),'b.-','MarkerSize',30);
%         xlim([reg_cr_RPE.regr_time(1) reg_cr_RPE.regr_time(end)]);
%         title('Coefficient of RPE and updated chosen value');
%         xlabel('Time from cue');
%         ylabel('Coefficients');
%         legend({'RPE','Updated Q'});
%         print(gcf,'-dpng','MLR-RL_RPE-updatedQ');    %png format
%         saveas(gcf, 'MLR-RL_RPE-updatdQ', 'fig');
%% running chosen Q controls       
         params.ifplot = 0;
         tlabel={'C(n)','QL-QR','C(n-1)','R(n-1)', 'updatedQ', 'Average reward', 'Cumulative reward'};
         reg_cr_chosenQ_ctrl = linear_regr_ctrl(pupil.dia, pupil.t, RL_event_chosenQ, params.trigTime, trialMask, params, tlabel);
         reg_cr_chosenQ_change_ctrl = linear_regr_ctrl(pupil.resp, pupil.respT, RL_event_chosenQ,  params.trigTime, trialMask, params, tlabel );

        %% save the regression results
         save(saveMLRmatpath, 'reg_cr_RPE','reg_cr_RPE_ctrl', 'reg_cr_pos', 'reg_cr_RPEpos_ctrl',...
            'reg_cr_neg','reg_cr_RPEneg_ctrl','reg_cr_chosenQ','reg_cr_chosenQ_ctrl');

        save(saveMLRmatpath_change, 'reg_cr_RPE_change', 'reg_cr_RPE_change_ctrl', 'reg_cr_RPEpos_change', ...
            'reg_cr_RPEpos_change_ctrl','reg_cr_RPEneg_change','reg_cr_RPEneg_change_ctrl', ...
            'reg_cr_chosenQ_change','reg_cr_chosenQ_change_ctrl');
        % clear params variable
        
        close all;
    end
end
end