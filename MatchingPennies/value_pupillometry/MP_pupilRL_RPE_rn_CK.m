function MP_pupilRL_RPE_rn_CK(dataIndex)

nFiles = size(dataIndex,1);
%% reference: Sul et al.
% s(t) = a0 + a1C(t) + a2R(t) + a3X(t) + a4deltaQ(t) + a5Qc(t)
% RPE + r_n regression
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
         saveMLRmatpath_change = fullfile(savematpath,[fn_beh.name(1:end-7),'regRL_RPE_rn_change_CK.mat']);  % regression for pupil change
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
        params.trigEvent3= reward;   % r(n)
        params.trigEvent4= [NaN;stats.r(1:end-1)];   % r(n-1)
        
        % dQ
        params.trigEvent5 = stats_new.ql-stats_new.qr; % delta q
        %RPE
        params.trigEvent6= stats_new.rpe;
        % dK
        params.trigEvent7 = stats_new.ckl-stats_new.ckr;
        % CKE
        % 1-choice choice
        params.trigEvent8 = NaN(size(trials.go));
        params.trigEvent8(stats_new.c==-1) = 1-stats_new.ckl(stats_new.c==-1); % choosing left
        params.trigEvent8(stats_new.c==1) = 1-stats_new.ckr(stats_new.c==1);
        
        params.trigEvent9 = NaN(size(trials.go));
        for kk = 1:length(trials.left)
            if kk <= 20
                params.trigEvent9(kk) = sum(trials.reward(1:kk))/kk;
            else
                params.trigEvent9(kk) = sum(trials.reward(kk-19:kk))/20;
            end
        end
        
        params.trigEvent10=NaN(size(trials.left));
        for kk = 1:length(trials.left)
            params.trigEvent10(kk) = sum(trials.reward(1:kk));
        end
        params.trigEvent10 = params.trigEvent9/sum(trials.reward);
        
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
        
         %tlabel={'C(n)','QL-QR','C(n-1)','R(n-1)', 'RPE', 'QLC-QRC', 'RPE_C', 'Average reward', 'Cumulative reward'};
        % reg_cr_RPE=linear_regr( pupil.dia, pupil.t, RL_event, params.trigTime, trialMask, params );
%         reg_cr_RPEctrl = linear_regr_ctrl( pupil.dia, pupil.t, RL_event, params.trigTime, trialMask, params, tlabel);
%         params.xtitle = {'Time from cue (s)'};
%         MP_plot_regrcoef_pupil(reg_cr_RPE,params.pvalThresh,tlabel,params.xtitle);
%         print(gcf,'-dpng','MLR-RL_RPE_fitall');    %png format
%         saveas(gcf, 'MLR-RL_RPE_fitall', 'fig');
        
        
        % pupil change
        reg_cr_RPE_change = linear_regr_PR(pupil.resp, pupil.respT, RL_event, params.trigTime, trialMask, params, trialData.cueTimes );
%% RPE controls    
%          params.ifplot = 0;
%          tlabel={'C(n)','QL-QR','C(n-1)','R(n-1)', 'RPE', 'QLC-QRC', 'RPE_C', 'Average reward', 'Cumulative reward'};
%          reg_cr_RPE_ctrl = linear_regr_ctrl(pupil.dia, pupil.t, RL_event, params.trigTime, trialMask, params, tlabel);
%          reg_cr_RPE_change_ctrl = linear_regr_ctrl(pupil.resp, pupil.respT, RL_event,  params.trigTime, trialMask, params, tlabel );

%% separate trials into pos/neg RPE trials 
%% examine the coefficients
% positive RPE trials
%         posIndex = stats.r > 0;
%         negIndex = stats.r == 0;
%         
%         % C(n)
%         paramsPos.trigEvent = params.trigEvent(posIndex);
%         % C(n-1)
%         paramsPos.trigEvent2= params.trigEvent2(posIndex);
%         %R(n-1)
%         paramsPos.trigEvent3= params.trigEvent3(posIndex);
%         %dQ
%         paramsPos.trigEvent4 = params.trigEvent4(posIndex); % delta q
%         %RPE
%         paramsPos.trigEvent5= params.trigEvent5(posIndex);
%         %dK
%         paramsPos.trigEvent6 = params.trigEvent6(posIndex);
%         % CKE
%         paramsPos.trigEvent7 = params.trigEvent7(posIndex);
%         %average reward
%         paramsPos.trigEvent8 = params.trigEvent8(posIndex);
%         %cumulative reward
%         paramsPos.trigEvent9=params.trigEvent9(posIndex);
%           % make matrix
%         RL_event_pos = concat_event(paramsPos);
%         
%         paramsPos.trigTime = trialData.cueTimes(posIndex);
%         paramsPos.window = [-1:0.1:5];
%         paramsPos.nback = 0;
%         paramsPos.pvalThresh = 0.01;   %p-value for coefficient be considered significant
%         paramsPos.interaction = false;
%         paramsPos.ifplot = 0;
%         % when align pupil signal to cue
% 
%         
%         %only perform analysis on this subset of trials
%         fieldname={'go'};
%         trialMask = getMask(trials,fieldname);
%         
%         tlabel={'C(n)','C(n-1)','R(n-1)','QL-QR', 'posRPE', 'Average reward', 'Cumulative reward'};
%         %reg_cr_pos=linear_regr( pupil.dia, pupil.t, RL_event_pos, params.trigTime, trialMask, params );
% %         % reg_cr_posctrl=linear_regr_ctrl( pupil.dia, pupil.t, RL_event_pos, params.trigTime, trialMask, params, tlabel );
% %         
% %         
% %         params.xtitle = {'Time from cue (s)'};
% %         MP_plot_regrcoef_pupil(reg_cr_pos,params.pvalThresh,tlabel,params.xtitle);
% %         
% %        
% %         print(gcf,'-dpng','MLR-RL_posRPE_fitall');    %png format
% %         saveas(gcf, 'MLR-RL_posRPE_fitall', 'fig');
%         
%         reg_cr_RPEpos_change = linear_regr_PR(pupil.resp, pupil.respT, RL_event_pos, paramsPos.trigTime, trialMask, paramsPos, trialData.cueTimes );
%         
%         %% pos RPE controls
% %          params.ifplot = 0;
% %          tlabel={'C(n)','QL-QR','C(n-1)','R(n-1)', 'Positive RPE', 'Average reward', 'Cumulative reward'};
% %          reg_cr_RPEpos_ctrl = linear_regr_ctrl(pupil.dia, pupil.t, RL_event_pos, params.trigTime, trialMask, params, tlabel);
% %          reg_cr_RPEpos_change_ctrl = linear_regr_ctrl(pupil.resp, pupil.respT, RL_event_pos,  params.trigTime, trialMask, params, tlabel );
% 
%         %% negative RPE trials
%        % C(n)
%         paramsNeg.trigEvent = params.trigEvent(negIndex);
%         % C(n-1)
%         paramsNeg.trigEvent2= params.trigEvent2(negIndex);
%         %R(n-1)
%         paramsNeg.trigEvent3= params.trigEvent3(negIndex);
%         %dQ
%         paramsNeg.trigEvent4 = params.trigEvent4(negIndex); % delta q
%         %RPE
%         paramsNeg.trigEvent5= params.trigEvent5(negIndex);
%         %dK
%         paramsNeg.trigEvent6 = params.trigEvent6(negIndex);
%         % CKE
%         paramsNeg.trigEvent7 = params.trigEvent7(negIndex);
%         %average reward
%         paramsNeg.trigEvent8 = params.trigEvent8(negIndex);
%         %cumulative reward
%         paramsNeg.trigEvent9=params.trigEvent9(negIndex);
%         
%          % make matrix
%         RL_event_neg = concat_event(paramsNeg);
%         paramsNeg.trigTime = trialData.cueTimes(negIndex);
%         paramsNeg.window = [-1:0.1:5];
%         paramsNeg.nback = 0;
%         paramsNeg.pvalThresh = 0.01;   %p-value for coefficient be considered significant
%         paramsNeg.interaction = false;
%         paramsNeg.ifplot = 0;
%         % when align pupil signal to cue
%         
%         %only perform analysis on this subset of trials
%         fieldname={'go'};
%         trialMask = getMask(trials,fieldname);
%         
%         tlabel={'C(n)','C(n-1)','R(n-1)','QL-QR', 'negRPE', 'Average reward', 'Cumulative reward'};
%          %reg_cr_neg=linear_regr( pupil.dia, pupil.t, RL_event_neg, params.trigTime, trialMask, params );
% %         % reg_cr_negctrl=linear_regr_ctrl( pupil.dia, pupil.t, RL_event_neg, params.trigTime, trialMask, params, tlabel );
% %         
% %         
% %         params.xtitle = {'Time from cue (s)'};
% %         MP_plot_regrcoef_pupil(reg_cr_neg,params.pvalThresh,tlabel,params.xtitle);
% %         
% %        
% %         print(gcf,'-dpng','MLR-RL_negRPE_fitall');    %png format
% %         saveas(gcf, 'MLR-RL_RPE_negfitall', 'fig');
%         
% 
%         reg_cr_RPEneg_change = linear_regr_PR(pupil.resp, pupil.respT, RL_event_neg, paramsNeg.trigTime, trialMask, paramsNeg, trialData.cueTimes );
% 
%         %% running control for neg RPE regression
% %          params.ifplot = 0;
% %          tlabel={'C(n)','QL-QR','C(n-1)','R(n-1)', 'Negative RPE', 'Average reward', 'Cumulative reward'};
% %          reg_cr_RPEneg_ctrl = linear_regr_ctrl(pupil.dia, pupil.t, RL_event_neg, params.trigTime, trialMask, params, tlabel);
% %          reg_cr_RPEneg_change_ctrl = linear_regr_ctrl(pupil.resp, pupil.respT, RL_event_neg,  params.trigTime, trialMask, params, tlabel );
% 
%         %% updated choosen value
% 
% %         params.trigEvent5= NaN(length(stats_new.rpe),1);
% %         params.trigEvent5(stats.c(:,1)==-1) = stats_new.ql(stats.c(:,1)==-1)+stats_new.alpha*stats_new.rpe(stats.c(:,1)==-1);
% %         params.trigEvent5(stats.c(:,1)==1) = stats_new.ql(stats.c(:,1)==1)+stats_new.alpha*stats_new.rpe(stats.c(:,1)==1);
% %         
% %         % updated choice kernel
% %         params.trigEvent7 = NaN(size(trials.go));
% %         params.trigEvent7(stats_new.c==-1) = stats_new.ckl(stats_new.c==-1)+stats_new.alphac*(1-stats_new.ckl(stats_new.c==-1)); % choosing left
% %         params.trigEvent7(stats_new.c==1) = stats_new.ckr(stats_new.c==1)+stats_new.alphac*(1-stats_new.ckr(stats_new.c==1));
% %         
% %         params.trigTime = trialData.cueTimes;
% %         RL_event_chosenQ = RL_event;
% %         RL_event_chosenQ(:,5) = params.trigEvent5; % change the RPE to chosen Q only from the last regression
% %         
% %          
% %        
% %         %only perform analysis on this subset of trials
% %         fieldname={'go'};
% %         trialMask = getMask(trials,fieldname);
% %         
% %         tlabel={'C(n)','QL-QR','C(n-1)','R(n-1)', 'updatedQ', 'QLC-QRC','updatedQC','Average reward', 'Cumulative reward'};
% %          reg_cr_chosenQ=linear_regr( pupil.dia, pupil.t, RL_event_chosenQ, params.trigTime, trialMask, params );
% %          reg_cr_updatedQctrl = linear_regr_ctrl( pupil.dia, pupil.t, RL_event_chosenQ, params.trigTime, trialMask, params, tlabel );
% % %         
% % %         
% % %         params.xtitle = {'Time from cue (s)'};
% % %         MP_plot_regrRL_pupil(reg_cr_updatedQ,params.pvalThresh,tlabel,params.xtitle);
% % %         
% % %        
% % %         print(gcf,'-dpng','MLR-RL_updatedQ_fitall');    %png format
% % %         saveas(gcf, 'MLR-RL_updatdQ_fitall', 'fig');
% % %         
% % 
% %          reg_cr_chosenQ_change = linear_regr(pupil.resp, pupil.respT, RL_event_chosenQ, params.trigTime, trialMask, params );
% %         
% %         % plot RPE and updated Q together
% % %         figure;
% % %         
% % %         plot(reg_cr_RPE.regr_time,reg_cr_RPE.coeff(:,end),'r.-','MarkerSize',30);
% % %         hold on;
% % %         plot(reg_cr_updatedQ.regr_time,reg_cr_updatedQ.coeff(:,end),'b.-','MarkerSize',30);
% % %         xlim([reg_cr_RPE.regr_time(1) reg_cr_RPE.regr_time(end)]);
% % %         title('Coefficient of RPE and updated chosen value');
% % %         xlabel('Time from cue');
% % %         ylabel('Coefficients');
% % %         legend({'RPE','Updated Q'});
% % %         print(gcf,'-dpng','MLR-RL_RPE-updatedQ');    %png format
% % %         saveas(gcf, 'MLR-RL_RPE-updatdQ', 'fig');
% %         
% % 
% % %% ruuning chosen Q regression control
% % 
% %          params.ifplot = 0;
% %          tlabel={'C(n)','QL-QR','C(n-1)','R(n-1)', 'updatedQ','QLC-QRC','updatedQC', 'Average reward', 'Cumulative reward'};
% %          reg_cr_chosenQ_ctrl = linear_regr_ctrl(pupil.dia, pupil.t, RL_event_chosenQ, params.trigTime, trialMask, params, tlabel);
% %          reg_cr_chosenQ_change_ctrl = linear_regr_ctrl(pupil.resp, pupil.respT, RL_event_chosenQ,  params.trigTime, trialMask, params, tlabel );
% 
% %% save the results 
% %         save(saveMLRmatpath, 'reg_cr_RPE','reg_cr_RPE_ctrl', 'reg_cr_pos', 'reg_cr_RPEpos_ctrl',...
% %             'reg_cr_neg','reg_cr_RPEneg_ctrl','reg_cr_chosenQ','reg_cr_chosenQ_ctrl');

        save(saveMLRmatpath_change, 'reg_cr_RPE_change');
        % clear params variable
        
        close all;
    end
end
end