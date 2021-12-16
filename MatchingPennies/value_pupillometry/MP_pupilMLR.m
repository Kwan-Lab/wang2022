function MP_pupilMLR(dataIndex)

%% reference to Sul et al.2011
% running multilinear regression 
% z-score = b0 + b1C(n+1) + b2C(n) + b3C(n-1) + b4C(n-2)+ b5R(n+1) + b6R(n) +
% b7R(n-1) + b8R(n-2) + b9X(n+1) + b10X(n) + b11X(n-1) + b12X(n-2) + averageReward + cumulativeReward;  
% C: choice, R: reward, X:interaction



nFiles = size(dataIndex,1);


%% load the behavior files first to get the ITI time
% no clear results

% iti_trueTime = cell(0);
% 
% for ii = 1:nFiles
%     fn_beh = dir(fullfile(dataIndex.BehPath{ii},[dataIndex.LogFileName{ii}(1:end-4),'_beh_cut.mat']));
%     load(fullfile(fn_beh.folder,fn_beh.name));
%     
% 
%     % load pupil files
%     date = num2str(dataIndex.DateNumber(ii));
%     pup_name = fullfile(dataIndex.BehPath{ii}, ['*',date(1:6),'*_pup.mat']);
%     fn_pup = dir(pup_name);
%     if length(fn_pup) == 1
%         
%         iti_time = zeros(1, length(trialData.cueTimes)-1);
%     
%         for tt=1:length(trialData.cueTimes)-1
%             iti_time(tt) = trialData.cueTimes(tt+1) - trialData.outcomeTimes(tt);
%         end
%         iti_trueTime{end+1} = iti_time;
%     end
% end


% check the 33, 66 percentile point of every session
% prc33 = cellfun(@(x) prctile(x,33),iti_trueTime);
% prc67 = cellfun(@(x) prctile(x,67), iti_trueTime);
% 
% % save it somewhere?
% iti_all = [iti_trueTime{:}];
% prc33_all = prctile(iti_all,33);
% prc67_all = prctile(iti_all,67);
% 
% % looks like the medians of every session is close to the percentile of all
% % sessions taken together, taken the median
% prc1 = median(prc33);
% prc2 = median(prc67);


%% go through every session, running MLR separately
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
        
      
        %saveRegName = fullfile(savematpath,[fn_beh.name(1:end-7),'regCR.mat']);
        saveRegName_change = fullfile(savematpath,[fn_beh.name(1:end-7),'regCR_change.mat']);  % regression for pupil change
        %saveRegName_ITI =  fullfile(savematpath,[fn_beh.name(1:end-7),'regCR_ITI.mat']);  
        
         %% linear regression with C(n+1), observable variables
  
   %  C(n-2) - C(n+1)
  %  R(n-2) - R(n+1)
  %   interaction term
        choice = NaN(size(trials.left));
        choice(trialData.response == 2) = 0;
        choice(trialData.response == 3) = 1;
              
        params_future.trigEvent1 = [choice(2:end);NaN];  % C(n+1)
        params_future.trigEvent2 = choice; %C(n)        
        params_future.trigEvent3 = [NaN;choice(1:end-1)];   % C(n-1)
        params_future.trigEvent4 = [NaN;NaN;choice(1:end-2);];  % C(n-2)
        
            %second predictor is outcome; dummy-code: reward=1, error=0, miss=NaN
            %params_future.trigEvent4 = [NaN;NaN;params_future.trigEvent(1:end-2)]; % C(n-2)
        
        reward = NaN(size(trials.left));
        reward(trialData.response ~= 0 & trials.reward == 0) = 0;
        reward(trialData.response ~= 0 & trials.reward == 1) = 1;
        
        params_future.trigEvent5 = [reward(2:end);NaN];  % R(n+1)
        params_future.trigEvent6 = reward; % R(n)
        params_future.trigEvent7 = [NaN;reward(1:end-1)]; % R(n-1)
        params_future.trigEvent8 = [NaN;NaN;reward(1:end-2)]; % R(n-2)
            
            % interaction
            params_future.trigEvent9 = params_future.trigEvent1 .* params_future.trigEvent5;
            params_future.trigEvent10 = params_future.trigEvent2 .* params_future.trigEvent6;
            params_future.trigEvent11 = params_future.trigEvent3 .* params_future.trigEvent7;
            params_future.trigEvent12 = params_future.trigEvent4 .* params_future.trigEvent8;
            
            % average reward rate on 20 trials window
            params_future.trigEvent13 = NaN(size(trials.go));
            for kk = 1:length(trials.left)
                if kk <= 20
                    params_future.trigEvent13(kk) = sum(trials.reward(1:kk))/kk;
                else
                    params_future.trigEvent13(kk) = sum(trials.reward(kk-19:kk))/20;
                end
            end
            
            % cumulative reward
            params_future.trigEvent14=NaN(size(trials.left));
            for kk = 1:length(trials.left)
                params_future.trigEvent14(kk) = sum(trials.reward(1:kk));
            end
            params_future.trigEvent14 = (params_future.trigEvent14)/sum(trials.reward);
        
        
        
            
            future_event = concat_event(params_future);
            % params.trigEvent2(trials.doublereward) = 2;

            fieldname={'go'};
            trialMask = getMask(trials,fieldname);
       
        
        params.xtitle = 'Time from cue (s)';
        params.window = [-1:0.1:5];
        params.nback = 0;       %how many trials back to regress against
        params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
        params.interaction = false;
        params.ifplot = 0;
        params.trigTime = trialData.cueTimes;
        %only perform analysis on this subset of trials
        
        tlabel={'C(n+1)','C(n)','C(n-1)','C(n-2)','R(n+1)','R(n)', 'R(n-1)','R(n-2)',...
                    'C(n+1)*R(n+1)','C(n)*R(n)','C(n-1)*R(n-1)','C(n-2)*R(n-2)','Reward Rate','Cumulative Reward'};
        
        % reg_cr_future=linear_regr( pupil.dia, pupil.t, future_event, params.trigTime, trialMask, params );
        
        reg_cr_future_change=linear_regr_PR( pupil.resp, pupil.respT, future_event,  params.trigTime,trialMask, params, trialData.cueTimes );
        % plot_regr_pupil(reg_cr,params.pvalThresh,tlabel,params.xtitle);
        % print(gcf,'-dpng','MLR-choiceoutcome_sig_choice0_1');    %png format
        % saveas(gcf, 'MLR-choiceoutcome_sigchoice0_1', 'fig');
        % MP_plot_regrcoef_pupil(reg_cr_future_change,params.pvalThresh,tlabel,params.xtitle);
       
%         print(gcf,'-dpng','MLR-choiceoutcome_cut_c(n+1)');    %png format
%         saveas(gcf, 'MLR-choiceoutcome_cut_c(n+1)', 'fig');
%         
           
        
        
        %% save all regression into one mat file
        % save all these in a structure
        %save(saveRegName, 'reg_cr', 'reg_cr1','reg_cr2','reg_cr3','reg_cr_future','reg_cr_future_ctrl','reg_cr_ctrl');
        save(saveRegName_change,'reg_cr_future_change');
        %save(saveRegName_ITI,'reg_cr1_change_1','reg_cr2_change_1','reg_cr3_change_1');
        close all;
    end
end
end
