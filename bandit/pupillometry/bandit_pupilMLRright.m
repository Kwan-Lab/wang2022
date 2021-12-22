function bandit_pupilMLRright(dataIndex)

%% reference to Sul et al.2011
% running multilinear regression 
% z-score = b0 + b1C(n) + b2C(n-1) + b3C(n-1) + b4R(n) + b5R(n-1) +
% b6R(n-2) + b7X(n) + b8X(n-1) + b9X(n-2);  
% C: choice, R: reward, X:interaction



nFiles = size(dataIndex,1);


%% load the behavior files first to get the ITI time
iti_trueTime = cell(0);

for ii = 1:nFiles
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},['bandit',dataIndex.LogFileName{ii}(end-30:end-4),'_beh.mat']));
    load(fullfile(fn_beh.folder,fn_beh.name));
    

    % load pupil files
    date = num2str(dataIndex.DateNumber(ii));
    pup_name = fullfile(dataIndex.BehPath{ii}, ['*',date(1:6),'*_pupright.mat']);
    fn_pup = dir(pup_name);
end


%% go through every session, running MLR separately
for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},['bandit',dataIndex.LogFileName{ii}(end-30:end-4),'_beh.mat']));
    load(fullfile(fn_beh.folder,fn_beh.name));
    

    % load pupil files
    date = num2str(dataIndex.DateNumber(ii));
    pup_name = fullfile(dataIndex.BehPath{ii}, ['*',date(1:6),'*_pupright.mat']);
    fn_pup = dir(pup_name);
    if length(fn_pup) == 1
        
        savepupilfigpath = fullfile(dataIndex.BehPath{ii},[date(1:6),'_figs-pupilright']);
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
        savepupilfigpath = fullfile(dataIndex.BehPath{ii},[date(1:6),'_figs-pupilright']);
        if ~exist(savepupilfigpath,'dir')
            mkdir(savepupilfigpath);
        end
        
        cd(savepupilfigpath);
    
        saveRegName_change = fullfile(savematpath,[fn_beh.name(1:end-7),'regCR_changeright.mat']);  % regression for pupil change
     
        %% linear regression with C(n+1) 
  % C(n-2) - C(n+1)
  % R(n-2) - R(n+1)
  % interaction term
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
            params_future.trigEvent14 = params_future.trigEvent14/sum(trials.reward);
        
          future_event = concat_event(params_future);

        
            params_future.xtitle = 'Time from cue (s)';
            params_future.window = [-1:0.1:5];
            params_future.nback = 0;       %how many trials back to regress against
            params_future.pvalThresh = 0.01;   %p-value for coefficient be considered significant
            params_future.interaction = false;
            params_future.ifplot = 0;
            params_future.trigTime = trialData.cueTimes;
            %only perform analysis on this subset of trials
            
            fieldname={'go'};
            trialMask = getMask(trials,fieldname);
            
            tlabel={'C(n+1)','C(n)','C(n-1)','C(n-2)','R(n+1)','R(n)', 'R(n-1)','R(n-2)',...
                    'C(n+1)*R(n+1)','C(n)*R(n)','C(n-1)*R(n-1)','C(n-2)*R(n-2)','Reward Rate','Cumulative Reward'};
        
            reg_cr_future_change=linear_regr_PR( pupil.resp, pupil.respT, future_event,  params_future.trigTime,trialMask, params_future, trialData.cueTimes);
 
      
        %% save all regression into one mat file
        % save all these in a structure
        save(saveRegName_change, 'reg_cr_future_change');
        close all;
    end
end
end
