function bandit_fittingPerSession(dataIndex,save_path)
% % bandit_fittingPerSession %
%PURPOSE:   Fit various learning algorithms to experimental data, on a
%           per-session basis
%AUTHORS:   H Atilgan and AC Kwan 191212
%
%INPUT ARGUMENTS
%   dataIndex:    a database index table for the sessions to analyze
%   save_path:    path for saving the analysis
%
%OUTPUT ARGUMENTS
%

model{1}.name = 'Win-stay lose-switch';   % text label to refer to the model
model{1}.fun = 'funWSLS'; % the corresponding .m code for the model
model{1}.initpar=0.8;     % initial parameter: [prob_WSLS]
model{1}.lb=0;            % upper bound of parameters
model{1}.ub=1;            % lower bound of parameters

model{2}.name = 'Q_RPE';            % text label to refer to the model
model{2}.fun = 'funQ_RPE';      % the corresponding .m code for the model
model{2}.initpar=[0.2 5];       % initial [alpha beta]
model{2}.lb=[0 0];              % upper bound of parameters
model{2}.ub=[1 inf];            % lower bound of parameters

model{3}.name = 'DQ_RPE';           % text label to refer to the model
model{3}.fun = 'funDQ_RPE';     % the corresponding .m code for the model
model{3}.initpar=[0.5 5 0.2];   % initial [alpha_reward beta alpha_noreward]
model{3}.lb=[0 0 0];            % upper bound of parameters
model{3}.ub=[1 inf 1];          % lower bound of parameters

model{4}.name = 'FQ_RPE';      % text label to refer to the model
model{4}.fun = 'funFQ_RPE';        % the corresponding .m code for the model
model{4}.initpar=[0.5 5];   % initial [alpha_reward beta lambda]
model{4}.lb=[0 0];            % upper bound of parameters
model{4}.ub=[1 inf];          % lower bound of parameters

%%
if ~exist(save_path,'dir')
    mkdir(save_path);
end

%% go through each session

disp('-----------------------------------------------------------');
disp(['--- Fitting models - summary of ', int2str(size(dataIndex,1)) ' sessions']);
disp('-----------------------------------------------------------');

for k=1:numel(model)
    
    disp(['Considering model "' model{k}.fun '"']);

    % Does the analysis file exist?
    fn = dir(fullfile(save_path, [model{k}.fun(4:end) '.mat']));
    if size(fn,1)>0
        answer = questdlg(['There is already .mat file for model fitting for ' model{k}.fun '. Run the analysis again? (may take time)'], ...
            'Run model fitting?', ...
            'Yes','No','Yes');
        if strcmp(answer,'Yes')
            runFit = true;
        else
            runFit = false;
        end
    else
        runFit = true;
    end
    
    if (runFit)
        for j = 1:size(dataIndex,1)
            
            disp(['Processing session # ' int2str(j) '...']);
            
            load(fullfile(dataIndex.BehPath{j},[dataIndex.LogFileName{j}(1:end-4),'_beh.mat']));
            
            trials = value_getTrialMasks(trialData);
            stats = value_getTrialStats(trials, sessionData.nRules);
            stats = value_getTrialStatsMore(stats);
            
            [fitpar{j}, ~, bic{j}, ~]=fit_fun(stats,model{k}.fun,model{k}.initpar,model{k}.lb,model{k}.ub);
            
            session{j} = dataIndex.LogFileName{j};
        end
        
        %save behavioral .mat file
        save(fullfile(save_path, [model{k}.fun(4:end) '.mat']),...
            'session','fitpar','bic');
    
    end
end

end

