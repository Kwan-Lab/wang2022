function bandit_fitOneSession(stats,modelName,save_path)
% % bandit_fittingPerSession %
%PURPOSE:   Fit one session with one learning algorithm
%AUTHORS:   AC Kwan 200210
%
%INPUT ARGUMENTS
%   stats:      the stats variable that contains choice/outcome information
%   modelName:  text that describe the name of the model
%   save_path:  path for saving the analysis
%
%OUTPUT ARGUMENTS
%

if strcmp(modelName,'WSLS')
    model.name = 'Win-stay lose-switch';   % text label to refer to the model
    model.fun = 'funWSLS'; % the corresponding .m code for the model
    model.initpar=0.8;     % initial parameter: [prob_WSLS]
    model.lb=0;            % upper bound of parameters
    model.ub=1;            % lower bound of parameters
    
elseif strcmp(modelName,'Q_RPE')
    model.name = 'Q_RPE';            % text label to refer to the model
    model.fun = 'funQ_RPE';      % the corresponding .m code for the model
    model.initpar=[0.2 5];       % initial [alpha beta]
    model.lb=[0 0];              % upper bound of parameters
    model.ub=[1 inf];            % lower bound of parameters
    
elseif strcmp(modelName,'DQ_RPE')
    model.name = 'DQ_RPE';           % text label to refer to the model
    model.fun = 'funDQ_RPE';     % the corresponding .m code for the model
    model.initpar=[0.5 5 0.2];   % initial [alpha_reward beta alpha_noreward]
    model.lb=[0 0 0];            % upper bound of parameters
    model.ub=[1 inf 1];          % lower bound of parameters
    
elseif strcmp(modelName,'FQ_RPE')
    model.name = 'FQ_RPE';      % text label to refer to the model
    model.fun = 'funFQ_RPE';        % the corresponding .m code for the model
    model.initpar=[0.5 5];   % initial [alpha_reward beta lambda]
    model.lb=[0 0];            % upper bound of parameters
    model.ub=[1 inf];          % lower bound of parameters
    
elseif strcmp(modelName,'bayesian_baseline')
    model.name = 'bayesian_baseline';      % text label to refer to the model
    model.fun = 'funbayesian_baseline';    % the corresponding .m code for the model
    model.initpar=[0.05];   % initial [alpha_reward beta lambda]
    model.lb=[0];          % upper bound of parameters
    model.ub=[1];          % lower bound of parameters

elseif strcmp(modelName,'bayesian_block')
    model.name = 'bayesian_block';      % text label to refer to the model
    model.fun = 'funbayesian_block';    % the corresponding .m code for the model
    model.initpar=[0.4];   % initial [alpha_reward beta lambda]
    model.lb=[0];          % upper bound of parameters
    model.ub=[1];          % lower bound of parameters
    
elseif strcmp(modelName,'bayesian_twoparams')
    model.name = 'bayesian_twoparams';      % text label to refer to the model
    model.fun = 'funbayesian_twoparams';    % the corresponding .m code for the model
    model.initpar=[0.05 0.05];   % initial [alpha_reward beta lambda]
    model.lb=[0.001 0.001];          % upper bound of parameters
    model.ub=[1 1];          % lower bound of parameters

end

%%
if ~exist(save_path,'dir')
    mkdir(save_path);
end

%% go through each session

disp(['Considering model "' model.fun '"']);

% Does the analysis file exist?
fn = dir(fullfile(save_path, [model.fun(4:end) '.mat']));
if size(fn,1)>0
    answer = questdlg(['There is already .mat file for model fitting for ' model.fun '. Run the analysis again? (may take time)'], ...
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
    
    [fitpar, ~, bic, ~]=fit_fun(stats,model.fun,model.initpar,model.lb,model.ub)
    
    %save behavioral .mat file
    save(fullfile(save_path, [model.fun(4:end) '.mat']),...
        'fitpar','bic');
    
end

end

