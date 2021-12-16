function MP_fittingDrift(dataIndex,save_path)
% %MP_fittingPerAnimal %
%PURPOSE:   Fit learning algorithms with drifting parameters to experimental data, on a
%           per-animal basis by merging sessions from same animal
%reference: J. Wang et al. Nat. Neurosci. 2018

%INPUT ARGUMENTS
%   dataIndex:    a database index table for the sessions to analyze
%   save_path:    path for saving the analysis
%
%OUTPUT ARGUMENTS
% some subject has multiple training periods, fit them separately
%

%% point:
% 1. negative correlation between alpha/beta

%% model{1}.name = 'WSLS';   % text label to refer to the model
% model{1}.fun = 'funWSLS'; % the corresponding .m code for the model
% model{1}.initpar=0.8;     % initial parameter: [prob_WSLS]
% model{1}.lb=0;            % upper bound of parameters
% model{1}.ub=1;            % lower bound of parameters
%
% model{2}.name = 'Q_RPE';            % text label to refer to the model
% model{2}.fun = 'funQ_RPE';      % the corresponding .m code for the model
% model{2}.initpar=[0.2 5];       % initial [alpha beta]
% model{2}.lb=[0 0];              % upper bound of parameters
% model{2}.ub=[1 inf];            % lower bound of parameters
%
% model{3}.name = 'DQ_RPE';           % text label to refer to the model
% model{3}.fun = 'funDQ_RPE';     % the corresponding .m code for the model
% model{3}.initpar=[0.5 5 0.2];   % initial [alpha_reward beta alpha_noreward]
% model{3}.lb=[0 0 0];            % upper bound of parameters
% model{3}.ub=[1 inf 1];          % lower bound of parameters
%
% model{4}.name = 'FQ_RPE';      % text label to refer to the model
% model{4}.fun = 'funFQ_RPE';    % the corresponding .m code for the model
% model{4}.initpar=[0.5 5];      % initial [alpha_reward beta]
% model{4}.lb=[0 0];             % upper bound of parameters
% model{4}.ub=[1 inf];           % lower bound of parameters

model{1}.name = 'FQ_RPE_CK_drift'; % with a choice autocorrelation term
model{1}.fun = 'funFQ_RPE_CK_drift';
model{1}.initpar = [0.1 0 0.1 0.1]; % initial [alpha beta alpha_K beta_K]
model{1}.lb = [0 0 0 0];
model{1}.ub = [1 inf 1 inf];

model{2}.name = 'FQ_RPE_drift';      % text label to refer to the model
model{2}.fun = 'funFQ_RPE_drift';    % the corresponding .m code for the model
model{2}.initpar=[0.5 5];      % initial [alpha_reward beta]
model{2}.lb=[0 0];             % upper bound of parameters
model{2}.ub=[1 4];           % lower bound of parameters

model{3}.name = 'Q_RPE_drift';            % text label to refer to the model
model{3}.fun = 'funQ_RPE_drift';      % the corresponding .m code for the model
model{3}.initpar=[0.2 5];       % initial [alpha beta]
model{3}.lb=[0 0];              % upper bound of parameters
model{3}.ub=[1 inf];
%%
if ~exist(save_path,'dir')
    mkdir(save_path);
end

%% go through each animal
animalFolder = unique(dataIndex.LogFilePath);
for ii = 1:length(animalFolder)
    Ind = strfind(animalFolder{ii},filesep);
    startInd = Ind(end);
    animalList{ii} = animalFolder{ii}(startInd+1:end);
end

disp('-----------------------------------------------------------');
disp(['--- Fitting models - summary of ', int2str(numel(animalList)) ' animals']);
disp('-----------------------------------------------------------');

for k=1:numel(model)
    
    disp(['Considering model "' model{k}.name '"']);
    
    % Does the analysis file exist?
    fn = dir(fullfile(save_path, [model{k}.fun(4:end) '_model_drift.mat']));
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
        for j = 1:numel(animalList)
            
            disp(['Processing animal # ' int2str(j) '...']);
            disp(['   ' num2str(sum(contains(dataIndex.LogFilePath,animalList(j)))) ' sessions associated with this animal']);
            
            %which session belong to this one animal
            currAnimalSessions = contains(dataIndex.LogFilePath,animalList(j));
            
            %concatenate the sessions for this one animal
            stats_all = MP_merge_sessions(dataIndex(currAnimalSessions,:));
            
            stats_fit.c = stats_all.c(:,1);
            stats_fit.r =stats_all.r;
            window  =20; % fit the parameter every 10 trials
            
            %             % find the window length with the smallest BIC
            %             window = 200:50:700;
            %
            %             for kk = 1:length(window)
            numWin = length(stats_fit.r)-window+1;
            
            % initiate data matrix
            if strcmp(model{k}.name, 'FQ_RPE_CK_drift')
                player1.label=['algo_',model{k}.name];
                player1.params_save.a=zeros(1, numWin);
                player1.params_save.b=zeros(1, numWin);
                player1.params_save.ac = zeros(1, numWin);
                player1.params_save.bc = zeros(1, numWin);
            elseif strcmp(model{k}.name, 'FQ_RPE_drift')
                player1.label=['algo_',model{k}.name];
                player1.params_save.a=zeros(1, numWin);
                player1.params_save.b=zeros(1, numWin);
            elseif strcmp(model{k}.name, 'Q_RPE_drift')
                player1.label=['algo_',model{k}.name];
                player1.params_save.a=zeros(1, numWin);
                player1.params_save.b=zeros(1, numWin);
            end
            
            for tt = 1:numWin
                stats_drift.c = stats_fit.c(tt:window+tt-1);
                stats_drift.r = stats_fit.r(tt:window+tt-1);
                if tt == 1
                    stats_drift.ql = 0; stats_drift.qr = 0;
                    stats_drift.pl = 0.5;  stats_drift.rpe = NaN; %stats_drift.pr = 0.5;
                    if strcmp(model{k}.name, 'FQ_RPE_CK_drift')
                        stats_drift.ckl = 0; stats_drift.ckr = 0;
                    end
                else
                    % update the value based on the last trial of the
                    % previous window
                    stats_drift.ql = stats_sim.ql(2); stats_drift.qr = stats_sim.qr(2);
                    stats_drift.pl = stats_sim.pl(2); %stats_drift.pr = stats_sim.pr(2);
                    stats_drift.rpe = stats_sim.rpe(2);
                    if strcmp(model{k}.name, 'FQ_RPE_CK_drift')
                        stats_drift.ckl = stats_sim.ckl(2); stats_drift.ckr = stats_sim.ckr(2);
                    end
                    Q = [stats_drift.pl 1-stats_drift.pl];
                    if strcmp(model{k}.name, 'FQ_RPE_CK_drift')
                        CK = [stats_drift.ckl stats_drift.ckr];
                        V = player1.params.b*Q + player1.params.bc*CK;
                    else
                        V = Q;
                    end
                    
                    p = exp(V)/sum(exp(V));
                    stats_drift.pl = p(1); %stats_drift.pr = p(2);
                    stats_save.rpe(end) = stats_drift.rpe;
                end
                if strcmp(model{k}.name, 'FQ_RPE_CK_drift')
                    stats_drift.latentV = [ stats_drift.qr,stats_drift.ql,stats_drift.ckr,stats_drift.ckl]; % action values & choice kernels from the previous trial
                else
                    stats_drift.latentV = [ stats_drift.qr,stats_drift.ql];
                end
                if isfield(model{k},'lb')
                    [fitpar{tt}, ~, bic{tt}, ~]=fit_fun(stats_drift,model{k}.fun,model{k}.initpar,model{k}.lb,model{k}.ub);
                else
                    [fitpar{tt}, ~, bic{tt}, ~]=fit_fun(stats_drift,model{k}.fun,model{k}.initpar);
                end
                
                % get latent variable estimation
                
                player1.params_save.a(tt) = fitpar{tt}(1);
                player1.params_save.b(tt) = fitpar{tt}(2);
                player1.params.a = fitpar{tt}(1);
                player1.params.b = fitpar{tt}(2);
                if strcmp(model{k}.name, 'FQ_RPE_CK_drift')
                    player1.params_save.ac(tt) = fitpar{tt}(3);
                    player1.params_save.bc(tt) = fitpar{tt}(4);
                    
                    player1.params.ac = fitpar{tt}(3);
                    player1.params.bc = fitpar{tt}(4);
                end
                stats_sim=predictAgent_drift(player1,stats_drift);
                
                if tt == 1
                    fieldname = fieldnames(stats_sim);
                    for uu = 1:length(fieldname)
                        if ~strcmp(fieldname{uu},'currTrial')
                            if strcmp(fieldname{uu},'playerlabel')
                                stats_save.(fieldname{uu})=  stats_sim.(fieldname{uu});
                            elseif strcmp(fieldname{uu},'playerparams')
                                stats_save.(fieldname{uu})=  stats_sim.(fieldname{uu});
                            else
                                stats_save.(fieldname{uu})(1)=  stats_sim.(fieldname{uu})(1);
                            end
                        end
                    end
                else
                    fieldname = fieldnames(stats_sim);
                    for uu = 1:length(fieldname)
                        if ~strcmp(fieldname{uu},'currTrial') && ~strcmp(fieldname{uu},'playerlabel') 
                            stats_save.(fieldname{uu}) = [stats_save.(fieldname{uu}); stats_sim.(fieldname{uu})(1)];
                        elseif  strcmp(fieldname{uu},'playerparams')
                            stats_save.(fieldname{uu}){end+1} = stats_sim.(fieldname{uu});
                        end
                        
                    end
                end
                
            end
            
            
            
            % only fit for FQ_RPE_CK model
            %get the latent variable, save them in sessions
            
            % compare the log likelihood of drifting model and stable
            % model
            
            al = zeros(1, length(stats_save.ql));
            be = zeros(1, length(stats_save.ql));
            
            for zz = 1:length(stats_save.ql)
                al(zz) = stats_save.playerparams{zz}.a;
                be(zz) = stats_save.playerparams{zz}.b;
            end
            stats_save.sLength = stats_all.sessionLength;
            stats_save.playerparams = player1.params_save;
            % save stats_sim
            if strcmp(model{k}.name, 'FQ_RPE_CK_drift')
                save(fullfile(save_path,[animalList{j},'_FQRPECK_drift.mat']),'stats_save');
            elseif strcmp(model{k}.name, 'FQ_RPE_CK_drift')
                save(fullfile(save_path,[animalList{j},'_FQRPE_drift.mat']),'stats_save');
            elseif strcmp(model{k}.name, 'Q_RPE_CK_drift')
                save(fullfile(save_path,[animalList{j},'_QRPE_drift.mat']),'stats_save');
            end
            
            animal{j} = animalList{j};
        end
        
        %save behavioral .mat file
        save(fullfile(save_path, [model{k}.fun(4:end) '_model_drift.mat']),...
            'animal','fitpar','bic');
        
        fitparMat = cell2mat(fitpar');
        disp('Median fitted parameters:');
        nanmedian(fitparMat,1)
    end
    
end

end

