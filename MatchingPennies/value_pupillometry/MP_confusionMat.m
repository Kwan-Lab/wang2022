function MP_confusionMat(save_path)

% calculate confusion matrix
model{1}.name = 'WSLS';   % text label to refer to the model
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
model{4}.fun = 'funFQ_RPE';    % the corresponding .m code for the model
model{4}.initpar=[0.5 5];      % initial [alpha_reward beta]
model{4}.lb=[0 0];             % upper bound of parameters
model{4}.ub=[1 inf];           % lower bound of parameters

% model{5}.name = 'FQ_CA'; % with a choice autocorrelation term
% model{5}.fun = 'funFQ_CA';
% model{5}.initpar = [0.1 0 0.1 0.1]; % initial [alpha beta tau phi]
% model{5}.lb = [0 0 -inf -inf];
% model{5}.ub = [1 inf inf inf];

numModel = length(model);

confusionMat = zeros(numModel);
gpuConf = gpuArray(confusionMat);
%%
if ~exist(save_path,'dir')
    mkdir(save_path);
end

n_stim=gpuArray(5000);    % number of trials to simulate
n_iter = 500;
for ii = 1:n_iter
    if mod(ii,10) == 0
        ii
    end
%% generate simulated data using different models
    for jj = 1:numModel
        player1.label=model{jj}.name;
        if jj == 1 % WSLS
            player1.params.p=rand(1);     % this should be generated from a random distribution     
        elseif jj == 2 % Q_RPE
            player1.params.a = rand(1);
            player1.params.b = exprnd(1);
        elseif jj == 3 % DQ_RPE
            player1.params.a = rand(1);
            player1.params.b = exprnd(1);
            player1.params.a2 = rand(1);
        elseif jj == 4 % FQ_RPE
            player1.params.a = rand(1);
            player1.params.b = exprnd(1);
%         elseif jj ==5
%             player1.params.a = rand(1);
%             player1.params.b = exprnd(1);
%             player1.params.tau = rand(1);
%             
        end
        player2.label='algo2';
        player2.params.trial_back=4;   % number of trial back to calculate conditional probabilities
        player2.params.trial_history=400;   %trials older than this number are not considered
        
        tlabel=['WSLS versus Algorithm 2, n=' num2str(n_stim) ' trials'];
    
    
        % generate data
        stats=simPennies(player1,player2,n_stim);
    
        % fit the data
        bic = zeros(1, numModel);
        stats_fit.c = stats.c(:,1);
        stats_fit.r =stats.r(:,1);
        for kk = 1:numModel
            

            if isfield(model{kk},'lb')
                [fitpar{kk}, ~, bic(1,kk), ~]=fit_fun(stats_fit,model{kk}.fun,model{kk}.initpar,model{kk}.lb,model{kk}.ub);
            else
                [fitpar{kk}, ~, bic(1,kk), ~]=fit_fun(stats_fit,model{kk}.fun,model{kk}.initpar);
            end  
            
        end
        
        % compare BIC, find the best fit model
        [~,I] = min(bic);
        
        % update confusion matrix
        gpuConf(jj,I) = gpuConf(jj,I) + 1;
        
    
    end    
end

gpuConf= gpuConf/ n_iter;
figure;
heatmap(gpuConf);
ax = gca;
ax.XData = ['WSLS' 'Q_RPE' 'DQ_RPE' 'FQ_RPE'];
ax.YData = ['WSLS' 'Q_RPE' 'DQ_RPE' 'FQ_RPE'];
save('confusionMat.mat','confusionMat');
x =1;
end
