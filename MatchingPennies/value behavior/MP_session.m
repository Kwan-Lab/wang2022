function MP_session(BehPath,LogFileName,savematpath)
% % bandit_session %
%PURPOSE:   Preparing to analyze a single session of mouse behavior%
%INPUT ARGUMENTS
%   BehPath:        path for the location of the analysis folder containing
%                   the behavioral .mat file
%   LogFileName:    name of the logfile
%
%OUTPUT ARGUMENTS
%
setup_figprop;
%% load the behavioral data

disp('-----------------------------------------------------------');
disp('--- Analyzing a single behavioral session ');
disp('-----------------------------------------------------------');
disp(['Loading ' LogFileName]);

load(fullfile(BehPath,[LogFileName(1:end-4),'_beh.mat']));

% Get trial information
if trialData.cutPoint ~= 0
    % only save the after cut data in beh_cut.mat
    trials = MP_getTrialMasks(trialData);
    fn = fieldnames(trials);
    for tt = 1:length(fn)
        trials.(fn{tt}) = trials.(fn{tt})(1:trialData.cutPoint);
    end
    stats = MP_getTrialStats(trials);
    
    fn2 = fieldnames(trialData);
    for tt = 1:length(fn2)
        if length(trialData.(fn2{tt})) > 1
            trialData.(fn2{tt}) = trialData.(fn2{tt})(1:trialData.cutPoint);
        end
    end
else
     trials = MP_getTrialMasks(trialData);
     stats = MP_getTrialStats(trials);
end
% What to put as title for some of the figures generated
tlabel=strcat('Subject=',sessionData.subject{1},', Time=',sessionData.dateTime(1),'-',sessionData.dateTime(2));

% Create a subfolder to save the images for this session
% folder named using year/month/day of file
yr=num2str(sessionData.dateTime{1}(9:10));
mo=num2str(sessionData.dateTime{1}(1:2));
day=num2str(sessionData.dateTime{1}(4:5));
savebehfigpath = fullfile(BehPath,[yr mo day]);

if ~exist(savebehfigpath)
    mkdir(savebehfigpath)
end
 %% analysis of behavioral performance
    
    % plot choice behavior - whole sessions
    cd(savebehfigpath);
    
    plot_session_game(stats,sessionData.nTrials,tlabel);
    
    %% plot lick rates
    trialType={'reward','noreward'};
    edges=[-2:0.1:5];   % edges to plot the lick rate histogram
    lick_trType=get_lickrate_byTrialType(trialData,trials,trialType,edges);
    MP_plot_lickrate_byTrialType(lick_trType);
    
    %% plot response times
    valLabel='Response time (s)';    
    trialType={'go','left','right'};
    edges=[-1:0.05:2];
    respTime_trType=get_val_byTrialType(trialData.rt,trials,trialType,edges,valLabel);
    plot_val_byTrialType(respTime_trType);
    
    print(gcf,'-dpng','rt_byTrialType');    %png format
    saveas(gcf, 'rt_byTrialType', 'fig');

    %% plot ITI
    valLabel='Inter-trial interval (s)';    
    trialType={'go','reward','noreward'};
    edges=[0:0.25:20];
    iti_dur = [trialData.cueTimes(2:end)-trialData.outcomeTimes(1:end-1); NaN];  %ill-defined ITI for last trial
    iti_trType=get_val_byTrialType(iti_dur,trials,trialType,edges,valLabel);
    plot_val_byTrialType(iti_trType);
    
    print(gcf,'-dpng','iti_byTrialType');    %png format
    saveas(gcf, 'iti_byTrialType', 'fig');

    disp(['Number of ITI > 10 s = ' int2str(sum(iti_dur>10))]);
    disp(['Number of ITI > 20 s = ' int2str(sum(iti_dur>20))]);
    disp(['Number of ITI > 30 s = ' int2str(sum(iti_dur>30))]);
    
    %% fit with logistic regression
    num_regressor=15;
    x=1;    %analyze player 1

    %regressors = rewarded choice, unrewarded choice
    [lregRCUC_output, ~, ~, ~]=logreg_RCUC(stats,x,num_regressor);
    plot_logreg(lregRCUC_output,tlabel);
    print(gcf,'-dpng','logregRCUC');    %png format
    saveas(gcf, 'logregRCUC', 'fig');
    
    %regressors = choice, choice x reward (equivalent to computer's choice)
    [lregCRInt_output, ~, ~, ~]=logreg_CRInt(stats,x,num_regressor);
    plot_logreg(lregCRInt_output,tlabel);
    print(gcf,'-dpng','logregCRInt');    %png format
    saveas(gcf, 'logregCRInt', 'fig');
  
    %% fit to WSLS and Q-learning algorithms
    % normalized likelihood (Ito and Doya, PLoS Comp Biol, 2015)
    
    model{1}.name = 'WSLS';
    model{1}.fun = 'funWSLS';
    model{1}.initpar=0.5; % initial [prob_WSLS]
    model{1}.lb=0;
    model{1}.ub=1;

    
    model{2}.name = 'Q_RPE';
    model{2}.fun = 'funQ_RPE';
    model{2}.initpar=[0.5 10]; % initial [alpha beta]
    model{2}.lb=[0 0];
    model{2}.ub=[1 inf];

    
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
    
    model{5}.name = 'FQ_RPE_CK'; % with a choice autocorrelation term
    model{5}.fun = 'funFQ_RPE_CK';
    model{5}.initpar = [0.1 1 0 0]; % initial [alpha beta tau phi]
    model{5}.lb = [0 0 0 0];
    model{5}.ub = [1 inf 1 inf];
    
    stats_fit.c = stats.c(:,1);
    stats_fit.r = stats.r;
    fitpar = struct;
    bic = struct;
    nlike = struct;
%     fitpar = cell(0);
%     bic = cell(0);
    for kk=1:5
        if isfield(model{kk},'lb')
            [fitpar.(model{kk}.name), ~, bic.(model{kk}.name), nlike.(model{kk}.name)]=fit_fun(stats_fit,model{kk}.fun,model{kk}.initpar,model{kk}.lb,model{kk}.ub);
        else
            [fitpar.(model{kk}.name), ~, bic.(model{kk}.name), nlike.(model{kk}.name)]=fit_fun(stats_fit,model{kk}.fun,model{kk}.initpar);
        end
    end
    
    
    [~, ~, bic.logregCRInt2, nlike.logregCRInt2]=logreg_CRInt(stats,1,2);
    [~, ~, bic.logregCRInt5, nlike.logregCRInt5]=logreg_CRInt(stats,1,5);
    [~, ~, bic.logregCRInt10, nlike.logregCRInt10]=logreg_CRInt(stats,1,10);
    
    %% compare simulated choice behavior to actual choice behavior
    % without cross-validation, this plot may suffer from over-fitting
    
    player1.label='algo_logreg_CRInt';
    player1.params.bias=lregCRInt_output.b_bias;
    player1.params.Ch=lregCRInt_output.b_coeff(:,1);
    player1.params.RC=lregCRInt_output.b_coeff(:,2);

    stats_sim=predictAgent(player1,stats);
    
    plot_session_PrL(stats,stats_sim,sessionData.nTrials,tlabel);
    
    %% calculate entropy
    
    % the possible choice combinations
    choiceBack = 3;
    combos = de2bi([0:2^choiceBack-1],choiceBack);
    combos = 2*(combos - 0.5);  %to make it -1 or 1
    
    % classify animal's choice sequence
    nTrial = size(stats.c,1);
    cumuoccur = zeros(2^choiceBack,nTrial);   %how many times each combo occurred, cumulatively
    for j = choiceBack:nTrial
        c = stats.c(j-choiceBack+1:j,1)';
        idx = ismember(combos,c,'rows');   %find if it matches one of the combos
        if sum(idx)==1
            cumuoccur(:,j) = cumuoccur(:,j-1) + idx;  % update the cumulative occurrence matrix
        else
            cumuoccur(:,j) = cumuoccur(:,j-1);
        end
    end
    
    % change the order of the sequence
    newCumu = cumuoccur;
    newCumu(2,:) = cumuoccur(5,:);
    newCumu(3,:) = cumuoccur(3,:);
    newCumu(4,:) = cumuoccur(2,:);
    newCumu(5,:) = cumuoccur(7,:);
    newCumu(6,:) = cumuoccur(6,:);
    newCumu(7,:) = cumuoccur(4,:);
    figure; hold on;
    colorCode = {[1 0 0 1],[1 0 0 0.7],[1 0 0 0.7],[1 0 0 0.7], [0 0 1 0.7], [0 0 1 0.7], [0 0 1 0.7], [0 0 1 1]};
    linestyle = {'-','--',':','-.','-.',':','--','-'};
    for k=1:2^choiceBack
        h=plot(newCumu(k,:),linestyle{k}, 'color',colorCode{k},'LineWidth',3);
    end
    ax.YAxis.FontSize = 30;
    ax.XAxis.FontSize = 30;
    set(gca,'LineWidth',3);
    xlabel('Trial','FontSize', 28);
    ylabel('Cumulative # of choice patterns','FontSize', 28); 
    Lgd = legend('LLL','LLR','LRL','RLL','LRR','RLR','RRL','RRR','Location','Northwest','FontSize', 28);
    set(Lgd,'EdgeColor','none');
    print(gcf,'-dpng','session-entropy');    %png format
    saveas(gcf, 'session-entropy', 'fig');
    saveas(gcf, 'sesson-entropy','svg');
    p = cumuoccur(:,end)/sum(cumuoccur(:,end));
    entro = -1*sum(p.*log2(p));
    
    
    %%
    save([savematpath,'_beh_cut.mat'],...
            'trialData','sessionData','logfileData','trials',...
            'lregRCUC_output','lregCRInt_output',...
            'lick_trType','iti_trType','respTime_trType',...
            'fitpar','bic','nlike',...
            'entro',...
            'stats');

    close all;
    clearvars -except i dirs expData;

end