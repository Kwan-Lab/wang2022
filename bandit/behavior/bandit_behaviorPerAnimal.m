function bandit_behaviorPerAnimal(dataIndex,save_path)
% % bandit_behaviorPerAnimal %
%PURPOSE:   Analyze bandit behavior averaged across animals
%AUTHORS:   H Atilgan and AC Kwan 191204
%
%INPUT ARGUMENTS
%   dataIndex:    a database index table for the sessions to analyze
%   save_path:    path for saving the plots
%
%OUTPUT ARGUMENTS
%
setup_figprop;


%%
if ~exist(save_path,'dir')
    mkdir(save_path);
end
cd(save_path)

%% go through each animal
animalList = unique(dataIndex.Animal);

disp('-----------------------------------------------------------');
disp(['--- Analyzing - summary of ', int2str(numel(animalList)) ' animals']);
disp('-----------------------------------------------------------');
nTrial_array = [];
rrate_array = [];
pFNL_array = [];

%for j = 1:numel(animalList)
 for j = 1:length(animalList)
    %which session belong to this one animal
    currAnimalSessions = ismember(dataIndex.Animal,animalList(j));
    
    %concatenate the sessions for this one animal
    [trialData, trials, nRules] = merge_sessions(dataIndex(currAnimalSessions,:));
    nTrial_array = [nTrial_array , trialData.sessionLength];
    rrate_array = [rrate_array, trialData.aveReward];
    stats = value_getTrialStats(trials, nRules);
    stats = value_getTrialStatsMore(stats);
    stats_all{j} = stats; % latent variable for all subjects
    %% plot choice behavior - around switches left to right
    trials_back=10;  % set number of previous trials
    
    sw_output{j}=choice_switch(stats,trials_back);
    
    %% plot choice behavior - around switch high-probability side to low-probability side
    sw_hrside_output{j}=choice_switch_hrside(stats,trials_back);

    %% plot choice behavior - around switches left to right, as a function of the statistics of the block preceding the switch
    L1_ranges=[10 20;10 20;10 20;10 20]; %consider only subset of blocks within the range, for trials to criterion
    L2_ranges=[0 4;5 9;10 14;15 30];      %consider only subset of blocks within the range, for random added number of trials
    sw_random_output{j}=choice_switch_random(stats,trials_back,L1_ranges,L2_ranges);

    sw_hrside_random_output{j}=choice_switch_hrside_random(stats,trials_back,L1_ranges,L2_ranges);    
    
    %% plot tendency to predict upcoming reversal
    %store the x and y variables to be made into a scatter plot
    sw_randomVsChoice{j}.dat=[stats.blockTrialRandomAdded stats.blockPreSwitchWorseChoiceAtSwitch]; 
    sw_randomVsChoice{j}.range{1}=[0 30];
    sw_randomVsChoice{j}.range{2}=[0 0.4];
    sw_randomVsChoice{j}.label{1}={'L_2 (random number of trials added)'};
    sw_randomVsChoice{j}.label{2}={'Fraction of trials';'selecting initial worse option'};    
    
    %% plot percentage of failed nolick periods
    start = 1;
    for kk = 1:length(trialData.sessionLength)
        pFNL = sum(trialData.numNolicks(start:start+trialData.sessionLength(kk)-1) == 5)/(trialData.sessionLength(kk));
        start = start+20+trialData.sessionLength(kk);
        pFNL_array = [pFNL_array, pFNL];
    end
    % get rid of these trials in the analysis?
    
    %% plot lick rates
%     trialType={{'left','reward'},{'left','noreward'},{'right','reward'},{'right','noreward'}};
%     edges=[-0.5:0.02:3];   % edges to plot the lick rate histogram
%     lick_trType{j}=get_lickrate_byTrialType(trialData,trials,trialType,edges);
%     
    %% plot response times
%     valLabel='Response time (s)';
%     trialType={'go','left','right'};
%     edges=[0:0.01:1];
%     respTime_trType{j}=get_val_byTrialType(trialData.rt,trials,trialType,edges,valLabel);
%     
    %% plot ITI
%     valLabel='Inter-trial interval (s)';
%     trialType={'go','reward','noreward'};
%     edges=[0:0.1:30];
%     iti_trType{j}=get_val_byTrialType(trialData.iti,trials,trialType,edges,valLabel);
    
    %% plot logistic regression analysis
%     num_regressor = 10;
%     [lreg{j}, ~, ~, ~]=logreg_RCUC(stats,num_regressor);
    
    %% plot logistic regression analysis - rewarded/unrewarded left/right choices
%    num_regressor = 10;
%    [lreg_LR{j}, ~, ~, ~]=logreg_RCUC_LR(stats,num_regressor);
    %% get the estimated parameter for every session
%     animalData = dataIndex(currAnimalSessions,:);
%     for jj=1:size(animalData,1)
%         load(fullfile(animalData.BehPath{jj}, ['bandit_',animalData.LogFileName{jj}(end-29:end-4),'_beh.mat']));
%        
%         %% fit the model to each session, not working for single session
%     % normalized likelihood (Ito and Doya, PLoS Comp Biol, 2015)
%     
%         choice = zeros(1,length(trials.left));
%         choice(trials.left==1) = -1;
%         choice(trials.right==1) = 1;
%         reward  =trials.reward;
%         stats_fit.c = choice;
%         stats_fit.r = double(reward');
%         fitpar = struct;
%         bic = struct;
%         nlike = struct;
% %     fitpar = cell(0);
% %     bic = cell(0);
%         for kk=1:5
%             if isfield(model{kk},'lb')
%                 [fitpar.(model{kk}.name), ~, bic.(model{kk}.name), nlike.(model{kk}.name)]=fit_fun(stats_fit,model{kk}.fun,model{kk}.initpar,model{kk}.lb,model{kk}.ub);
%             else
%                 [fitpar.(model{kk}.name), ~, bic.(model{kk}.name), nlike.(model{kk}.name)]=fit_fun(stats_fit,model{kk}.fun,model{kk}.initpar);
%             end
%         end
%         if exist('nlike_array','var')
%             fname = fieldnames(nlike);
%             for j=1:numel(fname)  %append for each field
%                 nlike_array.(fname{j})=[nlike_array.(fname{j}); nlike.(fname{j})];
%                 bic_array.(fname{j})=[bic_array.(fname{j}); bic.(fname{j})];
%             end
%             fname = fieldnames(fitpar);
%             for j=1:numel(fname)  %append for each field
%                 fitpar_array.(fname{j})=[fitpar_array.(fname{j}); fitpar.(fname{j})];
%             end
%         else
%             nlike_array = nlike;
%             bic_array = bic;
%             fitpar_array = fitpar;
%         end
%         
%         
%     end
end

tlabel = ['Summary of ' int2str(numel(animalList)) ' animals'];

% plot_switch(sw_output,tlabel,stats.rule_labels);
% print(gcf,'-dpng',fullfile(save_path,'switches'));
% saveas(gcf, fullfile(save_path,'switches'), 'fig');

% plot_switch_hrside(sw_hrside_output,tlabel);
% print(gcf,'-dpng',fullfile(save_path,'switches_hrside'));
% saveas(gcf, fullfile(save_path,'switches_hrside'), 'fig');

% plot_switch_random(sw_random_output,tlabel,stats.rule_labels);
% print(gcf,'-dpng',fullfile(save_path,'switches_random'));
% saveas(gcf, fullfile(save_path,'switches_random'), 'fig');
% plot_switch_hrside_random(sw_hrside_random_output,tlabel);
% print(gcf,'-dpng',fullfile(save_path,'switches_hrside_random'));
% saveas(gcf, fullfile(save_path,'switches_hrside_random'), 'fig');
% plot_switch_random_distn(sw_random_output,tlabel);
% print(gcf,'-dpng',fullfile(save_path,'switches_random_distn'));
% saveas(gcf, fullfile(save_path,'switches_random_distn'), 'fig');
% 
% plot_binxaveragey(sw_randomVsChoice,tlabel);
% print(gcf,'-dpng',fullfile(save_path,'switches_randomVsChoice'));
% saveas(gcf, fullfile(save_path,'switches_randomVsChoice'), 'fig');

% plot_lickrate_byTrialType(lick_trType);
% print(gcf,'-dpng',fullfile(save_path,'lickrates_byTrialType'));
% saveas(gcf, fullfile(save_path,'lickrates_byTrialType'), 'fig');
% 
% plot_val_byTrialType(respTime_trType);
% print(gcf,'-dpng',fullfile(save_path,'rt_byTrialType'));
% saveas(gcf, fullfile(save_path,'rt_byTrialType'), 'fig');
% 
% plot_val_byTrialType(iti_trType);
% print(gcf,'-dpng',fullfile(save_path,'iti_byTrialType'));
% saveas(gcf, fullfile(save_path,'iti_byTrialType'), 'fig');
% 
% plot_logreg(lreg,tlabel);
% print(gcf,'-dpng',fullfile(save_path,'logreg'));    
% saveas(gcf, fullfile(save_path,'logreg'), 'fig');

%plot_logreg(lreg_LR,tlabel);
%print(gcf,'-dpng',fullfile(save_path,'logreg_LR'));    
%saveas(gcf, fullfile(save_path,'logreg_LR'), 'fig');

%%  model fitting
disp(['Total number of trials: ' int2str(sum(~isnan(stats.c(:,1))))]);
disp(['Model fit, using concatenated choice behaviors:']);

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
    
    for hh = 1:length(animalList)
        stats_fit.c = stats_all{hh}.c(:,1);
        stats_fit.r = stats_all{hh}.r;
        fitpar = struct;
        bic = struct;
        nlike = struct;
        %     fitpar_array = cell(0);
        %     bic = cell(0);
        for kk=1:5
            if isfield(model{kk},'lb')
                [fitpar.(model{kk}.name), ~, bic.(model{kk}.name), nlike.(model{kk}.name)]=fit_fun(stats_fit,model{kk}.fun,model{kk}.initpar,model{kk}.lb,model{kk}.ub);
            else
                [fitpar.(model{kk}.name), ~, bic.(model{kk}.name), nlike.(model{kk}.name)]=fit_fun(stats_fit,model{kk}.fun,model{kk}.initpar);
            end
        end
    Fitpar{hh} = fitpar;
    Bic{hh} = bic;
    Nlike{hh} = nlike;
 end

    



%% simulate to get latent action value, compare with animal's choice
% player1.label='algo_FQ_RPE';   % change to CA later
% player1.params.a=fitpar_array{1}.FQ_RPE(1);
% player1.params.b=fitpar_array{1}.FQ_RPE(2);
% 
% stats_sim=predictAgent(player1,stats);
% 
% x = 1;  %player 1
% n_plot = 500;   %plot first 500 trials
% plot_session_qparam(stats_sim,x,n_plot);
% for jj = 1:length(animalList)
%     player1.params.a=fitpar_array{jj}.FQ_RPE(1);
%     player1.params.b=fitpar_array{jj}.FQ_RPE(2);
%     stats_simAll{jj} = predictAgent(player1,stats_all{jj});
% end
% plot_session_qhist(stats_simAll,x)


player1.label='algo_FQ_RPE_CK';   % change to CA later
player1.params.a=Fitpar{1}.FQ_RPE_CK(1);
player1.params.b=Fitpar{1}.FQ_RPE_CK(2);
player1.params.ac = Fitpar{1}.FQ_RPE_CK(3);
player1.params.bc = Fitpar{1}.FQ_RPE_CK(4);
stats_sim=predictAgent(player1,stats_all{1});

x = 1;  %player 1
n_plot = 500;   %plot first 500 trials
plot_session_qparam(stats_sim,x,n_plot);

%% plot the histogram of dQ and P_L
% get all latent variables

for jj = 1:length(animalList)
    player1.params.a=Fitpar{jj}.FQ_RPE_CK(1);
    player1.params.b=Fitpar{jj}.FQ_RPE_CK(2);
    player1.params.ac = Fitpar{jj}.FQ_RPE_CK(3);
    player1.params.bc = Fitpar{jj}.FQ_RPE_CK(4);
    stats_simAll{jj} = predictAgent(player1,stats_all{jj});
end
plot_session_qhist(stats_simAll,x)

%% whole session stats (reward rate, session length)
%%
figure;
subplot(2,3,1); hold on;
%plot(rand(1,numel(nTrial_array)),nTrial_array,'k^','MarkerSize',15);
boxplot(nTrial_array,'Colors','k','Notch','off','Labels',[char(956),'=',num2str(round(mean(nTrial_array)))]);
%boxplot(nTrial_array,'Colors','k','Notch','off');
ylim([0 1000]); 
ylabel(['Trials performed']);
set(gca,'box','off') 


subplot(2,3,2); hold on;
%plot(rand(1,numel(rrate_array)),100*rrate_array,'k^','MarkerSize',15);
boxplot(rrate_array*100,'Colors','k','Symbol','k+','Notch','off','Labels',{[char(956),'=',num2str((mean(rrate_array))*100),'%']});
%plot([-1 2],[50 50],'k--','LineWidth',2);
ylim([20 60]);
ylabel('Reward rate (%)');
set(gca,'box','off') 
print(gcf,'-dpng',['summary']);    %png format
saveas(gcf,['summary' ], 'fig');
saveas(gcf, ['summary' ],'svg');

% subplot(2,3,3); hold on;
% %plot(rand(1,numel(rrate_array)),100*rrate_array,'k^','MarkerSize',15);
% boxplot(pFNL_array*100,'Colors','k','Symbol','k+','Notch','off','Labels',{[char(956),'=',num2str((mean(rrate_array))*100),'%']});
% %plot([-1 2],[50 50],'k--','LineWidth',2);
% ylim([0 65]);
% ylabel('Failed nolick rate (%)');
% set(gca,'box','off') 
% print(gcf,'-dpng',['summary']);    %png format
% saveas(gcf,['summary' ], 'fig');
% saveas(gcf, ['summary' ],'svg');

%% plot the alphas and betas and relative value of beta/beta_k
beta = zeros(1,length(animalList));
beta_K = zeros(1,length(animalList));
for ii = 1:length(animalList)
    beta(ii) = Fitpar{ii}.FQ_RPE_CK(2);
    beta_K(ii) = Fitpar{ii}.FQ_RPE_CK(4);
end
rel_beta = beta_K./beta;
figure;
h = boxplot(rel_beta,'Colors','k','Symbol','k+','Notch','off');
ylabel('\beta_K/\beta','FontSize',50)
set(h,{'linew'},{4})
set(gca,'box','off');
set(gca,'linewidth',4);
%ylim([0 0.6]);
set(gca,'YTick',[0:0.2:1.8])
%ylim([-1,1]);

print(gcf,'-dpng',['betaK-beta' int2str(x)]);    %png format
saveas(gcf,['betaK-beta' int2str(x)], 'fig');
saveas(gcf, ['betaK-beta' int2str(x)],'svg');

% alpha_K/alpha
alpha = zeros(1,length(animalList));
alpha_K = zeros(1,length(animalList));
for ii = 1:length(animalList)
    alpha(ii) = Fitpar{ii}.FQ_RPE_CK(1);
    alpha_K(ii) = Fitpar{ii}.FQ_RPE_CK(3);
end
rel_alpha = alpha_K./alpha;


figure;
h = boxplot(rel_alpha,'Colors','k','Symbol','k+','Notch','off');
ylabel('\alpha_K/\alpha','FontSize',50)
ylim([0.8 1.6]);
set(gca,'YTick',[0:0.2:1.6]);
set(h,{'linew'},{4})
set(gca,'box','off');
set(gca,'linewidth',4);
print(gcf,'-dpng',['alphaK-alpha' int2str(x)]);    %png format
saveas(gcf,['alphaK-alpha' int2str(x)], 'fig');
saveas(gcf, ['alphaK-alpha' int2str(x)],'svg');

% beta_K + beta
sum_beta = beta_K + beta;

figure;
h = boxplot(sum_beta,'Colors','k','Symbol','k+','Notch','off');
ylabel('\beta_K+\beta','FontSize',50)
ylim([0 6]);
set(gca,'YTick',[0:2:6]);
set(h,{'linew'},{4})
set(gca,'box','off');
set(gca,'linewidth',4);
print(gcf,'-dpng',['beta_K+beta' int2str(x)]);    %png format
saveas(gcf,['beta_K+beta' int2str(x)], 'fig');
saveas(gcf, ['beta_K+beta' int2str(x)],'svg');

% beta / (beta_K + beta)
frac_beta = beta_K./(sum_beta);
figure;
h = boxplot(frac_beta,'Colors','k','Symbol','k+','Notch','off');
ylabel('\beta_K/(\beta_K+\beta)','FontSize',50)
ylim([0 1]);
set(gca,'YTick',[0:0.2:1]);
set(h,{'linew'},{4})
set(gca,'box','off');
set(gca,'linewidth',4);
print(gcf,'-dpng',['frac_beta' int2str(x)]);    %png format
saveas(gcf,['frac_beta' int2str(x)], 'fig');
saveas(gcf, ['frac_beta' int2str(x)],'svg');
% 
% figure;
% h = violinplot(rel_beta,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0],'ShowData',false);
% ylabel('\beta_K/\beta')
% ylim([-2,10]);
% set(h,{'linew'},{2})
% set(gca,'box','off');
% print(gcf,'-dpng',['betaK-beta' int2str(x)]);    %png format
% saveas(gcf,['betaK-beta' int2str(x)], 'fig');
% saveas(gcf, ['betaK-beta' int2str(x)],'svg');
