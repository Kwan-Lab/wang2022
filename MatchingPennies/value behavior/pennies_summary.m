function pennies_summary(dataIndex)
% summarize behavior of given dataset

% input: dataIndex


close all;

%% setup path and plotting formats

setup_figprop;  %set up default figure plotting parameters
%set(0, 'DefaultFigureRenderer', 'painters');
%% load data file list
animalFolder = unique(dataIndex.LogFilePath);
for ii = 1:length(animalFolder)
    Ind = strfind(animalFolder{ii},filesep);
    startInd = Ind(end);
    animalList{ii} = animalFolder{ii}(startInd+1:end);
end

%% load the data


% subMask = [];
logInd = 1;
for i = 1:length(animalList)   
    currAnimalSessions = ismember(dataIndex.Animal,animalList(i));
    animalData = dataIndex(currAnimalSessions,:);
    
    stats_sub.c=[];  %concatenate choices and outcomes across sessions
    stats_sub.r=[];
    for jj=1:size(animalData,1)
        load(fullfile(animalData.BehPath{jj},[animalData.LogFileName{jj}(1:end-4),'_beh_cut.mat']));


        % calculate ITI time for every session
        iti_time = zeros(1, length(trialData.cueTimes)-1);
        
        for tt=1:length(trialData.cueTimes)-1
            iti_time(tt) = trialData.cueTimes(tt+1) - trialData.outcomeTimes(tt);
        end
        iti_trueTime{logInd} = iti_time;
        iti_num(logInd) = length(trialData.itiTimes);
        trial_num(logInd) = length(trialData.cueTimes);
        
        lick_trType_array{logInd}=lick_trType;
    
        lregRCUC_array{logInd}=lregRCUC_output;
        lregCRInt_array{logInd}=lregCRInt_output;
    
        iti_array{logInd}=iti_trType;
        respTime_array{logInd}=respTime_trType;
        trueRespTime{logInd} = trialData.rt;
        choiceBySession{logInd} = stats;
        
    
        if exist('nlike_array','var')
            fname = fieldnames(nlike);
            for j=1:numel(fname)  %append for each field
                nlike_array.(fname{j})=[nlike_array.(fname{j}); nlike.(fname{j})];
                bic_array.(fname{j})=[bic_array.(fname{j}); bic.(fname{j})];
            end
            fname = fieldnames(fitpar);
            for j=1:numel(fname)  %append for each field
                fitpar_array.(fname{j})=[fitpar_array.(fname{j}); fitpar.(fname{j})];
            end
        else
            nlike_array = nlike;
            bic_array = bic;
            fitpar_array = fitpar;
        end
        
        nTrial_array(logInd)=sum(stats.c(:,1)==-1)+sum(stats.c(:,1)==1);
        entro_array(logInd)=entro;
        rrate_array(logInd)=sum(stats.r==1)/(sum(stats.r==1)+sum(stats.r==0));
        
        stats_sub.c=[stats_sub.c; stats.c];
        stats_sub.r=[stats_sub.r; stats.r];
        
        logInd = logInd + 1;
    end
    stats_all{i} = stats_sub;
    
    close all;
%     clearvars -except i dirs dataIndex ...
%         lick_trType_array lregRCUC_array lregCRInt_array iti_array respTime_array trueRespTime choiceBySession...
%         nlike_array bic_array fitpar_array iti_trueTime...
%         nTrial_array entro_array rrate_array subMask...
%         stats_all, animalList;
end


%%

%savebehfigpath = fullfile(dataIndex.BehPath{1}(1:22),'summary',dataIndex.Animal{1});
% savebehfigpath = fullfile(dataIndex.BehPath{1}(1:22),'summary');
% if ~exist(savebehfigpath,'dir')
%     mkdir(savebehfigpath);
% end
% 
% cd(savebehfigpath);
%save('stats_1.mat', 'choiceBySession','entro_array','respTime_array', 'rrate_array', 'subMask');
% tlabel=strcat('Group summary, n=',int2str(numel(iti_array)), ', subject=',dataIndex.Animal{1});
tlabel=strcat('Group summary, n=',int2str(numel(iti_array)), ', subject=', int2str(length(animalList)));
%% plot the previous ITI and response time
plot_rtITI(iti_trueTime, trueRespTime)
print(gcf,'-dpng','rt_ITI_cut');    %png format
saveas(gcf, 'rt_ITI_cut', 'fig');
saveas(gcf, 'rt_ITI_cut','svg');
%%
MP_plot_lickrate_byTrialType(lick_trType_array);

plot_val_byTrialType(respTime_array);
print(gcf,'-dpng','rt_byTrialType_cut');    %png format
saveas(gcf, 'rt_byTrialType_cut', 'fig');
saveas(gcf, 'rt_byTrialType_cut','svg');

plot_val_byTrialType(iti_array);
print(gcf,'-dpng','iti_byTrialType_cut');   %png format
saveas(gcf, 'iti_byTrialType_cut', 'fig');
saveas(gcf, 'iti_byTrialType_cut','svg');

%% plot iti_array in tertials
iti = [];
for ii = 1:length(iti_trueTime)
    iti = [iti,iti_trueTime{ii}];
end
shortEdge = 5.26;
mediumEdge = 6.50;

figure;
edges1 = 4:0.1:round(shortEdge*10)/10;
histogram(iti(iti<shortEdge),edges1,'EdgeColor','none','FaceColor',[0 0.4470 0.7410],'FaceAlpha',1);
hold on;
edges2 = round(shortEdge*10)/10:0.1:round(mediumEdge*10)/10;
histogram(iti(iti<mediumEdge & iti>=shortEdge),edges2,'EdgeColor','none','FaceColor',[0.8500 0.3250 0.0980],'FaceAlpha',1);
hold on;
edges3 = round(mediumEdge*10)/10:0.1:max(iti);
histogram(iti(iti>mediumEdge),edges3,'EdgeColor','none','FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',1);
set(gca,'box','off');
print(gcf,'-dpng','iti_byLength');   %png format
saveas(gcf, 'iti_byLength', 'fig');
saveas(gcf, 'iti_byLength','svg');

%% fit with logistic regression

%average of sessions
plot_logreg(lregRCUC_array,tlabel);
print(gcf,'-dpng','logregRCUC_cut');    %png format
saveas(gcf, 'logregRCUC_cut', 'fig');
saveas(gcf, 'logregRCUC_cut','svg');

plot_logreg(lregCRInt_array,tlabel);
print(gcf,'-dpng','logregCRInt_cut');    %png format
saveas(gcf, 'logregCRInt_cut', 'fig');
saveas(gcf,'logregCRInt_cut','svg');
%% concatenate sessions for each animal, fit with logistic regression 
% this part still needs revise: fit the logistic regression for different
% animals separately, then average them
num_regressor=15;

%regressors = rewarded choice, unrewarded choice
lregRCUC_output = cell(0);
for ii = 1:length(animalList)
    [lregRCUC_output{ii}, ~, ~, ~]=logreg_RCUC(stats_all{ii},1,num_regressor);
end
plot_logreg(lregRCUC_output{1},tlabel);
print(gcf,'-dpng','logregRCUC_concat_cut');    %png format
saveas(gcf, 'logregRCUC_concat_cut', 'fig');
saveas(gcf, 'logregRCUC_concat_cut','svg');

%regressors = choice, choice x reward (equivalent to computer's choice)
lregCRInt_output = cell(0);
for ii = 1:length(animalList)
    [lregCRInt_output{ii}, ~, ~, ~]=logreg_CRInt(stats_all{ii},1,num_regressor);
end
plot_logreg(lregCRInt_output{1},tlabel);
print(gcf,'-dpng','logregCRInt_concat_cut');    %png format
saveas(gcf, 'logregCRInt_concat_cut', 'fig');
saveas(gcf, 'logregCRInt_concat_cut','svg');
%%
%disp(['Total number of trials: ' int2str(sum(~isnan(stats_all.c(:,1))))]);
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
%     fitpar = cell(0);
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

%[~, ~, bic.logregCRInt2, nlike.logregCRInt2]=logreg_CRInt(stats_all,1,2);
%[~, ~, bic.logregCRInt5, nlike.logregCRInt5]=logreg_CRInt(stats_all,1,5);
%[~, ~, bic.logregCRInt10, nlike.logregCRInt10]=logreg_CRInt(stats_all,1,10);


%% simulate to get latent action value, compare with animal's choice

player1.label='algo_FQ_RPE_CK';   
player1.params.a=Fitpar{1}.FQ_RPE_CK(1);
player1.params.b=Fitpar{1}.FQ_RPE_CK(2);
player1.params.ac = Fitpar{1}.FQ_RPE_CK(3);
player1.params.bc = Fitpar{1}.FQ_RPE_CK(4);
stats_sim=predictAgent(player1,stats_all{1});

x = 1;  %player 1
n_plot = 500;   %plot first 500 trials
% Figure 2. C
plot_session_qparam(stats_sim,x,n_plot);
for jj = 1:length(animalList)
    player1.params.a=Fitpar{jj}.FQ_RPE_CK(1);
    player1.params.b=Fitpar{jj}.FQ_RPE_CK(2);
    player1.params.ac = Fitpar{jj}.FQ_RPE_CK(3);
    player1.params.bc = Fitpar{jj}.FQ_RPE_CK(4);
    stats_simAll{jj} = predictAgent(player1,stats_all{jj});
end
% Figure 3. E
plot_session_qhist(stats_simAll,x)


%%
disp('---Mean normalized likelihood for fits per session');
fname = fieldnames(nlike_array);
for j=1:numel(fname)  %append for each field
    disp([fname{j} ' - ' num2str(nanmean(nlike_array.(fname{j})))]);
end

disp('---Mean BIC for fits per session');
fname = fieldnames(bic_array);
for j=1:numel(fname)  %append for each field
    disp([fname{j} ' - ' num2str(nanmean(bic_array.(fname{j})))]);
end

%%
% Figure 1. F
figure;
subplot(2,3,1); hold on;
%plot(rand(1,numel(nTrial_array)),nTrial_array,'k^','MarkerSize',15);
%boxplot(nTrial_array,'Colors','k','Notch','off','Labels',[char(956),'=',num2str(round(mean(nTrial_array)))]);
boxplot(nTrial_array,'Colors','k','Notch','off');
%boxplot(nTrial_array,'Colors','k','Notch','off');
ylim([0 800]); 
ylabel(['Trials performed']);
set(gca,'box','off') 

subplot(2,3,2); hold on;
%plot(rand(1,numel(entro_array)),entro_array,'k^','MarkerSize',15);
%boxplot(entro_array,'Colors','k','Symbol','k+','Notch','off','Labels',{[char(956),'=',num2str((mean(entro_array)))]});
boxplot(entro_array,'Colors','k','Symbol','k+','Notch','off');

plot([-1 2],[3 3],'k--','LineWidth',2);
ylim([2 3.1]); 
ylabel('Entropy (bits)');
set(gca,'box','off') 

subplot(2,3,3); hold on;
%plot(rand(1,numel(rrate_array)),100*rrate_array,'k^','MarkerSize',15);
%boxplot(rrate_array*100,'Colors','k','Symbol','k+','Notch','off','Labels',{[char(956),'=',num2str((mean(rrate_array))*100),'%']});
boxplot(rrate_array*100,'Colors','k','Symbol','k+','Notch','off');
plot([-1 2],[50 50],'k--','LineWidth',2);
ylim([20 65]);
ylabel('Reward rate (%)');
set(gca,'box','off') 
print(gcf,'-dpng',['summary' int2str(x)]);    %png format
saveas(gcf,['summary' int2str(x)], 'fig');
saveas(gcf, ['summary' int2str(x)],'svg');

figure;
subplot(1,3,1); hold on;
plot([-0.1 1.1],[0 0],'k--','LineWidth',2);
plot([0 0],[-0.4 4],'k--','LineWidth',2);
plot(fitpar_array.FQ_RPE(:,1),fitpar_array.FQ_RPE(:,2),'k.','MarkerSize',30);
xlabel('\alpha, Learning rate');
ylabel('\beta, Inverse temperature');
xlim([-0.1 1.1]);
ylim([-0.4 4]);

subplot(1,3,2); hold on;
plot([-0.1 1.1],[0 0],'k--','LineWidth',2);
plot([0 0],[-0.4 4],'k--','LineWidth',2);
plot(fitpar_array.FQ_RPE_CK(:,1),fitpar_array.FQ_RPE_CK(:,2),'k.','MarkerSize',30);
xlabel('\alpha, Learning rate (CK)');
ylabel('\beta, Inverse temperature (CK)');
xlim([-0.1 1.1]);
ylim([-0.4 4]);

subplot(1,3,3); hold on;
plot([-0.1 1.1],[0 0],'k--','LineWidth',2);
plot([0 0],[-0.4 4],'k--','LineWidth',2);
plot(fitpar_array.FQ_RPE_CK(:,3),fitpar_array.FQ_RPE_CK(:,4),'k.','MarkerSize',30);
xlabel('\alpha_c, Learning rate (CK)');
ylabel('\beta_c, Inverse temperature (CK)');
xlim([-0.1 1.1]);
ylim([-0.4 4]);
print(gcf,'-dpng',['alpha-beta_cut' int2str(x)]);    %png format
saveas(gcf,['alpha-beta_cut' int2str(x)], 'fig');
saveas(gcf, ['alpha-beta_cut' int2str(x)],'svg');
% 
% alpha = fitpar_array.Q_RPE(:,1); beta = fitpar_array.Q_RPE(:,2);
% save('a_b_saline.mat','alpha','beta');

%% plot the alphas and betas and relative value of beta/beta_k

% Figure 3. D
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
ylim([0 inf]);
%set(gca,'YTick',[0:2:inf]);
set(h,{'linew'},{4})
set(gca,'box','off');
set(gca,'linewidth',4);
print(gcf,'-dpng',['betaK-beta' int2str(x)]);    %png format
saveas(gcf,['betaK-beta' int2str(x)], 'fig');
saveas(gcf, ['betaK-beta' int2str(x)],'svg');

alpha = zeros(1,length(animalList));
alpha_K = zeros(1,length(animalList));
for ii = 1:length(animalList)
    alpha(ii) = Fitpar{ii}.FQ_RPE_CK(1);
    alpha_K(ii) = Fitpar{ii}.FQ_RPE_CK(3);
end
rel_alpha = alpha_K./alpha;

% alpha/alpha_K
figure;
h = boxplot(rel_alpha,'Colors','k','Symbol','k+','Notch','off');
ylabel('\alpha_K/\alpha','FontSize',50)
%ax=gca;
%ax.YAxisLocation = 'origin';
ylim([0 1.6]);
%set(gca,'YTick',[0:0.2:1.2]);
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
%set(gca,'YTick',[0:1:4]);
set(h,{'linew'},{4})
set(gca,'box','off');
set(gca,'linewidth',4);
print(gcf,'-dpng',['beta_K+beta' int2str(x)]);    %png format
saveas(gcf,['beta_K+beta' int2str(x)], 'fig');
saveas(gcf, ['beta_K+beta' int2str(x)],'svg');

% beta_K / (beta_K + beta)
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

%% response time - total/different value correlation 
% no deterministic results

% plot_rt_latentV(stats_sim, trueRespTime);

close all;
