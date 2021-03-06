% Master file for analyzing pupillometry recordings of two-armed bandit task
% modified base on H Atilgan and AC Kwan, 191205
% H Wang, 2020
% To run this code:
% 1) Add all the subfolders to the Path in MATLAB
% 2) Change the variable 'root_path' below
%

clearvars;
close all;
setup_figprop;

root_path = 'J:\MP_bandit_pupil_data\bandit\pupilData';
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Preprocessing: parse each log file and tabulate the data set

disp('-----------------------------------------------------------');
disp('--- HERE WE GO ');
disp('-----------------------------------------------------------');

% Look for data files and create a database index
logfilepath = fullfile(root_path,'data');
analysispath = fullfile(root_path,'analysis');
dataIndex = makeDataIndex(logfilepath, analysispath);

% Parse and analyze each logfile, save as .mat files, as needed
dataIndex = createBehMatFiles(dataIndex);


% Add information about pupil
dataIndex = bandit_addIndexPupil(dataIndex);

dataIndex = sortdataIndex(dataIndex);

% Determine if each session fulfill performance criteria
[dataIndex, ~] = determineBehCriteria(dataIndex); 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Behavior - Example session

exampleLogName = '2019.10.13.SUB3.HO_phase3_R71NoCue_1910131512.log';

idx = find(ismember(dataIndex.LogFileName,exampleLogName)==1);
nFiles = size(dataIndex,1);
% Figure 3A

savematpath = [fullfile(dataIndex.BehPath{idx},dataIndex.LogFileName{idx}(1:end-4))];
bandit_session(dataIndex.BehPath{idx},dataIndex.LogFileName{idx},savematpath);

%bandit_session(dataIndex.BehPath{idx},dataIndex.LogFileName{idx});

%% In some sessions, the pupil recordings cannot be aligned with behavior logfiles
% use only matched trials
% for single-pupil recordings only
trigIndex = table(...
    NaN(length(dataIndex.Animal),1));

trigIndex.Properties.VariableNames = {...
    'match'};
mismatch = [46,48,53,56,57,60,89,93,95];
for tt = 1:length(dataIndex.Animal)
    if ismember(tt,mismatch)
        trigIndex.match(tt) = 0;
    else
        trigIndex.match(tt) = 1;
    end
end
dataIndex = [dataIndex,trigIndex];

%%
save_path = fullfile(root_path,'figs_summary');
%normalSubset = isnan(dataIndex.Lesioned);  %control or pre-lesion
%reversalSubset = (ismember(dataIndex.Phase,[3 8])==1);  %reversal

bandit_behaviorPerAnimal(dataIndex(dataIndex.pupil==1 & dataIndex.match==1,:),save_path);
%bandit_behaviorPerAnimal(dataIndex(dataIndex.pupil==1,:),save_path);
%bandit_behaviorPerSession(dataIndex(normalSubset & reversalSubset,:),save_path);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Behavior - Model fitting

model_path = fullfile(root_path,'mat_models');


bandit_fittingPerAnimal(dataIndex,model_path);

% save the predicted latent variable into each behavior mat files
bandit_saveLatent(dataIndex, model_path);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Behavior - Model comparison

model_path = fullfile(root_path,'mat_models');
%bandit_compareModels(model_path);   % model comparison
MP_compareModels(model_path);   % model comparison


print(gcf,'-dpng',fullfile(model_path,'modelcomp'));
saveas(gcf, fullfile(model_path,'modelcomp'), 'fig');
saveas(gcf, fullfile(model_path,'modelcomp'), 'svg');

%% pupil stuff
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  pupil analysis

% get the raw csv data, preprocess it

bandit_createPupilFiles(dataIndex(dataIndex.pupil==1,:));
% add a index for sessions that matlab and logfile did not match


%% pupil
%% save path
model_path = fullfile(root_path,'mat_models');
save_path_pupil = fullfile(root_path,'summary','figs_summary_pupil');

%% pupil simple plots
bandit_pupilSimpleplots(dataIndex(dataIndex.pupil==1 & dataIndex.match==1,:));

%% tonic activity
% no deterministic results
% bandit_tonic(dataIndex(dataIndex.pupil==1 & dataIndex.match==1,:), save_path_pupil);

%% calculate pupil change/pupil response
bandit_pupilChange(dataIndex(dataIndex.pupil==1& dataIndex.match==1,:), save_path_pupil);


%% running regression (choice and reward) on individual sessions
% Figure 7A
bandit_pupilMLR(dataIndex(dataIndex.pupil==1& dataIndex.match==1,:));
bandit_pupilMLR_change_acrossSessions(dataIndex(dataIndex.pupil==1& dataIndex.match==1,:), save_path_pupil);


%% linear regression with latent variables
% Figure 8B
bandit_pupilRL_MLR_CK(dataIndex(dataIndex.pupil==1& dataIndex.match==1,:));
bandit_pupilRL_change_acrossSessions_CK(dataIndex(dataIndex.pupil==1& dataIndex.match==1,:), save_path_pupil);

%% linear regression to get the RPE or Q
bandit_pupilRL_RPE_CK(dataIndex(dataIndex.pupil==1& dataIndex.match==1,:));
%bandit_pupilRLRPE_acrossSessions_CK(dataIndex(dataIndex.pupil==1& dataIndex.match==1,:), save_path_pupil);
bandit_pupilRLRPE_change_acrossSessions_CK(dataIndex(dataIndex.pupil==1& dataIndex.match==1,:), save_path_pupil);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Simulation - using algorithm and actual choices and outcomes from session

model_path = fullfile(root_path,'mat_models/');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Simulation - using algorithm and task design
saveSimPath = fullfile(root_path,'summary','sim');
if ~exist(saveSimPath)
    mkdir(saveSimPath)
end

n_stim=100000;    % number of trials to simulate

taskparams.p_pairs=[0.7 0.1; 0.1 0.7];
taskparams.rule_labels={'0.7:0.1','0.1:0.7'};
taskparams.crit_hit=10;       % switching criterion: number of times picked the high prob. side
taskparams.crit_geo=1/11;     % switching criterion: after hit, random number of trials defined by geometric dist. p=1/mean

player_sim.model_type='FQ_RPE_CK';

% using median value from fit with animal data
model_path = fullfile(root_path,'mat_models');
load(fullfile(model_path, [player_sim.model_type '_model.mat']));
fitparMat = cell2mat(fitpar');
params = nanmedian(fitparMat,1);
player.params=params;
player_sim.params.a=player.params(1);    % learning rate (also = 1 minus the forgetting rate)
player_sim.params.b = player.params(2);
player_sim.params.ac = player.params(3);
player_sim.params.bc = player.params(4);
player_sim.label = 'algo_FQ_RPE_CK';

tlabel=[player_sim.model_type ', n=' num2str(n_stim) ' trials'];
save_path = fullfile(root_path, ['figs_simulation ' player_sim.model_type]);

% simulate for animal's "true" parameters
stats_sim=simBandit(player_sim,taskparams,n_stim);

pStayFit = sum(stats_sim.c(1:end-1,1) == stats_sim.c(2:end,1))/n_stim;
rewardFit = sum(stats_sim.r(:,1))/n_stim;
stats_sim.rule_labels = taskparams.rule_labels;
%analyze_session(stats_sim,tlabel,save_path);

%% alter the beta_K/beta_sum
bkOverSum = 0:0.05:1;
pStayList = zeros(1, length(bkOverSum));
rewardList = zeros(1, length(bkOverSum));
entropyList = zeros(1,length(bkOverSum));
bSum = player.params(2) + player.params(4);

for ii = 1:length(bkOverSum)
    player_sim.label='algo_FQ_RPE_CK';
    player_sim.params.a=player.params(1);    % learning rate (also = 1 minus the forgetting rate)
    player_sim.params.b = bSum*(1-bkOverSum(ii));
    player_sim.params.ac = player.params(3);
    player_sim.params.bc = bkOverSum(ii)*bSum;
    %player1.bias = 1;  % if it is a biased payoff matrix
    
    stats_sim=simBandit(player_sim,taskparams,n_stim);
    
    % calculate probability of stay, reward rate and entropy
    
    pStay = sum(stats_sim.c(1:end-1,1) == stats_sim.c(2:end,1))/n_stim;
    
    pStayList(ii) = pStay;
    rewardList(ii) = sum(stats_sim.r(:,1))/n_stim;
    
end

 %% plot the results
 
fig = figure;
left_color = [0 0 0];
right_color = [1 0 1];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);


plot(bkOverSum,pStayList,'k');
hold on;
plot(bkOverSum,rewardList,'k--');
hold on;
scatter(player.params(4)/(player.params(2)+player.params(4)), pStayFit,120,'k','filled','HandleVisibility','off');
hold on;
scatter(player.params(4)/(player.params(2)+player.params(4)), rewardFit,120,'k','HandleVisibility','off');
ylim([0 1]);

xlabel('\beta_K/(\beta+\beta_K)');
lgd = legend('pStay','Reward');
set(lgd,'color','none','box','off');

title(['\alpha = ',num2str(player.params(1)),' Sum(\beta) = ',num2str(player.params(2)+player.params(4)),' \alpha_K = ',num2str(player.params(3))]);
set(gca,'box','off');

savepath = fullfile(saveSimPath, 'alter_beta');
print(gcf,'-dpng',savepath );   %png format
saveas(gcf, savepath , 'fig');
saveas(gcf, savepath , 'svg');

%% drift parameters
% bandit_fittingDrift(dataIndex,model_path);
% MP_plotDrift(dataIndex,model_path);

%% for two-pupil analysis

clearvars;
close all;
setup_figprop;

root_path = 'J:\MP_bandit_pupil_data\bandit\twopupildata';

disp('-----------------------------------------------------------');
disp('--- HERE WE GO ');
disp('-----------------------------------------------------------');

% Look for data files and create a database index
logfilepath = fullfile(root_path,'data');
analysispath = fullfile(root_path,'analysis');
dataIndex = makeDataIndex(logfilepath, analysispath);

% Parse and analyze each logfile, save as .mat files, as needed
dataIndex = createBehMatFiles(dataIndex);

% for twp-pupil recordings
dataIndex = bandit_addIndextwoPupil(dataIndex);

dataIndex = sortdataIndex(dataIndex);

% Determine if each session fulfill performance criteria
[dataIndex, ~] = determineBehCriteria(dataIndex); 

% for two pupil
save_path = fullfile(root_path,'figs_summary');
bandit_behaviorPerAnimal(dataIndex(dataIndex.pupilSide==3,:),save_path);

%% pupil
model_path = fullfile(root_path,'mat_models');
save_path_pupil = fullfile(root_path,'summary','figs_summary_pupil');

bandit_createtwoPupilFiles(dataIndex(dataIndex.pupil==1,:));
bandit_pupilChange_twopupil(dataIndex(dataIndex.pupil==1,:), save_path_pupil);

% for two-pupil analysis, Figure 8B, C
bandit_pupilMLRleft(dataIndex(dataIndex.pupilSide==3,:));
bandit_pupilMLRright(dataIndex(dataIndex.pupilSide==3,:));
bandit_pupilMLR_change_acrossSessions_twopupil(dataIndex(dataIndex.pupilSide==3,:), save_path_pupil);
