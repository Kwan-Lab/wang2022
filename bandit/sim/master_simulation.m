% Master file for simulating different model options for two-choice tasks
% This model simulation uses BADS optimizer and requires;
%           Optimization Toolbox
%           Global Optimization Toolbox
%           Installed BADS 
%               to install: (bandit-master\simulation\bads-master\install)
%               to test: bads('test')    
%   
% H Atilgan, H Wang and AC Kwan, 03/11/20

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Simulation - task design
n_stim=100000;    % number of trials to simulate

taskparams.p_pairs=[0.7 0.1; 0.1 0.7];
taskparams.rule_labels={'0.7:0.1','0.1:0.7'};
taskparams.crit_hit=10;       % switching criterion: number of times picked the high prob. side
taskparams.crit_geo=1/11;     % switching criterion: after hit, random number of trials defined by geometric dist. p=1/mean 

%% Different optimizer comparison on simulated data figure 
player_sim.model_type='FQ_RPE';
player_sim.params=[0.3 5];
stats_sim = simBandit(player_sim,taskparams,n_stim);
model.fun = @funFQ_RPE;

filename = ['figs/optimviz-init_0.3_5.tiff'];
optimviz_bandit(model.fun,[0.3 5],0,100,stats_sim)
print(gcf,'-dpng',filename);
saveas(gcf, [filename, '.fig']);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Model free learning: Q learning model stimulation 

player_sim.model_type='algo_FQ_RPE';
% avaliable options: 'Q_RPE';

% using median value from fit with animal data
model_path = fullfile(root_path,'mat_models');
load(fullfile(model_path, [player_sim.model_type '.mat']));
fitparMat = cell2mat(fitpar');
player_sim.params=nanmedian(fitparMat,1) 
%player_sim.params=[0.3 5];

tlabel=[player_sim.model_type ', n=' num2str(n_stim) ' trials'];
save_path = fullfile(root_path, ['figs_simulation ' player_sim.model_type]);

stats_sim=simBandit(player_sim,taskparams,n_stim);

analyze_session(stats_sim,tlabel,save_path);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Model based learning: Bayesian model 

% --- agent has no knowledge of blocks
%           prior to criterion: p_trans = params(1); expect some baseline probability of switch
%           after meeting criterion: p_crit = params(1); expect some baseline probability of switch

player_sim.model_type='bayesian_baseline'; 
player_sim.params=[0.03]; %transition probability, same for all trials

% --- agent has knowledge of blocks, 
%           prior to criterion: p_trans = 0; expect no switch prior to criterion,
%           after meeting criterion: p_crit = params(1); expect some probability of switch

player_sim.model_type='bayesian_block'; 
player_sim.params=[1/11]; %p_crit transition probability after criterion is met

% --- agent has knowledge of blocks, 
%           using observable: running # of rewards from same option to predict upcoming switch
%           could have some baseline expectation of probability switch

player_sim.model_type='bayesian_twoparams';
player_sim.params=[0.04 1/5]; %[p_trans p_crit] transition probability before and after criterion is met

% Stimulate accordingly: 
tlabel=[player_sim.model_type ', n=' num2str(n_stim) ' trials'];
save_path = fullfile(root_path, ['test_simulation ' player_sim.model_type]);

stats_sim=simBandit(player_sim,taskparams,n_stim);

% Plots for stimulation
analyze_session(stats_sim,tlabel,save_path);

bandit_fitOneSession(stats_sim,'bayesian_twoparams',save_path);

