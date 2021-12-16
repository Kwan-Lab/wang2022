function bandit_pupilReward_acrossSessions(dataIndex, savefigpath)

%% reference to Sul et al.2011
% averaged within subject
nFiles = size(dataIndex,1);

subject_mask = [];
animalList = unique(dataIndex.Animal);
% load the first file to get some parameters

all_coeff_reward = [];
all_pval_reward = [];
all_mse_reward = [];
all_sumsq_reward = [];

all_coeff_change_reward = [];
all_pval_change_reward = [];
all_mse_change_reward = [];
all_sumsq_change_reward = [];


for ii = 1:nFiles
    savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},['bandit',dataIndex.LogFileName{ii}(end-30:end-4),'_beh.mat']));
    
    saveRegName = fullfile(savematpath,[fn_beh.name(1:end-7),'regReward.mat']);
    saveRegName_change = fullfile(savematpath,[fn_beh.name(1:end-7),'regReward_change.mat']);
    
    if exist(saveRegName) & exist(saveRegName_change)
        load(saveRegName)
        load(saveRegName_change)
        
        % get subject mask
        subject_mask(end+1) = find(strcmp(animalList, dataIndex.Animal{ii}));
        
        % load choice and reward regression
        
        reg_all.regr_time = reg_crR1.regr_time;
        reg_all.numPredictor = reg_crR1.numPredictor;
        reg_all.nback = reg_crR1.nback;
        reg_all.interaction = reg_crR1.interaction;
        
     % load the regression with past rewards
     % 1:cumulative reward; 2:cumulative average reward; 3:running ave
     % reward 10; 4:running average reward 20; 5:running average reward 40;
        all_coeff_reward = cat(3, all_coeff_reward, [reg_crR1.coeff(:,2),reg_crR2.coeff(:,2),reg_crR3.coeff(:,2),reg_crR4.coeff(:,2),reg_crR5.coeff(:,2)]);
        all_pval_reward = cat(3, all_pval_reward, [reg_crR1.pval(:,2),reg_crR2.pval(:,2),reg_crR3.pval(:,2),reg_crR4.pval(:,2),reg_crR5.pval(:,2)]);
        all_mse_reward = cat(2, all_mse_reward, [reg_crR1.mse;reg_crR2.mse;reg_crR3.mse;reg_crR4.mse;reg_crR5.mse]);
        % normalize the sum of squared variance
        all_sumsq_reward = cat(3, all_sumsq_reward, [reg_crR1.SumSq(:,2)./sum(reg_crR1.SumSq,2),...
                                                     reg_crR2.SumSq(:,2)./sum(reg_crR2.SumSq,2),...
                                                     reg_crR3.SumSq(:,2)./sum(reg_crR3.SumSq,2),...
                                                     reg_crR4.SumSq(:,2)./sum(reg_crR4.SumSq,2),...
                                                     reg_crR5.SumSq(:,2)./sum(reg_crR5.SumSq,2)]);
        
        % pupil change
        all_coeff_change_reward = cat(3, all_coeff_change_reward, [reg_crR1_change.coeff(:,2),...
                                                     reg_crR2_change.coeff(:,2),...
                                                     reg_crR3_change.coeff(:,2),...
                                                     reg_crR4_change.coeff(:,2),...
                                                     reg_crR5_change.coeff(:,2)]);
        all_pval_change_reward = cat(3, all_pval_change_reward, [reg_crR1_change.pval(:,2),...
                                                   reg_crR2_change.pval(:,2),...
                                                   reg_crR3_change.pval(:,2),...
                                                   reg_crR4_change.pval(:,2),...
                                                   reg_crR5_change.pval(:,2)]);
        all_mse_change_reward = cat(2, all_mse_change_reward, [reg_crR1_change.mse;...
                                                 reg_crR2_change.mse;...
                                                 reg_crR3_change.mse;...
                                                 reg_crR4_change.mse;...
                                                 reg_crR5_change.mse]);
        % fraction of the sum of squared variance
        all_sumsq_change_reward = cat(3, all_sumsq_change_reward, [reg_crR1_change.SumSq(:,2)./sum(reg_crR1_change.SumSq,2),...
                                                     reg_crR2_change.SumSq(:,2)./sum(reg_crR2_change.SumSq,2),...
                                                     reg_crR3_change.SumSq(:,2)./sum(reg_crR3_change.SumSq,2),...
                                                     reg_crR4_change.SumSq(:,2)./sum(reg_crR4_change.SumSq,2),...
                                                     reg_crR5_change.SumSq(:,2)./sum(reg_crR5_change.SumSq,2)]);
    end
end

if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);

%% check the 5 reward regression
reg_reward.coeff= all_coeff_reward;
% fill the first "bias" term
sizeCoeff = size(all_coeff_reward);
bias = NaN(sizeCoeff(1),1, sizeCoeff(3));
reg_reward.coeff = cat(2, bias,reg_reward.coeff);
% use bootstrp to get coefficient
reg_reward = getBootstrp(reg_reward, 0);
reg_reward.regr_time = reg_all.regr_time;
reg_reward.numPredictor = 5;
reg_reward.nback = 0;
reg_reward.interaction = 0;

pvalThresh = NaN;
xtitle = 'Time from cue (s)';
tlabel={'cumu-r','ave-cumu-r','run-r-10','run-r-20','run-r-40'};
% tlabel={'C(n)','R(n)','C(n)xR(n)'};
if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);

MP_plot_regrcoef_pupil(reg_reward,pvalThresh,tlabel,xtitle);

print(gcf,'-dpng','MLR-differentR_averageSession');    %png format
saveas(gcf, 'MLR-MLR-differentR_averageSession', 'fig');

% check the amount of variance explained by different regression factor
sumsq_reward.coeff= all_sumsq_reward;
% fill the first "bias" term
sizeCoeff = size(all_sumsq_reward);
bias = NaN(sizeCoeff(1),1, sizeCoeff(3));
sumsq_reward.coeff = cat(2, bias,sumsq_reward.coeff);
% use bootstrp to get coefficient
sumsq_reward = getBootstrp(sumsq_reward, 0);
sumsq_reward.regr_time = reg_all.regr_time;
sumsq_reward.numPredictor = 5;
sumsq_reward.nback = 0;
sumsq_reward.interaction = 0;

pvalThresh = NaN;
xtitle = 'Frac of var explained';
tlabel={'cumu-r','ave-cumu-r','run-r-10','run-r-20','run-r-40'};
% tlabel={'C(n)','R(n)','C(n)xR(n)'};


MP_plot_regrcoef_pupil(sumsq_reward,pvalThresh,tlabel,xtitle);
print(gcf,'-dpng','MLR-RewardVar_averageSession');    %png format
saveas(gcf, 'MLR-MLR-RewardVar_averageSession', 'fig');

%% check the 5 reward regression for pupil change
reg_reward_change.coeff= all_coeff_change_reward;
% fill the first "bias" term
sizeCoeff_change = size(all_coeff_change_reward);
bias = NaN(sizeCoeff_change(1),1, sizeCoeff_change(3));
reg_reward_change.coeff = cat(2, bias,reg_reward_change.coeff);
% use bootstrp to get coefficient
reg_reward_change = getBootstrp(reg_reward_change, 0);
reg_reward_change.regr_time = reg_all.regr_time;
reg_reward_change.numPredictor = 5;
reg_reward_change.nback = 0;
reg_reward_change.interaction = 0;

pvalThresh = NaN;
xtitle = 'Time from cue (s)';
tlabel={'cumu-r','ave-cumu-r','run-r-10','run-r-20','run-r-40'};
% tlabel={'C(n)','R(n)','C(n)xR(n)'};


MP_plot_regrcoef_pupil(reg_reward_change,pvalThresh,tlabel,xtitle);

print(gcf,'-dpng','MLR-differentR-change_averageSession');    %png format
saveas(gcf, 'MLR-MLR-differentR-change_averageSession', 'fig');

% check the amount of variance explained by different regression factor
sumsq_reward_change.coeff= all_sumsq_change_reward;
% fill the first "bias" term
sizeCoeff_change = size(all_sumsq_change_reward);
bias = NaN(sizeCoeff_change(1),1, sizeCoeff_change(3));
sumsq_reward_change.coeff = cat(2, bias,sumsq_reward_change.coeff);
% use bootstrp to get coefficient
sumsq_reward_change = getBootstrp(sumsq_reward_change, 0);
sumsq_reward_change.regr_time = reg_all.regr_time;
sumsq_reward_change.numPredictor = 5;
sumsq_reward_change.nback = 0;
sumsq_reward_change.interaction = 0;

pvalThresh = NaN;
xtitle = 'Frac of var explained';
tlabel={'cumu-r','ave-cumu-r','run-r-10','run-r-20','run-r-40'};
% tlabel={'C(n)','R(n)','C(n)xR(n)'};
if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);

MP_plot_regrcoef_pupil(sumsq_reward_change,pvalThresh,tlabel,xtitle);
print(gcf,'-dpng','MLR-RewardVar-change_averageSession');    %png format
saveas(gcf, 'MLR-MLR-RewardVar-change_averageSession', 'fig');


end
        
