function bandit_pupilMLR_change_acrossSessions(dataIndex, savefigpath)

%% reference to Sul et al.2011
set(0, 'DefaultFigureRenderer', 'painters');
% averaged within subject
nFiles = size(dataIndex,1);

subject_mask = [];
animalList = unique(dataIndex.Animal);
% load the first file to get some parameters
all_coeff_future = [];
all_pval_future = [];

for ii = 1:nFiles
    savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},['bandit',dataIndex.LogFileName{ii}(end-30:end-4),'_beh.mat']));
    
    saveRegName = fullfile(savematpath,[fn_beh.name(1:end-7),'regCR_change.mat']);

    if exist(saveRegName)
        load(saveRegName)
        
        % get subject mask
        subject_mask(end+1) = find(strcmp(animalList, dataIndex.Animal{ii}));

    % load the MLR with C(n+1)

        all_coeff_future = cat(3,all_coeff_future, reg_cr_future_change.coeff);
        all_pval_future = cat(3, all_pval_future, reg_cr_future_change.pval);
     
    end
end
        



%% original linear regression
% 1. use bootstrap to get the average and 95% CI for each factor, plot the
% bar plot


% go to the save path
if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);

%%  linear regression with C(n+1)
reg_cr_future_change.coeff= all_coeff_future;

% use bootstrp to get coefficient
reg_cr_future_change = getBootstrp(reg_cr_future_change, 0,0.05);

if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);
pvalThresh = 0.01;
xtitle='Time from cue (s)';
tlabel={'C(n+1)','C(n)','C(n-1)','C(n-2)','R(n+1)','R(n)', 'R(n-1)','R(n-2)',...
        'C(n+1)*R(n+1)','C(n)*R(n)','C(n-1)*R(n-1)','C(n-2)*R(n-2)','Reward Rate','Cumulative Reward'};
   
MP_plot_regrcoef_pupil(reg_cr_future_change,pvalThresh,tlabel,xtitle);
print(gcf,'-dpng','MLR-change-choiceoutcome-future-averageSession');    %png format
saveas(gcf, 'MLR-change-choiceoutcome-future-averageSession', 'fig');
saveas(gcf, 'MLR-change-choiceoutcome-future-averageSession', 'svg');

% plot the figure as number of session that is significant
reg_sig.coeff = all_coeff_future;
reg_sig.pval = all_pval_future;
reg_sig.regr_time = reg_cr_future_change.regr_time;
reg_sig.numPredictor = reg_cr_future_change.numPredictor;
reg_sig.nback = reg_cr_future_change.nback;
reg_sig.interaction = reg_cr_future_change.interaction;
reg_sig.pvalThresh= 0.01;

% preprocessing the control coefficient and pval: bootstrap the pval to
% determine a baseline percentage of significant session
%reg_pval_futurectrl = getBootstrp(all_pval_futurectrl, 0.01);

if exist('all_pval_future_ctrl')
    reg_pval_future_ctrl = getBootstrp(all_pval_future_ctrl, 0.01,0.01);

% MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
    MP_plot_regr(reg_sig,reg_pval_future_ctrl, reg_sig.pvalThresh,tlabel,xtitle);
else
    MP_plot_regr(reg_sig,[], reg_sig.pvalThresh,tlabel,xtitle);
end
print(gcf,'-dpng','MLR-change-choiceoutcome-future-sigSession');    %png format
saveas(gcf, 'MLR-change-choiceoutcome-future-sigSession', 'fig');
saveas(gcf, 'MLR-change-choiceoutcome-future-sigSession', 'svg');

close all

end