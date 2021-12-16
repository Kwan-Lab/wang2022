function bandit_pupilRL_change_acrossSessions_CK(dataIndex, savefigpath);

%% reference to Sul et al.2011
% averaged within subject
nFiles = size(dataIndex,1);

% load the first file to get some parameters
all_coeff = [];
all_pval = [];

for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},['bandit',dataIndex.LogFileName{ii}(end-30:end-4),'_beh.mat']));
    savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
    saveRegName = fullfile(savematpath,[fn_beh.name(1:end-7),'regRL_change_CK.mat']);
    
    if exist(saveRegName)
        load(saveRegName);
        reg_all.regr_time = reg_cr_change.regr_time;
        reg_all.numPredictor = reg_cr_change.numPredictor;
        reg_all.interaction = reg_cr_change.interaction;

        all_coeff = cat(3,all_coeff, reg_cr_change.coeff);
        all_pval = cat(3,all_pval, reg_cr_change.pval);
    end
    
    if ~ exist('all_coeff_ctrl') && exist('reg_cr_change_ctrl')
        % initialize the control matrix
        fieldsName = fieldnames(reg_cr_change_ctrl);
        for tt = 1:length(fieldsName)
            all_coeff_ctrl.(fieldsName{tt}) = [];
            all_pval_ctrl.(fieldsName{tt}) = [];
        end
        
    end
    % load the regression results
    if exist('reg_cr_change_ctrl')
        fieldsName = fieldnames(reg_cr_change_ctrl);
        for uu = 1:length(fieldsName)
            all_coeff_ctrl.(fieldsName{uu}) = cat(3, all_coeff_ctrl.(fieldsName{uu}), reg_cr_change_ctrl.(fieldsName{uu}).coeff);
            all_pval_ctrl.(fieldsName{uu}) = cat(3, all_pval_ctrl.(fieldsName{uu}), reg_cr_change_ctrl.(fieldsName{uu}).pval);
        end
    end

end

%% get the average
reg_all.coeff= all_coeff;

% using bootstrap to get p value
reg_all = getBootstrp(reg_all, 0,0.05);

if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);

pvalThresh = 0.01;
tlabel={'c(n)','r(n)','c(n)xr(n)','c(n-1)','r(n-1)', 'c(n-1)xr(n-1)','Q_L-Q_R', 'ChosenQ','K_L-K_R', 'ChosenK', 'Reward Rate','Cumulative reward'};
  
xtitle = {'Time from cue (s)'};

reg_all.nback = 0;
MP_plot_regrcoef_pupil(reg_all,pvalThresh,tlabel,xtitle);
print(gcf,'-dpng','MLR-RL_change_coeff_chosenQ');    %png format
saveas(gcf, 'MLR-RL_change_coeff_chosenQ', 'fig');
saveas(gcf, 'MLR-RL_change_coeff_chosenQ', 'svg');

% plot the figure as number of session that is significant
reg_sig.coeff = all_coeff;
reg_sig.pval = all_pval;
reg_sig.regr_time = reg_all.regr_time;
reg_sig.numPredictor = reg_all.numPredictor;
reg_sig.interaction = reg_all.interaction;
reg_sig.pvalThresh= 0.01;
reg_sig.nback = 0;
xtitle = {'Time from cue (s)'};


% preprocessing the control coefficient and pval: bootstrap the pval to
% determine a baseline percentage of significant session
%reg_pval_ctrl = getBootstrp(all_pval_ctrl, 0.01);

if exist('all_coeff_ctrl')
    reg_pval_ctrl = getBootstrp(all_pval_ctrl, 0.01,0.01);

% MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
    MP_plot_regr(reg_sig,reg_pval_ctrl, reg_sig.pvalThresh,tlabel,xtitle);
else
    MP_plot_regr(reg_sig,[], reg_sig.pvalThresh,tlabel,xtitle);
end
print(gcf,'-dpng','MLR-RL_change_chosenQ');    %png format
saveas(gcf, 'MLR-RL_change_chosenQ', 'fig');
saveas(gcf, 'MLR-RL_change_chosenQ', 'svg');


close all

end