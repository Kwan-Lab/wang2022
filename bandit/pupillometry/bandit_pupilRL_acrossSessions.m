function bandit_pupilRL_acrossSessions(dataIndex, savefigpath);

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
    saveRegName = fullfile(savematpath,[fn_beh.name(1:end-7),'regRL.mat']);
    
    if exist(saveRegName)
        load(saveRegName);
        reg_all.regr_time = reg_cr.regr_time;
        reg_all.numPredictor = reg_cr.numPredictor;
        reg_all.interaction = reg_cr.interaction;

        all_coeff = cat(3,all_coeff, reg_cr.coeff);
        all_pval = cat(3,all_pval, reg_cr.pval);
    end
    
    if ~ exist('all_coeff_ctrl') && exist('reg_cr_ctrl')
        % initialize the control matrix
        fieldsName = fieldnames(reg_cr_ctrl);
        for tt = 1:length(fieldsName)
            all_coeff_ctrl.(fieldsName{tt}) = [];
            all_pval_ctrl.(fieldsName{tt}) = [];
        end
        
    end
    % load the regression results
    if exist('reg_cr_ctrl')
        fieldsName = fieldnames(reg_cr_ctrl);
        for uu = 1:length(fieldsName)
            all_coeff_ctrl.(fieldsName{uu}) = cat(3, all_coeff_ctrl.(fieldsName{uu}), reg_cr_ctrl.(fieldsName{uu}).coeff);
            all_pval_ctrl.(fieldsName{uu}) = cat(3, all_pval_ctrl.(fieldsName{uu}), reg_cr_ctrl.(fieldsName{uu}).pval);
        end
    end

end

% get the average
reg_all.coeff= all_coeff;

% using bootstrap to get p value
reg_all = getBootstrp(reg_all, 0);

pvalThresh = 0.01;
tlabel={'C(n)','R(n)','C(n)xR(n)','C(n-1)','R(n-1)', 'C(n-1)xR(n-1)','QL-QR', 'ChosenQ', 'Reward Rate','Cumulative reward'};
        
if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);
xtitle = {'Time from cue (s)'};

MP_plot_regrcoef_RLpupil(reg_all,pvalThresh,tlabel,xtitle);
print(gcf,'-dpng','MLR-RL_coeff_chosenQ');    %png format
saveas(gcf, 'MLR-RL_coeff_chosenQ', 'fig');



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
% reg_pval_ctrl = getBootstrp(all_pval_ctrl, 0.01);

if exist('all_coeff_ctrl')
    reg_pval_ctrl = getBootstrp(all_pval_ctrl, 0.01);

% MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
    MP_plot_regr(reg_sig,reg_pval_ctrl, reg_sig.pvalThresh,tlabel,xtitle);
else
    MP_plot_regr(reg_sig,[], reg_sig.pvalThresh,tlabel,xtitle);
end

print(gcf,'-dpng','MLR-RL_sigSession');    %png format
saveas(gcf, 'MLR-RL_sigSession', 'fig');

close all

end