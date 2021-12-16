function MP_pupilRL_acrossSessions_delta(dataIndex, savefigpath)

%% reference to Sul et al.2011
% averaged within subject
setup_figprop
nFiles = size(dataIndex,1);

% load the first file to get some parameters
all_coeff_high = [];
all_pval_high = [];
all_coeff_low = [];
all_pval_low = [];

for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},[dataIndex.LogFileName{ii}(1:end-4),'_beh_cut.mat']));
    savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
    saveRegName = fullfile(savematpath,[fn_beh.name(1:end-7),'regRL_change_CK_value.mat']);
    
    if exist(saveRegName)
        load(saveRegName);
        %if reg_cr_change_high.mse < 
            reg_all_high.regr_time = reg_cr_change_high.regr_time;
            reg_all_high.numPredictor = reg_cr_change_high.numPredictor;
            reg_all_high.interaction = reg_cr_change_high.interaction;

            all_coeff_high = cat(3,all_coeff_high, reg_cr_change_high.coeff);
            all_pval_high = cat(3,all_pval_high, reg_cr_change_high.pval);
        
            reg_all_low.regr_time = reg_cr_change_low.regr_time;
            reg_all_low.numPredictor = reg_cr_change_low.numPredictor;
            reg_all_low.interaction = reg_cr_change_low.interaction;

            all_coeff_low = cat(3,all_coeff_low, reg_cr_change_low.coeff);
            all_pval_low = cat(3,all_pval_low, reg_cr_change_low.pval);
        end
    
            if ~ exist('all_coeff_ctrl_high') && exist('reg_cr_change_high_ctrl')
        % initialize the control matrix
        fieldsName = fieldnames(reg_cr_change_high_ctrl);
        for tt = 1:length(fieldsName)
            all_coeff_ctrl_high.(fieldsName{tt}) = [];
            all_pval_ctrl_high.(fieldsName{tt}) = [];
        end
        
    end
    % load the regression results
    if exist('reg_cr_change_high_ctrl')
        fieldsName = fieldnames(reg_cr_change_high_ctrl);
        for uu = 1:length(fieldsName)
            all_coeff_ctrl_high.(fieldsName{uu}) = cat(3, all_coeff_ctrl_high.(fieldsName{uu}), reg_cr_change_high_ctrl.(fieldsName{uu}).coeff);
            all_pval_ctrl_high.(fieldsName{uu}) = cat(3,all_pval_ctrl_high.(fieldsName{uu}), reg_cr_change_high_ctrl.(fieldsName{uu}).pval);
        end
    end
    
    if ~ exist('all_coeff_ctrl_low') && exist('reg_cr_change_low_ctrl')
        % initialize the control matrix
        fieldsName = fieldnames(reg_cr_change_low_ctrl);
        for tt = 1:length(fieldsName)
            all_coeff_ctrl_low.(fieldsName{tt}) = [];
            all_pval_ctrl_low.(fieldsName{tt}) = [];
        end
        
    end
    % load the regression results
    if exist('reg_cr_change_low_ctrl')
        fieldsName = fieldnames(reg_cr_change_low_ctrl);
        for uu = 1:length(fieldsName)
            all_coeff_ctrl_low.(fieldsName{uu}) = cat(3, all_coeff_ctrl_low.(fieldsName{uu}), reg_cr_change_low_ctrl.(fieldsName{uu}).coeff);
            all_pval_ctrl_low.(fieldsName{uu}) = cat(3, all_pval_ctrl_low.(fieldsName{uu}), reg_cr_change_low_ctrl.(fieldsName{uu}).pval);
        end
    end
end

% get the average
reg_all_high.coeff= all_coeff_high;
reg_all_low.coeff = all_coeff_low;
% using bootstrap to get p value
reg_all_high = getBootstrp(reg_all_low, 0);
reg_all_low = getBootstrp(reg_all_low,0);

pvalThresh = 0.01;
tlabel={'C(n)','C(n-1)','C(n-2)','R(n)','R(n-1)','R(n-2)','C(n)xR(n)','C(n-1)xR(n-1)','C(n-2)xR(n-2)'};

if ~exist(savefigpath)
    mkdir(savefigpath)
end

cd(savefigpath);
reg_all_high.numPredictor = 9;
xtitle = {'Time from cue (s)'};
reg_all_high.nback = 0;
MP_plot_regrcoef_pupil(reg_all_high,pvalThresh,tlabel,xtitle);
print(gcf,'-dpng','MLR-RL-coeff-chosenQ');    %png format
saveas(gcf, 'MLR-RL-coeff-chosenQ', 'fig');



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
if exist('all_coeff_ctrl')
    reg_pval_ctrl = getBootstrp(all_pval_ctrl, 0.01);

% MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
    MP_plot_regr(reg_sig,reg_pval_ctrl, reg_sig.pvalThresh,tlabel,xtitle);
else
    MP_plot_regr(reg_sig,[], reg_sig.pvalThresh,tlabel,xtitle);
end

print(gcf,'-dpng','MLR-RL-chosenQ');    %png format
saveas(gcf, 'MLR-RL-chosenQ', 'fig');



close all

end