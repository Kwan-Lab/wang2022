function MP_pupilMLR_change_acrossSessions(dataIndex, savefigpath)

%% reference to Sul et al.2011
% read in the individual sessions of linear regression results, average
% them to get average coefficient and significant fraction of sessions

setup_figprop
set(0, 'DefaultFigureRenderer', 'painters');
nFiles = size(dataIndex,1);

subject_mask = [];
animalList = unique(dataIndex.Animal);

% load the first file to get some parameters
all_coeff_future = [];
all_pval_future = [];


for ii = 1:nFiles
    savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},[dataIndex.LogFileName{ii}(1:end-4),'_beh_cut.mat']));
    
    saveRegName = fullfile(savematpath,[fn_beh.name(1:end-7),'regCR_change.mat']);
   %saveRegName_ITI = fullfile(savematpath,[fn_beh.name(1:end-7),'regCR_ITI.mat']);

    if exist(saveRegName)
        load(saveRegName)
        load(saveRegName_ITI)
        % get subject mask
        subject_mask(end+1) = find(strcmp(animalList, dataIndex.Animal{ii}));
        
        % load choice and reward regression
%          reg_all.regr_time = reg_cr_change.regr_time;
%         reg_all.numPredictor = reg_cr_change.numPredictor;
%         reg_all.nback = reg_cr_change.nback;
%         reg_all.interaction = reg_cr_change.interaction;
%           all_coeff = cat(3,all_coeff, reg_cr_change.coeff);
%         all_pval = cat(3,all_pval, reg_cr_change.pval);
        
    
    % load the MLR with C(n+1)

        all_coeff_future = cat(3,all_coeff_future, reg_cr_future_change.coeff);
        all_pval_future = cat(3, all_pval_future, reg_cr_future_change.pval);
        
        reg_all.regr_time = reg_cr_future_change.regr_time;
        reg_all.numPredictor = reg_cr_future_change.numPredictor;
        reg_all.nback = reg_cr_future_change.nback;
        reg_all.interaction = reg_cr_future_change.interaction;
     % load the ITI regression (n+1 and n)
        
       
%         all_coeff_iti1 = cat(3,all_coeff_iti1, reg_cr1_change.coeff);
%         all_pval_iti1 = cat(3,all_pval_iti1, reg_cr1_change.pval);
%         
%         all_coeff_iti2 = cat(3,all_coeff_iti2, reg_cr2_change.coeff);
%         all_pval_iti2 = cat(3,all_pval_iti2, reg_cr2_change.pval);
%        
%         all_coeff_iti3 = cat(3,all_coeff_iti3, reg_cr3_change.coeff);
%         all_pval_iti3 = cat(3,all_pval_iti3, reg_cr3_change.pval);
%         
         % load the ITI regression (n-1 and n)
        
        
        
    end
end
        

%% other things can be done for pupil response:

% go to the save path
if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);

%%  linear regression with C(n+1)
reg_cr_future_change.coeff= all_coeff_future;

% use bootstrp to get coefficient
reg_cr_future_change = getBootstrp(reg_cr_future_change, 0, 0.05);

if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);

xtitle='Time from cue (s)';
tlabel={'c(n+1)','c(n)','c(n-1)','c(n-2)','r(n+1)','r(n)', 'r(n-1)','r(n-2)',...
        'c(n+1)*r(n+1)','c(n)*r(n)','c(n-1)*r(n-1)','c(n-2)*r(n-2)','Reward Rate','Cumulative Reward'};
pvalThresh=0.01;
MP_plot_regrcoef_pupil(reg_cr_future_change,pvalThresh,tlabel,xtitle);
print(gcf,'-dpng','MLR-change-choiceoutcome_future_averageSession');    %png format
saveas(gcf, 'MLR-change-choiceoutcome_future_averageSession', 'fig');
saveas(gcf, 'MLR-change-choiceoutcome_future_averageSession','svg');

% plot the figure as number of session that is significant
reg_sig.coeff = all_coeff_future;
reg_sig.pval = all_pval_future;
reg_sig.regr_time = reg_cr_future_change.regr_time;
reg_sig.numPredictor = reg_cr_future_change.numPredictor;
reg_sig.nback = reg_cr_future_change.nback;
reg_sig.interaction = reg_cr_future_change.interaction;
reg_sig.pvalThresh= 0.01;


if exist('all_pval_future_ctrl')
    reg_pval_future_ctrl = getBootstrp(all_pval_future_ctrl, 0.01, 0.05);

% MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
    MP_plot_regr(reg_sig,reg_pval_future_ctrl, reg_sig.pvalThresh,tlabel,xtitle);
else
    MP_plot_regr(reg_sig,[], reg_sig.pvalThresh,tlabel,xtitle);
end
print(gcf,'-dpng','MLR-change-choiceoutcome_future_sigSession');    %png format
saveas(gcf, 'MLR-change-choiceoutcome_future_sigSession', 'fig');
saveas(gcf, 'MLR-change-choiceoutcome_future_sigSession','svg');

%%
close all

end