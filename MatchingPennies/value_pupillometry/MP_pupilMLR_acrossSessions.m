function MP_pupilMLR_acrossSessions(dataIndex, savefigpath)

%% reference to Sul et al.2011
% averaged within subject
nFiles = size(dataIndex,1);

subject_mask = [];
animalList = unique(dataIndex.Animal);
% load the first file to get some parameters
all_coeff = [];
all_pval = [];
    
all_coeff_future = [];
all_pval_future = [];


all_coeff_iti1 = [];
all_pval_iti1 = [];
all_coeff_iti2 = [];
all_pval_iti2 = [];
all_coeff_iti3 = [];
all_pval_iti3 = [];



for ii = 1:nFiles
    savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},[dataIndex.LogFileName{ii}(1:end-4),'_beh_cut.mat']));
    
    saveRegName = fullfile(savematpath,[fn_beh.name(1:end-7),'regCR.mat']);
    saveRegName_future = fullfile(savematpath,[fn_beh.name(1:end-7),'regCR_cut_future.mat']);

    
    if exist(saveRegName)
        load(saveRegName)
        
        % get subject mask
        subject_mask(end+1) = find(strcmp(animalList, dataIndex.Animal{ii}));
        
        % load choice and reward regression
        
        reg_all.regr_time = reg_cr.regr_time;
        reg_all.numPredictor = reg_cr.numPredictor;
        reg_all.nback = reg_cr.nback;
        reg_all.interaction = reg_cr.interaction;
        all_coeff = cat(3,all_coeff, reg_cr.coeff);
        all_pval = cat(3,all_pval, reg_cr.pval);
        
    
    % load the MLR with C(n+1)

        all_coeff_future = cat(3,all_coeff_future, reg_cr_future.coeff);
        all_pval_future = cat(3, all_pval_future, reg_cr_future.pval);
    
     % load the ITI regression
        
       
        all_coeff_iti1 = cat(3,all_coeff_iti1, reg_cr1.coeff);
        all_pval_iti1 = cat(3,all_pval_iti1, reg_cr1.pval);
        
        all_coeff_iti2 = cat(3,all_coeff_iti2, reg_cr2.coeff);
        all_pval_iti2 = cat(3,all_pval_iti2, reg_cr2.pval);
       
        all_coeff_iti3 = cat(3,all_coeff_iti3, reg_cr3.coeff);
        all_pval_iti3 = cat(3,all_pval_iti3, reg_cr3.pval);

     % load the control
        if ~exist('all_coeff_future_ctrl')  && exist('reg_cr_future_ctrl')
            % initialize the control matrix
            fieldsName = fieldnames(reg_cr_future_ctrl);
            for tt = 1:length(fieldsName)
                all_coeff_future_ctrl.(fieldsName{tt}) = [];
                all_pval_future_ctrl.(fieldsName{tt}) = [];
            end
            
        end   
        
        if exist('reg_cr_future_ctrl')
        % load the regression results
            fieldsName = fieldnames(reg_cr_future_ctrl);
            for uu = 1:length(fieldsName)
                all_coeff_future_ctrl.(fieldsName{uu}) = cat(3, all_coeff_future_ctrl.(fieldsName{uu}), reg_cr_future_ctrl.(fieldsName{uu}).coeff);
                all_pval_future_ctrl.(fieldsName{uu}) = cat(3, all_pval_future_ctrl.(fieldsName{uu}), reg_cr_future_ctrl.(fieldsName{uu}).pval);
            end
        end  

        % load control regression for choice and reward regression
        if ~exist('all_coeff_ctrl')  && exist('reg_cr_ctrl')
            % initialize the control matrix
            fieldsName = fieldnames(reg_cr_ctrl);
            for tt = 1:length(fieldsName)
                all_coeff_ctrl.(fieldsName{tt}) = [];
                all_pval_ctrl.(fieldsName{tt}) = [];
            end
        end
        
        if exist('reg_cr_ctrl')
        % load the regression results
            fieldsName = fieldnames(reg_cr_ctrl);
            for uu = 1:length(fieldsName)
                all_coeff_ctrl.(fieldsName{uu}) = cat(3, all_coeff_ctrl.(fieldsName{uu}), reg_cr_ctrl.(fieldsName{uu}).coeff);
                all_pval_ctrl.(fieldsName{uu}) = cat(3, all_pval_ctrl.(fieldsName{uu}), reg_cr_ctrl.(fieldsName{uu}).pval);
            end
        end
    end
end
        



%% original linear regression
% get the average
reg_all.coeff= all_coeff;

% use bootstrp to get coefficient
reg_all = getBootstrp(reg_all, 0);


pvalThresh = NaN;
xtitle = 'Time from cue (s)';
tlabel={'c(n)','c(n-1)','c(n-2)','r(n)','r(n-1)','r(n-2)','c(n)xr(n)','c(n-1)xr(n-1)','c(n-2)xr(n-2)'};

if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);

MP_plot_regrcoef_pupil(reg_all,pvalThresh,tlabel,xtitle);
print(gcf,'-dpng','MLR-choiceoutcome_averageSession');    %png format
saveas(gcf, 'MLR-choiceoutcome_averageSession', 'fig');
saveas(gcf, 'MLR-choiceoutcome_averageSession','svg');
% plot the figure as number of session that is significant
reg_sig.coeff = all_coeff;
reg_sig.pval = all_pval;
reg_sig.regr_time = reg_all.regr_time;
reg_sig.numPredictor = reg_all.numPredictor;
reg_sig.nback = reg_all.nback;
reg_sig.interaction = reg_all.interaction;
reg_sig.pvalThresh= 0.01;

% preprocessing the control coefficient and pval: bootstrap the pval to
% determine a baseline percentage of significant session

% MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
if exist('all_coeff_ctrl')
    reg_pval_ctrl = getBootstrp(all_pval_ctrl, 0.01);

% MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
    MP_plot_regr(reg_sig,reg_pval_ctrl, reg_sig.pvalThresh,tlabel,xtitle);
else
    MP_plot_regr(reg_sig,[], reg_sig.pvalThresh,tlabel,xtitle);
end

print(gcf,'-dpng','MLR-choiceoutcome_sigSession');    %png format
saveas(gcf, 'MLR-choiceoutcome_sigSession', 'fig');
saveas(gcf, 'MLR-choiceoutcome_sigSession','svg');
%% separate the choice coefficient into pos/neg 

tlabel={'C(n)-pos','C(n)-neg','C(n)-|pos-neg|'};
MP_plot_regr_PN(reg_sig,reg_sig.pvalThresh,tlabel,xtitle);
print(gcf,'-dpng','MLR-choiceoutcome_sigSession_abs_pos-neg');    %png format
saveas(gcf, 'MLR-choiceoutcome_sigSession_abs_pos-neg', 'fig');

% choice effect should be more consistent within subject, not so for motor
% effects?

% loop through all animals, plot then separatly
for tt = 1:length(animalList)
    tlabel={'C(n)-pos','C(n)-neg','C(n)-|pos-neg|'};
    xtitle = ['Time from cue (s)',animalList{tt}] ;
    if sum(subject_mask == tt) > 0   % if pupil data for certain animal exists
        reg_sub = reg_sig;
        reg_sub.coeff = reg_sig.coeff(:,:,subject_mask == tt);
        reg_sub.pval = reg_sig.pval(:,:, subject_mask == tt);
        MP_plot_regr_PN(reg_sub,reg_sub.pvalThresh,tlabel,xtitle);
        print(gcf,'-dpng',['MLR-choiceoutcome_sigSession_abs_pos-neg',animalList{tt}]);    %png format
        saveas(gcf, ['MLR-choiceoutcome_sigSession_abs_pos-neg',animalList{tt}], 'fig');
    end
end


%%  linear regression with C(n+1)
reg_cr_future.coeff= all_coeff_future;

% use bootstrp to get coefficient
reg_cr_future = getBootstrp(reg_cr_future, 0);
pvalThresh = NaN;
xtitle = 'Time from cue (s)';
tlabel={'c(n+1)','c(n)','c(n-1)','c(n-2)','r(n+1)','r(n)', 'r(n-1)','r(n-2)',...
                    'c(n+1)*r(n+1)','c(n)*r(n)','c(n-1)*r(n-1)','c(n-2)*r(n-2)','Reward Rate','Cumulative Reward'};

if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);

MP_plot_regrcoef_pupil(reg_cr_future,pvalThresh,tlabel,xtitle);
print(gcf,'-dpng','MLR-choiceoutcome_future_averageSession');    %png format
saveas(gcf, 'MLR-choiceoutcome_future_averageSession', 'fig');
saveas(gcf, 'MLR-choiceoutcome_future_averageSession','svg');
% plot the figure as number of session that is significant
reg_sig.coeff = all_coeff_future;
reg_sig.pval = all_pval_future;
reg_sig.regr_time = reg_cr_future.regr_time;
reg_sig.numPredictor = reg_cr_future.numPredictor;
reg_sig.nback = reg_cr_future.nback;
reg_sig.interaction = reg_cr_future.interaction;
reg_sig.pvalThresh= 0.01;

% preprocessing the control coefficient and pval: bootstrap the pval to
% determine a baseline percentage of significant session
if exist('all_pval_future_ctrl')
    reg_pval_future_ctrl = getBootstrp(all_pval_future_ctrl, 0.01);

% MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
    MP_plot_regr(reg_sig,reg_pval_future_ctrl, reg_sig.pvalThresh,tlabel,xtitle);
else
    MP_plot_regr(reg_sig,[], reg_sig.pvalThresh,tlabel,xtitle);
end

MP_plot_regr(reg_sig,reg_pval_future_ctrl, reg_sig.pvalThresh,tlabel,xtitle);
print(gcf,'-dpng','MLR-choiceoutcome_future_sigSession');    %png format
saveas(gcf, 'MLR-choiceoutcome_future_sigSession', 'fig');
saveas(gcf, 'MLR-choiceoutcome_future_sigSession','svg');
%% plot the coefficient in different ITIs
% plot the figure as number of session that is significant
reg_sig_iti1.coeff = all_coeff_iti1;
reg_sig_iti1.pval = all_pval_iti1;
reg_sig_iti1.regr_time = reg_cr1.regr_time;
reg_sig_iti1.numPredictor = reg_cr1.numPredictor;
reg_sig_iti1.nback = reg_cr1.nback;
reg_sig_iti1.interaction = reg_cr1.interaction;
reg_sig_iti1.pvalThresh= 0.01;

reg_sig_iti2.coeff = all_coeff_iti2;
reg_sig_iti2.pval = all_pval_iti2;
reg_sig_iti2.regr_time = reg_cr2.regr_time;
reg_sig_iti2.numPredictor = reg_cr2.numPredictor;
reg_sig_iti2.nback = reg_cr2.nback;
reg_sig_iti2.interaction = reg_cr2.interaction;
reg_sig_iti2.pvalThresh= 0.01;

reg_sig_iti3.coeff = all_coeff_iti3;
reg_sig_iti3.pval = all_pval_iti3;
reg_sig_iti3.regr_time = reg_cr3.regr_time;
reg_sig_iti3.numPredictor = reg_cr3.numPredictor;
reg_sig_iti3.nback = reg_cr3.nback;
reg_sig_iti3.interaction = reg_cr3.interaction;
reg_sig_iti3.pvalThresh= 0.01;

tlabel={'R(n)','R(n-1)','R(n-2)'};
MP_plot_regr_iti(reg_sig_iti1, reg_sig_iti2, reg_sig_iti3, reg_sig.pvalThresh,tlabel,xtitle);
print(gcf,'-dpng','MLR-choiceoutcome_ITI');    %png format
saveas(gcf, 'MLR-choiceoutcome_ITI', 'fig');
close all

end