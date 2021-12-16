function MP_pupil_AS_arousal_acrossSessions(dataIndex, savefigpath)

%% reference to Sul et al.2011
% averaged within subject
setup_figprop
set(0, 'DefaultFigureRenderer', 'painters');
nFiles = size(dataIndex,1);

subject_mask = [];
animalList = unique(dataIndex.Animal);

% load the first file to get some parameters
all_coeff_low = [];
all_pval_low = [];

all_coeff_high = [];
all_pval_high = [];

for ii = 1:nFiles
    savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},[dataIndex.LogFileName{ii}(1:end-4),'_beh_cut.mat']));
    
    saveRegName = fullfile(savematpath,[fn_beh.name(1:end-7),'regRL_AS_arousal.mat']);
    %saveRegName_ITI = fullfile(savematpath,[fn_beh.name(1:end-7),'regCR_ITI.mat']);

    if exist(saveRegName)
        load(saveRegName)
        
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

        all_coeff_low = cat(3,all_coeff_low, reg_cr_low.coeff);
        all_pval_low = cat(3, all_pval_low, reg_cr_low.pval);
        
        all_coeff_high = cat(3,all_coeff_high, reg_cr_high.coeff);
        all_pval_high = cat(3, all_pval_high, reg_cr_high.pval);
        
        reg_all.regr_time = reg_cr_low.regr_time;
        reg_all.numPredictor = reg_cr_low.numPredictor;
        reg_all.nback = reg_cr_low.nback;
        reg_all.interaction = reg_cr_low.interaction;
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
        
%         all_coeff_iti_1 = cat(3,all_coeff_iti_1, reg_cr1_change_1.coeff);
%         all_pval_iti_1 = cat(3,all_pval_iti_1, reg_cr1_change_1.pval);
%         
%         all_coeff_iti_2 = cat(3,all_coeff_iti_2, reg_cr2_change_1.coeff);
%         all_pval_iti_2 = cat(3,all_pval_iti_2, reg_cr2_change_1.pval);
%        
%         all_coeff_iti_3 = cat(3,all_coeff_iti_3, reg_cr3_change_1.coeff);
%         all_pval_iti_3 = cat(3,all_pval_iti_3, reg_cr3_change_1.pval);

        %%  load control regression for choice and reward regression
        %no control for now
%         
%         if ~exist('all_coeff_ctrl') && exist('reg_cr_change_ctrl')
%             % initialize the control matrix
%             fieldsName = fieldnames(reg_cr_change_ctrl);
%             for tt = 1:length(fieldsName)
%                 all_coeff_ctrl.(fieldsName{tt}) = [];
%                 all_pval_ctrl.(fieldsName{tt}) = [];
%             end
%             
%         end
%         if exist('reg_cr_change_ctrl')
%         % load the regression results
%             fieldsName = fieldnames(reg_cr_change_ctrl);
%             for uu = 1:length(fieldsName)
%                 all_coeff_ctrl.(fieldsName{uu}) = cat(3, all_coeff_ctrl.(fieldsName{uu}), reg_cr_change_ctrl.(fieldsName{uu}).coeff);
%                 all_pval_ctrl.(fieldsName{uu}) = cat(3, all_pval_ctrl.(fieldsName{uu}), reg_cr_change_ctrl.(fieldsName{uu}).pval);
%             end
%             
%         end   
%       
%         % load future control
%         if ~exist('all_coeff_future_ctrl')  && exist('reg_cr_future_change_ctrl')
%             % initialize the control matrix
%             fieldsName_future = fieldnames(reg_cr_future_change_ctrl);
%             for tt = 1:length(fieldsName_future)
%                 all_coeff_future_ctrl.(fieldsName_future{tt}) = [];
%                 all_pval_future_ctrl.(fieldsName_future{tt}) = [];
%             end
%         end
%         % load the regression results
%         
%         if exist('reg_cr_future_change_ctrl')
%             fieldsName_future = fieldnames(reg_cr_future_change_ctrl);
%             for uu = 1:length(fieldsName_future)
%                 all_coeff_future_ctrl.(fieldsName_future{uu}) = cat(3, all_coeff_future_ctrl.(fieldsName_future{uu}), reg_cr_future_change_ctrl.(fieldsName_future{uu}).coeff);
%                 all_pval_future_ctrl.(fieldsName_future{uu}) = cat(3, all_pval_future_ctrl.(fieldsName_future{uu}), reg_cr_future_change_ctrl.(fieldsName_future{uu}).pval);
%             end
%         end
    end
end
        



%% original linear regression
% 1. use bootstrap to get the average and 95% CI for each factor, plot the
% bar plot

%% other things can be done for pupil response:
% correlation: pupil response - latent variable
%                                                     - response time

% reg_all.coeff= all_coeff;
% 
% % use bootstrp to get coefficient
% reg_all = getBootstrp(reg_all, 0, 0.05);
% 

% go to the save path
if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);


%%  linear regression with C(n+1)
reg_cr_low.coeff= all_coeff_low;
reg_cr_high.coeff= all_coeff_high;
% use bootstrp to get coefficient
reg_cr_low = getBootstrp(reg_cr_low, 0, 0.05);
reg_cr_high = getBootstrp(reg_cr_high, 0, 0.05);

if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);
% 
% figure;
% xAxis = 1:8;
% bar(xAxis(2:8),reg_cr_future.coeff_bootave(2:8),'FaceColor',[0.7,0.7,0.7])                
% 
% hold on
% erneg = reg_cr_future.coeff_bootave(2:8)-reg_cr_future.bootlow(2:8);
% erpos =reg_cr_future.boothigh(2:8)-reg_cr_future.coeff_bootave(2:8);
% er = errorbar(xAxis(2:8),reg_cr_future.coeff_bootave(2:8),erneg,erpos);    
% 
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% set(gca,'xticklabel',{'C(n)','C(n+1)','C(n-1)','C(n-2)','R(n)','R(n-1)','R(n-2)'})
% hold off
% 
% xtickangle(45)
% ylabel('Coefficients (a.u.)');
% title('Coefficient for pupil change - choice and reward');
xtitle='Time from cue (s)';
 tlabel={'C(n)','R(n)','C(n)xR(n)','C(n-1)','R(n-1)', 'C(n-1)xR(n-1)','QL-QR', 'ChosenQ', 'QLC-QRC', 'ChosenQC','Reward Rate', 'Cumulavtive reward'};
pvalThresh=0.01;
MP_plot_regrcoef_pupil_two(reg_cr_low,reg_cr_high,{'low','high'},pvalThresh,tlabel,xtitle);
print(gcf,'-dpng','MLR-AS_arousal_averageSession');    %png format
saveas(gcf, 'MLR-AS_arousal_averageSession', 'fig');
saveas(gcf, 'MLR-AS_arousal_averageSession','svg');

% plot the figure as number of session that is significant
reg_sig1.coeff = all_coeff_low;
reg_sig1.pval = all_pval_low;
reg_sig2.coeff = all_coeff_high;
reg_sig2.pval = all_pval_high;
reg_sig1.regr_time = reg_cr_low.regr_time;
reg_sig1.numPredictor = reg_cr_low.numPredictor;
reg_sig1.nback = reg_cr_low.nback;
reg_sig1.interaction = reg_cr_low.interaction;
reg_sig1.pvalThresh= 0.01;

% preprocessing the control coefficient and pval: bootstrap the pval to
% determine a baseline percentage of significant session
%reg_pval_futurectrl = getBootstrp(all_pval_futurectrl, 0.01);

if exist('all_pval_future_ctrl')
    reg_pval_future_ctrl = getBootstrp(all_pval_future_ctrl, 0.01, 0.05);

% MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
    MP_plot_regr(reg_sig1,reg_pval_future_ctrl, reg_sig1.pvalThresh,tlabel,xtitle);
else
    MP_plot_regr_two(reg_sig1,reg_sig2,[], reg_sig1.pvalThresh,tlabel,xtitle);
end
print(gcf,'-dpng','MLR-AS_future_sigSession');    %png format
saveas(gcf, 'MLR-AS_future_sigSession', 'fig');
saveas(gcf, 'MLR-AS_future_sigSession','svg');


%%
close all

end