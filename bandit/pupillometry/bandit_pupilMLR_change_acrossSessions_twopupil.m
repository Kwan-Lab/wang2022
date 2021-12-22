function bandit_pupilMLR_change_acrossSessions_twopupil(dataIndex, savefigpath)

%% reference to Sul et al.2011
set(0, 'DefaultFigureRenderer', 'painters');
% averaged within subject
nFiles = size(dataIndex,1);

subject_mask = [];
animalList = unique(dataIndex.Animal);
% load the first file to get some parameters

all_coeff_future.left = [];all_coeff_future.right = [];
all_pval_future.left = [];all_pval_future.right = [];

for ii = 1:nFiles
    savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},['bandit',dataIndex.LogFileName{ii}(end-30:end-4),'_beh.mat']));
    
    saveRegNameleft = fullfile(savematpath,[fn_beh.name(1:end-7),'regCR_changeleft.mat']);
    saveRegNameright = fullfile(savematpath,[fn_beh.name(1:end-7),'regCR_changeright.mat']);
    if exist(saveRegNameleft)
        changeleft = load(saveRegNameleft);
        changeright = load(saveRegNameright);
        % get subject mask
        subject_mask(end+1) = find(strcmp(animalList, dataIndex.Animal{ii}));
        
        % load choice and reward regression
         
        all_coeff_future.left = cat(3,all_coeff_future.left, changeleft.reg_cr_future_change.coeff);
        all_pval_future.left = cat(3, all_pval_future.left, changeleft.reg_cr_future_change.pval);
    
        all_coeff_future.right = cat(3,all_coeff_future.right, changeright.reg_cr_future_change.coeff);
        all_pval_future.right = cat(3, all_pval_future.right, changeright.reg_cr_future_change.pval);
     
        %%  load control regression for choice and reward regression
       
        % load future control
        if ~exist('all_coeff_future_ctrl')  && isfield(changeleft,'reg_cr_change_ctrl')
            % initialize the control matrix
            fieldsName_future = fieldnames(changeleft.reg_cr_future_change_ctrl);
            for tt = 1:length(fieldsName_future)
                all_coeff_future_ctrl.left.(fieldsName_future{tt}) = [];
                all_pval_future_ctrl.left.(fieldsName_future{tt}) = [];
                all_coeff_future_ctrl.right.(fieldsName_future{tt}) = [];
                all_pval_future_ctrl.right.(fieldsName_future{tt}) = [];
            end
        end
        % load the regression results
        
        if isfield(changeleft,'reg_cr_future_change_ctrl')
            fieldsName_future = fieldnames(changeleft.reg_cr_future_change_ctrl);
            for uu = 1:length(fieldsName_future)
                all_coeff_future_ctrl.left.(fieldsName_future{uu}) = cat(3, all_coeff_future_ctrl.left.(fieldsName_future{uu}), changeleft.reg_cr_future_change_ctrl.(fieldsName_future{uu}).coeff);
                all_pval_future_ctrl.left.(fieldsName_future{uu}) = cat(3, all_pval_future_ctrl.left.(fieldsName_future{uu}), changeleft.reg_cr_future_change_ctrl.(fieldsName_future{uu}).pval);
            end
        end
        
        if isfield(changeright,'reg_cr_change_ctrl')
            fieldsName_future = fieldnames(changeright.reg_cr_future_change_ctrl);
            for uu = 1:length(fieldsName_future)
                all_coeff_future_ctrl.right.(fieldsName_future{uu}) = cat(3, all_coeff_future_ctrl.right.(fieldsName_future{uu}), changeright.reg_cr_future_change_ctrl.(fieldsName_future{uu}).coeff);
                all_pval_future_ctrl.right.(fieldsName_future{uu}) = cat(3, all_pval_future_ctrl.right.(fieldsName_future{uu}), changeright.reg_cr_future_change_ctrl.(fieldsName_future{uu}).pval);
            end
        end
    end
end
        



%% original linear regression
% 1. use bootstrap to get the average and 95% CI for each factor, plot the
% bar plot

%% other things can be done for pupil response:
% correlation: pupil response - latent variable
%                                                     - response time

reg_all.left.coeff= all_coeff_future.left(:,:,[1:11,13]);
reg_all.right.coeff= all_coeff_future.right(:,:,[1:11,13]);


% use bootstrp to get coefficient
reg_all.left = getBootstrp(reg_all.left, 0,0.05);
reg_all.right = getBootstrp(reg_all.right, 0,0.05);

reg_all.left.regr_time = changeleft.reg_cr_future_change.regr_time;
reg_all.left.numPredictor = changeleft.reg_cr_future_change.numPredictor;
reg_all.left.nback = changeleft.reg_cr_future_change.nback;
reg_all.left.interaction = changeleft.reg_cr_future_change.interaction;
reg_all.left.pvalThresh= 0.01;

reg_all.right.regr_time = changeright.reg_cr_future_change.regr_time;
reg_all.right.numPredictor = changeright.reg_cr_future_change.numPredictor;
reg_all.right.nback = changeright.reg_cr_future_change.nback;
reg_all.right.interaction = changeright.reg_cr_future_change.interaction;
reg_all.right.pvalThresh= 0.01;

% go to the save path
if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);

% %figure 1
% figure;
% xAxis = 1:10;
% bar(xAxis(2:10),reg_all.coeff_bootave(2:10),'FaceColor',[0.7,0.7,0.7])                
% 
% hold on
% erneg = reg_all.coeff_bootave(2:10)-reg_all.bootlow(2:10);
% erpos = reg_all.boothigh(2:10)-reg_all.coeff_bootave(2:10);
% er = errorbar(xAxis(2:10),reg_all.coeff_bootave(2:10),erneg,erpos);    
% 
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% set(gca,'xticklabel',{'C(n)','C(n-1)','C(n-2)','R(n)','R(n-1)','R(n-2)','C(n)xR(n)','C(n-1)xR(n-1)','C(n-2)xR(n-2)'})
% hold off
% pvalThresh = NaN;
% xtickangle(45)
% ylabel('Coefficients (a.u.)');
% title('Coefficient for pupil change - choice and reward');
pvalThresh = NaN;
xtitle = 'Time from cue (s)';
tlabel={'C(n+1)','C(n)','C(n-1)','C(n-2)','R(n+1)','R(n)', 'R(n-1)','R(n-2)',...
                     'C(n+1)*R(n+1)','C(n)*R(n)','C(n-1)*R(n-1)','C(n-2)*R(n-2)','Reward Rate','Cumulative Reward'};
%             
MP_plot_regrcoef_pupil_both(reg_all,pvalThresh,tlabel,xtitle);
MP_plot_regrcoef_pupil(reg_all.left,pvalThresh,tlabel,xtitle);
print(gcf,'-dpng','MLR-change-choiceoutcome-averageSession-left');    %png format
saveas(gcf, 'MLR-change-choiceoutcome-averageSession-left', 'fig');
saveas(gcf, 'MLR-change-choiceoutcome-averageSession-left', 'svg');
MP_plot_regrcoef_pupil(reg_all.right,pvalThresh,tlabel,xtitle);
print(gcf,'-dpng','MLR-change-choiceoutcome-averageSession-right');    %png format
saveas(gcf, 'MLR-change-choiceoutcome-averageSession-right', 'fig');
saveas(gcf, 'MLR-change-choiceoutcome-averageSession-right', 'svg');

% scatter plot

figure;
cmap = linspace(1,10,8);
reg_leftChoice = [];
reg_rightChoice = [];
for tt = 1:8
    
    sz  =40;
    reg_leftChoice = [reg_leftChoice;reg_all.left.coeff(41:60,3,tt)];
    reg_rightChoice = [reg_rightChoice; reg_all.right.coeff(41:60,3,tt)];
    for ii = 41:60
        
        scatter(reg_all.left.coeff(ii,3,tt), reg_all.right.coeff(ii,3,tt),sz, cmap(tt),'filled');
        hold on;
    end
end
plot([-2 2],[-2 2],'k','LineWidth',1);
xlabel('Coefficients of c(n), left pupil');
ylabel('Coefficients of c(n), right pupil');
axis square
print(gcf,'-dpng','MLR-change-choice-averageSession-scatter');    %png format
saveas(gcf, 'MLR-change-choice-averageSession-scatter', 'fig');
saveas(gcf, 'MLR-change-choice-averageSession-scatter', 'svg');
[rho_choice,pval_choice] = corr(reg_leftChoice, reg_rightChoice);


figure;
cmap = linspace(1,10,8);
reg_leftReward = [];
reg_rightReward = [];

for tt = 1:8
    
    sz  =40;
    for ii = 41:60
        reg_leftReward = [reg_leftReward;reg_all.left.coeff(41:60,7,tt)];
        reg_rightReward = [reg_rightReward; reg_all.right.coeff(41:60,7,tt)];
        scatter(reg_all.left.coeff(ii,7,tt), reg_all.right.coeff(ii,7,tt),sz, cmap(tt),'filled');
        hold on;
    end
end
plot([-2 2],[-2 2],'k','LineWidth',1);
axis square
xlabel('Coefficients of r(n), left pupil');
ylabel('Coefficients of r(n), right pupil');
print(gcf,'-dpng','MLR-change-outcome-averageSession-scatter');    %png format
saveas(gcf, 'MLR-change-outcome-averageSession-scatter', 'fig');
saveas(gcf, 'MLR-change-outcome-averageSession-scatter', 'svg');
[rho_reward,pval_reward] = corr(reg_leftReward, reg_rightReward);

%% plot the figure as number of session that is significant
reg_sig.left.coeff = all_coeff_future.left(:,:,[1:11,13]);
reg_sig.left.pval = all_pval_future.left(:,:,[1:11,13]);
reg_sig.left.regr_time = reg_all.left.regr_time;
reg_sig.left.numPredictor = reg_all.left.numPredictor;
reg_sig.left.nback = reg_all.left.nback;
reg_sig.left.interaction = reg_all.left.interaction;
reg_sig.left.pvalThresh= 0.01;

reg_sig.right.coeff = all_coeff_future.right(:,:,[1:11,13]);
reg_sig.right.pval = all_pval_future.right(:,:,[1:11,13]);
reg_sig.right.regr_time = reg_all.right.regr_time;
reg_sig.right.numPredictor = reg_all.right.numPredictor;
reg_sig.right.nback = reg_all.right.nback;
reg_sig.right.interaction = reg_all.right.interaction;
reg_sig.right.pvalThresh= 0.01;


if exist('all_coeff_ctrl')
    reg_pval_ctrl.left = getBootstrp(all_pval_ctrl.left, 0.01);
    reg_pval_ctrl.right = getBootstrp(all_pval_ctrl.right, 0.01);
% MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
    MP_plot_regr(reg_sig.left,reg_pval_ctrl.left, reg_sig.left.pvalThresh,tlabel,xtitle);
    print(gcf,'-dpng','MLR-change-choiceoutcome-sigSession-left');    %png format
    saveas(gcf, 'MLR-change-choiceoutcome-sigSession-left', 'fig');
    saveas(gcf, 'MLR-change-choiceoutcome-averageSession-left', 'svg');
    MP_plot_regr(reg_sig.right,reg_pval_ctrl.right, reg_sig.right.pvalThresh,tlabel,xtitle);
    print(gcf,'-dpng','MLR-change-choiceoutcome-sigSession-right');    %png format
    saveas(gcf, 'MLR-change-choiceoutcome-sigSession-right', 'fig');
    saveas(gcf, 'MLR-change-choiceoutcome-averageSession-right', 'svg');
else
    MP_plot_regr(reg_sig.left,[], reg_sig.left.pvalThresh,tlabel,xtitle);
    print(gcf,'-dpng','MLR-change-choiceoutcome-sigSession-left');    %png format
    saveas(gcf, 'MLR-change-choiceoutcome-sigSession-left', 'fig');
    saveas(gcf, 'MLR-change-choiceoutcome-averageSession-left', 'svg');
    MP_plot_regr(reg_sig.right,[], reg_sig.right.pvalThresh,tlabel,xtitle);
    print(gcf,'-dpng','MLR-change-choiceoutcome-sigSession-right');    %png format
    saveas(gcf, 'MLR-change-choiceoutcome-sigSession-right', 'fig');
    saveas(gcf, 'MLR-change-choiceoutcome-averageSession-right', 'svg');
end




close all

end