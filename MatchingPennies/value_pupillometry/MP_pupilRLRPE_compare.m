function MP_pupilRLRPE_compare(dataIndex, savefigpath);

%% reference to Sul et al.2011

nFiles = size(dataIndex,1);

% load the first file to get some parameters

all_coeff_pos_new = [];
all_coeff_neg_new = [];
all_pval_pos_new = [];
all_pval_neg_new = [];

all_coeff_pos = [];
all_coeff_neg = [];
all_pval_pos = [];
all_pval_neg = [];

all_coeff = [];
all_coeff_new = [];
% all_coeff_updatedQ = [];
% all_pval_updatedQ = [];

for ii = 1:nFiles
    
    % load behavior files
    savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
    savematpath_new = fullfile('J:\MatchingPennies\pupilData\analysis_new',dataIndex.Animal{ii},'analysis-pupil');
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},[dataIndex.LogFileName{ii}(1:end-4),'_beh_cut.mat']));
    saveRegName = fullfile(savematpath,[fn_beh.name(1:end-7),'regRL_RPE_change_CK.mat']);
    saveRegName_new = fullfile(savematpath_new,[fn_beh.name(1:end-7),'regRL_RPE_change_CK.mat']);
    if exist(saveRegName)
        load(saveRegName);  
        
         all_coeff= cat(3,all_coeff, reg_cr_RPE_change.coeff);
        
        
        all_coeff_pos = cat(3,all_coeff_pos, reg_cr_RPEpos_change.coeff);
        all_pval_pos = cat(3,all_pval_pos, reg_cr_RPEpos_change.pval);
        
        all_coeff_neg = cat(3,all_coeff_neg, reg_cr_RPEneg_change.coeff);
        all_pval_neg = cat(3,all_pval_neg, reg_cr_RPEneg_change.pval);
    end
    
    if exist(saveRegName_new)
        load(saveRegName_new);      
         all_coeff_new= cat(3,all_coeff_new, reg_cr_RPE_change.coeff);
        
        all_coeff_pos_new = cat(3,all_coeff_pos_new, reg_cr_RPEpos_change.coeff);
        all_pval_pos_new = cat(3,all_pval_pos_new, reg_cr_RPEpos_change.pval);
        
        all_coeff_neg_new = cat(3,all_coeff_neg_new, reg_cr_RPEneg_change.coeff);
        all_pval_neg_new = cat(3,all_pval_neg_new, reg_cr_RPEneg_change.pval);
    end
      reg_all.regr_time = reg_cr_RPEneg_change.regr_time;
        reg_all.numPredictor = reg_cr_RPEneg_change.numPredictor;
        reg_all.interaction = reg_cr_RPEneg_change.interaction;
        
end



% % plot pos/neg RPE

reg_all_pos = reg_all;  % get other parameters from reg_all
reg_all_pos.coeff= all_coeff_pos;
reg_all_pos = getBootstrp(reg_all_pos,0,0.05);

reg_all_neg = reg_all;
reg_all_neg.coeff= all_coeff_neg;
reg_all_neg = getBootstrp(reg_all_neg, 0,0.05);

reg_all_pos_new = reg_all;  % get other parameters from reg_all
reg_all_pos_new.coeff= all_coeff_pos_new(:,:,1:55);
reg_all_pos_new = getBootstrp(reg_all_pos_new,0,0.05);

reg_all_neg_new = reg_all;
reg_all_neg_new.coeff= all_coeff_neg_new;
reg_all_neg_new = getBootstrp(reg_all_neg_new, 0,0.05);

pvalThresh = 0.01;
xtitle = 'Time from cue (s)';
MP_plot_regrcoef_pupil(reg_all_pos,pvalThresh,tlabel3,xtitle);
MP_plot_regrcoef_pupil(reg_all_pos_new,pvalThresh,tlabel3,xtitle);
tlabel3={'C(n)','C(n-1)','R(n-1)', 'dQ','posRPE','dK','CKE','Average reward', 'Cumulative reward'};
tlabel4={'C(n)','C(n-1)','R(n-1)', 'dQ','negRPE','dK','CKE','Average reward', 'Cumulative reward'};
% MP_plot_regrcoef_pupil(reg_all_pos,pvalThresh,tlabel3,xtitle);
% print(gcf,'-dpng','MLR-RL_averageSession_posRPE_change');    %png format
% saveas(gcf, 'MLR-RL_averageSession_posRPE_change', 'fig');
% 
% MP_plot_regrcoef_pupil(reg_all_neg,pvalThresh,tlabel4,xtitle);
% print(gcf,'-dpng','MLR-RL_averageSession_negRPE_change');    %png format
% saveas(gcf, 'MLR-RL_averageSession_negRPE_change', 'fig');

%% plot pos/neg RPE only
reg_all_sign = reg_all_pos;
reg_all_sign.numPredictor = 2;
% the first one doesnt matter
reg_all_sign.coeff = [reg_all_pos.coeff(:,1,:),reg_all_pos.coeff(:,6,:),reg_all_neg.coeff(:,6,:)];
reg_all_sign.bootSig = [reg_all_pos.bootSig(:,1,:),reg_all_pos.bootSig(:,6,:),reg_all_neg.bootSig(:,6,:)];
reg_all_sign.coeff_bootave = [reg_all_pos.coeff_bootave(:,1),reg_all_pos.coeff_bootave(:,6),reg_all_neg.coeff_bootave(:,6)];
reg_all_sign.bootlow = [reg_all_pos.bootlow(:,1),reg_all_pos.bootlow(:,6), reg_all_neg.bootlow(:,6)];
reg_all_sign.boothigh = [reg_all_pos.boothigh(:,1),reg_all_pos.boothigh(:,6), reg_all_neg.boothigh(:,6)];
tlabel = {'posRPE', 'negRPE'};

MP_plot_regrcoef_pupil(reg_all_sign,pvalThresh,tlabel,xtitle);
print(gcf,'-dpng','MLR-RL_averageSession_allRPE_change');    %png format
saveas(gcf, 'MLR-RL_averageSession_allRPE_change', 'fig');

%% plot the figure as number of session that is significant

% RPE
reg_sig.coeff = all_coeff;
reg_sig.pval = all_pval;
reg_sig.regr_time = reg_all.regr_time;
reg_sig.numPredictor = reg_all.numPredictor;
reg_sig.interaction = reg_all.interaction;
reg_sig.pvalThresh= 0.01;
reg_sig.nback = 0;
xtitle = {'Time from cue (s)'};


%reg_pval_ctrl = getBootstrp(all_pval_ctrl, 0.01);
if exist('all_coeff_RPE_change_ctrl')
    reg_pval_ctrl = getBootstrp(all_pval_RPE_change_ctrl, 0.01,0.01);

% MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
    MP_plot_regr(reg_sig,reg_pval_ctrl, reg_sig.pvalThresh,tlabel1,xtitle);
else
    MP_plot_regr(reg_sig,[], reg_sig.pvalThresh,tlabel1,xtitle);
end

print(gcf,'-dpng','MLR-RL_RPE_change_CK');    %png format
saveas(gcf, 'MLR-RL_RPE_change_CK', 'fig');
saveas(gcf, 'MLR-RL_RPE_change_CK','svg');

% updatedQ
% reg_sig_updatedQ.coeff = all_coeff_updatedQ;
% reg_sig_updatedQ.pval = all_pval_updatedQ;
% reg_sig_updatedQ.regr_time = reg_all_updatedQ.regr_time;
% reg_sig_updatedQ.numPredictor = reg_all_updatedQ.numPredictor;
% reg_sig_updatedQ.interaction = reg_all_updatedQ.interaction;
% reg_sig_updatedQ.pvalThresh= 0.01;
% reg_sig_updatedQ.nback = 0;
% xtitle = {'Time from cue (s)'};
% 
% %reg_pval_updatedQctrl = getBootstrp(all_pval_updatedQctrl, 0.01);
% 
% if exist('all_coeff_updatedQctrl')
%     reg_pval_ctrl = getBootstrp(all_pval_updatedQctrl, 0.01);
% 
% % MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
%     MP_plot_regr(reg_sig_updatedQ,reg_pval_ctrl, reg_sig.pvalThresh,tlabel2,xtitle);
% else
%     MP_plot_regr(reg_sig_updatedQ,[], reg_sig.pvalThresh,tlabel2,xtitle);
% end
% print(gcf,'-dpng','MLR-RL_updatedQ_change_CK');    %png format
% saveas(gcf, 'MLR-RL_updatedQ_change_CK', 'fig');
% saveas(gcf, 'MLR-RL_updatedQ_change_CK','svg');


%% plot the positive and negative RPE
reg_sig_pos.coeff = all_coeff_pos;
reg_sig_pos.pval = all_pval_pos;
reg_sig_pos.regr_time = reg_all.regr_time;
reg_sig_pos.numPredictor = reg_all.numPredictor;
reg_sig_pos.interaction = reg_all.interaction;
reg_sig_pos.pvalThresh= 0.01;
reg_sig_pos.nback = 0;
xtitle = {'Time from cue (s)'};

%reg_pval_posctrl = getBootstrp(all_pval_ctrl, 0.01);
if exist('all_coeff_pos_ctrl')
    reg_pval_ctrl = getBootstrp(all_pval_pos_ctrl, 0.01);

% MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
    MP_plot_regr(reg_sig_pos,reg_pval_ctrl, reg_sig.pvalThresh,tlabel3,xtitle);
else
    MP_plot_regr(reg_sig_pos,[], reg_sig.pvalThresh,tlabel3,xtitle);
end

print(gcf,'-dpng','MLR-RL_posRPE_change');    %png format
saveas(gcf, 'MLR-RL_posRPE_change', 'fig');

reg_sig_neg.coeff = all_coeff_neg;
reg_sig_neg.pval = all_pval_neg;
reg_sig_neg.regr_time = reg_all.regr_time;
reg_sig_neg.numPredictor = reg_all.numPredictor;
reg_sig_neg.interaction = reg_all.interaction;
reg_sig_neg.pvalThresh= 0.01;
reg_sig_neg.nback = 0;
xtitle = {'Time from cue (s)'};

if exist('all_coeff_neg_ctrl')
    reg_pval_ctrl = getBootstrp(all_pval_neg_ctrl, 0.01);

% MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
    MP_plot_regr(reg_sig_neg,reg_pval_ctrl, reg_sig.pvalThresh,tlabel4,xtitle);
else
    MP_plot_regr(reg_sig_neg,[], reg_sig.pvalThresh,tlabel4,xtitle);
end
print(gcf,'-dpng','MLR-RL_negRPE_change');    %png format
saveas(gcf, 'MLR-RL_negRPE_change', 'fig');
% 
% % plot pos/neg together
% reg_all.coeff = [all_coeff_pos(:,1,:),all_coeff_pos(:,end,:),all_coeff_neg(:,end,:)];
% reg_all.pvalThresh = 0.01;
% reg_all.nback = 0;
% reg_all.numPredictor = 2;
% reg_all.pval = [reg_sig_pos.pval(:,1,:),reg_sig_pos.pval(:,end,:),reg_sig_neg.pval(:,end,:)];
% MP_plot_regr(reg_all,[],reg_all.pvalThresh,tlabel,xtitle);
% print(gcf,'-dpng','MLR-RL_posnegRPE_change');    %png format
% saveas(gcf, 'MLR-RL_posnegRPE_change', 'fig');


close all

end