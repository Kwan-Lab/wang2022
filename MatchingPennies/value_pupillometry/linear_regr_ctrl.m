function reg_cr_ctrl = linear_regr_ctrl(signal, t, event, eventTime, trialSubset, params, tlabel);


%% shuffle every factors to get a ctrl regression results
% input: event: matrix
%         params.ifplot: true for plot control regression

numFactors = size(event,2);

 for kk = 1:numFactors
     shuffled_mat = event;
     shuffled_mat(:,kk) = shuffled_mat(randperm(size(shuffled_mat,1)),kk);
     reg_cr_ctrl.(['trigEvent',num2str(kk)]) = linear_regr( signal, t, shuffled_mat, eventTime, trialSubset, params );
     % plot them to make sure shuffle is right
     if params.ifplot == 1
        MP_plot_regrcoef_pupil(reg_cr_ctrl.(['trigEvent',num2str(kk)]),params.pvalThresh,tlabel,params.xtitle);
        print(gcf,'-dpng',['MLR-choiceoutcome_ctrl_-',['trigEvent',num2str(kk)]]);    %png format
        saveas(gcf, ['MLR-choiceoutcome_ctrl_-',['trigEvent',num2str(kk)]], 'fig');
     end
 end