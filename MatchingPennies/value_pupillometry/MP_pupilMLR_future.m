function MP_pupilMLR_future(dataIndex)

%% reference to Sul et al.2011

% now try to get rid of the fatigue trials
% add future one trial into regression

nFiles = size(dataIndex,1);

for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},[dataIndex.LogFileName{ii}(1:end-4),'_beh_cut.mat']));
    load(fullfile(fn_beh.folder,fn_beh.name));
    

    % load pupil files
    date = num2str(dataIndex.DateNumber(ii));
    pup_name = fullfile(dataIndex.BehPath{ii}, ['*',date(1:6),'*_pup.mat']);
    fn_pup = dir(pup_name);
    if length(fn_pup) == 1
        
        savepupilfigpath = fullfile(dataIndex.BehPath{ii},[date(1:6),'_figs-pupil']);
        if ~exist(savepupilfigpath,'dir')
            mkdir(savepupilfigpath);
        end
        
        cd(savepupilfigpath);
        
    
        load(fullfile(fn_pup.folder,fn_pup.name));
        
        % make folders to save analysis and plots
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
        if ~exist(savematpath,'dir')
            mkdir(savematpath);
        end
        savepupilfigpath = fullfile(dataIndex.BehPath{ii},[date(1:6),'_figs-pupil']);
        if ~exist(savepupilfigpath,'dir')
            mkdir(savepupilfigpath);
        end
        
        cd(savepupilfigpath);
        
        params=[];
        
        if trialData.cutPoint ~= 0
            % cut the trial
            fn = fieldnames(trials);
            for tt = 1:length(fn)
                trials.(fn{tt}) = trials.(fn{tt})(1:trialData.cutPoint);
            end
            params.trigTime = trialData.cueTimes(1:trialData.cutPoint);
        else
            params.trigTime = trialData.cueTimes;
        end
        %first predictor is choice; dummy-code: left=0, right=1, miss=NaN
        
            params.trigEvent=NaN(size(trials.left));
            params.trigEvent(trials.left) = 0;
            params.trigEvent(trials.right) = 1;
            params.trigEvent2 = [params.trigEvent(2:end);NaN];  % C(n+1)
            
            params.trigEvent3 = [NaN;params.trigEvent(1:end-1)];   % C(n-1)
            %second predictor is outcome; dummy-code: reward=1, error=0, miss=NaN
            params.trigEvent4 = [NaN;NaN;params.trigEvent(1:end-2)]; % C(n-2)
            
            params.trigEvent5=NaN(size(trials.go)); % R(n)
            
            params.trigEvent5(trials.reward) = 1;
            params.trigEvent5(trials.noreward) = 0;
            % params.trigEvent2(trials.omitreward) = 0;
            params.trigEvent6 = [NaN;params.trigEvent5(1:end-1)]; % R(n-1)
            params.trigEvent7 = [NaN;NaN;params.trigEvent5(1:end-2)];
            
            % params.trigEvent2(trials.doublereward) = 2;

            fieldname={'go'};
            trialMask = getMask(trials,fieldname);
       
        
        params.xtitle = 'Time from cue (s)';
        params.window = [-3:0.1:5];
        params.nback = 0;       %how many trials back to regress against
        params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
        params.interaction = false;
        %only perform analysis on this subset of trials
        
        reg_cr=linear_regr( pupil.dia, pupil.t, [params.trigEvent2 params.trigEvent params.trigEvent3 params.trigEvent4 params.trigEvent5 params.trigEvent6 params.trigEvent7], params.trigTime, trialMask, params );
        
        tlabel={'C(n+1)','C(n)','C(n-1)','C(n-2)','R(n)','R(n-1)','R(n-2)'};
        % plot_regr_pupil(reg_cr,params.pvalThresh,tlabel,params.xtitle);
        %print(gcf,'-dpng','MLR-choiceoutcome_sig_choice0_1');    %png format
        % saveas(gcf, 'MLR-choiceoutcome_sigchoice0_1', 'fig');
        MP_plot_regrcoef_pupil(reg_cr,params.pvalThresh,tlabel,params.xtitle);
        saveRegName = fullfile(savematpath,[fn_beh.name(1:end-7),'regCR_cut_future.mat']);
        
        save(saveRegName, 'reg_cr');
        print(gcf,'-dpng','MLR-choiceoutcome_cut_c(n+1)');    %png format
        saveas(gcf, 'MLR-choiceoutcome_cut_c(n+1)', 'fig');
        
        close all;
    end
end
end
