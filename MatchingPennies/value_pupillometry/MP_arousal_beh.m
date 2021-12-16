function MP_arousal_beh(dataIndex, savefigpath)

%% analyze some behavior for low/high baseline pupil
% run separate linear regressions\

nFiles = size(dataIndex,1);


%% load the behavior files first to get the ITI time

%% go through every session, running MLR separately
lowTrialList = zeros(1, nFiles);
highTrialList = zeros(1, nFiles);
lowRewardList = zeros(1, nFiles);
highRewardList =zeros(1, nFiles);
lowWSLSList = zeros(1, nFiles);
highWSLSList = zeros(1, nFiles);
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
        savepupilfigpath = fullfile(dataIndex.BehPath{ii},[date(1:6),'figs_summary']);
        if ~exist(savepupilfigpath,'dir')
            mkdir(savepupilfigpath);
        end
        
        cd(savepupilfigpath);
        
        choice = NaN(size(trials.left));
        choice(trialData.response == 2) = 0;
        choice(trialData.response == 3) = 1;
    
        reward = NaN(size(trials.left));
        reward(trialData.response ~= 0 & trials.reward == 0) = 0;
        reward(trialData.response ~= 0 & trials.reward == 1) = 1;
        
         trialInd = 1:length(pupil.base);
      midPoint = prctile(pupil.base,50);
       lowInd = trialInd(pupil.base<midPoint); highInd = trialInd(pupil.base>=midPoint);
       
        % trial number
        lowTrialList(ii) = sum(pupil.base<midPoint);
        highTrialList(ii) = sum(pupil.base>=midPoint);
      % average reward rate
        lowRewardList(ii) = nanmean(reward(pupil.base<midPoint));
        highRewardList(ii) = nanmean(reward(pupil.base>=midPoint));
        
      % probability of win-stay-lose-switch
     
     
      lowWSLS = 0; highWSLS = 0; totalLow = 0; totalHigh = 0;
      for ww = 1:length(lowInd)
          if lowInd(ww) - 1 > 0 & ~isnan(choice(lowInd(ww))) & ~isnan(choice(lowInd(ww)-1))
              totalLow = totalLow+1;
            if reward(lowInd(ww)-1) == 1 % win
                if choice(lowInd(ww)) == choice(lowInd(ww) -1)
                    lowWSLS = lowWSLS + 1; % win-stay
                end
            else % lose
                if choice(lowInd(ww)) ~= choice(lowInd(ww) - 1)
                    lowWSLS = lowWSLS + 1; % lose-switch
                end
            end
          end
      end
      lowWSLSList(ii) = lowWSLS/totalLow;
      for ww = 1:length(highInd)
          if highInd(ww) - 1 > 0 & ~isnan(choice(highInd(ww))) & ~isnan(choice(highInd(ww)-1))
              totalHigh = totalHigh+1;
            if reward(highInd(ww)-1) == 1 % win
                if choice(highInd(ww)) == choice(highInd(ww) -1)
                    highWSLS = highWSLS + 1; % win-stay
                end
            else % lose
                if choice(highInd(ww)) ~= choice(highInd(ww) - 1)
                    highWSLS = highWSLS + 1; % lose-switch
                end
            end
          end
      end
        highWSLSList(ii) = highWSLS/totalHigh;
     
        %% save all regression into one mat file
        % save all these in a structure
        %save(saveRegName, 'reg_cr', 'reg_cr1','reg_cr2','reg_cr3','reg_cr_future','reg_cr_future_ctrl','reg_cr_ctrl');
        
        %save(saveRegName_ITI,'reg_cr1_change_1','reg_cr2_change_1','reg_cr3_change_1');
        close all;
    end
end

cd(savefigpath)

totalInd = 1:nFiles;
nanZero = totalInd(highWSLSList~=0);


tempTrial=[lowTrialList(nanZero);highTrialList(nanZero)];

figure;
boxplot(tempTrial','Labels',{'Low','High'});
set(gca,'box','off');
print(gcf,'-dpng','TrialNum-lowhigh');    %png format
saveas(gcf, 'TrialNum-lowhigh', 'fig');
saveas(gcf, 'TrialNum-lowhigh','svg');

[h,p] = ttest(lowTrialList(nanZero),highTrialList(nanZero))

tempReward = [lowRewardList(nanZero);highRewardList(nanZero)];

figure;
boxplot(tempReward','Labels',{'Low','High'});
set(gca,'box','off');
print(gcf,'-dpng','Reward-lowhigh');    %png format
saveas(gcf, 'Reward-lowhigh', 'fig');
saveas(gcf, 'Reward-lowhigh','svg');
[h,p] = ttest(lowRewardList(nanZero),highRewardList(nanZero))

tempWSLS = [lowWSLSList(nanZero);highWSLSList(nanZero)];

figure;
boxplot(tempWSLS','Labels',{'Low','High'});
set(gca,'box','off');
print(gcf,'-dpng','WSLS-lowhigh');    %png format
saveas(gcf, 'WSLS-lowhigh', 'fig');
saveas(gcf, 'WSLS-lowhigh','svg');
[h,p] = ttest(lowWSLSList(nanZero),highWSLSList(nanZero))

end
