clear all;

setup_figprop;
%add a plot showing how many trials in each session

%combine the whole sessions (including algorithm1 and algorithm2)
base_dir='J:\MatchingPennies\pupilData\learningData';
savebehfigpath=[base_dir(1:end-12),'figs_summary_learning\'];


cd(base_dir);
subjects_dir = dir();

%A1_dir = 'model_A1';
%A2_dir = 'model_A2';

for ii = 3:length(subjects_dir)
%need in algorithm 2
    sub_dir = fullfile(subjects_dir(ii).folder, subjects_dir(ii).name);


    phase='phase2';
    %phase1: lick suppression analysis
    %algorithm 2

    %bandit_setPathList(base_dir);
    rewardCode=[10, 100, 111];
    incorrectCode=[110, 101];
%data=
%'/Users/phoenix/Documents/Kwanlab/reinforcement_learning/logfile/human/170511/';\


        %get successful training subject
%         subList=cell(0);
%         for i=1:length(subjects)
%         if subjects(i).name(end-1:end)~='NA'
%             subList{end+1}=subjects(i).name;
%         end
%         end

    if phase=='phase1' %lick suppression analysis
        noLick=cell(1, length(subList));
        for ii=1:length(subList)
            path=[sub_dir,subList{ii},'\phase1'];
            cd(path);
            logfiles=dir('*.log');
            nolickList=zeros(1, length(logfiles));
            for x =1:length(logfiles)
                disp(x);
                [ logData ] = parseLogfileHW(path, logfiles(x).name);
                
                [ sessionData, trialData] = MP_getSessionData( logData );
                episodeStart=-1;
                hit_miss = zeros(1, sessionData.nTrials); %record the hit/miss of every trial for later
                nolick_list = zeros(1, sessionData.nTrials);
                start_nolick = 1;
                
                for i=1:length(trialData.startTimes)-1
                    episodeEnd = trialData.startTimes(i+1) - trialData.startTimes(i) - 0.1;
                    if i>1      %if it's not first episode, also look at negative time as specified by episodeStart
                        episodeIndex = sessionData.time > (trialData.startTimes(i)+episodeStart) & sessionData.time < (trialData.startTimes(i)+episodeEnd);
                    else
                        episodeIndex = sessionData.time >= trialData.startTimes(i) & sessionData.time < (trialData.startTimes(i)+episodeEnd);
                    end
                    episodeIndexStart=find(episodeIndex,1,'first');     %getting the events that within [episodeStart, episodeEnd]
                    episodeIndexEnd=find(episodeIndex,1,'last');     %getting the events that within [episodeStart, episodeEnd]
                    if i < numel(trialData.startTimes)
                        endtrialIdx = find(sessionData.time == trialData.startTimes(i+1), 1);
                    end
                    clear episodeIndex;
                    
                    num_nolick = 0;
                    %for ii =1:numel(sessionData.nTrials)
                    for j=episodeIndexStart:episodeIndexEnd
                        if sessionData.code(j) == 19 & j < endtrialIdx
                            num_nolick = num_nolick+1;
                        end
                    end
                    %end
                    nolick_list(i) = num_nolick;
                    
                    ifReward = 0;
                    for j=episodeIndexStart:episodeIndexEnd
                        if sessionData.code(j) == 10 & j < endtrialIdx
                            ifReward = 1;
                            rewardtime = sessionData.time(j);
                            hit_miss(i)=1;
                            break;
                        end
                    end
                    ave_nolick = mean(nolick_list(logical(hit_miss)));
                    nolickList(x)=ave_nolick;
                end
            end
            nolick{ii}=nolickList;
        end
        
        figure;
        
        for jj=1:length(nolick)
            plot(nolick{jj}, 'Linewidth',2)
            hold on;
        end
        legend(subList);
        xlabel('Sessions');
        ylabel('# of ');
        
        
    
    elseif phase=='phase2'
        cd(sub_dir);
        logfiles = dir('*.log');
        
        % sort the files
        [~,idx] = sort([logfiles.datenum]);
        logfiles = logfiles(idx);
        c_nomiss1=[];r_nomiss1=[]; com_nomiss1=[];
        num_session1=[];
        cutTrial = [];
        
        %num_session2=[];
        %get concatenated choice and reward history
        for i =1:length(logfiles)
            [ logData ] = MP_parseLogfile(sub_dir, logfiles(i).name);
        
            [ sessionData, trialData] = MP_getSessionData( logData );
            
            % check for cut point
            trials = MP_getTrialMasks(trialData);
            % calculate the running entropy
            running_window = 30;
            running_entropy = cal_runningEntropy(trialData, running_window);
            
            % get the cut point
            [cutPoint, ave_entropy] = cutoff(running_entropy);
            
            % save the information in .mat file
            trialData.cutPoint = cutPoint;
            trialData.aveEntropy = ave_entropy;
            trialData.runningEntropy = running_entropy;
            
            % save the cutoff figure
                % save the figure for later check
            saveFigPath = [sub_dir, '\cut\'];
            if ~exist(saveFigPath)
                mkdir(saveFigPath)
            end
            dataIndex.Animal = {subjects_dir(ii).name};
            dataIndex.DateNumber = logfiles.datenum;
            plot_cutoff(trialData, saveFigPath, dataIndex);
            if trialData.cutPoint == 0
                cutTrial = [cutTrial,0];
            else
                cutTrial = [cutTrial, length(trialData.cue)-trialData.cutPoint];
            end
            
            c=double(trialData.response(1:trialData.cutPoint));
                for i=1:length(c)
                if c(i)~=0
                    c(i)=(c(i)-2.5)*2;
                end
            
                end
            r=trialData.outcome(1:trialData.cutPoint);
            com=trialData.rule(1:trialData.cutPoint);
            n_missed=sum(r(1:trialData.cutPoint)==77);
            r(~ismember(r, rewardCode))=0;
            r(ismember(r, rewardCode))=1;
            num_session1=[num_session1, length(c(c~=-1))];
            %concatenate
            c_nomiss1=[c_nomiss1;c(c~=0)]; r_nomiss1=[r_nomiss1;r(c~=0)]; com_nomiss1=[com_nomiss1;com(c~=0)];
        
        end
    
        % the following par includes algorithm 1 and 2
%         cd(base_dir2);
%         logfiles = dir('*.log');
%         c_nomiss2=[];r_nomiss2=[]; com_nomiss2=[];
%         [~,idx] = sort([logfiles.datenum]);
%         logfiles = logfiles(idx);
%         %get concatenated choice and reward history
%         for i =1:length(logfiles)
%             [ logData ] = parseLogfileHW(base_dir2, logfiles(i).name);
%         
%             [ sessionData, trialData] = MP_getSessionData( logData );
%             c=double(trialData.response);
%             for i=1:length(c)
%                 if c(i)~=0
%                     c(i)=(c(i)-2.5)*2;
%                 end
%             end
%             r=trialData.outcome;
%             com=trialData.comChoiceCode;
%             n_missed=sum(r==77);
%             r(~ismember(r, rewardCode))=0;
%             r(ismember(r, rewardCode))=1;
%             num_session2=[num_session2, length(c(c~=0))];
%             %concatenate
%             c_nomiss2=[c_nomiss2;c(c~=0)]; r_nomiss2=[r_nomiss2;r(c~=0)]; com_nomiss2=[com_nomiss2; com(c~=0)];
%         
%         end
    
    
        %exclude the missed trials
%         c_nomiss1=c_all1(c_all1~=0);
%         r_nomiss1=r_nomiss1(c_all1~=0);
%         com_nomiss1=com_nomiss1(c_all1~=0);
%     
%         c_nomiss2=c_nomiss2(c_nomiss2~=0);
%         r_nomiss2=r_nomiss2(c_nomiss2~=0);
%         com_nomiss2=com_all2(c_nomiss2~=0);
        
        choice{ii-2}.a1 = c_nomiss1; %choice{ii}.a2 = c_nomiss2;
        reward{ii-2}.a1 = r_nomiss1; %reward{ii}.a2 = r_nomiss2;
        comChoice{ii-2}.a1 = com_nomiss1; %comChoice{ii}.a2 = com_nomiss2;
        cutTrials{ii-2}.a1 = cutTrial;

    end
end
    %---------figure 0: number of trials per session
    
%     transit_session=length(num_session1)+0.5;
%     figure;
%     plot([num_session1,num_session2], 'k.', 'MarkerSize',35);
%     ylim([0 500]);
%     xlabel('Session');
%     ylabel('Number of trials');
%     title('Number of trials per session');
%     hold on;plot([transit_session transit_session],[0 500], 'k-','LineWidth',1);
%     
%     print(gcf, '-r0', [savebehfigpath sessionData.subject{1}, '_n_trials'], '-dpng'); %png format
%     saveas(gcf, [savebehfigpath sessionData.subject{1}, '_n_trials'], 'fig'); %fig format
%     


    %%--------figure 1: P_right per 200 trials  no missed trials
    
    
%% get some average
if ~exist(savebehfigpath)
    mkdir(savebehfigpath)
end
numSubjects = length(subjects_dir)-2;
figure;
%Transit_Point = zeros(1, 12);
for jj = 1:length(subjects_dir)-2
    block=200;
    len1=length(choice{jj}.a1);
    transit_point=200*(floor(len1/block)+0.5);
    
    pRight1=zeros(1, floor(len1/block));
    
    
    for i =1:length(pRight1)
        pRight1(i)=sum(choice{jj}.a1((i-1)*block+1:i*block)==1)/block;
    end
    
%     Transit_Point(jj-2) = length(pRight1) + 1;
%     
%     len2=length(choice{jj}.a2);
%     pRight2=zeros(1, floor(len2/block));
%     for i =1:length(pRight2)
%         pRight2(i)=sum(choice{jj}.a2((i-1)*block+1:i*block)==1)/block;
%     end
%     
%     
%     pRight=[pRight1, pRight2];
%     PRight{jj-2} = pRight;
    x_axis=[1:length(pRight1)]*block;
    sz=80;

    % set transit point to zero
    hold on;
    plot(x_axis,pRight1,'LineWidth',1);
    
    pRight_cell{jj} = pRight1;
end
hold on;plot([0 0],[0 1], 'k-','LineWidth',1);
hold on;plot([0 40000], [0.5,0.5],'k--','Linewidth',1);
ylim([0 1]);
xlabel('Trial');
ylabel('P(right)');
title('P(right)');
       print(gcf, '-r0', [savebehfigpath , 'p_right'], '-dpng'); %png format 
    saveas(gcf, [savebehfigpath , 'p_right'], 'fig'); %fig format
    
 % average pRight, align to the start
 len = cellfun('length', pRight_cell);
 maxLength = max(len);
 xAxis = [1:maxLength]*block;
 pRightMat = zeros(numSubjects, maxLength);
 for ii = 1:length(pRight_cell)
     pRightMat(ii,1:length(pRight_cell{ii})) = pRight_cell{ii};
 end
 
 % set empty entries to NaN
 pRightMat(pRightMat == 0) = NaN;
 
 % get mean and std
 meanpRight = nanmean(pRightMat,1);
 stdpRight = nanstd(pRightMat,1);
 figure;
s = shadedErrorBar(xAxis(1:maxLength),pRightMat(:,1:maxLength),{@nanmean,@nanstd},'lineprops','-black','patchSaturation',0.2);
set(s.edge,'LineWidth',1,'LineStyle',':')
 xlabel('block of 200 trials');
ylabel('probability to choose right');
title('Averaged probability to choose right');
print(gcf, '-r0', [savebehfigpath , 'p_right_ave_start'], '-dpng'); %png format
saveas(gcf, [savebehfigpath , 'p_right_ave_start'], 'fig'); %fig format
    
% align to the transition point (only get 5 blocks before and 5 blocks
% after
% numBlocks = 5;
% pRightMatTrans = zeros(numBlocks*2, numSubjects);
% for ii = 1:numSubjects
%     pRightMatTrans(:,ii) = PRight{ii}(Transit_Point(ii) - numBlocks: Transit_Point(ii)+numBlocks-1);
% end
% meanpRightTrans = mean(pRightMatTrans, 2);
% stdpRightTrans = std(pRightMatTrans,0,2);
% 
% figure; 
% errorbar(meanpRightTrans, stdpRightTrans, 'black', 'LineWidth',2);
% hold on; plot([5.5 5.5],[0 1],'k--', 'LineWidth',1);
% xlim([-1,11]);
% xlabel('block of 200 trials');
% ylabel('probability to choose right');
% title('Averaged probability to choose right');
% print(gcf, '-r0', [savebehfigpath , 'p_right_ave_transit'], '-dpng'); %png format
% saveas(gcf, [savebehfigpath , 'p_right_ave_transit'], 'fig'); %fig format
%     
%%---------figure 2: P_reward per 200 trials

figure;
for jj = 1:length(subjects_dir)-2
    block=200;
    len1=length(choice{jj}.a1);
    
    pReward1=zeros(1, floor(len1/block));
    for i =1:length(pReward1)
        pReward1(i)=sum(reward{jj}.a1((i-1)*block+1:i*block)==1)/block;
    end
    
%     pReward2=zeros(1, floor(len2/block));
%     for i =1:length(pReward2)
%         pReward2(i)=sum(reward{jj}.a2((i-1)*block+1:i*block)==1)/block;
%     end
    
%     pReward=[pReward1, pReward2];
    pReward_cell{jj} = pReward1;
    x_axis=[1:length(pReward1)]*block;
    hold on;
    plot(x_axis,pReward1,'LineWidth',1);
end
hold on;plot([0 0],[0 1], 'k-','LineWidth',1);
hold on;plot([0 40000], [0.5,0.5],'k--','Linewidth',1);
ylim([0 1]);
xlabel('Trial'); ylabel('P(reward)');
title('P(reward), n=5');
    %add the transit point
    
    print(gcf, '-r0', [savebehfigpath , 'p_reward'], '-dpng'); %png format
    saveas(gcf, [savebehfigpath , 'p_reward'], 'fig'); %fig format
 
    % average pRight, align to the start
 len = cellfun('length', pReward_cell);
 maxLength = max(len);
 pRewardMat = zeros(numSubjects, maxLength);
 for ii = 1:length(pRight_cell)
     pRewardMat(ii,1:length(pReward_cell{ii})) = pReward_cell{ii};
 end
 
 % set empty entries to NaN
 pRewardMat(pRewardMat == 0) = NaN;
 
 % get mean and std
meanpReward = nanmean(pRewardMat,1);
stdpReward = nanstd(pRewardMat,1);
figure;

s = shadedErrorBar(xAxis(1:maxLength),pRewardMat(:,1:maxLength),{@nanmean,@nanstd},'lineprops','-black','patchSaturation',0.2);
set(s.edge,'LineWidth',1,'LineStyle',':')
hold on;
hold on;plot([0 150], [0.5,0.5],'k--','Linewidth',1);
xlabel('block of 200 trials');
ylabel('probability of reward');
title('Averaged probability of reward, n=5');
print(gcf, '-r0', [savebehfigpath , 'p_reward_ave_start'], '-dpng'); %png format
saveas(gcf, [savebehfigpath , 'p_reward_ave_start'], 'fig'); %fig format
    
% align to the transition point (only get 5 blocks before and 5 blocks
% after
% numBlocks = 20;
% pRewardMatTrans = zeros(numBlocks*2, numSubjects);
% for ii = 1:numSubjects
%     pRewardMatTrans(:,ii) = PReward{ii}(Transit_Point(ii) - numBlocks: Transit_Point(ii)+numBlocks-1);
% end

% figure; 
% x_trans = -numBlocks:numBlocks-1;
% s = shadedErrorBar(x_trans,pRewardMatTrans',{@mean,@std},'lineprops','-black','patchSaturation',0.2);
% set(s.edge,'LineWidth',1,'LineStyle',':')
% hold on; plot([-0.5 -0.5],[0 1],'k--', 'LineWidth',1);
% hold on;plot([-numBlocks-1,numBlocks+1], [0.5,0.5],'k--','Linewidth',1);
% xlim([-numBlocks-1,numBlocks+1]);
% xlabel('Block of 200 trials from transition');
% ylabel('Probability of reward');
% title('Averaged probability of reward, n=12');
% print(gcf, '-r0', [savebehfigpath , 'p_reward_ave_transit'], '-dpng'); %png format
% saveas(gcf, [savebehfigpath , 'p_reward_ave_transit'], 'fig'); %fig format

%%------------figure 3: P_stay per block trials
figure;
for jj = 1:length(subjects_dir)-2
    block=200;
    len1=length(choice{jj}.a1);
    %transit_point=200*(floor(len1/block)+0.5);

    pStay1=zeros(1, floor(len1/block));
    for i =1:length(pStay1)
        sumStay=0;
        if i==1
            for k=2:block
                if choice{jj}.a1(k)==choice{jj}.a1(k-1)
                    sumStay=sumStay+1;
                end
            end
        else
            for k=(i-1)*block+1:i*block
                if choice{jj}.a1(k)==choice{jj}.a1(k-1)
                    sumStay=sumStay+1;
                end
            end
        end
        pStay1(i)=sumStay/block;
    end
    
%     pStay2=zeros(1, floor(len2/block));
%     for i =1:length(pStay2)
%         sumStay=0;
%         if i==1
%             for k=2:block
%                 if choice{jj}.a2(k)==choice{jj}.a2(k-1)
%                     sumStay=sumStay+1;
%                 end
%             end
%         else
%             for k=(i-1)*block+1:i*block
%                 if choice{jj}.a2(k)==choice{jj}.a2(k-1)
%                     sumStay=sumStay+1;
%                 end
%             end
%         end
%         pStay2(i)=sumStay/block;
%     end
    
    %pStay=[pStay1, pStay2];
    pStay_cell{jj} = pStay1;
    x_axis=[1:length(pStay1)]*block;
    hold on;
    plot(x_axis,pStay1,'LineWidth',1);
end
hold on;plot([0 0],[0 1], 'k-','LineWidth',1);
hold on;plot([0, 40000], [0.5,0.5],'k--','Linewidth',1);
ylim([0 1]);
xlabel('Trial'); ylabel('P(stay)');
title('P(stay)');
    print(gcf, '-r0', [savebehfigpath , 'p_stay'], '-dpng'); %png format
    saveas(gcf, [savebehfigpath , 'p_stay'], 'fig'); %fig format

    % average pRight, align to the start
 len = cellfun('length', pStay_cell);
 maxLength = max(len);
 pStayMat = zeros(numSubjects, maxLength);
 for ii = 1:numSubjects
     pStayMat(ii,1:length(pStay_cell{ii})) = pStay_cell{ii};
 end
 
 % set empty entries to NaN
 pStayMat(pStayMat == 0) = NaN;
 
 % get mean and std
meanpStay = nanmean(pStayMat,1);
stdpStay = nanstd(pStayMat,1);
figure;
xAxis = 1:maxLength;
s = shadedErrorBar(xAxis(1:maxLength),pStayMat(:,1:maxLength),{@nanmean,@nanstd},'lineprops','-black','patchSaturation',0.2);
set(s.edge,'LineWidth',1,'LineStyle',':')
hold on;plot([0 150], [0.5,0.5],'k--','Linewidth',1);
xlabel('Block of 200 trials');
ylabel('Probability to stay');
title('Averaged probability to stay, n=5');
print(gcf, '-r0', [savebehfigpath , 'p_stay_ave_start'], '-dpng'); %png format
saveas(gcf, [savebehfigpath , 'p_stay_ave_start'], 'fig'); %fig format
    
% align to the transition point (only get 5 blocks before and 5 blocks
% after
% numBlocks = 20;
% pStayMatTrans = zeros(numBlocks*2, numSubjects);
% for ii = 1:numSubjects
%     pStayMatTrans(:,ii) = PStay{ii}(Transit_Point(ii) - numBlocks: Transit_Point(ii)+numBlocks-1);
% end
% meanpStayTrans = mean(pStayMatTrans, 2);
% stdpStayTrans = std(pStayMatTrans,0,2);
% 
% figure; 
% x_trans = -numBlocks:numBlocks-1;
% s = shadedErrorBar(x_trans,pStayMatTrans',{@mean,@std},'lineprops','-black','patchSaturation',0.2);
% set(s.edge,'LineWidth',1,'LineStyle',':')
% hold on; plot([-0.5 -0.5],[0 1],'k--', 'LineWidth',1);
% hold on;plot([-numBlocks-1,numBlocks+1], [0.5,0.5],'k--','Linewidth',1);
% xlim([-numBlocks-1,numBlocks+1]);
% xlabel('Block of 200 trials');
% ylabel('Probability to stay');
% title('Averaged probability to stay, n=12');
% print(gcf, '-r0', [savebehfigpath , 'p_stay_ave_transit'], '-dpng'); %png format
% saveas(gcf, [savebehfigpath , 'p_stay_ave_transit'], 'fig'); %fig format
% 


    %%--------figure 4: P_wsls per block trials
figure;
for jj = 1:length(subjects_dir)-2
    block=200;
    len1=length(choice{jj}.a1);
    %len2=length(choice{jj}.a2);
    %transit_point=200*(floor(len1/block)+0.5);
    pWSLS1=zeros(1, floor(len1/block));
    for i =1:length(pWSLS1)
        sumWSLS=0;
        
        for k=(i-1)*block+1:i*block
            if reward{jj}.a1(k)==1
                if choice{jj}.a1(k+1)==choice{jj}.a1(k)
                    sumWSLS=sumWSLS+1;
                end
            else
                if choice{jj}.a1(k+1)~=choice{jj}.a1(k)
                    sumWSLS=sumWSLS+1;
                end
            end
        end
        
        pWSLS1(i)=sumWSLS/block;
    end
    
%     pWSLS2=zeros(1, floor(len2/block));
%     for i =1:length(pWSLS2)
%         sumWSLS=0;
%         
%         for k=(i-1)*block+1:i*block
%             if reward{jj}.a2(k)==1
%                 if choice{jj}.a2(k+1)==choice{jj}.a2(k)
%                     sumWSLS=sumWSLS+1;
%                 end
%             else
%                 if choice{jj}.a2(k+1)~=choice{jj}.a2(k)
%                     sumWSLS=sumWSLS+1;
%                 end
%             end
%         end
%         
%         pWSLS2(i)=sumWSLS/block;
%     end
    
    %pWSLS=[pWSLS1,pWSLS2];
    pWSLS_cell{jj} = pWSLS1;
    x_axis=[1:length(pWSLS1)]*block;
    hold on;
    
    plot(x_axis,pWSLS1,'LineWidth',1);
end

    hold on;plot([0 0],[0 1], 'k-','LineWidth',1);
    hold on;plot([0 40000], [0.5,0.5],'k--','Linewidth',1);
    ylim([0 1]);
    xlabel('Trial'); ylabel('P(WSLS)');
    title('P(WSLS)');
    print(gcf, '-r0', [savebehfigpath , 'p_WSLS'], '-dpng'); %png format
    saveas(gcf, [savebehfigpath , 'p_WSLS'], 'fig'); %fig format
 
        % average pRight, align to the start
 len = cellfun('length', pWSLS_cell);
 maxLength = max(len);
 pWSLSMat = zeros(numSubjects, maxLength);
 for ii = 1:length(pWSLS_cell)
     pWSLSMat(ii,1:length(pWSLS_cell{ii})) = pWSLS_cell{ii};
 end
 
 % set empty entries to NaN
 pWSLSMat(pWSLSMat == 0) = NaN;
 
 % get mean and std

figure;
xAxis = 1:maxLength;
s = shadedErrorBar(xAxis(1:maxLength),pStayMat(:,1:maxLength),{@nanmean,@nanstd},'lineprops','-black','patchSaturation',0.2);
set(s.edge,'LineWidth',1,'LineStyle',':')
hold on;plot([0 150], [0.5,0.5],'k--','Linewidth',1);
xlabel('Block of 200 trials');
ylabel('Probability to WSLS');
title('Averaged probability to WSLS, n=5');
print(gcf, '-r0', [savebehfigpath , 'p_WSLS_ave_start'], '-dpng'); %png format
saveas(gcf, [savebehfigpath , 'p_WSLS_ave_start'], 'fig'); %fig format
    
% align to the transition point (only get 5 blocks before and 5 blocks
% after
% numBlocks = 20;
% pWSLSMatTrans = zeros(numBlocks*2, numSubjects);
% for ii = 1:numSubjects
%     pWSLSMatTrans(:,ii) = pWSLS_cell{ii}(Transit_Point(ii) - numBlocks: Transit_Point(ii)+numBlocks-1);
% end
% 
% 
% figure; 
% x_trans = -numBlocks:numBlocks-1;
% s = shadedErrorBar(x_trans,pStayMatTrans',{@mean,@std},'lineprops','-black','patchSaturation',0.2);
% set(s.edge,'LineWidth',1,'LineStyle',':')
% hold on; plot([-0.5 -0.5],[0 1],'k--', 'LineWidth',1);
% hold on;plot([-numBlocks-1,numBlocks+1], [0.5,0.5],'k--','Linewidth',1);
% xlim([-numBlocks-1,numBlocks+1]);
% xlabel('Block of 200 trials');
% ylabel('Probability to WSLS');
% title('Averaged probability to WSLS, n=12');
% print(gcf, '-r0', [savebehfigpath , 'p_WSLS_ave_transit'], '-dpng'); %png format
% saveas(gcf, [savebehfigpath , 'p_WSLS_ave_transit'], 'fig'); %fig format
    
    %%-------figure 5: entropy (combinations of 3)
figure;
for jj = 1:length(subjects_dir)-2
    block=200;
    len1=length(choice{jj}.a1);
    %len2=length(choice{jj}.a2);
    %transit_point=200*(floor(len1/block)+0.5);
    
    enA13=zeros(1, floor(len1/block));
    %enA23=zeros(1, floor(len2/block));
    
    enA15=zeros(1, floor(len1/block));
    %enA25=zeros(1, floor(len2/block));
    
    %change left representation from  -1 to 0
    choice{jj}.a1(choice{jj}.a1==-1)=0;
    %choice{jj}.a2(choice{jj}.a2==-1)=0;
    choice{jj}.a1(choice{jj}.a1==-3)=0;  % there are several manual response
    %choice{jj}.a2(choice{jj}.a2==-3)=0;
    for i =1:length(enA13)
        entropyCount3=zeros(1,8);
        entropyCount5=zeros(1,32);
        for k =(i-1)*block+1:i*block
            ind3=bin2dec(num2str(choice{jj}.a1(k:k+2)'))+1;
            entropyCount3(ind3)=entropyCount3(ind3)+1;
            ind5=bin2dec(num2str([choice{jj}.a1(k:k+2); reward{jj}.a1(k:k+1)]'))+1;
            entropyCount5(ind5)=entropyCount5(ind5)+1;
        end
        %calculate the entropy
        p3=entropyCount3/sum(entropyCount3);
        if ismember(0,p3)
            warning('certain patterns did not occur, zero probability generated');
            IndexZero=find(~p3);
            for x=1:length(IndexZero)
                p3(IndexZero(x))=realmin;
            end
        end
        en1=-sum(p3.*log2(p3))+(length(p3)-1)/(1.3863*sum(entropyCount3));
        enA13(i)=en1;
        p5=entropyCount5/sum(entropyCount5);
        if ismember(0,p5)
            warning('certain patterns did not occur, zero probability generated');
            IndexZero=find(~p5);
            for x=1:length(IndexZero)
                p5(IndexZero(x))=realmin;
            end
        end
        en2=-sum(p5.*log2(p5))+(length(p5)-1)/(1.3863*sum(entropyCount5));
        enA15(i)=en2;
    end
    
%     for i =1:length(enA23)
%         entropyCount3=zeros(1,8);
%         entropyCount5=zeros(1,32);
%         for k =(i-1)*block+1:i*block
%             ind3=bin2dec(num2str(choice{jj}.a2(k:k+2)'))+1;
%             entropyCount3(ind3)=entropyCount3(ind3)+1;
%             ind5=bin2dec(num2str([choice{jj}.a2(k:k+2);reward{jj}.a2(k:k+1)]'))+1;
%             entropyCount5(ind5)=entropyCount5(ind5)+1;
%         end
%         %calculate the entropy
%         p3=entropyCount3/sum(entropyCount3);
%         if ismember(0,p3)
%             warning('certain patterns did not occur, zero probability generated');
%             IndexZero=find(~p3);
%             for x=1:length(IndexZero)
%                 p3(IndexZero(x))=realmin;
%             end
%         end
%         en1=-sum(p3.*log2(p3))+(length(p3)-1)/(1.3863*sum(entropyCount3));
%         enA23(i)=en1;
%         p5=entropyCount5/sum(entropyCount5);
%         if ismember(0,p5)
%             warning('certain patterns did not occur, zero probability generated');
%             IndexZero=find(~p5);
%             for x=1:length(IndexZero)
%                 p5(IndexZero(x))=realmin;
%             end
%         end
%         en2=-sum(p5.*log2(p5))+(length(p5)-1)/(1.3863*sum(entropyCount5));
%         enA25(i)=en2;
%     end
    
    %enA3=[enA13, enA23];
    %enA5=[enA15, enA25];
    ENA5{jj} = enA15;
    ENA3{jj} = enA13;
    x_axis=[1:length(enA13)]*block;
    hold on;
    
    subplot(2,1,1);
    hold on;
    plot(x_axis,enA13,'LineWidth',1);
   
    
    subplot(2,1,2); 
    hold on;
    plot(x_axis,enA15,'LineWidth',1);
end
subplot(2,1,1)
hold on;plot([0 0],[0 4], 'k-','LineWidth',1);
hold on;plot([0 40000], [3,3],'k--','Linewidth',1);
ylim([0 4]);
xlabel('Trial'); ylabel('Entropy (bits)');
title('Entropy of 3 consecutive choices');
subplot(2,1,2)
hold on;plot([0 0],[0 6], 'k-','LineWidth',1);
hold on;plot([0 40000], [5,5],'k--','Linewidth',1);
ylim([0 6]);
xlabel('Trial'); ylabel('Entropy (bits)');
title('Entropy of 3 choices and 2 rewards, n=12');
print(gcf, '-r0', [savebehfigpath , 'entropy_5'], '-dpng'); %png format
saveas(gcf, [savebehfigpath , 'entropy_5'], 'fig'); %fig format

 len = cellfun('length', ENA5);
 maxLength = max(len);

 ENA5Mat = zeros(numSubjects, maxLength);
 for ii = 1:numSubjects
     ENA5Mat(ii,1:length(ENA5{ii})) = ENA5{ii};
 end
 
 % set empty entries to NaN
 ENA5Mat(ENA5Mat == 0) = NaN;
 
 ENA3Mat = zeros(numSubjects, maxLength);
 for ii = 1:numSubjects
     ENA3Mat(ii,1:length(ENA3{ii})) = ENA3{ii};
 end
 
 % set empty entries to NaN
 ENA3Mat(ENA3Mat == 0) = NaN;
 % get mean and std

figure;
subplot(2,1,1)
xAxis = 1:maxLength;
s = shadedErrorBar(xAxis(1:maxLength),ENA5Mat(:,1:maxLength),{@nanmean,@nanstd},'lineprops','-black','patchSaturation',0.2);
set(s.edge,'LineWidth',1,'LineStyle',':')
hold on;plot([0 150], [5,5],'k--','Linewidth',1);
xlabel('Block of 200 trials');
ylabel('Entropy of 3 choices and 2 outcomes');
title('Averaged entropy of 3 choices and 2 outcomes, n=5');

subplot(2,1,2)
xAxis = 1:maxLength;
s = shadedErrorBar(xAxis(1:maxLength),ENA3Mat(:,1:maxLength),{@nanmean,@nanstd},'lineprops','-black','patchSaturation',0.2);
set(s.edge,'LineWidth',1,'LineStyle',':')
hold on;plot([0 150], [3,3],'k--','Linewidth',1);
xlabel('Block of 200 trials');
ylabel('Entropy of 3 choices');
title('Averaged entropy of 3 choices, n=5');
print(gcf, '-r0', [savebehfigpath , 'entropy5_ave_start'], '-dpng'); %png format
saveas(gcf, [savebehfigpath , 'entropy5_ave_start'], 'fig'); %fig format
    
%% plot the number of trials that was cut off

len = zeros(1,length(cutTrials));
for kk = 1:length(cutTrials)
    len(kk) = length(cutTrials{kk}.a1)
end
 maxLength = max(len);

 cutTrialsMat = NaN(numSubjects, maxLength);
 %cutTrialsMat = NaN;
 for ii = 1:numSubjects
     cutTrialsMat(ii,1:length(cutTrials{ii}.a1)) = cutTrials{ii}.a1;
 end
 
 % set empty entries to NaN
 figure;
xAxis = 1:maxLength;
s = shadedErrorBar(xAxis(1:maxLength),cutTrialsMat(:,1:maxLength),{@nanmean,@nanstd},'lineprops','-black','patchSaturation',0.2);
set(s.edge,'LineWidth',1,'LineStyle',':')
hold on;plot([0 50], [3,3],'k--','Linewidth',1);
xlabel('Session');
ylabel('Trials');
title('Averaged cutoff trials, n=5');
print(gcf, '-r0', [savebehfigpath , 'cutoff_trials'], '-dpng'); %png format
saveas(gcf, [savebehfigpath , 'cutoff_trials'], 'fig'); %fig format
% align to the transition point (only get 5 blocks before and 5 blocks
% after
% numBlocks = 20;
% ENA5MatTrans = zeros(numBlocks*2, numSubjects);
% for ii = 1:numSubjects
%     ENA5MatTrans(:,ii) = ENA5{ii}(Transit_Point(ii) - numBlocks: Transit_Point(ii)+numBlocks-1);
% end
% 
% 
% figure; 
% x_trans = -numBlocks:numBlocks-1;
% s = shadedErrorBar(x_trans,ENA5MatTrans',{@mean,@std},'lineprops','-black','patchSaturation',0.2);
% set(s.edge,'LineWidth',1,'LineStyle',':')
% hold on; plot([-0.5 -0.5],[0 6],'k--', 'LineWidth',1);
% hold on;plot([-numBlocks-1,numBlocks+1], [5,5],'k--','Linewidth',1);
% xlim([-numBlocks-1,numBlocks+1]);
% xlabel('Block of 200 trials');
% ylabel('Entropy of 3 choices and 2 rewards');
% title('Entropy of 3 choices and 2 rewards, n=12');
% print(gcf, '-r0', [savebehfigpath , 'ENA5_ave_transit'], '-dpng'); %png format
% saveas(gcf, [savebehfigpath , 'ENA5_ave_transit'], 'fig'); %fig format
% 
% % ena3
%  len = cellfun('length', ENA3);
%  maxLength = max(len);
%  numSubjects = 12;
%  ENA3Mat = zeros(numSubjects, maxLength);
%  for ii = 1:length(ENA3)
%      ENA3Mat(ii,1:length(ENA3{ii})) = ENA3{ii};
%  end
%  
%  % set empty entries to NaN
%  ENA3Mat(ENA3Mat == 0) = NaN;
%  
%  % get mean and std
% 
%  figure;
% xAxis = 1:maxLength;
% s = shadedErrorBar(xAxis(1:150),ENA3Mat(:,1:150),{@nanmean,@nanstd},'lineprops','-black','patchSaturation',0.2);
% set(s.edge,'LineWidth',1,'LineStyle',':')
% hold on;plot([0 150], [3,3],'k--','Linewidth',1);
% 
%  xlabel('Block of 200 trials');
% ylabel('Entropy of 3 choices');
% title('Averaged entropy of 3 choices,n=12');
% print(gcf, '-r0', [savebehfigpath , 'entropy3_ave_start'], '-dpng'); %png format
% saveas(gcf, [savebehfigpath , 'entropy3_ave_start'], 'fig'); %fig format
%     
% % align to the transition point (only get 5 blocks before and 5 blocks
% % after
% numBlocks = 20;
% ENA3MatTrans = zeros(numBlocks*2, numSubjects);
% for ii = 1:numSubjects
%     ENA3MatTrans(:,ii) = ENA3{ii}(Transit_Point(ii) - numBlocks: Transit_Point(ii)+numBlocks-1);
% end
% 
% figure; 
% x_trans = -numBlocks:numBlocks-1;
% s = shadedErrorBar(x_trans,ENA3MatTrans',{@mean,@std},'lineprops','-black','patchSaturation',0.2);
% set(s.edge,'LineWidth',1,'LineStyle',':')
% hold on; plot([-0.5 -0.5],[0 4],'k--', 'LineWidth',1);
% hold on;plot([-numBlocks-1,numBlocks+1], [3,3],'k--','Linewidth',1);
% xlim([-numBlocks-1,numBlocks+1]);
% xlabel('block of 200 trials');
% ylabel('Entropy of 3 choices');
% title('Entropy of 3 choices, n=12');
% print(gcf, '-r0', [savebehfigpath , 'ENA3_ave_transit'], '-dpng'); %png format
% saveas(gcf, [savebehfigpath , 'ENA3_ave_transit'], 'fig'); %fig format

close all
    %% -- find mutual information
    MI1=zeros(1, floor(len1/block));
    MI2=zeros(1, floor(len2/block));
    com_nomiss1(com_nomiss1==-1)=0;
    com_nomiss2(com_nomiss2==-1)=0;
    for i=1:length(MI1)
        choiceList=zeros(1,64);
        cList=zeros(1,2);
        C_RList=zeros(1,64,2);
        %for k =1:floor(block/3)
        for k=(i-1)*block+1:i*block
            %cPattern=[c_nomiss1(3*k+block*(i-1)-2:3*k+block*(i-1)); com_nomiss1(3*k+block*(i-1)-2:3*k+block*(i-1))];
            cPattern=[c_nomiss1(k:k+2);com_nomiss1(k:k+2)];
            ind=bin2dec(num2str(cPattern'))+1;
            choiceList(ind)=choiceList(ind)+1;
            if c_nomiss1(k+3)==1
                cList(1)=cList(1)+1;
                C_RList(1,ind, 1)=C_RList(1,ind,1)+1;
            else
                cList(2)=cList(2)+1;
                C_RList(1,ind, 2)=C_RList(1,ind,2)+1;
            end
        end
        %calculate the mutual information
        MI_sum=0;
        p1=cList(1)/sum(cList);
        p2=cList(2)/sum(cList);
        for h=1:64
            p_j=choiceList(h)/sum(choiceList);
            p_1j=C_RList(1,h, 1)/sum(sum(C_RList));
            p_2j=C_RList(1,h,2)/sum(sum(C_RList));
            if p_1j==0
                p_1j=realmin;
            end
            if p_2j==0
                p_2j=realmin;
            end
            if p_j==0
                p_j=realmin;
            end
            %         if p_1j==0
            %             warning('certain patterns did not occur, zero probability generated');
            %            p_whole1=realmin;
            %         end
            %         if p_2j==0
            %             warning('certain patterns did not occur, zero probability generated');
            %            p_whole2=realmin;
            %         end
            %         if (p_1j~=0 && p_j==0)
            %         warning('certain patterns did not occur, zero probability generated');
            %            p_j=realmin;
            %            p_whole1=p_1j/(p_j*p1);
            %         end
            %         if (p_2j~=0 && p_j==0)
            %         warning('certain patterns did not occur, zero probability generated');
            %            p_j=realmin;
            %            p_whole2=p_2j/(p_j*p2);
            %         end
            %         if (p_1j~=0 && p_j~=0)
            %             p_whole1=p_1j/(p_j*p1);
            %         end
            %         if (p_2j~=0 && p_j~=0)
            %             p_whole2=p_2j/(p_j*p2);
            %         end
            MI=p_1j*log2(p_1j/(p_j*p1))+p_2j*log2(p_2j/(p_j*p2));
            %MI=p_1j*log2(p_whole1)+p_2j*log2(p_whole2);
            MI_sum=MI_sum+MI;
        end
        MI_sum1=MI_sum%-63/(1.3863*block);
        MI1(i)=MI_sum1;
    end
    for i=1:length(MI2)
        choiceList=zeros(1,64);
        cList=zeros(1,2);
        C_RList=zeros(1,64,2);
        %for k =1:floor(block/3)
        for k=(i-1)*block+1:i*block
            %cPattern=[c_nomiss2(3*k+block*(i-1)-2:3*k+block*(i-1)); com_nomiss2(3*k+block*(i-1)-2:3*k+block*(i-1))];
            cPattern=[c_nomiss2(k:k+2);com_nomiss2(k:k+2)];
            ind=bin2dec(num2str(cPattern'))+1;
            choiceList(ind)=choiceList(ind)+1;
            if c_nomiss2(k+3)==1
                cList(1)=cList(1)+1;
                C_RList(1,ind, 1)=C_RList(1,ind,1)+1;
            else
                cList(2)=cList(2)+1;
                C_RList(1,ind, 2)=C_RList(1,ind,2)+1;
            end
        end
        %calculate the mutual information
        MI_sum=0;
        p1=cList(1)/sum(cList);
        p2=cList(2)/sum(cList);
        for h=1:64
            p_j=choiceList(h)/sum(choiceList);
            p_1j=C_RList(1,h, 1)/sum(sum(C_RList));
            p_2j=C_RList(1,h,2)/sum(sum(C_RList));
            if p_1j==0
                p_1j=realmin;
            end
            if p_2j==0
                p_2j=realmin;
            end
            if p_j==0
                p_j=realmin;
            end
            
            %         if p_1j==0
            %             warning('certain patterns did not occur, zero probability generated');
            %            p_whole1=realmin;
            %         end
            %         if p_2j==0
            %             warning('certain patterns did not occur, zero probability generated');
            %            p_whole2=realmin;
            %         end
            %         if (p_1j~=0 && p_j==0)
            %         warning('certain patterns did not occur, zero probability generated');
            %            p_j=realmin;
            %            p_whole1=p_1j/(p_j*p1);
            %         end
            %         if (p_2j~=0 && p_j==0)
            %         warning('certain patterns did not occur, zero probability generated');
            %            p_j=realmin;
            %            p_whole2=p_2j/(p_j*p2);
            %         end
            %         if (p_1j~=0 && p_j~=0)
            %             p_whole1=p_1j/(p_j*p1);
            %         end
            %         if (p_2j~=0 && p_j~=0)
            %             p_whole2=p_2j/(p_j*p2);
            %         end
            %MI=p_1j*log2(p_whole1)+p_2j*log2(p_whole2);
            MI=p_1j*log2(p_1j/(p_j*p1))+p_2j*log2(p_2j/(p_j*p2));
            MI_sum=MI_sum+MI;
            
        end
        MI_sum1=MI_sum%-63/(1.3863*block);
        MI2(i)=MI_sum1;
    end
    
    MItotal=[MI1,MI2];
    
    figure;
    scatter(x_axis,MItotal,sz,'k','filled');
    hold on;plot([0 x_axis(end)], [0,0],'k--','Linewidth',1);
    ylim([-0.2 0.8]);
    hold on;plot([transit_point transit_point],[-0.2 0.8], 'k-','LineWidth',1);
    %hold on;plot([0 x_axis(end)], [3,3],'k--','Linewidth',1);
    %ylim([0 4]);
    xlabel('Trial'); ylabel('Mutual information (bits)');
    title('Mutual information');
    print(gcf, '-r0', [savebehfigpath sessionData.subject{1}, '_mutualinfo'], '-dpng'); %png format
    saveas(gcf, [savebehfigpath sessionData.subject{1}, '_mutualinfo'], 'fig'); %fig format
    
    %significance test
    %probability right
    [h1,p1]=ttest2(pRight1,pRight2);
    [h2,p2]=ttest2(pReward1,pReward2);
    [h3,p3]=ttest2(pStay1,pStay2);
    [h4,p4]=ttest2(pWSLS1,pWSLS2);
    [h5,p5]=ttest2(enA13,enA23);
    [h6,p6]=ttest2(enA15,enA25);