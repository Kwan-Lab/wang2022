function MP_pupilSimpleplots(dataIndex)

% simple plots

nFiles = size(dataIndex,1);

for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},[dataIndex.LogFileName{ii}(1:end-4),'_beh_cut.mat']));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    % load pupil files
    date = num2str(dataIndex.DateNumber(ii));
    pup_name = fullfile(dataIndex.BehPath{ii},['*',date(1:6),'*_pup.mat']);
    fn_pup = dir(pup_name);
    if length(fn_pup) == 1
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
        
        %% Plot pupil for all trials
        
        % MP_plot_pupil( pupil, trialData );
        
        %% Plot cue-aligned pupil
        
        % plot a single trace, create a gif
        
       
        % tCue = trialData.cueTimes(1):trialData.cueTimes(2);
        
        %% an animation of one example trial
%         tPupil = pupil.t(pupil.t<trialData.cueTimes(4) & pupil.t>trialData.cueTimes(3)-2);
%         rawPupil = pupil.dia(pupil.t<trialData.cueTimes(4) & pupil.t>trialData.cueTimes(3)-2);
%         cueInd = find(tPupil<trialData.cueTimes(2),1,'last');
%         responseInd = find(tPupil<trialData.rt(2),1,'last');
%         outcomeInd = find(tPupil<trialData.outcomeTimes(1),1,'last');
%         h=figure;
%         filename = 'gif_animaNoreward.gif';
%         for n = 1:length(tPupil)
%             yA = NaN(1, length(tPupil));
%             yA(1:n) = rawPupil(1:n);
%             plot(tPupil,yA,'black')
%             ylabel('z-score');
%             ylim([-1.5,2.5]);
%             xlim([-1,6]);
%             set(gca,'XColor', 'none', 'YColor','none');
%             
%             if n==cueInd
%                 hold on;plot([trialData.cueTimes(1),trialData.cueTimes(1)],[-1, 6]);
%             elseif n==responseInd
%                 hold on; plot([trialData.rt(1),trialData.rt(1)],[-1,6]);
%             end
%         Capture the plot as an image 
%             frame = getframe(h); 
%             im = frame2im(frame); 
%             [imind,cm] = rgb2ind(im,256); 
%             Write to the GIF File 
%             if n == 1 
%                 imwrite(imind,cm,filename,'gif','DelayTime',0.05, 'Loopcount',1); 
%             else 
%                 imwrite(imind,cm,filename,'gif','DelayTime',0.05,'WriteMode','append'); 
%             end 
%         end
%         
%         figure;
%         plot(tPupil-trialData.cueTimes(3),rawPupil,'black');
%         
%         plot the arrow indicating cue and outcome
%         
%         
%         x1 = trialData.cueTimes(3);[minValue,closestIndex] = min(abs(x1-tPupil));
%         y1 = rawPupil(closestIndex)-0.2;
%         y12 = rawPupil(closestIndex)-1.2;
%         ylim([-3,3]);
%         xlim([tPupil(1)-trialData.cueTimes(3),tPupil(end)-trialData.cueTimes(3)]);
%         hold on;
%         plot cue time
%         plot([x1,x1],[-3,3],'k--');
%         plot scale bar
%         hold on;
%         plot([4,4],[3,2],'k');
%         q1 = quiver(x1,y12,x1-x1,y1-y12,0);
%         q1.Color = 'black';
%         q1.LineWidth = 2;
%         x2 = trialData.rt(1);
%         [minValue,closestIndex] = min(abs(x2-tPupil));
%         y2 = rawPupil(closestIndex)-0.2;
%         y22 = y12;
%         hold on;
%         q2 = quiver(x2,y22,x1-x1,y2-y22,0);
%         q2.Color = 'black';
%         q2.LineWidth = 2;
%         xlabel('Time (s)','FontSize',38);
%         set(gca,'YTick',[-2:2:2]);
%         ylabel('Pupil diameter (z)','FontSize',38);
%         set(gca,'box','off');
%         set(gca,'LineWidth',4);
%         set(gca,'XColor', 'none', 'YColor','none');
%         print(gcf,'-dpng',['sample-trace1']);    %png format
%         saveas(gcf,['sample-trace1'], 'fig');
%         saveas(gcf, 'sample-trace1','svg');
% 
%         figure;
%         numTrial = 80;
%         tPupil = pupil.t(pupil.t<trialData.cueTimes(numTrial+1) & pupil.t>trialData.cueTimes(numTrial)-1);
%         rawPupil = pupil.dia(pupil.t<trialData.cueTimes(numTrial+1) & pupil.t>trialData.cueTimes(numTrial)-1);
%         plot(tPupil,rawPupil,'black');
%         
%         plot the arrow indicating cue and outcome
%         
%        
%         x1 = trialData.cueTimes(numTrial);
%         ylim([-3,3]);
%         xlim([tPupil(1),tPupil(end)]);
%         hold on;
%         plot cue time
%         plot([x1,x1],[-3,3],'k--');
%         plot scale bar
%         hold on;
%         plot([x1+4,x1+4],[3,2],'k');
%         
%         xlabel('Time (s)','FontSize',38);
%         ylabel('Pupil diameter (z)','FontSize',38);
%         set(gca,'box','off');
%         set(gca,'LineWidth',4); set(gca,'YTick',[-2:2:2]);
%         set(gca,'XColor', 'none', 'YColor','none');
%         print(gcf,'-dpng',['sample-trace2']);    %png format
%         saveas(gcf,['sample-trace2'], 'fig');
%         saveas(gcf, 'sample-trace2','svg');
%         
%         figure;
%         numTrial = 200;
%         tPupil = pupil.t(pupil.t<trialData.cueTimes(numTrial+1) & pupil.t>trialData.cueTimes(numTrial)-1);
%         rawPupil = pupil.dia(pupil.t<trialData.cueTimes(numTrial+1) & pupil.t>trialData.cueTimes(numTrial)-1);
%         plot(tPupil,rawPupil,'black');
%         
%         plot the arrow indicating cue and outcome
%         
%        
%         x1 = trialData.cueTimes(numTrial);
%         ylim([-3,3]);
%         hold on;
%         plot cue time
%         plot([x1,x1],[-3,3],'k--');
%         plot scale bar
%         hold on;
%         plot([x1+4,x1+4],[3,2],'k');
%        
%         xlabel('Time (s)','FontSize',38);
%         ylabel('Pupil diameter (z)','FontSize',38);
%         set(gca,'box','off');
%         set(gca,'LineWidth',4); set(gca,'YTick',[-2:2:2]);
%         set(gca,'XColor', 'none', 'YColor','none');
%         print(gcf,'-dpng',['sample-trace3']);    %png format
%         saveas(gcf,['sample-trace3'], 'fig');
%         saveas(gcf, 'sample-trace3','svg');
%         
        % plot the whole session
        figure;

        plot(pupil.t,pupil.dia,'black','LineWidth',1);
        set(gca,'box','off');
        %set(gca,'XColor', 'none', 'YColor','none');
        hold on;
        xlim([-100 5000]);
        %plot([5000,5000],[3,2],'k');
        xlabel('Time (s)');
        ylabel('Pupil diameter (z)');
        print(gcf,'-dpng',['whole-trace']);    %png format
        saveas(gcf,['whole-trace'], 'fig');
        saveas(gcf, 'whole-trace','svg');

        %% plot trial averaged pupil trace
        
        % left reward trials
        LRCue = trialData.cueTimes(trials.left == 1 & trials.reward == 1);
        window = [-3:0.05:5-0.05];
        LRTrace = zeros(length(window), length(LRCue));
        for ii = 1:length(LRCue)
            t1 = LRCue(ii)-3;
            t2 = LRCue(ii) + 5;
            if t2<pupil.t(end)
                LRTrace(:,ii) = pupil.dia(pupil.t>= t1 & pupil.t<= t2)';
            end
        end
        gray=[0.7 0.7 0.7];
        figure;plot(window,nanmean(LRTrace,2),'k');
        hold on;
        errorshade(window,nanmean(LRTrace,2)'+nanstd(LRTrace,0,2)',nanmean(LRTrace,2)'-nanstd(LRTrace,0,2)',gray);
        hold on;
        plot(window,nanmean(LRTrace,2),'k');
        hold on;
        plot([0 0], [-2 2],'k--');
        set(gca, 'box','off');
        title('Left reward trials');
        ylabel('Z-score');
        xlabel('Time from cue (s)');
         print(gcf,'-dpng',['left reward']);    %png format
        saveas(gcf,['left reward'], 'fig');
        saveas(gcf, 'left reward','svg');
        
          RRCue = trialData.cueTimes(trials.right == 1 & trials.reward == 1);
        window = [-3:0.05:5-0.05];
        LRTrace = zeros(length(window), length(RRCue));
        for ii = 1:length(RRCue)
            t1 = RRCue(ii)-3;
            t2 = RRCue(ii) + 5;
            if t2<pupil.t(end)
                RRTrace(:,ii) = pupil.dia(pupil.t>= t1 & pupil.t<= t2)';
            end
        end
        gray=[0.7 0.7 0.7];
        figure;plot(window,nanmean(RRTrace,2),'k');
        hold on;
        errorshade(window,nanmean(RRTrace,2)'+nanstd(RRTrace,0,2)',nanmean(RRTrace,2)'-nanstd(RRTrace,0,2)',gray);
        hold on;
        plot(window,nanmean(RRTrace,2),'k');
        hold on;
        plot([0 0], [-2 2],'k--');
        set(gca, 'box','off');
        title('Right reward trials');
        ylabel('Z-score');
        xlabel('Time from cue (s)');
         print(gcf,'-dpng',['right reward']);    %png format
        saveas(gcf,['right reward'], 'fig');
        saveas(gcf, 'right reward','svg');
        
          LNRCue = trialData.cueTimes(trials.left == 1 & trials.noreward == 1);
        window = [-3:0.05:5-0.05];
        LNRTrace = zeros(length(window), length(LNRCue));
        for ii = 1:length(LNRCue)
            t1 = LNRCue(ii)-3;
            t2 = LNRCue(ii) + 5;
            if t2<pupil.t(end)
                LNRTrace(:,ii) = pupil.dia(pupil.t>= t1 & pupil.t<= t2)';
            end
        end
        gray=[0.7 0.7 0.7];
        figure;plot(window,nanmean(LNRTrace,2),'k');
        hold on;
        errorshade(window,nanmean(LNRTrace,2)'+nanstd(LNRTrace,0,2)',nanmean(LNRTrace,2)'-nanstd(LNRTrace,0,2)',gray);
        hold on;
        plot(window,nanmean(LNRTrace,2),'k');
        hold on;
        plot([0 0], [-2 2],'k--');
        set(gca, 'box','off');
        title('Left noreward trials');
        ylabel('Z-score');
        xlabel('Time from cue (s)');
         print(gcf,'-dpng',['left noreward']);    %png format
        saveas(gcf,['left noreward'], 'fig');
        saveas(gcf, 'left noreward','svg');
        
         RNRCue = trialData.cueTimes(trials.right == 1 & trials.noreward == 1);
        window = [-3:0.05:5-0.05];
        RNRTrace = zeros(length(window), length(RNRCue));
        for ii = 1:length(RNRCue)
            t1 = RNRCue(ii)-3;
            t2 = RNRCue(ii) + 5;
            if t2<pupil.t(end)
                RNRTrace(:,ii) = pupil.dia(pupil.t>= t1 & pupil.t<= t2)';
            end
        end
        gray=[0.7 0.7 0.7];
        figure;plot(window,nanmean(RNRTrace,2),'k');
        hold on;
        errorshade(window,nanmean(RNRTrace,2)'+nanstd(RNRTrace,0,2)',nanmean(RNRTrace,2)'-nanstd(RNRTrace,0,2)',gray);
        hold on;
        plot(window,nanmean(RNRTrace,2),'k');
        hold on;
        plot([0 0], [-2 2],'k--');
        set(gca, 'box','off');
        title('Right noreward trials');
        ylabel('Z-score');
        xlabel('Time from cue (s)');
         print(gcf,'-dpng',['right noreward']);    %png format
        saveas(gcf,['right noreward'], 'fig');
        saveas(gcf, 'right noreward','svg');
        % left-reward; right-reward; left-noreward; right-noreward
%         tPupil = pupil.t(pupil.t<trialData.cueTimes(4) & pupil.t>trialData.cueTimes(3)-1);
%         rawPupil = pupil.dia(pupil.t<trialData.cueTimes(4) & pupil.t>trialData.cueTimes(3)-1);
%         cueInd = find(tPupil<trialData.cueTimes(3),1,'last');
%         responseInd = find(tPupil<(trialData.rt(3)+trialData.cueTimes(3)),1,'last');
%         outcomeInd = find(tPupil<trialData.outcomeTimes(3),1,'last');
%         h=figure;
%         filename = 'gif_animareward.gif';
%         for n = 1:length(tPupil)
%             yA = NaN(1, length(tPupil));
%             yA(1:n) = rawPupil(1:n);
%             plot(tPupil,yA,'black')
%             ylabel('z-score');
%             ylim([-1.5,2.5]);
%             xlim([15,24]);
%            % set(gca,'XColor', 'none', 'YColor','none');
%             
%             if n==cueInd
%                 hold on;plot([trialData.cueTimes(3),trialData.cueTimes(3)],[-2, 6]);
%             elseif n==outcomeInd
%                 hold on; plot([trialData.outcomeTimes(3),trialData.outcomeTimes(3)],[-2,6]);
%             end
%         % Capture the plot as an image 
%             frame = getframe(h); 
%             im = frame2im(frame); 
%             [imind,cm] = rgb2ind(im,256); 
%             % Write to the GIF File 
%             if n == 1 
%                 imwrite(imind,cm,filename,'gif','DelayTime',0.05, 'Loopcount',inf); 
%             else 
%                 imwrite(imind,cm,filename,'gif','DelayTime',0.05,'WriteMode','append'); 
%             end 
%         end
        
        % plot the curve
        
        %% plot pupil change against cue (psth plot)
       
%         params=[];
%         params.trigTime = trialData.cueTimes;
%         params.xtitle = 'Time from cue (s)';
%         params.window = [-3:0.1:8];
%         params.minNumTrial = 5;
%         params.numBootstrapRepeat = 1000;   %number of repeats for bootstrap (for estimating CI)
%         params.CI = 0.95;  %confidence interval
%         psth_output=[];
%         
%         for k=1:2
%             fieldname=[];
%             if k==1     %panel 1
%                 fieldname{1}={'left','reward'}; col{1} = 'r';
%                 fieldname{2}={'right','reward'}; col{2} = 'b';
%             elseif k==2 %panel 2
%                 fieldname{1}={'left','noreward'}; col{1} = 'r';
%                 fieldname{2}={'right','noreward'}; col{2} = 'b';
%             end
%             for kk=1:numel(fieldname)
%                 trialMask = getMask(trials,fieldname{kk});
%                 % find the shortest timeline
%                 
%                 
% 
%                 trialMask = trialMask(trialData.cueTimes<(pupil.t(end)-8));
%                 
%                 psth_panel(k).sig{kk} = MP_get_psth(pupil.dia, pupil.t, params.trigTime(trialMask), strjoin(fieldname{kk}), params );
%                 psth_panel(k).col{kk} = col{kk};
%             end
%         end
%         
%         
%         
%         tlabel = 'pupil';
%         
%         MP_plot_psth(psth_panel,tlabel,params.xtitle);
%         print(gcf,'-dpng',['cell_choice' int2str(j)]);
%         
%         % plot reward aligned to rewardtime
%         params=[];
%         params.trigTime = trialData.outcomeTimes;
%         params.xtitle = 'Time from outcome (s)';
%         params.window = [-3:0.1:5];
%         params.minNumTrial = 5;
%         params.numBootstrapRepeat = 1000;   %number of repeats for bootstrap (for estimating CI)
%         params.CI = 0.95;  %confidence interval
%         psth_output=[];
%         for k=1:2
%         if k==1 % panel 3
%                 fieldname{1}={'left','reward'}; col{1} = 'r';
%                 fieldname{2}={'left','noreward'}; col{2} = 'b';
%             elseif k==2
%                 fieldname{1}={'right','reward'}; col{1} = 'r';
%                 fieldname{2}={'right','noreward'}; col{2} = 'b';
%             end
%             for kk=1:numel(fieldname)
%                 trialMask = getMask(trials,fieldname{kk});
%                 trialMask = trialMask(trialData.cueTimes<(pupil.t(end)-8));
%                 psth_panel(k).sig{kk} = MP_get_psth(pupil.dia, pupil.t, params.trigTime(trialMask), strjoin(fieldname{kk}), params );
%                 psth_panel(k).col{kk} = col{kk};
%             end
%         end
%         
%         
%         
%         tlabel = 'pupil';
%         
%         MP_plot_psth(psth_panel,tlabel,params.xtitle);
%         print(gcf,'-dpng',['cell_reward' int2str(j)]);
        
        % plot a figure to show the extremties influence the error bar
        % (bootstrap)
        % MP_plot_extreme(psth_panel, tlabel, params.xtitle);
        close all;
    end

end

end