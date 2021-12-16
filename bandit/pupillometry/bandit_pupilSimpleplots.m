function bandit_pupilSimpleplots(dataIndex)

% simple plots

nFiles = size(dataIndex,1);

for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},['bandit',dataIndex.LogFileName{ii}(end-30:end-4),'_beh.mat']));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    % load pupil files
    date = num2str(dataIndex.DateNumber(ii));
    pup_name = fullfile(dataIndex.BehPath{ii}, ['*',date(1:6),'*_pup.mat']);
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
%         tPupil = pupil.t(pupil.t<trialData.cueTimes(2) & pupil.t>trialData.cueTimes(1)-1);
%         rawPupil = pupil.dia(pupil.t<trialData.cueTimes(2) & pupil.t>trialData.cueTimes(1)-1);
%         cueInd = find(tPupil<trialData.cueTimes(1),1,'last');
%         responseInd = find(tPupil<trialData.rt(1),1,'last');
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
%         % Capture the plot as an image 
%             frame = getframe(h); 
%             im = frame2im(frame); 
%             [imind,cm] = rgb2ind(im,256); 
%             % Write to the GIF File 
%             if n == 1 
%                 imwrite(imind,cm,filename,'gif','DelayTime',0.05, 'Loopcount',1); 
%             else 
%                 imwrite(imind,cm,filename,'gif','DelayTime',0.05,'WriteMode','append'); 
%             end 
%         end
%         
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
        
        %% plot pupil change against cue
        % Fig. 3d in paper came from 140605 data set, cell 8 10 37 74
        params=[];
        params.trigTime = trialData.cueTimes;
        params.xtitle = 'Time from cue (s)';
        params.window = [-3:0.1:8];
        params.minNumTrial = 5;
        params.numBootstrapRepeat = 1000;   %number of repeats for bootstrap (for estimating CI)
        params.CI = 0.95;  %confidence interval
        psth_output=[];
        
        for k=1:2
            fieldname=[];
            if k==1     %panel 1
                fieldname{1}={'left','reward'}; col{1} = 'r';
                fieldname{2}={'right','reward'}; col{2} = 'b';
            elseif k==2 %panel 2
                fieldname{1}={'left','noreward'}; col{1} = 'r';
                fieldname{2}={'right','noreward'}; col{2} = 'b';
            end
            for kk=1:numel(fieldname)
                trialMask = getMask(trials,fieldname{kk});
                % find the shortest timeline
                
                

                trialMask = trialMask(trialData.cueTimes<(pupil.t(end)-8));
                
                psth_panel(k).sig{kk} = MP_get_psth(pupil.dia, pupil.t, params.trigTime(trialMask), strjoin(fieldname{kk}), params );
                psth_panel(k).col{kk} = col{kk};
            end
        end
        
        
        
        tlabel = 'pupil';
        
        MP_plot_psth(psth_panel,tlabel,params.xtitle);
        print(gcf,'-dpng',['cell_choice' int2str(j)]);
        
        % plot reward aligned to rewardtime
        params=[];
        params.trigTime = trialData.outcomeTimes;
        params.xtitle = 'Time from outcome (s)';
        params.window = [-3:0.1:5];
        params.minNumTrial = 5;
        params.numBootstrapRepeat = 1000;   %number of repeats for bootstrap (for estimating CI)
        params.CI = 0.95;  %confidence interval
        psth_output=[];
        for k=1:2
        if k==1 % panel 3
                fieldname{1}={'left','reward'}; col{1} = 'r';
                fieldname{2}={'left','noreward'}; col{2} = 'b';
            elseif k==2
                fieldname{1}={'right','reward'}; col{1} = 'r';
                fieldname{2}={'right','noreward'}; col{2} = 'b';
            end
            for kk=1:numel(fieldname)
                trialMask = getMask(trials,fieldname{kk});
                trialMask = trialMask(trialData.cueTimes<(pupil.t(end)-8));
                psth_panel(k).sig{kk} = MP_get_psth(pupil.dia, pupil.t, params.trigTime(trialMask), strjoin(fieldname{kk}), params );
                psth_panel(k).col{kk} = col{kk};
            end
        end
        
        
        
        tlabel = 'pupil';
        
        MP_plot_psth(psth_panel,tlabel,params.xtitle);
        print(gcf,'-dpng',['cell_reward' int2str(j)]);
        
        % plot a figure to show the extremties influence the error bar
        % (bootstrap)
        % MP_plot_extreme(psth_panel, tlabel, params.xtitle);
        close all;
    end

end

end