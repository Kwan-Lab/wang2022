function  pennies_pupilSimpleplots_perAnimal(dataIndex, modelpath, savefigpath, animalInd);

nFiles = size(dataIndex,1);

currAnimal = dataIndex.Animal{1};
% subMask = [];
savepath = fullfile(savefigpath,[currAnimal,'_fig-pupil']);

for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},[dataIndex.LogFileName{ii}(1:end-4),'_beh.mat']));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    % load pupil files
    date = num2str(dataIndex.DateNumber(ii));
    pup_name = [dataIndex.BehPath{ii},'/', '*',date(1:6),'*_pup.mat'];
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