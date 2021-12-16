%% Make simple plots of pupillometry data (no need to save anything)
% run this after running start_beh and start_dff_compute

clearvars;
close all;

tic;    %set clock to estimate how long this takes

%% setup path and plotting formats

value_setPathList;

setup_figprop;  %set up default figure plotting parameters

%% which set of data to load?
m=7;

switch m
    case 1  % Load the data set from the NN paper - M2 imaging
        data_subdir = fullfile(data_dir,'M2 NN');
        [ dirs, expData ] = expData_M2_NN(data_subdir);
        interpulse_present=false;

    case 4  % Load the data set from the M2 SST data set
        data_subdir = fullfile(data_dir,'M2 SST');
        [ dirs, expData ] = expData_M2_SST(data_subdir);
        interpulse_present=false;
        
    case 7 % load the data set from pupillometry
        data_subdir = fullfile(data_dir);
        [ dirs, expData ] = expData_pupil(data_subdir);
        
end

%% process data files
for i = 1:numel(expData)
    disp(['Processing file ' int2str(i) ' out of ' int2str(numel(expData)) '.']);

    % setup/create subdirectories to save analysis and figures
    savematpath = fullfile(dirs.analysis,expData(i).sub_dir);
    if ~exist(savematpath,'dir')
        mkdir(savematpath);
    end
    savefluofigpath = fullfile(dirs.analysis,expData(i).sub_dir,'figs-fluo');
    if ~exist(savefluofigpath,'dir')
        mkdir(savefluofigpath);
    end

    % load the saved behavioral analysis (from start_beh.m)
    % load the saved dF/F (from start_dff_compute.m)
    cd(fullfile(dirs.analysis,expData(i).sub_dir));
    load('beh.mat');
    load('pupZ.mat');

    cd(savefluofigpath);

    %% Plot pupil for all trials
    MP_plot_pupil( pupil, trialData );
    
    %% Plot cue-aligned pupil
    % Fig. 3d in paper came from 140605 data set, cell 8 10 37 74
    params=[];
    params.trigTime = trialData.cueTimes;
    params.xtitle = 'Time from cue (s)';
    params.window = [-3:0.1:8];
    params.numBootstrapRepeat = 1000;   %number of repeats for bootstrap (for estimating CI)
    params.CI = 0.95;  %confidence interval
    psth_output=[];
    
    for k=1:2
        fieldname=[];
        if k==1     %panel 1
            fieldname{1}={'left','reward'};
            fieldname{2}={'right','reward'};
        elseif k==2 %panel 2
            fieldname{1}={'left','noreward'};
            fieldname{2}={'right','noreward'};
        end
        for kk=1:numel(fieldname)
            trialMask = getMask(trials,fieldname{kk});
            psth_panel(k).sig{kk} = MP_get_psth(pupil.dia, pupil.t, params.trigTime(trialMask), strjoin(fieldname{kk}), params );
        end
    end
    
 

    tlabel = 'pupil';
        
    MP_plot_psth(psth_panel,tlabel,params.xtitle);
    print(gcf,'-dpng',['cell' int2str(j)]);
    
    % plot a figure to show the extremties influence the error bar
    % (bootstrap)
    % MP_plot_extreme(psth_panel, tlabel, params.xtitle);
    close;

    %%
    clearvars -except i dirs expData;
    close all;
end

% plays sound when done
load train;
sound(y,Fs);

toc