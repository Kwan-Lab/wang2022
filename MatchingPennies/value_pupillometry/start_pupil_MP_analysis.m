%% Analyze fluorescence data, for discrim tasks with no flexibility
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
    case 5  % Load the data set from the M2 Discrim data set
        data_subdir = fullfile(data_dir,'');
        [ dirs, expData ] = expData_M2_Discrim(data_subdir);
        interpulse_present=false;
    case 6  % Load the data set from the M2 Discrim data set
        data_subdir = fullfile(data_dir,'M2 omission');
        [ dirs, expData ] = expData_M2_Omission(data_subdir);
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

    %% Choice selectivity during sound trials
%     params=[];
%     params.trigTime = trialData.cueTimes;
%     params.xtitle = 'Time from cue (s)';
%     params.window = [-3:0.1:8];
% 
%     psth_soundleft=[]; psth_soundright=[];
%         fieldname={'sound','upsweep','left','hit'};
%         trialMask = getMask(trials,fieldname);
%         psth_soundleft{j} = MP_get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname),params);
%         
%         fieldname={'sound','downsweep','right','hit'};
%         trialMask = getMask(trials,fieldname);
%         psth_soundright{j} = MP_get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname),params);
%     end
%     
%     tlabel = 'Choice selectivity';
%     plot_selectivity(psth_soundleft,psth_soundright,tlabel,params.xtitle);
%     print(gcf,'-dpng','choice-selectivity');    %png format
%     saveas(gcf, 'choice-selectivity', 'fig');
%     
    %% Multiple linear regression  - choice and reward for sound trials
    % Fig. 8a
    params=[];
    
    %first predictor is choice; dummy-code: left=-1, right=1, miss=NaN
    params.trigEvent=NaN(size(trials.left)); 
    params.trigEvent(trials.left) = -1;  
    params.trigEvent(trials.right) = 1;
    %second predictor is outcome; dummy-code: reward=1, error=0, miss=NaN
    params.trigEvent2=NaN(size(trials.go)); 
    params.trigEvent2(trials.reward) = 1;  
    params.trigEvent2(trials.noreward) = 0;
    % params.trigEvent2(trials.omitreward) = 0;
    % params.trigEvent2(trials.doublereward) = 2;
    
    params.trigTime = trialData.cueTimes;
    params.xtitle = 'Time from cue (s)';
    params.window = [-3:0.1:8];
    params.nback = 3;       %how many trials back to regress against
    params.pvalThresh = 0.01;   %p-value for coefficient be considered significant

    %only perform analysis on this subset of trials
    fieldname={'go'};
    trialMask = getMask(trials,fieldname);
    reg_cr=linear_regr( pupil.dia, pupil.t, [params.trigEvent params.trigEvent2], params.trigTime, trialMask, params );

    tlabel={'C(n)','C(n-1)','C(n-2)','C(n-3)','R(n)','R(n-1)','R(n-2)','R(n-3)'};
    plot_regr_pupil(reg_cr,params.pvalThresh,tlabel,params.xtitle);
    print(gcf,'-dpng','MLR-choiceoutcome');    %png format
    saveas(gcf, 'MLR-choiceoutcome', 'fig');

    
 
    
    %% save the dF/F so don't have to re-compute every time
    save(fullfile(savematpath,'dff_and_beh.mat'),...
        'reg_cr')
        

    clearvars -except i dirs expData;
    close all;
end

% plays sound when done
load train;
sound(y,Fs);

toc