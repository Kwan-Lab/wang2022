%% Analyze pupillometry data -- compute running average z-score of pupil diamter
% run this after running start_beh.m

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

    case 4  % Load the data set from the M2 SST data set
        data_subdir = fullfile(data_dir,'M2 SST');
        [ dirs, expData ] = expData_M2_SST(data_subdir);
        
    case 5  % Load the data set from the M2 Discrim data set
        data_subdir = fullfile(data_dir,'M2 discrim');
        [ dirs, expData ] = expData_M2_Discrim(data_subdir);
    case 6  % Load the data set from the M2 Discrim data set
        data_subdir = fullfile(data_dir,'M2 omission');
        [ dirs, expData ] = expData_M2_Omission(data_subdir);
    
    case 7 % load the data set from pupillometry
        data_subdir = fullfile(data_dir);
        [ dirs, expData ] = expData_pupil(data_subdir);
        
    case 20  % Load the data set from the NN paper - M2 imaging
        data_subdir = fullfile(data_dir,'V1 NN');
        [ dirs, expData ] = expData_V1_NN(data_subdir);
   
    case 99  % Load a test data set
        data_subdir = fullfile(data_dir,'test');
        [ dirs, expData ] = expData_test(data_subdir);

end

%% process data files
for i = 1:numel(expData)
    disp(['Processing file ' int2str(i) ' out of ' int2str(numel(expData)) '.']);

    % setup/create subdirectories to save analysis and figures
    if isfield(expData(i),'onefolder')   %all the log files saved in one folder
        temp = sprintf('%04d',i);     %then create subfolder with names based on the file orders
        savematpath = fullfile(dirs.analysis,expData(i).sub_dir,temp);
        savefluofigpath = fullfile(dirs.analysis,expData(i).sub_dir,temp,'figs-fluo');
    else                              %otherwise, one directory = one data set, create subfolder using the same directory name
        savematpath = fullfile(dirs.analysis,expData(i).sub_dir);
        savefluofigpath = fullfile(dirs.analysis,expData(i).sub_dir,'figs-fluo');
    end
    if ~exist(savefluofigpath,'dir')
        mkdir(savefluofigpath);
    end
    
    % load the saved behavioral analysis (from start_beh.m)
    cd(savematpath);
    load('beh.mat');
    
    % load the stackinfo
    cd(fullfile(dirs.data,expData(i).sub_dir))
    stackFile = dir('*.mat');
    stackInfo = load(stackFile.name);
    
    % load the csv file
    pupilFile = dir('*.csv');
    pupilDia = loadPupil(pupilFile.name);
    
    %% make an initial plot of the imaging and behavior trigger times to see
    %if there are obvious problems
    cd(savefluofigpath);    
    % tlabel = char(stackInfo.savFile_name);
    check_pupilTriggerTimes ( pupilDia, stackInfo, trialData);
    
    %% calculate the pupil diameter z-score (running window 100 trials)
    % set the blink to NaN, looks like the normal pupil diameter won't exceed
% 70 (horizontal)

    dia_noBlink = pupilDia;
    dia_noBlink(dia_noBlink > 70 | dia_noBlink < 25) = NaN;

% get z-score
    dia_z = (dia_noBlink - nanmean(dia_noBlink)) / nanstd(dia_noBlink);


% use running average to get z-score / 10 mins:
    windowSize = 10*60*20;
    dia_runningz = zeros(1, length(dia_z));
    for jj = 1:windowSize:length(dia_noBlink)
        if jj+windowSize-1<length(dia_z)
            dia_runningz(jj:jj+windowSize-1) = (dia_noBlink(jj:jj+windowSize-1) - nanmean(dia_noBlink(jj:jj+windowSize-1))) / nanstd(dia_noBlink(jj:jj+windowSize-1));
        else
            dia_runningz(jj:length(dia_z)) = (dia_noBlink(jj:length(dia_z)) - nanmean(dia_noBlink(jj:length(dia_z)))) / nanstd(dia_noBlink(jj:length(dia_z)));
        end
    end

    pupil.dia = dia_runningz;
    
    % get frame time
    frameTime = zeros(1, length(pupilDia));
    for kk = 1:length(pupilDia)
        frameTime(kk) = trialData.triggerTimes(1)-stackInfo.delayTime + 0.05*(kk-1);
    end
    pupil.t = frameTime;
    %% save the dF/F so don't have to re-compute every time
    save(fullfile(savematpath,'pupZ.mat'),...
        'pupil');

    clearvars -except i dirs expData;
end

% plays sound when done
load train;
sound(y,Fs);

toc