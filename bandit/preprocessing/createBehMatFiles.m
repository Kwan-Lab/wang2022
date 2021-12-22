function newDataIndex = createBehMatFiles(dataIndex)
% % createBehMatFiles %
%PURPOSE:   Analyze each logfile specified in dataIndex and save the
%           results in a behavioral .mat file
%AUTHORS:   H Atilgan and AC Kwan 191127
%
%INPUT ARGUMENTS
%   dataIndex:  a table of the data files
%
%OUTPUT ARGUMENTS
%   newDataIndex:  a table of the data files, now with the BehFileCreated =
%                  true for those files that are processed
%

%% Create logfile-info related table

nFile = size(dataIndex,1);

behIndex = table(...
    cell(nFile,1),...
    cell(nFile,1),...
    NaN(nFile,1),...
    NaN(nFile,1),...
    NaN(nFile,1)...
    );

behIndex.Properties.VariableNames = {...
    'Animal',...
    'Experiment',... Name of the scenario file ran by NBS Presentation
    'Phase',...      Phase3=Reversal, Phase6=6 sets of reward prob, Phase8=same as phase 3, with pupil recording
    'DateNumber',... Date/Time e.g, 1806291321 = 2018 June 29 13:21
    'BehCreated'...  Has the behavioral .mat file been created
    };

%% parse and plot the logfiles specified in dataIndex

disp(['-----------------------------------------------------------']);
disp(['--- Detected: ' int2str(nFile) ' behavioral logfiles.']);
disp(['-----------------------------------------------------------']);

for i = 1:nFile
    
    % Is the analysis file already created?
    fn = dir(fullfile(dataIndex.BehPath{i},[dataIndex.LogFileName{i}(1:end-4),'_beh.mat']));
    if size(fn,1)>0
        behIndex.BehCreated(i) = 1;
    end
    
    if isnan(behIndex.BehCreated(i))
        disp(['Parsing file #' int2str(i) '.']);
        disp(['    ' dataIndex.LogFileName{i}]);
        
        % extract data for each trial from raw logfile
        [ logData ] = parseLogfile(dataIndex.LogFilePath{i}, dataIndex.LogFileName{i});
        
        logfileData.Animal = logData.subject;
        logfileData.Experiment = logData.scenario;
        yr=num2str(logData.dateTime{1}(9:10));
        mo=num2str(logData.dateTime{1}(1:2));
        day=num2str(logData.dateTime{1}(4:5));
        hr=num2str(logData.dateTime{2}(1:2));
        min=num2str(logData.dateTime{2}(4:5));
        logfileData.DateNumber=str2num([yr mo day hr min]);
        
        logfileData.Experiment
        
        if strcmp(logfileData.Experiment,'Phase3R_71_NoCue')...
                || strcmp(logfileData.Experiment,'Phase3_R71_NoCue')...
                || strcmp(logfileData.Experiment,'Phase3_R71NoCue')
            logfileData.Phase=3;
        end
        if strcmp(logfileData.Experiment,'Phase6_Value')
            logfileData.Phase=6;
        end
        if strcmp(logfileData.Experiment,'Phase8_R71NoCueWithPupil')
            logfileData.Phase=8;
        end
        if strcmp(logfileData.Experiment,'Phase3_R71NoCue_Opto')
            logfileData.Phase=21;
        end
        
        [ sessionData, trialData ] = value_getSessionData(logData,logfileData.Phase);
        
        trials = bandit_getTrialMasks(trialData);
        
        % 
        %save behavioral .mat file
        save(fullfile(dataIndex.BehPath{i}, ['bandit',dataIndex.LogFileName{i}(end-30:end-4),'_beh.mat'])',...
            'sessionData','trialData','logfileData','trials');
        
        behIndex.BehCreated(i) = 1;
        
    else  %else if the behavioral .mat file already created
        disp(['Loading file #' int2str(i) '.']);
        disp(['    ' dataIndex.LogFileName{i}]);

        load(fullfile(dataIndex.BehPath{i},['bandit_',dataIndex.LogFileName{i}(end-29:end-4),'_beh.mat']));
    end
    
    behIndex.Animal(i) = logfileData.Animal;
    behIndex.Experiment(i) = logfileData.Experiment;
    behIndex.Phase(i) = logfileData.Phase;
    behIndex.DateNumber(i) = logfileData.DateNumber;
        
    clearvars -except i dataIndex behIndex
    
end

%% Add the logfile-extracted information into the database index

newDataIndex = [dataIndex behIndex];

end
