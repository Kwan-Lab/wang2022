function bandit_createPupilFiles(dataIndex)

% read in csv and .mat files, save the pupil diameter and center position into a new .mat file
% also plot a control figure to show the recording did not miss any frames
% apply low-pass filter on raw data with a 4 Hz cutoff frequency

nFiles = size(dataIndex,1);

for ii = 1:nFiles
    pupCreated = 0;
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},['bandit',dataIndex.LogFileName{ii}(end-30:end-4),'_beh.mat']));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    
    
    % get the pupillometry files
    pupilFolder = dataIndex.LogFilePath{ii};
    date = num2str(dataIndex.DateNumber(ii));
    cd(pupilFolder)
    dateStrCSV = ['*', date(1:6),'*.csv'];
    dateStrMat = ['*', date(1:6),'*.mat'];
    pupilCSV = dir(dateStrCSV);
    pupilMat = dir(dateStrMat);
    
    
    % check if only one file is found
    
    % check whether pup mat already exist
    if length(pupilCSV) == 1
        disp(['--- Detected: ' pupilCSV.name]);
        
        fn = dir(fullfile(fn_beh.folder,[pupilCSV.name(1:end-9),'_pup.mat']));
        if size(fn,1)>0
            pupCreated = 1;
        end
        
        if pupCreated == 0
            
            
            % if only one match file found
            
            savematpath = fn_beh.folder;
            savefigpath = fullfile(fn_beh.folder, date(1:6));
            % load the csv file
            [stackInfo, pupilData] = loadPupilFiles(pupilCSV, pupilMat);
            
            frameTime = zeros(1, length(pupilData));
            for kk = 1:length(pupilData)
                frameTime(kk) = trialData.triggerTimes(1)-stackInfo.delayTime + 0.05*(kk-1);
            end
            pupil.t = frameTime;
            %  pupil.center = pupilCenter;
            % use time from logfile to get rid of the invalid
            % recordings
            latestTime = trialData.outcomeTimes(end)+6;  %this is as far as logfile goes
            % check whether frames are equal in log file and
            pupilData = pupilData(pupil.t<latestTime);
            pupil.t = pupil.t(pupil.t<latestTime);
            % pupil.center = pupil.center(pupil.t<latestTime);
            
            check_pupilTriggerTimes (stackInfo, trialData, savefigpath);
            
            
            % try isoutlier, looks good
            pupil_temp = pupilData;
            %center_temp = pupil.center;
            
            % low-pass filter the data
            cutF = 4;  % cutoff frequency 4 Hz
            fs = 20;  % recorded 20 Hz
            pupil_filtered = lowpass(pupil_temp, cutF, fs);
            %center_filtered = lowpass(center_temp,cutF, fs);
            
            % get tonic pupil activity
            % Tomas Knapen 2016 Plos One
            % low-pass filter, cut off frequency 0.02 Hz
            cutTonicF = 0.02;
            pupil_tonic = lowpass(pupil_temp, cutTonicF, fs);
            pupil.tonic = pupil_tonic;
            
            % get rid of blinks by set outlier to NaN
            TF = isoutlier(pupil_filtered);
            pupil_filtered(TF) = NaN;
            
            % TF_center = isoutlier(center_filtered);
            % center_filtered(TF_center) = NaN;
            
            
            % plot the time around scenario start (start  to 30 second)
            figure;
            plot(pupil.t(pupil.t<30), pupil_filtered(pupil.t<30));
            %hold on; plot(pupil.t(pupil.t<30), center_filtered(pupil.t<30));
            %legend('Diameter','center');
            title('Raw pupil size around scenario start');
            xlabel('Time from scenario start (s)');
            ylabel('Raw pupil size (px)');
            print(gcf,'-dpng',[savefigpath,'\pupilsize_aroundstart']);    %png format
            saveas(gcf, [savefigpath,'\pupilsize_aroundstart'], 'fig');
            
            close;
            
            % get running z-score in 10 min running windows
            windowSize = 10*60*20;
            dia_runningz = zeros(1, length(pupil_filtered));
            for jj = 1:windowSize:length(dia_runningz)
                if jj+windowSize-1<length(dia_runningz)
                    dia_runningz(jj:jj+windowSize-1) = (pupil_filtered(jj:jj+windowSize-1) - nanmean(pupil_filtered(jj:jj+windowSize-1))) / nanstd(pupil_filtered(jj:jj+windowSize-1));
                else
                    dia_runningz(jj:length(dia_runningz)) = (pupil_filtered(jj:length(dia_runningz)) - nanmean(pupil_filtered(jj:length(dia_runningz)))) / nanstd(pupil_filtered(jj:length(dia_runningz)));
                end
            end
            
            pupil.dia = dia_runningz;
            
            % get frame time
            
            
            % check the time of trialDat
            save(fullfile(savematpath,[pupilMat.name(1:end-9),'_pup.mat']),...
                'pupil');
            
            clearvars -except ii dirs dataIndex ;
        end
        
    elseif length(pupilCSV) == 0
        disp('None csv files found');
    else
        disp('Multiple csv files found');
    end
    
end

end
