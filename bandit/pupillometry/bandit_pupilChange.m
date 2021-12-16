function bandit_pupilChange(dataIndex, savefigpath)

%%  analyze pupil change
% reference Nassar, 2012, NAt. Neurosci.

%% to do:
% 1. calculate pupil change in: ITI - outcome; ITI - ITI?
% 2. taking derivative

% then run multilinear regression on the pupil change and derivative

%nFiles = size(dataIndex,1);


if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);

nFiles = size(dataIndex,1);

for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},['bandit',dataIndex.LogFileName{ii}(end-30:end-4),'_beh.mat']));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    % load pupil files
    date = num2str(dataIndex.DateNumber(ii));
    pup_name = [dataIndex.BehPath{ii},'/', '*',date(1:6),'*_pup.mat'];
    fn_pup = dir(pup_name);
    
    % pupil response with 20 Hz signal after cue
    % baseline: -1-0s; reponse: -3-5s
    % timeStep = 0.1;
    % startTime = 0; endTime = 5;
    % numDP = (endTime - startTime)/timeStep;
    
    if length(fn_pup) == 1
        load(fullfile(fn_pup.folder,fn_pup.name));
        
        % calculate the pupil response between cue and outcome
        % average Z-score of 0-1s after outcome - average z-score of 0-1s
        % before cue?
        pupilResp = [];
        pupilRespT = [];
        
        % find the maximum length of the pupil recording/behavior file
        % 4 second is the minimum trial duration
        maxTime = min(pupil.t(end), trialData.cueTimes(end)+4);
        maxTrial = find(trialData.cueTimes+4 <= maxTime, 1, 'last');
        for ll = 1:maxTrial
            
            % baseline
            timeCueStart = trialData.cueTimes(ll)-2; timeCueEnd = trialData.cueTimes(ll)-1;
            pupilBefore = pupil.dia(pupil.t >= timeCueStart & pupil.t < timeCueEnd);
            
            pupilAfter = pupil.dia(pupil.t > trialData.cueTimes(ll)-2 & pupil.t < trialData.cueTimes(ll)+6);
            tAfter = pupil.t(pupil.t > trialData.cueTimes(ll)-2 & pupil.t < trialData.cueTimes(ll)+6);
            tempResp = [pupilAfter - nanmean(pupilBefore)];
            
             if length(tempResp) < 160
                tempResp = [tempResp, NaN(1, (160-length(tempResp)))];
                tempt = tAfter;
                tAfter = [tempt, tempt(end)+0.05:0.05:tempt(end)+0.05*(160-length(tempt))];
            end
            pupilResp = [pupilResp; tempResp];
            pupilRespT = [pupilRespT; tAfter];

        end
        
        % reshape 
        %[dimx,dimy] = size(pupilResp);
        pupil.resp = pupilResp;
         pupil.respT = pupilRespT;

        % get derivate of pupil zscore
%         pupil.derivative = diff(pupil.dia)/mean(diff(pupil.t));
%         pupil.derivative = [pupil.derivative,NaN];  
        
        % save the pupil data
        save(fullfile(fn_pup.folder,fn_pup.name), 'pupil');
    end
    
    
end

