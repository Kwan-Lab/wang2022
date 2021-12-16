function MP_pupilChange(dataIndex, savefigpath)

%%  analyze pupil change
% reference Nassar, 2012, NAt. Neurosci.
% calculate transient pupil change
% pupilChange = pupilz-score - baseline
% baseline = mean(z-score(-2 - -1s)
% grouped into trials

% then run multilinear regression on the pupil change

%nFiles = size(dataIndex,1);


if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);

nFiles = size(dataIndex,1);

for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},[dataIndex.LogFileName{ii}(1:end-4),'_beh_cut.mat']));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    % load pupil files
    date = num2str(dataIndex.DateNumber(ii));
    pup_name = [dataIndex.BehPath{ii},'/', '*',date(1:6),'*_pup.mat'];
    fn_pup = dir(pup_name);
    
    % pupil response in different time window
    % baseline: -3--2s; response: timewindow 0.5s; -2 - 6s;
    % will only use -2-6 in the regression, however calculate the response
    % in a larger time window: -3-7
    if length(fn_pup) == 1
        load(fullfile(fn_pup.folder,fn_pup.name));
        
        % calculate the pupil response
        % average z-score of 0-1s before cue?
        pupilResp = [];
        pupilRespT = [];
        pupilBase = [];
        % find the maximum length of the pupil recording/behavior file
        maxTime = min(pupil.t(end), trialData.cueTimes(end)+6);
        maxTrial = find(trialData.cueTimes+6 <= maxTime, 1, 'last');
        for ll = 1:maxTrial
            % baseline
            timeCueStart = trialData.cueTimes(ll)-2; timeCueEnd = trialData.cueTimes(ll)-1;
            pupilBefore = pupil.dia(pupil.t >= timeCueStart & pupil.t < timeCueEnd);
            
            pupilAfter = pupil.dia(pupil.t > trialData.cueTimes(ll)-2 & pupil.t < trialData.cueTimes(ll)+6);
            tAfter = pupil.t(pupil.t > trialData.cueTimes(ll)-2 & pupil.t < trialData.cueTimes(ll)+6);
            tempResp = [pupilAfter - nanmean(pupilBefore)];
            tempBase = nanmean(pupilBefore);
            % the length should be 160 from -2s to 6s, if not long enough,
            % use NaN 
            if length(tempResp) < 160
                tempResp = [tempResp, NaN(1, (160-length(tempResp)))];
                tempt = tAfter;
                tAfter = [tempt, tempt(end)+0.05:0.05:tempt(end)+0.05*(160-length(tempt))];
            end
            pupilResp = [pupilResp; tempResp];
            pupilRespT = [pupilRespT; tAfter];
            pupilBase = [pupilBase;tempBase];
        end
        
        % reshape 
%         [dimx,dimy] = size(pupilResp);
         pupil.resp = pupilResp;
         pupil.respT = pupilRespT;
         pupil.base = pupilBase;
        % get derivate of pupil zscore
%         pupil.derivative = diff(pupil.dia)/mean(diff(pupil.t));
%         pupil.derivative = [pupil.derivative,NaN];  
%         
        % save the pupil data
        save(fullfile(fn_pup.folder,fn_pup.name), 'pupil');
    end
    
    
end

