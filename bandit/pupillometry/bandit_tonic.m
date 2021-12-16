function MP_tonic(dataIndex, savefigpath)

% analyze the topic pupil dynamics

nFiles = size(dataIndex,1);
rawTonic = cell(0);
time = cell(0);
numTrials = [];
cueTimes = cell(0);
outcomeTimes = cell(0);
cutPoints = [];
files = cell(0);

if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);

for ii = 1:nFiles
    
    % load behavior files
    fn_beh = fullfile(dataIndex.BehPath{ii},['bandit',dataIndex.LogFileName{ii}(end-30:end-4),'_beh.mat']);  % for tonic activity, no need to do cut
    date = num2str(dataIndex.DateNumber(ii));
    
    % load behavior data
    load(fn_beh);
    
    %load pupil data
    pup_name = fullfile(dataIndex.BehPath{ii}, ['*',date(1:6),'*_pup.mat']);
    fn_pup = dir(pup_name);
    
    if length(fn_pup) == 1
        load(fullfile(fn_pup.folder,fn_pup.name));
        rawTonic{end+1} = pupil.tonic;
        time{end+1} = pupil.t;
        
        % load behavior
        cueTimes{end+1} = trialData.cueTimes;
        outcomeTimes{end+1} = trialData.outcomeTimes;
        numTrials(end+1) = length(trialData.cue);
        
        
        files{end+1} = fn_pup.name;
    end
        
end

% check the end time first
endTime_pupil = zeros(1, length(cueTimes));
endTime_beh = zeros(1, length(cueTimes));

for ii = 1:length(cueTimes)
    endTime_pupil(ii) = time{ii}(end)/(10^3);
    endTime_beh(ii) = outcomeTimes{ii}(end)/(10^3);
end
figure;
scatter(endTime_beh, endTime_pupil, 'black', 'filled');
hold on; plot([0 10], [0 10], 'lineWidth',1);
axis square;
xlabel('Behavior end time (10^3 s)');
ylabel('Video end time (10^3 s)');
title('Value should not deviate from diagonal');

%% calculate the z-score (baseline: pupil size before session start)
tonicpup_zscore = cell(0);
baseLine = [];

for jj = 1:length(time)
    baseMean = nanmean(rawTonic{jj}(time{jj}<0));
    baseStd = nanstd(rawTonic{jj}(time{jj}<0));
    baseLine(end+1) = baseMean;
    %calculate z-score
    tonicpup_zscore{end+1} = (rawTonic{jj} - baseMean) / baseStd;
end

% calculate average z-score (before session, first 100 trials, 200 trials
% ....)

maxTrials = max(numTrials);
aveZscore = zeros(floor(maxTrials/100),length(numTrials));


for jj = 1:length(numTrials)
    for kk = 1:size(aveZscore,1)
        % find the corresponding time interval
        if kk == 1
            timeStart = -100;
            timeEnd = 0;
        else
            if length(cueTimes{jj}) > kk*100
                timeStart = cueTimes{jj}((kk-1)*100+1);
                timeEnd = cueTimes{jj}(kk*100);
            else
                timeStart = NaN;
                timeEnd = NaN;
            end
        end
        if ~isnan(timeStart)
            aveZscore(kk,jj) = nanmean( tonicpup_zscore{jj}(time{jj} > timeStart & time{jj} < timeEnd));
        else
            aveZscore(kk,jj) = NaN;
        end
        
    end
end

% bootstrap to get median and CI, get 25%, 50%, 75% percentile
Prc1 = prctile(numTrials, 25);
Prc2 = prctile(numTrials, 50);
Prc3 = prctile(numTrials, 75);

ave1= aveZscore(:,numTrials<Prc1);
ave2 = aveZscore(:,numTrials>Prc1&numTrials<Prc2);
ave3 = aveZscore(:,numTrials>Prc2&numTrials<Prc3);
ave4 = aveZscore(:,numTrials>Prc3);

% get bootstrp
bootstat1 = bootstrp(1000, @nanmedian,ave1');
bootstat2 = bootstrp(1000, @nanmedian,ave2');
bootstat3 = bootstrp(1000, @nanmedian,ave3');
bootstat4 = bootstrp(1000, @nanmedian,ave4');

lower1 = prctile(bootstat1,5);
high1 = prctile(bootstat1, 95);
lower2 = prctile(bootstat2,5);
high2 = prctile(bootstat2,95);
lower3 = prctile(bootstat3,5);
high3 = prctile(bootstat3,95);
lower4 = prctile(bootstat4,5);
high4 = prctile(bootstat4,95);

xaxis = [-100:100:700];
figure;
% errorshade(xaxis(~isnan(lower1)),lower1(~isnan(lower1)),high1(~isnan(lower1)),c);
% hold on;
% errorshade(xaxis(~isnan(lower2)),lower2(~isnan(lower2)),high2(~isnan(lower2)),c);
% hold on;
% errorshade(xaxis(~isnan(lower3)),lower3(~isnan(lower3)),high3(~isnan(lower3)),c);
% hold on;
% errorshade(xaxis(~isnan(lower4)),lower4(~isnan(lower4)),high4(~isnan(lower4)),c);
% hold on;
median1 = nanmedian(ave1,2);
len1 = length(median1(~isnan(median1)));
errorbar(xaxis(1:len1),median1(1:len1),lower1(1:len1),high1(1:len1), 'LineWidth',2); hold on;
median2 = nanmedian(ave2,2);
len2 = length(median2(~isnan(median2)));
errorbar(xaxis(1:len2),median2(1:len2),lower2(1:len2),high2(1:len2), 'LineWidth',2);hold on;
median3 = nanmedian(ave3,2);
len3 = length(median3(~isnan(median3)));
errorbar(xaxis(1:len3),median3(1:len3),lower3(1:len3), high3(1:len3), 'LineWidth',2);hold on;
median4 = nanmedian(ave4,2);
len4 = length(median4(~isnan(median4)));
errorbar(xaxis(1:len4),median4(1:len4), lower4(1:len4), high4(1:len4), 'LineWidth',2);
xlim([-100,700]);
legend({'25 percentile','50 percentile','75 percentile','100 percentile'});
xlabel('Trials');
ylabel('Median average pupil Z-score');
print(gcf,'-dpng','tonic_pupil_by_session_length');    %png format
saveas(gcf, 'tonic_pupil_by_session_length', 'fig');



%% session length and absolute tonic pupil size

% 1) raw baseline and session length
% baseLine1 = baseLine(numTrials<Prc1);
% baseLine2 = baseLine(numTrials >= Prc1 & numTrials < Prc2);
% baseLine3 = baseLine(numTrials >= Prc2 & numTrials < Prc3);
% baseLine4 = baseLine(numTrials >= Prc3);

SessionLengthMasks= zeros(1, length(numTrials));
SessionLengthMasks(numTrials<Prc1) = 1;
SessionLengthMasks(numTrials >= Prc1 & numTrials < Prc2) = 2;
SessionLengthMasks(numTrials >= Prc2 & numTrials < Prc3) = 3;
SessionLengthMasks(numTrials >= Prc3) = 4;
%SessionLengthMasks(numTrials >= Prc2) = 3;
p = kruskalwallis(baseLine, SessionLengthMasks,'on')

% 2) raw average pupil size and session length
aveRaw = cellfun(@nanmean,rawTonic);
p = kruskalwallis(aveRaw, SessionLengthMasks,'on')
medianRaw = cellfun(@nanmedian,rawTonic);
p = kruskalwallis(medianRaw, SessionLengthMasks,'on')

%% baseline pupil size before and after session start  ~ 70s each
Index = 1:length(cueTimes);

zscore_before = zeros(1, length(time));
zscore_after = zeros(1, length(time));
for u = 1:length(time)
    zscore_before(u) = nanmean(tonicpup_zscore{u}(time{u} <0 ));
    zscore_after(u) = nanmean(tonicpup_zscore{u}(time{u} > 0 & time{u} < 70 ));
end

figure;
for kk = 1:length(zscore_before)
    plot([1,2],[zscore_before(kk), zscore_after(kk)],'-', 'Color', [0.7 0.7 0.7]);
    hold on;
end
hold on; plot([1,2], [nanmean(zscore_before),nanmean(zscore_after)],'+-black');
p = signrank(zscore_before, zscore_after);

xlim([0 3]);
ylim([-1 4]);
sigString = ['p = ',num2str(p)];
text(1.5, 3.5,sigString);
xticks([0 1 2 3]);
xticklabels({' ','Before','After',' '});
ylabel('Average z-score');
title('Average z-score before and after the session start');
print(gcf,'-dpng','tonic_pupil_by_sessionstart');    %png format
saveas(gcf, 'tonic_pupil_by_sessionstart', 'fig');
