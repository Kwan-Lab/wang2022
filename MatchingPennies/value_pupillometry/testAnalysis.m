% read in the pupil mark data and calculate the std
clear all
DLCFile = "F:\pupildata\pupilDLCTest\testVideoDeepCut_resnet50_PupillometryMay31shuffle1_200000.csv";
setup_figprop;

% skip first 4 rows 
M = csvread(DLCFile,3);
upX_DLC = M(42:192, 2); upY_DLC = M(42:192,3); 
downX_DLC = M(42:192,5); downY_DLC = M(42:192,6);
leftX_DLC = M(42:192,8); leftY_DLC = M(42:192,9);
rightX_DLC = M(42:192,11); rightY_DLC = M(42:192,12);
centerX_DLC = M(42:192,14); centerY_DLC = M(42:192,15);


%% read the pupil coordinate marked by experimenter

rootPath = 'F:\pupildata\pupilDLCTest\matData\markResults\';
cd(rootPath);
matFiles = dir('*.mat');

% the fifth one is marked by another experimenter
% read the mat files
upX = zeros(5, 151); upX_2 = zeros(1, 151);
downX = zeros(5,151); downX_2 = zeros(1,151);
leftX = zeros(5, 151); leftX_2 = zeros(1,151);
rightX = zeros(5,151); rightX_2 = zeros(1, 151);
centerX = zeros(5, 151); centerX_2 = zeros(1, 151);
upY= zeros(5, 151); upY_2 = zeros(1, 151);
downY = zeros(5,151); downY_2 = zeros(1,151);
leftY = zeros(5, 151); leftY_2 = zeros(1,151);
rightY = zeros(5,151); rightY_2 = zeros(1, 151);
centerY = zeros(5, 151); centerY_2 = zeros(1, 151);

for ii = 1:length(matFiles)
    load(matFiles(ii).name);
    if ii == 5
        upX_2(1,:) = data(1, 1, :); upY_2(1,:) = data(1, 2, :);
        downX_2(1,:) = data(2, 1, :); downY_2(1,:) = data(2, 2, :);
        leftX_2(1,:) = data(3, 1, :); leftY_2(1,:) = data(3, 2, :);
        rightX_2(1,:) = data(4, 1, :); rightY_2(1,:) = data(4, 2, :);
        centerX_2(1,:) = data(5, 1, :); centerY_2(1,:) = data(5, 2, :);
    else
        upX(ii,:) = data(1, 1, :); upY(ii,:) = data(1, 2, :);
        downX(ii,:) = data(2, 1, :); downY(ii,:) = data(2, 2, :);
        leftX(ii,:) = data(3, 1, :); leftY(ii,:) = data(3, 2, :);
        rightX(ii,:) = data(4, 1, :); rightY(ii,:) = data(4, 2, :);
        centerX(ii,:) = data(5, 1, :); centerY(ii,:) = data(5, 2, :);
    end
        
end

% get rid of the emtpy one
upX(5,:) = []; upY(5,:) = [];
downX(5,:) = []; downY(5,:) = [];
leftX(5,:) = []; leftY(5,:) = [];
rightX(5,:) = []; rightY(5,:) = [];
centerX(5,:) = []; centerY(5,:) = [];

%% plot some basic stats of the data
figure;
hold on;scatter(mean(upX),mean(upY), 'b');
hold on;scatter(upX_2,upY_2, 'r');
hold on;scatter(mean(downX),mean(downY), 'b');
hold on;scatter(downX_2,downY_2, 'r');
hold on;scatter(mean(leftX),mean(leftY), 'b');
hold on;scatter(leftX_2,leftY_2, 'r');
hold on;scatter(mean(rightX),mean(rightY), 'b');
hold on;scatter(rightX_2,rightY_2, 'r');
hold on;scatter(mean(centerX),mean(centerY), 'b');
hold on;scatter(centerX_2,centerY_2, 'r');
hold on; scatter(leftX_DLC, leftY_DLC, 'y');
hold on; scatter(rightX_DLC, rightY_DLC, 'y');
hold on; scatter(upX_DLC, upY_DLC, 'y');
hold on; scatter(downX_DLC, downY_DLC, 'y' );
hold on; scatter(centerX_DLC, centerY_DLC, 'y');
    
title('Estimated coordinates');
print(gcf,'-dpng','estimated coordinates');    %png format
saveas(gcf, 'estimated coordinates', 'fig');


%% load a example frame to plot the stats
dataPath = 'F:\pupildata\pupilDLCTest';
cd(dataPath);
matFiles = dir('*.mat');
frames = cell(0);
for ii = 2:length(matFiles)
    load(matFiles(ii).name);
    for jj = 1:length(data)
        frames{end+1} = data{jj};
    end
end
exampleFrame = frames{1};

%% get the std within subjects
meanLeftX = mean(leftX); meanLeftY = mean(leftY);
meanRightX = mean(rightX); meanRightY = mean(rightY);
meanUpX = mean(upX); meanUpY = mean(upY);
meanDownX = mean(downX); meanDownY = mean(downY);
meanCenterX = mean(centerX); meanCenterY = mean(centerY);

% calculate the distance from mean
leftDis = mean((sqrt((leftX - meanLeftX).^2 + (leftY - meanLeftY).^2)));
rightDis = mean((sqrt((rightX - meanRightX).^2 + (rightY - meanRightY).^2)));
upDis = mean((sqrt((upX - meanUpX).^2 + (upY - meanUpY).^2)));
downDis = mean((sqrt((downX - meanDownX).^2 + (downY - meanDownY).^2)));
centerDis = mean((sqrt((centerX - meanCenterX).^2 + (centerY - meanCenterY).^2)));
% figure;plot(leftDis')
% hold on; plot(rightDis');
% hold on; plot(upDis');
% hold on; plot(downDis');
% hold on; plot(centerDis');

figure;boxplot([upDis; downDis; leftDis; rightDis; centerDis]'); 
title('Distance from mean, manual');
print(gcf,'-dpng','distance from mean manual');    %png format
saveas(gcf, 'distance from mean manual', 'fig');

% bootstrap to estimate the 95% distance
% resample the sample, calculate the mean, get 2.5% and 97.5 quantile
maxIter = 10000;

% left
[leftStat, ~] = bootstrp(maxIter, @mean, leftDis);
% get percentile
leftStatSort = sort(leftStat);
% leftLow = leftStatSort(maxIter * 0.025); 
leftHigh = leftStatSort(maxIter * 0.95);

%right
[rightStat, ~] = bootstrp(maxIter, @mean, rightDis);
% get percentile
rightStatSort = sort(rightStat);
%rightLow = rightStatSort(maxIter * 0.025); 
rightHigh = rightStatSort(maxIter * 0.95);

%up
[upStat, ~] = bootstrp(maxIter, @mean, upDis);
% get percentile
upStatSort = sort(upStat);
% upLow = upStatSort(maxIter * 0.025); 
upHigh = upStatSort(maxIter * 0.95);

%down
[downStat, ~] = bootstrp(maxIter, @mean, downDis);
% get percentile
downStatSort = sort(downStat);
% downLow = downStatSort(maxIter * 0.025); 
downHigh = downStatSort(maxIter * 0.95);

%center
[centerStat, ~] = bootstrp(maxIter, @mean, centerDis);
% get percentile
centerStatSort = sort(centerStat);
% centerLow = centerStatSort(maxIter * 0.025); 
centerHigh = centerStatSort(maxIter * 0.95);


%% get the bewteen experimenters variance
leftDis_2 = sqrt((leftX_2 - meanLeftX).^2 + (leftY_2 - meanLeftY).^2);
rightDis_2 = sqrt((rightX_2 - meanRightX).^2 + (rightY_2 - meanRightY).^2);
upDis_2 = sqrt((upX_2 - meanUpX).^2 + (upY_2 - meanUpY).^2);
downDis_2 = sqrt((downX_2 - meanDownX).^2 + (downY_2 - meanDownY).^2);
centerDis_2 = sqrt((centerX_2 - meanCenterX).^2 + (centerY_2 - meanCenterY).^2);
% figure;plot(leftDis_2')
% hold on; plot(rightDis_2');
% hold on; plot(upDis_2');
% hold on; plot(downDis_2');
% hold on; plot(centerDis_2');

figure;boxplot([upDis_2; downDis_2; leftDis_2; rightDis_2; centerDis_2]'); 
title('Distance from mean, manual2');
print(gcf,'-dpng','distance from mean manual 2');    %png format
saveas(gcf, 'distance from mean manual 2', 'fig');

% bootstrap
% left
[leftStat_2, ~] = bootstrp(maxIter, @mean, leftDis_2);
% get percentile
leftStatSort_2 = sort(leftStat_2);
% leftLow_2 = leftStatSort_2(maxIter * 0.025); 
leftHigh_2 = leftStatSort_2(maxIter * 0.95);

%right
[rightStat_2, ~] = bootstrp(maxIter, @mean, rightDis_2);
% get percentile
rightStatSort_2 = sort(rightStat_2);
% rightLow_2 = rightStatSort_2(maxIter * 0.025); 
rightHigh_2 = rightStatSort_2(maxIter * 0.95);

%up
[upStat_2, ~] = bootstrp(maxIter, @mean, upDis_2);
% get percentile
upStatSort_2 = sort(upStat_2);
% upLow_2 = upStatSort_2(maxIter * 0.025); 
upHigh_2 = upStatSort_2(maxIter * 0.95);

%down
[downStat_2, ~] = bootstrp(maxIter, @mean, downDis_2);
% get percentile
downStatSort_2 = sort(downStat_2);
% downLow_2 = downStatSort_2(maxIter * 0.025); 
downHigh_2 = downStatSort_2(maxIter * 0.95);

%center
[centerStat_2, ~] = bootstrp(maxIter, @mean, centerDis_2);
% get percentile
centerStatSort_2 = sort(centerStat_2);
% centerLow_2 = centerStatSort_2(maxIter * 0.025); 
centerHigh_2 = centerStatSort_2(maxIter * 0.95);

%% get variance from DLC
leftDis_DLC = sqrt((leftX_DLC' - meanLeftX).^2 + (leftY_DLC' - meanLeftY).^2);
rightDis_DLC = sqrt((rightX_DLC' - meanRightX).^2 + (rightY_DLC' - meanRightY).^2);
upDis_DLC = sqrt((upX_DLC' - meanUpX).^2 + (upY_DLC' - meanUpY).^2);
downDis_DLC = sqrt((downX_DLC' - meanDownX).^2 + (downY_DLC' - meanDownY).^2);
centerDis_DLC = sqrt((centerX_DLC' - meanCenterX).^2 + (centerY_DLC' - meanCenterY).^2);
% figure;plot(leftDis_DLC')
% hold on; plot(rightDis_DLC');
% hold on; plot(upDis_DLC');
% hold on; plot(downDis_DLC');
% hold on; plot(centerDis_DLC');

figure;boxplot([upDis_DLC; downDis_DLC; leftDis_DLC; rightDis_DLC; centerDis_DLC]'); 
title('Distance from mean, DLC');
print(gcf,'-dpng','distance from mean DLC');    %png format
saveas(gcf, 'distance from mean DLC', 'fig');

% bootstrap
% left
[leftStat_DLC, ~] = bootstrp(maxIter, @mean, leftDis_DLC);
% get percentile
leftStatSort_DLC = sort(leftStat_DLC);
% leftLow_DLC = leftStatSort_DLC(maxIter * 0.025); 
leftHigh_DLC = leftStatSort_DLC(maxIter * 0.95);

%right
[rightStat_DLC, ~] = bootstrp(maxIter, @mean, rightDis_DLC);
% get percentile
rightStatSort_DLC = sort(rightStat_DLC);
% rightLow_DLC = rightStatSort_DLC(maxIter * 0.025); 
rightHigh_DLC = rightStatSort_DLC(maxIter * 0.95);

%up
[upStat_DLC, ~] = bootstrp(maxIter, @mean, upDis_DLC);
% get percentile
upStatSort_DLC = sort(upStat_DLC);
% upLow_DLC = upStatSort_DLC(maxIter * 0.025); 
upHigh_DLC = upStatSort_DLC(maxIter * 0.95);

%down
[downStat_DLC, ~] = bootstrp(maxIter, @mean, downDis_DLC);
% get percentile
downStatSort_DLC = sort(downStat_DLC);
% downLow_DLC = downStatSort_DLC(maxIter * 0.025); 
downHigh_DLC = downStatSort_DLC(maxIter * 0.95);

%center
[centerStat_DLC, ~] = bootstrp(maxIter, @mean, centerDis_DLC);
% get percentile
centerStatSort_DLC = sort(centerStat_DLC);
% centerLow_DLC = centerStatSort_DLC(maxIter * 0.025); 
centerHigh_DLC = centerStatSort_DLC(maxIter * 0.95);


%% down estimation is the worst (influenced by whisker)
figure;
plot(downX_DLC); hold on; plot(mean(downX));hold on;plot(downX_2);
legend('DLC estimation','manual marks', 'manual marks 2'); 
title('Down X');
print(gcf,'-dpng','Down X');    %png format
saveas(gcf, 'Down X', 'fig');

figure;
plot(downY_DLC); hold on; plot(mean(downY));
legend('DLC estimation','manual marks'); 
title('Down Y');
print(gcf,'-dpng','Down Y');    %png format
saveas(gcf, 'Down Y', 'fig');

%% plot the estimated position
figure;imshow(exampleFrame,'InitialMagnification',400);
hold on; scatter(175,95, 30,'filled');
hold on; viscircles([175, 95], upHigh, 'Color','r', 'LineWidth',1);
hold on; viscircles([175, 95], upHigh_2, 'Color','b','LineWidth',1);
hold on; viscircles([175, 95], upHigh_DLC, 'Color','g','LineWidth',1);
hold on; scatter(175, 140, 30, 'filled');
hold on; viscircles([175, 140], downHigh, 'Color','r','LineWidth',1);
hold on; viscircles([175, 140], downHigh_2, 'Color','b','LineWidth',1);
hold on; viscircles([175, 140], downHigh_DLC, 'Color','g','LineWidth',1);
hold on; scatter(155, 115, 30, 'filled');
hold on; viscircles([155, 115], leftHigh, 'Color','r','LineWidth',1);
hold on; viscircles([155, 115], leftHigh_2, 'Color','b','LineWidth',1);
hold on; viscircles([155, 115], leftHigh_DLC, 'Color','g','LineWidth',1);
hold on; scatter(199, 115, 30, 'filled');
hold on; viscircles([199, 115], rightHigh, 'Color','r','LineWidth',1);
hold on; viscircles([199, 115], rightHigh_2, 'Color','b','LineWidth',1);
hold on; viscircles([199, 115], rightHigh_DLC, 'Color','g','LineWidth',1);
hold on; scatter(175, 115, 30, 'filled');
hold on; viscircles([175, 115], centerHigh, 'Color','r','LineWidth',1);
hold on; viscircles([175, 115], centerHigh_2, 'Color','b','LineWidth',1);
hold on; viscircles([175, 115], centerHigh_DLC, 'Color','g','LineWidth',1);
title('95% percentile of the estimations');
print(gcf,'-dpng','95% CI');    %png format
saveas(gcf, '95 CI', 'fig');

%% plot the pupil diameter estimation
diaVer = sqrt((mean(upX) - mean(downX)).^2 + (mean(upY) - mean(downY)).^2);
diaHor = sqrt((mean(leftX) - mean(rightX)).^2 + (mean(leftY) - mean(rightY)).^2);
diaVer_DLC = sqrt((upX_DLC - downX_DLC).^2 + (upY_DLC - downY_DLC).^2);
diaHor_DLC = sqrt((leftX_DLC - rightX_DLC).^2 + (leftY_DLC - rightY_DLC).^2);
diaVer_2 = sqrt((upX_2 - downX_2).^2 + (upY_2 - downY_2).^2);
diaHor_2 = sqrt((leftX_2 - rightX_2).^2 + (leftY_2 - rightY_2).^2);

figure;
subplot(2,1,1)
plot(diaVer, 'LineWidth', 2);
hold on; plot(diaVer_2, 'LineWidth', 2);
hold on; plot(diaVer_DLC, 'LineWidth', 2);
legend('manual', 'manual2', 'DLC');
title('Vertical diameter');
subplot(2, 1, 2)
plot(diaHor, 'LineWidth', 2);
hold on; plot(diaHor_2, 'LineWidth', 2);
hold on; plot(diaHor_DLC, 'LineWidth', 2);
legend('manual', 'manual2', 'DLC');
title('Horizontal diameter');
print(gcf,'-dpng','hor_ver diameter');    %png format
saveas(gcf, 'hor_ver diameter', 'fig');

figure;
plot((diaVer+diaHor)/2, 'LineWidth', 2);
hold on;plot((diaVer_2 + diaHor_2)/2, 'LineWidth', 2);
hold on; plot((diaVer_DLC + diaHor_DLC)/2, 'LineWidth', 2);
legend('manual', 'manual2', 'DLC');
title('Averaged diameter');
print(gcf,'-dpng','averaged diameter');    %png format
saveas(gcf, 'averaged diameter', 'fig');

%% plot the 78th mark to actually see the results
figure; imshow(frames{78},'InitialMagnification',400);
hold on; scatter([mean(upX(:,78)),mean(downX(:,78))],[mean(upY(:,78)),mean(downY(:,78))],'filled');
hold on; scatter([upX_2(78), downX_2(78)], [upY_2(78), downY_2(78)], 'filled');
hold on; scatter([upX_DLC(78), downX_DLC(78)],[upY_DLC(78), downY_DLC(78)], 'filled');
legend('manual', 'manual2', 'DLC');
title('The coordinates of frame #78');
print(gcf,'-dpng','frame #78');    %png format
saveas(gcf, 'frame #78', 'fig');

% indeed, the DLC mark is a little off, but does it influence the result?
figure; imshow(frames{63},'InitialMagnification',400);
hold on; scatter([mean(upX(:,63)),mean(downX(:,63))],[mean(upY(:,63)),mean(downY(:,63))],'filled');
hold on; scatter([upX_2(63), downX_2(63)], [upY_2(63), downY_2(63)], 'filled');
hold on; scatter([upX_DLC(63), downX_DLC(63)],[upY_DLC(63), downY_DLC(63)], 'filled');
legend('manual', 'manual2', 'DLC');
title('The coordinates the frame # 63');
print(gcf,'-dpng','frame #63');    %png format
saveas(gcf, 'frame #63', 'fig');
