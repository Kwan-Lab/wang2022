function [distance, xoff, yoff] = loadLabel(humanLabel, DLCLabel)

% modified from function loadPupil
% to load csv file with training dataset labeled by human (without
% likelihood) and DLC

% then calculate the distance between human label and DLC label

M1 = readmatrix(humanLabel);
%M1 = csvread(filename,3);
upX1 = M1(:, 2); upY1 = M1(:,3);
downX1 = M1(:,4); downY1 = M1(:,5);
leftX1 = M1(:,6); leftY1 = M1(:,7);
rightX1 = M1(:,8); rightY1 = M1(:,9);
centerX1 = M1(:,10); centerY1 = M1(:,11);


M2 = readmatrix(DLCLabel);
%M1 = csvread(filename,3);
upX2 = M2(:, 2); upY2 = M2(:,3);
downX2 = M2(:,5); downY2 = M2(:,6);
leftX2 = M2(:,8); leftY2 = M2(:,9);
rightX2 = M2(:,11); rightY2 = M2(:,12);
centerX2 = M2(:,14); centerY2 = M2(:,15);

distance = zeros(5, size(M1,1));
distance(1,:) = sqrt((upX1-upX2).^2 + (upY1-upY2).^2);
distance(2,:) = sqrt((downX1-downX2).^2 + (downY1-downY2).^2);
distance(3,:) = sqrt((leftX1-leftX2).^2 + (leftY1-leftY2).^2);
distance(4,:) = sqrt((rightX1-rightX2).^2 + (rightY1-rightY2).^2);
distance(5,:) = sqrt((centerX1-centerX2).^2 + (centerY1-centerY2).^2);

%% get x and y offset separately
xoff = zeros(5, size(M1,1));
xoff(1,:) = upX1-upX2; xoff(2,:) = downX1-downX2;
xoff(3,:) = leftX1-leftX2; xoff(4,:) = rightX1-rightX2;
xoff(5,:) = centerX1- centerX2;

yoff = zeros(5, size(M1,1));
yoff(1,:) = upY1-upY2; yoff(2,:) = downY1-downY2;
yoff(3,:) = leftY1-leftY2; yoff(4,:) = rightY1-rightY2;
yoff(5,:) = centerY1- centerY2;

%% calculate the pupil diameter from the dlc label
diaVer = sqrt((upX2-downX2).^2 + (upY2 - downY2).^2);
diaHor = sqrt((leftX2-rightX2).^2 + (leftY2 - rightY2).^2);
verOut = isoutlier(diaVer);
figure; plot(diaVer)
hold on; plot(verOut*20);
hold on; plot(diaHor);
end