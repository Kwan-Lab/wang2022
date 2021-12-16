% randomly pick one frame from each subject, concatenate them in .avi
% files and .mat files, save

overAllPath = 'F:\pupildata\pupilDLCTest';
dataPath = 'matData';
videoPath = 'aviData';

cd(overAllPath);
if ~exist(videoPath,'dir')
    mkdir(videoPath);
end

if ~exist(dataPath, 'dir')
    mkdir(dataPath)
end

matFiles = dir('*.mat');



%% this is a bad attempt
% % randomly picking frames
% testFrames = cell(1);
% for ii =1:length(matFiles)
%     matName = matFiles(ii).name;
%     load(matName);
%     numFrames = length(data);
%     framesPick = randsample(numFrames, 1);
%     testFrames{ii} = data{framesPick};
% end
% 
% filename = [overAllPath, '\', videoPath,'\testVideo.avi'];
% v = VideoWriter(filename, 'Uncompressed AVI');
%     
% open(v);
% videoIter = 1000;
% for ll = 1:videoIter
%     for jj = 1:length(testFrames)
%         writeVideo(v,mat2gray(data{jj})); 
%     end    
% end
% close(v);
% 
% % write down the .mat file for manual mark
% matData = cell(0);
% matIter = 100;
% matfilename = [overAllPath, '\', dataPath,'\testMat.mat'];
% for uu = 1:matIter
%     for kk = 1:length(testFrames)
%         matData{end+1} = data{kk}; 
%     end
% end
% save(matfilename, 'matData');

%% try another way this time
filename = [overAllPath, '\', videoPath,'\testVideo.avi'];
v = VideoWriter(filename, 'Uncompressed AVI');
    
open(v);
for ll = 1:length(matFiles)
    load(matFiles(ll).name);
    for jj = 1:length(data)
        writeVideo(v,mat2gray(data{jj})); 
    end    
end
close(v);

% save .mat data
matData = cell(0);
matfilename = [overAllPath, '\', dataPath,'\testMat.mat'];
for ll = 1:length(matFiles)
    load(matFiles(ll).name);
    for jj = 1:length(data)
        matData{end+1} = data{jj}; 
    end    
end

save(matfilename, 'matData');