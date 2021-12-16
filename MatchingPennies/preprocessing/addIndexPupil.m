function newDataIndex = addIndexPupil(dataIndex)
% % addIndexLesion %
%PURPOSE:   Append to the DataIndex information about the pupillometry
%AUTHORS:   Hongli Wang
%
%INPUT ARGUMENTS
%   dataIndex:  a table of the data files
%
%OUTPUT ARGUMENTS
%   newDataIndex:  a table of the data files, now including information about
%                  lesions
%

%% Create lesion-related table

nFile = size(dataIndex,1);

pupilIndex = table(...
    NaN(nFile,1),...
    NaN(nFile,1)... pupil recording Side
    );

pupilIndex.Properties.VariableNames = {...
    'pupil',...  NaN=no pupillometry recording, 1=pupillometry
    'pupilSide'... 1=Left, 2=Right
    };


%%  Get info about each logfile
for b = 1:nFile
    
    % Info about pupil recordings
    if contains(dataIndex.LogFileName{b}, 'pupil')
        pupilIndex.pupil(b)=1;
        pupilIndex.pupilSide(b) = 2;        % pupil recording always on the right for single-pupil recordings
    else
        pupilIndex.pupil(b) = 0;
    end
    
end

%% Add the lesion information into the database index

newDataIndex = [dataIndex pupilIndex];

end
