function [ logData ] = parseLogfile( data_dir, logfile )
% % parseLogfile %
%
%PURPOSE: To read and parse Presentation logfile for further analysis.
%AUTHORS: MJ Siniscalchi and AC Kwan, 161209.
%         modified H Atilgan and AC Kwan, 191210.
%           modified H Wang to analyze H. Ortega's files
%INPUT ARGUMENTS
%   data_dir:    Path for logfile.   
%   logfile:     Filename for logfile.
%
%OUTPUT VARIABLES
%   logData:     Structure containing these fields:
%                {subject, dateTime, header, values}.

%% read header information and find line corresponding to Trial 1
fID=fopen(fullfile(data_dir,logfile));  %use fullfile to avoid backslash/slash difference in mac/pc
header{1} = textscan(fID,'%s %s %s %*[^\n]',1);
header{2} = textscan(fID,'%s %s %c %s %s',1);
header{3} = textscan(fID,'%s',5,'delimiter','\t');  %Labels for first 5 columns 
header{4} = textscan(fID,'%s %s %s %*[^\n]',1);     %Labels for the remainder of the columns (skip)
firstLine = 3;
while ~ismember(header{4}{2},'1')             %Find the first line that has Trial = 1
    header{4} = textscan(fID,'%s %s %*[^\n]',1);
    firstLine = firstLine + 1;
end
fclose(fID);

logData.scenario = header{1}{3};

%% modified to find the subject number, for Heather's data
pattern = 'SUB*';
startInd = regexp(header{4}{1},pattern);
if isempty(startInd{1})
    subString = header{4}{1}{1};
else
    subString = header{4}{1}{1}(startInd{1}+3:startInd{1}+4);
end
% there are two possible patterns of the subString: '3.' or '03', so
% lets check which one is this, the first two conditions are coded
% specifically to find the subject name for Heather's data
% the third one is for Huriye/Cayla's data
if ismember('.',subString)
    logData.subject{1} = subString(1);
elseif length(subString) == 2 & ~ismember('.',subString)
    logData.subject{1} = subString(2);
else
    logData.subject{1} = subString;
end


%% end of the modification
logData.dateTime = [header{2}{4:5}];
logData.header = header{3}{1}(1:5)'; %'Subject' 'Trial' 'Event Type' 'Code' 'Time'

%%  read the rest of the log file
fID=fopen(fullfile(data_dir,logfile));
if firstLine >= 3
    for j = 1:firstLine-1                        %Skip all lines before line with Trial = 1
        temp = textscan(fID,'%s %*[^\n]',1);     %Skip line
    end
else
    error('Error parsing log file');
end

% instead of reading in the typical format %s %u %s %u %d
% here going for more flexibility in the logfile structure by reading in
% everything as strings, then later convert certain columns to numbers
data = textscan(fID,'%s %s %s %s %s %*[^\n]','delimiter','\t','EmptyValue',-1);

fclose(fID);

logData.values = data;

end

