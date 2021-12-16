% % bandit_setPathList %
%
%PURPOSE: To set up the paths to run the code.
%AUTHORS: AC Kwan, 170518.
%

data_dir = 'E:\data\matching_pennies\pupil';
%data_dir = '/Volumes/haha/MatchingPennies/pupilData'
code_dir = 'J:\MatchingPennies';
%code_dir = '/Volumes/haha/MatchingPennies/pupilData'
% add the paths needed for this code
path_list = {...
    code_dir;...
    fullfile(code_dir,'common functions');...
    fullfile(code_dir,'common functions','cbrewer');...
    fullfile(code_dir,'value behavior');...
    fullfile(code_dir,'value sim');...
    fullfile(code_dir,'value exp lists');...
    };
addpath(path_list{:});