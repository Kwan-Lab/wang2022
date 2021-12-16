function [ dirs, expData ] = expData_pupil(data_dir)
% % expData_pupil %
%
%PURPOSE: Create data structure for pupillometry files and behavioral log files
%AUTHORS: AC Kwan, 170518.
%
%INPUT ARGUMENTS
%   data_dir:    The base directory to which the raw data are stored.  
%
%OUTPUT VARIABLES
%   dirs:        The subfolder structure within data_dir to work with
%   expData:     Info regarding each experiment

dirs.data = fullfile(data_dir,'data');
dirs.analysis = fullfile(data_dir,'analysis');
dirs.summary = fullfile(data_dir,'summary');

i=1;
        expData(i).sub_dir = '874_1022';
        expData(i).logfile = '874-phase2_MP_2A_pupil11.log';
% i=i+1;
%         expData(i).sub_dir = '140530 M8';
%         expData(i).logfile = 'M8-phase2_auditory_discrim4.log';
% i=i+1;
%         expData(i).sub_dir = '140531 L4';
%         expData(i).logfile = 'L4-phase2_auditory_discrim11.log';