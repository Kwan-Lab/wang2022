function newDataIndex = addIndexLesion(dataIndex)
% % addIndexLesion %
%PURPOSE:   Append to the DataIndex information about the lesions
%AUTHORS:   H Atilgan and AC Kwan 191127
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

lesionIndex = table(...
    NaN(nFile,1),...
    NaN(nFile,1)... LesionSide
    );

lesionIndex.Properties.VariableNames = {...
    'Lesioned',...  NaN=no lesion or pre-lesion, 1=post-lesion
    'LesionSide'... 1=Left, 2=Right, 3=Bilateral, 4=Saline
    };

%% Animals with lesions, and the date of lesion

% Animal, Surgery date: 18 06 20 0000, 1:Left, 2:Right, 3:Bilateral, 4:Saline
lesionRecord.name{1}='1806'; lesionRecord.date(1)=1806290000; lesionRecord.lesionside(1)=1;
lesionRecord.name{2}='1807'; lesionRecord.date(2)=1806290000; lesionRecord.lesionside(2)=2;

lesionRecord.name{3}='1808'; lesionRecord.date(3)=1810030000; lesionRecord.lesionside(3)=1;
lesionRecord.name{4}='1810'; lesionRecord.date(4)=1810030000; lesionRecord.lesionside(4)=1;

lesionRecord.name{5}='18102'; lesionRecord.date(5)=1902240000; lesionRecord.lesionside(5)=1;
lesionRecord.name{6}='18103'; lesionRecord.date(6)=1902240000; lesionRecord.lesionside(6)=4;
lesionRecord.name{7}='18104'; lesionRecord.date(7)=1902240000; lesionRecord.lesionside(7)=2;

lesionRecord.name{8}='18106'; lesionRecord.date(8)=1902140000; lesionRecord.lesionside(8)=3;
lesionRecord.name{9}='18107'; lesionRecord.date(9)=1902140000; lesionRecord.lesionside(9)=1;
lesionRecord.name{10}='18109'; lesionRecord.date(10)=1902140000; lesionRecord.lesionside(10)=2;

lesionRecord.name{11}='19102'; lesionRecord.date(11)=1905200000; lesionRecord.lesionside(11)=3;
lesionRecord.name{12}='19106'; lesionRecord.date(12)=1905200000; lesionRecord.lesionside(12)=3;
lesionRecord.name{13}='19107'; lesionRecord.date(13)=1905200000; lesionRecord.lesionside(13)=1;
lesionRecord.name{14}='19109'; lesionRecord.date(14)=1905200000; lesionRecord.lesionside(14)=2;

lesionRecord.name{15}='19114'; lesionRecord.date(15)=1908260000; lesionRecord.lesionside(15)=4;

lesionRecord.name{16}='19116'; lesionRecord.date(16)=1907220000; lesionRecord.lesionside(16)=4;
lesionRecord.name{17}='19117'; lesionRecord.date(17)=1907220000; lesionRecord.lesionside(17)=3;
lesionRecord.name{18}='19118'; lesionRecord.date(18)=1907220000; lesionRecord.lesionside(18)=4;

%%  Get info about each logfile
for b = 1:nFile
    
    % Info about lesion
    if ismember(dataIndex.Animal{b},lesionRecord.name) %was it a lesioned animal?
        lesionIndex.LesionSide(b) = lesionRecord.lesionside(strcmp(lesionRecord.name,dataIndex.Animal{b}));        
        
        if  dataIndex.DateNumber(b)>lesionRecord.date(strcmp(lesionRecord.name,dataIndex.Animal{b})) %after the date of lesion?
            lesionIndex.Lesioned(b)  = 1;
        end
    end
    
end

%% Add the lesion information into the database index

newDataIndex = [dataIndex lesionIndex];

end
