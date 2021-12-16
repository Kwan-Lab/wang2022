function newdataIndex = sortdataIndex(dataIndex)

% sort the dataindex by acquiring date within subject
% Author: Hongli Wang

animal = unique(dataIndex.Animal);
startInd = 0;
for ii = 1:length(animal)
    sortIndex = find(strcmp(dataIndex.Animal, animal(ii)));
    dataInd = sortrows(dataIndex(sortIndex, :),'DateNumber');
    
    % get new Index
    for jj = 1:height(dataInd)
        dataInd.Index(jj)=startInd + jj;
    end
    startInd = startInd + height(dataInd);
    if ii ==1
        newdataIndex = dataInd;
    else
        newdataIndex = [newdataIndex;dataInd];
    end
end

