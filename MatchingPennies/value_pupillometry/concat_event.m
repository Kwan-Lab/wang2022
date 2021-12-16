function triggerMat = concat_event(params)

% take a structure of trigger event, concatenate them into a matrix

triggerMat = [];
fields = fieldnames(params);
for jj = 1:length(fields)
    if contains(fields{jj},'trigEvent')
        triggerMat = [triggerMat params.(fields{jj})];
    end
end