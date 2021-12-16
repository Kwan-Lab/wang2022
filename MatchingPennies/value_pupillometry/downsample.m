dataPath = 'E:\data\pupildata\training';
cd(dataPath);

videoFile= dir('*.avi');



for ii = 1:length(videoFile)
    frames = cell(0);
    newName = [videoFile(ii).name(1:end-4),'ds.avi'];
    v = VideoReader(videoFile(ii).name);
    numFrames = v.Duration*v.FrameRate;
    
    k = 1;
    vWrite = VideoWriter(newName, 'Uncompressed AVI');
    %open(v);
    open(vWrite);
    for jj = 1:numFrames
        if mod(jj,10000) == 1
            frames{k} = read(v, jj);
            writeVideo(vWrite,mat2gray(frames{k})); 
            k = k + 1;
        end
    end
    %close(v);
    close(vWrite);
    
end
    