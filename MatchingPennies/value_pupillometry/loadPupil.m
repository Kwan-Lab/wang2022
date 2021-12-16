function pupilDia = loadPupil(filename)

%load csv file labeled by DeepLabCut

M1 = readmatrix(filename);
%M1 = csvread(filename,3);
upX = M1(:, 2); upY = M1(:,3);
downX = M1(:,5); downY = M1(:,6);
leftX = M1(:,8); leftY = M1(:,9);
rightX = M1(:,11); rightY = M1(:,12);
if size(M1,2) > 13
    centerX = M1(:,14); centerY = M1(:,15);
end

diaVer = sqrt((upX-downX).^2+(upY-downY).^2);
diaHor = sqrt((leftX-rightX).^2+(leftY-rightY).^2);

pupilDia = diaHor;

end

