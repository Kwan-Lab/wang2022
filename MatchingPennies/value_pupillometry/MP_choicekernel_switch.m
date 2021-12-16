function MP_choicekernel_switch(dataIndex,model_path,save_path)

%check how many times the choice kernal value changes in a session

nFiles = size(dataIndex,1);
numSwitch = zeros(1, nFiles);

animalIndex = zeros(1,nFiles);
animalFolder = unique(dataIndex.LogFilePath);
for ii = 1:length(animalFolder)
    Ind = strfind(animalFolder{ii},filesep);
    startInd = Ind(end);
    animalList{ii} = animalFolder{ii}(startInd+1:end);
end


for ii = 1:nFiles
    Ind = strfind(dataIndex.LogFilePath{ii},filesep);
    startInd = Ind(end);
    currAni = dataIndex.LogFilePath{ii}(startInd+1:end);
    Index = find(contains(animalList,currAni));
    animalIndex(ii) = Index;
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},[dataIndex.LogFileName{ii}(1:end-4),'_beh_cut.mat']));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    fn_latent = fullfile(dataIndex.BehPath{ii}, [dataIndex.LogFileName{ii}(1:end-4),'_FQRPECKlatentV.mat']);
    load(fn_latent);
    % load pupil files
    date = num2str(dataIndex.DateNumber(ii));
    pup_name = fullfile(dataIndex.BehPath{ii}, ['*',date(1:6),'*_pup.mat']);
    
    diffCK = stats_new.ckl-stats_new.ckr;
    
    % find the point where diffCK cross 0
    numX = 0;
    for tt = 1:length(diffCK)-1
        if diffCK(tt)*diffCK(tt+1) < 0
            numX = numX+1;
        end
    end
    numSwitch(ii) = numX;
end

nSwitchGroup = [];
for aa = 1:length(animalFolder)
    nSwitchGroup = [nSwitchGroup;numSwitch(animalIndex==aa)'];
end
figure;boxplot(numSwitch,animalIndex)


% reviewer figure 1
figure;
v=violinplot(numSwitch,animalIndex,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0],'ShowData',false);
v(1).ViolinColor = [6 6 122]/255;
v(2).ViolinColor = [2 124 246]/255;
v(3).ViolinColor = [108 206 108]/255;
v(4).ViolinColor = [245 145 1]/255;
v(5).ViolinColor = [122 1 1]/255;

set(gca, 'box','off');
ylabel('# of choice kernel switches','FontSize',28);
print(gcf,'-dpng',fullfile(save_path,'CKSwitches-violin'));    %png format
saveas(gcf, fullfile(save_path,'CKSwitches-violin'), 'fig');
saveas(gcf, fullfile(save_path,'CKSwitches-violin'),'svg');


% reviewer figure 1
% histogram of choice kernel difference (per animal)
edges=-1:0.02:1;
for ii = 1:length(animalList)
     filename2 = fullfile(model_path,[animalList{ii},'_FQRPECK.mat']);
     load(filename2);
     figure;histogram(stats_sim.ckl-stats_sim.ckr,edges,'EdgeColor','black','FaceColor','black','FaceAlpha',0.7);
     xlabel(['Subject ',animalList{ii}]);
     set(gca,'box','off')
     print(gcf,'-dpng',['CK-distribution-',animalList{ii}]);    %png format
    saveas(gcf, ['CK-distribution-',animalList{ii}], 'fig');
    saveas(gcf, ['CK-distribution-',animalList{ii}],'svg');
end

close all