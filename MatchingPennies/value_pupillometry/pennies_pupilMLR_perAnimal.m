function pennies_pupilMLR_perAnimal(dataIndex, model_path, savefigpath, animalInd)

%% load the data

stats_all.c=[];  %concatenate choices and outcomes across sessions
stats_all.r=[];

if animalInd ~= 0
    currAnimal = dataIndex.Animal{1};
else
    currAnimal = 'all';
end
% subMask = [];
savepath = fullfile(savefigpath,[currAnimal,'_fig-pupil']);
savematpath = fullfile(savepath, 'regcr_all.mat');

for i=1:size(dataIndex,1)
    
    load(fullfile(dataIndex.BehPath{i},[dataIndex.LogFileName{i}(1:end-4),'_beh.mat']));
    
    cueTimes{i} = trialData.cueTimes;
    % outcomeTimes{i} = trialData.outcomeTimes;
    % Trials{i} = trials;
    % respTime_array{i}=respTime_trType;
    % trueRespTime{i} = trialData.rt;
    %  hoiceBySession{i} = stats;
    
    % runningEntropy{i} = trialData.runningEntropy;
    
    stats_all.c=[stats_all.c; stats.c];
    stats_all.r=[stats_all.r; stats.r];
    
    % load pupildata
    if dataIndex.pupil(i) == 1
        % if pupil data exist
        date = num2str(dataIndex.DateNumber(i));
        pupFile = fullfile(dataIndex.BehPath{i}, ['*',date(1:6),'*_pup.mat']);
        fn_pup = dir(pupFile);
        if ~isempty(fn_pup)
            load(fullfile(fn_pup.folder,fn_pup.name));
            pupilAnimal.dia{i} = pupil.dia;
            pupilAnimal.t{i} = pupil.t;
        else
            pupilAnimal.dia{i} = [];
            pupilAnimal.t{i} = [];
        end
    else
        pupilAnimal.dia{i} = [];
        pupilAnimal.t{i} = [];
    end
    
    close all;
    clearvars -except i dirs dataIndex model_path savefigpath animalInd pupilAnimal...
        lick_trType_array cueTimes outcomeTimes iti_array respTime_array trueRespTime choiceBySession...
        nlike_array Trials fitpar_array...
        nTrial_array currAnimal rrate_array subMask...
        stats_all runningEntropy savepath;
end

if sum(cellfun(@isempty,pupilAnimal.dia))< size(dataIndex,1)
    % load the model-fitting results
   
    % linear regression
    params=[];
    
    %first predictor is choice. dummy code: left:0, right:1
    stats_all.c(stats_all.c(:,1)==-1,1) = 0;  % change left from -1 to 0
    params.trigEvent=stats_all.c(:,1);
    
    %second predictor is action value right
    params.trigEvent2=stats_all.r;

    
    % when align pupil signal to cue
    params.trigTime = cueTimes;
    params.xtitle = 'Time from cue (s)';
    params.window = [-3:0.1:8];
    params.nback = 2;       %how many trials back to regress against
    params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
    params.interaction = true;
    %only perform analysis on this subset of trials

    reg_cr=linear_regr_perAnimal( pupilAnimal.dia, pupilAnimal.t, [params.trigEvent params.trigEvent2], params.trigTime,  params );
    
    save(savematpath,'reg_cr');
    
    reg_cr.interaction = 1;
    if ~exist(savepath)
        mkdir(savepath)
    end
    cd(savepath)
    tlabel={'C(n)','C(n-1)','C(n-2)','R(n)','R(n-1)','R(n-2)','C(n)xR(n)','C(n-1)xR(n-1)','C(n-2)xR(n-2)'};
    MP_plot_regrcoef_pupil(reg_cr,params.pvalThresh,tlabel,params.xtitle);
    if currAnimal == 'all'
        print(gcf,'-dpng','MLR-interaction-all');    %png format
        saveas(gcf, 'MLR-interaction-all', 'fig');
    else
        print(gcf,'-dpng','MLR-interaction');    %png format
        saveas(gcf, 'MLR-interaction', 'fig');
    end
    close all;
    
end



end

        


