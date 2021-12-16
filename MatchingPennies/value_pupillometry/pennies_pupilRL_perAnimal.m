function pennies_pupilRL_perAnimal(dataIndex, model_path, savefigpath, animalInd)

%% load the data

if animalInd ~= 0
    currAnimal = dataIndex.Animal{1};
else
    currAnimal = 'all';
end
% subMask = [];
savepath = fullfile(savefigpath,[currAnimal,'_fig-pupil']);
savematpath = fullfile(savepath, 'regcr_all.mat');

stats_all.c=[];  %concatenate choices and outcomes across sessions
stats_all.r=[];
% subMask = [];

for i=1:size(dataIndex,1)
    
    load(fullfile(dataIndex.BehPath{i},[dataIndex.LogFileName{i}(1:end-4),'_beh_cut.mat']));
    
    cueTimes{i} = trialData.cueTimes;
    outcomeTimes{i} = trialData.outcomeTimes;
    
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
        nlike_array trials fitpar_array...
        nTrial_array entro_array rrate_array subMask...
        stats_all runningEntropy savepath;
end

if sum(cellfun(@isempty,pupilAnimal.dia))< size(dataIndex,1)
    % load the model-fitting results
    modelFile = fullfile(model_path, 'FQ_RPE.mat');
    load(modelFile);
    
    % get simulated latent variables using fitpar
    FQpar = fitpar{animalInd};
    
    player1.label='algo_FQ_RPE';   % change to CA later
    player1.params.a=FQpar(1);
    player1.params.b=FQpar(2);
    
    stats_sim=predictAgent(player1,stats_all);
    
    % linear regression
    params=[];
    
    %first predictor is action value left
    params.trigEvent=stats_sim.ql;
    
    %second predictor is action value right
    params.trigEvent2=stats_sim.qr;
    
    %third predictor is rpe
    params.trigEvent3 = stats_sim.rpe;
    
    % when align pupil signal to cue
    params.trigTime = cueTimes;
    params.xtitle = 'Time from cue (s)';
    params.window = [-3:0.1:5];
    params.nback = 2;       %how many trials back to regress against
    params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
    params.interaction = true;
    %only perform analysis on this subset of trials
    fieldname={'go'};
    trialMask = getMask(trials,fieldname);
    reg_cr=linear_regr_perAnimal( pupilAnimal.dia, pupilAnimal.t, [params.trigEvent params.trigEvent2 params.trigEvent3], params.trigTime,  params );
    
    
    params.trigTime = outcomeTimes;
    params.xtitle= 'Time from outcome (s)';
    params.window = [-3:0.1:5];
    params.nback = 3;       %how many trials back to regress against
    params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
    params.interaction = true;
    %only perform analysis on this subset of trials
    fieldname={'go'};
    trialMask = getMask(trials,fieldname);
    reg_cr_outcome=linear_regr_perAnimal( pupilAnimal.dia, pupilAnimal.t, [params.trigEvent params.trigEvent2 params.trigEvent3], params.trigTime, params );
    
    tlabel={'Q_L','Q_R','RPE'};
    params.xtitle = {'Time from cue (s)','Time from outcome (s)'};
    
    if ~exist(savepath)
        mkdir(savepath)
    end
    cd(savepath)
    MP_plot_regrRL_pupil([reg_cr,reg_cr_outcome],params.pvalThresh,tlabel,params.xtitle);
    print(gcf,'-dpng','MLR-RL-perAnimal_cut');    %png format
    saveas(gcf, 'MLR-RL-perAnimal_cut', 'fig');
    close all;
    
end



end

        


