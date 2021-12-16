function MP_pupilResp_latent(dataIndex,save_path_pupil,model_path)


nFiles = size(dataIndex,1);
Animal = unique(dataIndex.Animal);
  disp('-----------------------------------------------------------');
    disp(['--- Fitting models - summary of ', int2str(numel(Animal)) ' animals']);
    disp('-----------------------------------------------------------');
    
    fileInd = 1:nFiles;
for aa = 1:length(Animal)
  
    currAnimalSessions = contains(dataIndex.LogFilePath,Animal{aa});
    
    animalpupilResp = [];
    animaldQ = []; animaldK = [];
    for ii = fileInd(currAnimalSessions==1)
    % load behavior files
        fn_beh = dir(fullfile(dataIndex.BehPath{ii},[dataIndex.LogFileName{ii}(1:end-4),'_beh_cut.mat']));
        load(fullfile(fn_beh.folder,fn_beh.name));
        
        fn_latent = fullfile(dataIndex.BehPath{ii}, [dataIndex.LogFileName{ii}(1:end-4),'_FQRPECKlatentV.mat']);
        load(fn_latent);

    % load pupil files
        date = num2str(dataIndex.DateNumber(ii));
        pup_name = fullfile(dataIndex.BehPath{ii}, ['*',date(1:6),'*_pup.mat']);
        fn_pup = dir(pup_name);
        if length(fn_pup) == 1
        
        savepupilfigpath = fullfile(dataIndex.BehPath{ii},[date(1:6),'_figs-pupil']);
        if ~exist(savepupilfigpath,'dir')
            mkdir(savepupilfigpath);
        end
        
        cd(savepupilfigpath);
        
    
        load(fullfile(fn_pup.folder,fn_pup.name));
        
        % calculate the average pupil response
        pResp = zeros(1, size(pupil.resp,1));
        for tt = 1:size(pupil.resp,1)
            respInd = pupil.respT(tt,:)<trialData.cueTimes(tt)+3 & pupil.respT(tt,:)>trialData.cueTimes(tt)+2;
            pResp(tt) = nanmean(pupil.resp(respInd));
        end
        
        if size(pupil.resp,1) == length(stats_new.ql)
            animalpupilResp = [animalpupilResp,pResp];
            animaldQ = [animaldQ, stats_new.ql'-stats_new.qr'];
            animaldK = [animaldK, stats_new.ckl'-stats_new.ckr'];
        elseif size(pupil.resp,1) < length(stats_new.ql)
            animalpupilResp = [animalpupilResp,pResp];
            animaldQ = [animaldQ, stats_new.ql(1:size(pupil.resp,1))'-stats_new.qr(1:size(pupil.resp,1))'];
            animaldK = [animaldK, stats_new.ckl(1:size(pupil.resp,1))'-stats_new.ckr(1:size(pupil.resp,1))'];
        else
            animalpupilResp = [animalpupilResp,pResp(1:length(stats_new.ql))];
            animaldQ = [animaldQ, stats_new.ql'-stats_new.qr'];
            animaldK = [animaldK, stats_new.ckl'-stats_new.ckr'];
        end
        end
    end
     
    % minimize the mean square error
    % softmax
    maxit=1e6;
    maxeval=1e6;
    op=optimset('fminsearch');
    op.MaxIter=maxit;
    op.MaxFunEvals=maxeval;
    
    initpar = [0.1 0.1]; %  beta  beta_K
    lb = [0 0];
    ub = [inf inf];
    
    dat = [animaldQ; animaldK; animalpupilResp];
    func_handle = @softmaxLatent;
    [qpar{aa}, mse{aa}, exitflag]=fmincon(func_handle, initpar, [], [], [], [], lb, ub, [], op, dat);
     
end
save(fullfile(model_path, ['pupilResp_model.mat']),...
           'Animal', 'qpar','mse');