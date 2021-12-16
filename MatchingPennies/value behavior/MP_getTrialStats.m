function stats = MP_getTrialStats(trials);

%% setting up data for plotting
    
    stats.playerlabel{1} = 'Mouse';
    stats.playerlabel{2} = 'Algo2';
    
    %choice: left=-1; right=1; miss=NaN
    stats.c=nan(length(trials.go),2);
    stats.c(trials.left==1,1)=-1;
    stats.c(trials.right==1,1)=1;
    stats.c(trials.compleft==1,2)=-1;
    stats.c(trials.compright==1,2)=1;
 
    %outcome: reward=1; no reward:0; miss=NaN
    stats.r=nan(length(trials.go),1);
    stats.r(trials.reward==1)=1;
    stats.r(trials.noreward==1)=0;
    
end
    