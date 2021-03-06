function MP_plot_lickrate_byTrialType(input)
% % plot_lickrate_byTrialType %
%PURPOSE:   Plot lick rates for different trial types
%AUTHORS:   AC Kwan 170518
%
%INPUT ARGUMENTS
%   input:        Structure generated by get_lickrate_byTrialType().
%   tlabel:       Text to put as title of the plot.

%%
% if called from single-session analysis, input is a struct
% if called from summary analysis, input is a cell array
% here convert everything to a cell array first
if ~iscell(input)
    temp = input;
    clear input;
    input{1} = temp;
end

% load from the cell array
edges=input{1}.edges;
trialType=input{1}.trialType;

for j=1:numel(input)
    for l=1:numel(trialType)
        if j==1     %first time, load the array
            temp_lTimes{l}=input{j}.leftTimes{l};
            temp_rTimes{l}=input{j}.rightTimes{l};
        else        %otherwise, append
            temp_lTimes{l}=[(leftTimes{l}) input{j}.leftTimes{l}];
            temp_rTimes{l}=[(rightTimes{l}) input{j}.rightTimes{l}];
        end
    end
    leftTimes=temp_lTimes; clear temp_lTimes;
    rightTimes=temp_rTimes; clear temp_rTimes;
end

%% calculate mean and sem
edges=edges(1:end-1)+nanmean(diff(edges))/2;   %plot using the center of the histogram bins

lTimes=nan(numel(edges),numel(trialType)); % of choosing left
rTimes=nan(numel(edges),numel(trialType)); % of choosing right
for j=1:numel(trialType)
    lTimes(:,j)=nanmean(leftTimes{j},2);
    rTimes(:,j)=nanmean(rightTimes{j},2);
    
    lTimes_sem(:,j)=nanstd(leftTimes{j},[],2)./sqrt(numel(input));
    rTimes_sem(:,j)=nanstd(rightTimes{j},[],2)./sqrt(numel(input));
end

%% plot
figure;

panel_h=numel(trialType);

for j=1:numel(trialType)
    
    subplot(2,panel_h,j); hold on;
    plot(edges,lTimes(:,j),'-','Linewidth',3,'Color','r');
    plot([0 0],[0 10],'k--','LineWidth',2);
    axis([edges(1) edges(end) 0 10]);
    title({trialType{j};'Left lick'},'interpreter','none');
    if j==1 
        ylabel('Lick density (Hz)');
    end
    
    subplot(2,panel_h,j+numel(trialType)); hold on;
    plot(edges,rTimes(:,j),'-','Linewidth',3,'Color','b');
    plot([0 0],[0 10],'k--','LineWidth',2);
    title({trialType{j};'Right lick'},'interpreter','none');
    axis([edges(1) edges(end) 0 10]);
    if j==1
        ylabel('Lick density (Hz)');
    end
    ylabel('Lick density (Hz)');
    xlabel('Time from Cue (s)')
end

print(gcf,'-dpng','lickrates_byTrialType');    %png format
saveas(gcf, 'lickrates_byTrialType', 'fig');
saveas(gcf, 'lickrates_byTrialType','svg');
end


