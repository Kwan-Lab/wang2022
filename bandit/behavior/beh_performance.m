function output=beh_performance(stats)
% % beh_performance %
%PURPOSE:   Calculate basic performance metrics for a session
%AUTHORS:   AC Kwan 191210
%
%INPUT ARGUMENTS
%   stats:  stats of the task
%
%OUTPUT ARGUMENTS
%   output:     numbers used to plot figure
%
% To plot the output, use plot_behperf().

%%

output.nTrial = sum(stats.c==-1 | stats.c==1);
output.nLeftResp = sum(stats.c==-1);
output.nRightResp = sum(stats.c==1);

output.meanRewardRate = sum((stats.r>0) & (stats.c==-1 | stats.c==1)) / sum(stats.c==-1 | stats.c==1);

output.nSwitch = numel(stats.blockLength) - 1;

output.rule_labels = stats.rule_labels;

for k=1:numel(stats.rule_labels)
    output.meanTrialtoCrit(k) = nanmean(stats.blockTrialtoCrit(stats.blockRule(1:end-1) == k));
end

end


