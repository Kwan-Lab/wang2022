function plot_binxaveragey(input,tlabel)
% % plot_binxaveragey %
%PURPOSE:   Bin x and average y
%AUTHORS:   AC Kwan 191211
%
%INPUT ARGUMENTS
%   input:      Paired data

%%
% if called from single-session analysis, input is a struct
% if called from summary analysis, input is a cell array
% here convert everything to a cell array first
if ~iscell(input)
    temp = input;
    clear input;
    input{1} = temp;
end

label=input{1}.label;
xrange=input{1}.range{1};
yrange=input{1}.range{2};

val=[];
for j=1:numel(input)
    val=[val; input{j}.dat];
end

%% plot

figure;
    
subplot(2,2,1); hold on;
for j=xrange(1):1:xrange(2)
    y=val(val(:,1)==j,2);
    plot(j,nanmean(y),'k.','MarkerSize',30);
    plot([j j],nanmean(y)+nanstd(y)/sqrt(numel(y))*[-1 1],'k-','LineWidth',3);
end
xlim([xrange(1)-1 xrange(2)+1]);
ylim([yrange(1) yrange(2)]);
xlabel(label{1});
ylabel(label{2});
title(tlabel,'interpreter','none');

end


