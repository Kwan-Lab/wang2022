function plot_switch_random_distn(input,tlabel)
% % plot_switch_random_distn %
%PURPOSE:   A plot that would accompanies plot_switch_random(), to show the
%           distributions of L1 and L2
%AUTHORS:   H Atilgan, AC Kwan 191205
%
%INPUT ARGUMENTS
%   input:  Structure generated by choice_switch().
%   tlabel: Text string that will be put on top of the figure

%%
if numel(input)>1  %more than 1 data set, plot mean+-sem
    L1_ranges=input{1}.L1_ranges;
    L2_ranges=input{1}.L2_ranges;
    L1=[];
    L2=[];
    for j=1:numel(input)
        L1=[L1; input{j}.L1];
        L2=[L2; input{j}.L2];
    end
else                %plot the 1 data set
    L1_ranges=input.L1_ranges;
    L2_ranges=input.L2_ranges;
    L1=input.L1;
    L2=input.L2;
end

numRange = size(L2_ranges,1);

figure;

% plot distribution of L1
edges=[-0.5:1:100.5];

subplot(2,2,1); hold on;
x=edges+nanmean(diff(edges))/2; %plot at the center of the bins
n=histc(L1,edges);
bar(x,n,1,'FaceColor','k');
for k=1:numRange
    nn=histc(L1((L1 >= L1_ranges(k,1)) & (L1 <= L1_ranges(k,2))),edges);
    bar(x,nn,1,'FaceColor',[1 k/(numRange+1) 0]);
end
    
xlim([edges(1) edges(end)+1]);
ylim([0 1.1*nanmax(n)]);
title(tlabel,'interpreter','none');
ylabel('Number of blocks');
xlabel('L_1 (trial to meet criterion)');

% plot distribution of L2
edges=[-0.5:1:30.5];

subplot(2,2,2); hold on;
x=edges+nanmean(diff(edges))/2; %plot at the center of the bins
n=histc(L2,edges);
bar(x,n,1,'FaceColor','k');
for k=1:numRange
    nn=histc(L2((L2 >= L2_ranges(k,1)) & (L2 <= L2_ranges(k,2))),edges);
    bar(x,nn,1,'FaceColor',[1 k/(numRange+1) 0]);
end
xlim([edges(1) edges(end)+1]);
ylim([0 1.1*nanmax(n(:))]);
title(tlabel,'interpreter','none');
ylabel('Number of blocks');
xlabel('L_2 (random number of trials added)');

end


