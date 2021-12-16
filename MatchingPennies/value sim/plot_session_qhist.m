function plot_session_qhist(stats,x)

%action value difference
qdiff = [];
stats_all.c = [];

for ii = 1:length(stats)
    if contains(stats{ii}.playerlabel{x},'CK')
        % if there is choice autocorrelation, plot the histogram with beta
        qdiff_t = stats{ii}.playerparams{x}.b*(stats{ii}.ql(:,x)-stats{ii}.qr(:,x))+stats{ii}.playerparams{1}.bc*(stats{ii}.ckl(:,x)-stats{ii}.ckr(:,x));
    else
        qdiff_t=stats{ii}.ql(:,x)-stats{ii}.qr(:,x);
    end
    qdiff = [qdiff;qdiff_t];
    stats_all.c = [stats_all.c; stats{ii}.c];
    
end
%plot histogram of action value differences
maxVal = ceil(nanmax(qdiff));
minVal = floor(nanmin(qdiff));
% find the one with larger absolute value
limit = max(abs(maxVal),abs(minVal));
edges=[-limit:0.05:limit];
edges_center=edges(1:end-1)+mean(diff(edges))/2; %center of the bins, for plotting

n=histcounts(qdiff,edges);
%n=n(1:end-1);

probL=nan(numel(edges_center),1);
for j=1:numel(edges_center)
    numL=sum(stats_all.c(qdiff>=edges(j) & qdiff<edges(j+1))==-1);
    numR=sum(stats_all.c(qdiff>=edges(j) & qdiff<edges(j+1))==1);
    probL(j)=numL/(numL+numR);
end


fig=figure;
handle = subplot(2,2,1); hold on;
left_color = [0 0 0];
right_color = [1 0 1];

yyaxis left


bar(edges_center,n,'k');
if contains(stats{1}.playerlabel{1},'CK')
    xlabel('\beta\DeltaQ + \beta_K\DeltaK');
else
    xlabel('Q_L - Q_R');
end
ylabel('Occurrence');
ylim([0 1.1*max(n)]);
xlim([edges(1) edges(end)]);
title(['Player ' int2str(x) ': ' stats{1}.playerlabel{x}],'interpreter','none');

%plot probability of choosing left as a function of action value diff
yyaxis right


plot(edges_center,probL,'m.','MarkerSize',35);
if contains(stats{1}.playerlabel{1},'RPE') %if softmax rule was used, what is the expected curve?
    if contains(stats{1}.playerlabel{1},'CK')
        exp_probL=1./(1+exp(-(edges_center)));
    else
        exp_probL=1./(1+exp(-stats{1}.playerparams{x}.b*(edges_center)));    
    end
    plot(edges_center,exp_probL,'m--','LineWidth',3);
elseif strcmp(stats{1}.playerlabel{1},'algo_DFQ') %if softmax rule was used, what is the expected curve?
    exp_probL=1./(1+exp(-1*edges_center));    
    plot(edges_center,exp_probL,'m--','LineWidth',3);
end
ylabel('P_L');
ylim([-0.02 1.02]);
handle.YAxis(1).Color = 'k';
handle.YAxis(2).Color = 'm';


if contains(stats{1}.playerlabel{1},'CK')
    print(gcf,'-dpng',['session-qhist-player-CK' int2str(x)]);    %png format
    saveas(gcf,['session-qhist-player-CK' int2str(x)], 'fig');
    saveas(gcf, ['session-qhist-player-CK' int2str(x)],'svg');
else
    print(gcf,'-dpng',['session-qhist-player-FQ' int2str(x)]);    %png format
    saveas(gcf,['session-qhist-player-FQ' int2str(x)], 'fig');
    saveas(gcf, ['session-qhist-player-FQ' int2str(x)],'svg');
end