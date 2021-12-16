function plot_rt_latentV(stats_sim, trueRespTime)

% plot the correlation between response time and latent variables, by
% sessions

qSum = stats_sim.ql + stats_sim.qr;
qDif = stats_sim.ql - stats_sim.qr;

% seperate the latent variable into sessions
Ind = 1;
for ii = 1:length(trueRespTime)
    qSum_session{ii} = qSum(Ind:Ind+length(trueRespTime{ii})-1);
    qDif_session{ii} = qDif(Ind:Ind+length(trueRespTime{ii})-1);
    Ind = Ind + length(trueRespTime{ii});
end






%% plot the sum of action-value and the response time correlation
corrSum = zeros(1, length(qSum_session));
pCorrSum = zeros(1, length(qSum_session));

for tt = 1:length(qSum_session)
    noNan=~isnan(trueRespTime{tt});
    [R,P] = corrcoef(abs(qSum_session{tt}(noNan)), trueRespTime{tt}(noNan)); 
    corrSum(tt) = R(2,1);
    pCorrSum(tt) = P(2,1);
end

binSize = 0.1;
corrBin_total = zeros(1, length([-1:0.1:1]));
corrBin_sigtotal = zeros(1, length([-1:0.1:1]));
for ii = 1:length(corrSum)
    corrBin_total(ceil((corrSum(ii)+1)/0.1)) = corrBin_total(ceil((corrSum(ii)+1)/0.1)) + 1;
    if pCorrSum(ii) < 0.05
        corrBin_sigtotal(ceil((corrSum(ii)+1)/0.1)) =  corrBin_sigtotal(ceil((corrSum(ii)+1)/0.1)) + 1;
    end
end

figure;
bar([-1:0.1:1], corrBin_total, 'FaceColor', 'black', 'FaceAlpha',0.3);
hold on; bar(-1:0.1:1, corrBin_sigtotal, 'FaceColor','Black','FaceAlpha',0.8);
title('Correlations between total value and response time');
xlabel('Correlation coefficients');
ylabel('# of sessions')
print(gcf,'-dpng','pos_neg_sessions_total_cut');    %png format
saveas(gcf, 'pos_neg_sessions_total_cut', 'fig');


%% plot the difference of action-value and the response time correlation
corrDif = zeros(1, length(qDif_session));
pCorrDif = zeros(1, length(qDif_session));
for tt = 1:length(qDif_session)
    noNan=~isnan(trueRespTime{tt});
    [R,P] = corrcoef(abs(qDif_session{tt}(noNan)), trueRespTime{tt}(noNan)); 
    corrDif(tt) = R(2,1);
    pCorrDif(tt) = P(2,1);
end

sessionPos_dif = corrDif > 0;
sessionNeg_dif = corrDif < -0;
sessionPos_difsig = corrDif > 0 & pCorrDif < 0.05;
sessionNeg_difsig = corrDif < 0 & pCorrDif < 0.05;


corrBin_total = zeros(1, length([-1:0.1:1]));
corrBin_sigtotal = zeros(1, length([-1:0.1:1]));
for ii = 1:length(corrDif)
    corrBin_total(ceil((corrDif(ii)+1)/0.1)) = corrBin_total(ceil((corrDif(ii)+1)/0.1)) + 1;
    if pCorrDif(ii) < 0.05
        corrBin_sigtotal(ceil((corrDif(ii)+1)/0.1)) =  corrBin_sigtotal(ceil((corrDif(ii)+1)/0.1)) + 1;
    end
end

figure;
bar([-1:0.1:1], corrBin_total, 'FaceColor', 'black', 'FaceAlpha',0.3);
hold on; bar(-1:0.1:1, corrBin_sigtotal, 'FaceColor','Black','FaceAlpha',0.8);
title('Correlations between absolute value difference and response time');
xlabel('Correlation coefficients');
ylabel('# of sessions')
print(gcf,'-dpng','pos_neg_sessions_dif_cut');    %png format
saveas(gcf, 'pos_neg_sessions_total_cut', 'fig');


%% compare average response time in first 50 trials and last 50 trials in a session
sessionLength = cellfun(@length,trueRespTime);
meanRTStart = zeros(1, sum(sessionLength>300));
meanRTEnd = zeros(1, sum(sessionLength>300));
oInd = 1:length(trueRespTime);
Ind = oInd(sessionLength>300);
for ii=1:length(Ind)
    rt = trueRespTime{Ind(ii)}(~isnan(trueRespTime{Ind(ii)}));
    meanRTStart(ii) = mean(rt(1:50));
    meanRTEnd(ii) = mean(rt(end-49:end));
end

[h1,p1] = ttest(meanRTStart, meanRTEnd, 'tail','left')
meanStart = mean(meanRTStart);
meanEnd = mean(meanRTEnd);
stdStart =std(meanRTStart);
stdEnd = std(meanRTEnd);

diff = meanRTEnd - meanRTStart;

figure;
subplot(1,2,1) % box plot
boxplot(diff, 'Notch','on');
xlabel('End of session - start of session');
ylabel('Response time difference (s)');
title('Response time difference between session end and session start');
subplot(1,2,2) % bar plot
bar([meanStart, meanEnd],0.6,'black');
set(gca,'xticklabel', {'First 50 trials', 'Last 50 trials'});
hold on;
er = errorbar([meanStart,meanEnd], [stdStart,stdEnd]);
er.Color = [0 0 0];
er.LineStyle = 'none';
ylim([0 0.5])
title('Averaged response time');
print(gcf,'-dpng','rt_stage_cut');    %png format
saveas(gcf, 'st_stage', 'fig');
print(gcf,'-dpng','rt_box_cut');    %png format
saveas(gcf, 'rt_box', 'fig');


end