function MP_pupilCoef_change(dataIndex, savefigpath)

% plot the coefficient of chosen q and R to determine whether the session
% encodes RPE or updated chosen Q
setup_figprop;
nFiles = size(dataIndex,1);

if ~exist([savefigpath,filesep,'coeff'])
    mkdir([savefigpath,filesep,'coeff'])
end
cd([savefigpath,filesep,'coeff'])
% try time = 2.45 second first
coeff_chosenQ = [];
coeff_R = [];
coeff_chosenQP = [];
coeff_RP = [];
subject_mask = [];
animalList = unique(dataIndex.Animal);

t_seq = 2:0.05:4;
%RPE_coeff = zeros(1, length(t_seq));
for ll = 1:length(t_seq)
    t_coeff = t_seq(ll);
    for ii = 1:nFiles
        fn_beh = dir(fullfile(dataIndex.BehPath{ii},[dataIndex.LogFileName{ii}(1:end-4),'_beh_cut.mat']));
        saveRegName_cue = fullfile(fn_beh.folder,'analysis-pupil',[fn_beh.name(1:end-7),'regRL_CK.mat']);
        
        animal = dataIndex.Animal{ii};
        
        
        if exist(saveRegName_cue)
            subject_mask = [subject_mask, find(contains(animalList, animal))];
            load(saveRegName_cue)
            % the R is the second factor, chosen q is the 5th factor
            coeff_chosenQ = [coeff_chosenQ;reg_cr_change.coeff(:,5)'];
            coeff_R =[coeff_R; reg_cr_change.coeff(:,2)'];
            coeff_chosenQP =[coeff_chosenQP; reg_cr_change.pval(:,5)'];
            coeff_RP = [coeff_RP,; reg_cr_change.pval(:,2)'];
        end
        
    end
    
    
    plot1 = coeff_chosenQ(:,(find((reg_cr_change.regr_time<t_coeff),1,'last')));
    pval1 = coeff_chosenQP(:,(find((reg_cr_change.regr_time<t_coeff),1,'last')));
    plot2 = coeff_R(:,(find((reg_cr_change.regr_time<t_coeff),1,'last')));
    pval2 = coeff_RP(:,(find((reg_cr_change.regr_time<t_coeff),1,'last')));
    figure;
    cmap = colormap('prism');
    for ii = 1:length(animalList)  % plot different animals in different color
        index_sig = (subject_mask==ii)' & pval1<0.05 & pval2<0.05;
        index_unsig = (subject_mask==ii)' & (pval1>0.05 | pval2>0.05);
        scatter(plot1(index_sig),plot2(index_sig),80,cmap(ii*2,:),'filled');
        hold on;
        scatter(plot1(index_unsig),plot2(index_unsig),80,cmap(ii*2,:));
        hold on;
    end
    ax = gca;
    xlabel('Coefficient of chosen value');
    ylabel('Coefficient of R(n)');
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    %legend(animalList)
    title(['Coefficient of R(n) and chosen Q at t = ',num2str(t_coeff)]);
    print(gcf,'-dpng',['coeff_r_chosenq_change_t=',num2str(t_coeff*100)]);    %png format
    saveas(gcf, ['coeff_r_chosenq_change_t=',num2str(t_coeff*100)], 'fig');
end

% calculate the percentage of the data points that fall in the second and
% fourth qua
t = -2.9:0.1:5;
figure;
for ll = 1:length(animalList)
    p24(:,ll) = mean((coeff_chosenQ(subject_mask==ll,:) >0 & coeff_R(subject_mask==ll,:)< 0) | (coeff_chosenQ(subject_mask==ll,:) < 0 & coeff_R(subject_mask==ll,:) > 0), 1);
    subplot(2,4,ll)
    plot(t, p24(:,ll))
    title(animalList{ll})
end
print(gcf,'-dpng','coeff_r_chosenq_change_time');    %png format
    saveas(gcf, 'coeff_r_chosenq_change_time', 'fig');

    % plot the total percentage
p24_whole = mean((coeff_chosenQ >=0 & coeff_R<= 0) | (coeff_chosenQ <= 0 & coeff_R >= 0), 1);
figure;
plot(t, p24_whole*100, 'k');
hold on;
plot([0 0],[0 100],'k--');
set(gca,'box','off') 
%title('Percentage of coefficients fall in 2 & 4')
xlabel('Time from cue (s)');
ylabel('Fraction of RPE-coding sessions (%)');
ylim([0,100]);
xlim([-3,5]);set(gca,'linewidth',5)
print(gcf,'-dpng','coeff_r_chosenq_change_time_total');    %png format
saveas(gcf, 'coeff_r_chosenq_change_time_total', 'fig');
saveas(gcf, 'coeff_r_chosenq_change_time_total', 'svg');


