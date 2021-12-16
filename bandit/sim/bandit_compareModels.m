function bandit_compareModels(model_path)
% % bandit_compareModels %
%PURPOSE:   Compare BIC values for all the models fitted
%AUTHORS:   H Atilgan and AC Kwan 191212
%
%INPUT ARGUMENTS
%   model_path:    path where all the model analyses are saved
%
%OUTPUT ARGUMENTS
%

%% find all the models and BIC values
AllModelfiles = dir(fullfile(model_path,'*_model.mat'));

nFile = numel(AllModelfiles);

for k = 1:nFile
    
    load(fullfile(AllModelfiles(k).folder,AllModelfiles(k).name));
    
    models(k).name = AllModelfiles(k).name(1:end-4);
    models(k).bic = cell2mat(bic');
    
end

%% Plot
figure; hold on;

subplot(2,1,1); hold on;

x=[];
y=[];

%if Q_RPE exist, use it as the reference
idx = ismember({models.name},'Q_RPE');
if any(idx) == 0
    idx = 1;
end

for k = 1:nFile
    x = [x; k*ones(size(models(k).bic))];
    y = [y; models(k).bic - models(idx).bic];  %calculate everything relative to reference model (default Q_RPE)
    label{k} = models(k).name;
end

if nFile > 50 %violin plot if a lot of data points
    violinplot(y,x,'ViolinColor',[0 0 0],'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0],'ShowData',false);
else          %beeswarm plot works well for fewer number of points
    % beeswarm(x,y,'rand','none',3); hold on;
    
    % box plot
    h = boxplot(y,x,'Colors','k','Notch','on','Symbol',' ');
    set(h,{'linew'},{2})
end

plot([0 nFile+1],[0 0],'k--','LineWidth',3);

set(gca,'TickLabelInterpreter','none');
ylabel(['BIC (relative to ' models(idx).name ')'],'interpreter','none');
set(gca,'xtick',[1:nFile]);
set(gca,'xticklabel',label);
xlim([0 nFile+1]);
ylim([nanmin(y)-100 nanmax(y)+100]);

print(gcf,'-dpng',fullfile(model_path,'modelcomparison-bic'));
print(gcf,'-dsvg',fullfile(model_path,'modelcomparison-bic'));
saveas(gcf, fullfile(model_path,'modelcomparison-bic'), 'fig');