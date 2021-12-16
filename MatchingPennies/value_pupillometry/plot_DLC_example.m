% plot_DLC_example

figurePath = 'J:\figures\Training-M873_phase2_algo2_WithPupil_1911100627-img01694.png';
cd('J:\figures');
fig = imread(figurePath);

figure
image(fig)
set(gca,'XColor','none','YColor','none')
hold on;
plot([10,60],[10,10],'LineWidth',2,'color','black');
print(gcf,'-dpng','DLC_example');    %png format
saveas(gcf,'DLC_example', 'fig');
saveas(gcf, 'DLC_example','svg');
