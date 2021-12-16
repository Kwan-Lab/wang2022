x = linspace(0, 100, 100)';
M = zeros(100,3);
for tt = 1:100
    M(tt,:) = [1, 1-x(tt)/100, 1-x(tt)/100];
end

% Define the vertices: the points at (x, f(x)) and (x, 0)
N = length(x);
verts = [x(:), zeros(N,1); x(:) zeros(N,1)+1];
% Define the faces to connect each adjacent f(x) and the corresponding points at y = 0.
q = (1:N-1)';
faces = [q, q+1, q+N+1, q+N];
figure;
p = patch('Faces', faces, 'Vertices', verts, 'FaceVertexCData', [M;M], 'FaceColor', 'interp', 'EdgeColor', 'none');
set(gca,'visible','off')

% calculate three p values: p = 0.001, p = 0.01, p = 0.1;
pvalue = [1e-20, 1e-10, 0.001, 1];
pstring = {'1e-20';'1e-10';'0.001';'0.05'};
xvalue = 100*(log(pvalue)/log(1e-20));
for tt = 1:length(pvalue)
    hold on;plot([xvalue(tt),xvalue(tt)],[-0.05,0],'k','LineWidth',2);
    text(xvalue(tt),-0.1+(tt-1)*(-0.05), num2str(pstring{tt}));
end

savefigpath = 'J:\MatchingPennies\pupilData\summary';
cd(savefigpath);
print(gcf,'-dpng','color_map_log');    %png format
saveas(gcf, 'color_map_log', 'fig');
saveas(gcf, 'color_map_log','svg');
