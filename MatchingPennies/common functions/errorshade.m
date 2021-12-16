function errorshade(x,h,l,c,alpha)
% % errorshade %
%PURPOSE:   Plot shaded area, e.g. for confidence intervals
%AUTHORS:   AC Kwan 170515
%
%INPUT ARGUMENTS
%   x:  the independent variable
%   h:  the upper bound for the dependent variable
%   l:  the lower bound for the dependent variable
%   c:  the color of the shading, e.g., [0.7 0.7 0.7] for gray
%   alpha:  the transparency of the shading (optional, default = 1)
%%

if nargin > 4
  facealpha = alpha;
else
  facealpha = 1;
end

x=x(:)';
h=h(:)';
l=l(:)';

h=fill([x fliplr(x)],[h fliplr(l)],c,'LineStyle','none','FaceAlpha',facealpha);