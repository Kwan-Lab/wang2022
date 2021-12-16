function mse_pupilResp = softmaxLatent(xpar,dat)

% xpar: beta; beta_K
% dat: pupilResp, dQ, dK

beta = xpar(1); beta_K = xpar(2);
predPupilResp = 1./(1+exp(-(beta*dat(1,:)+beta_K*dat(2,:))));

mse_pupilResp = sum((dat(3,:)-predPupilResp).^2);
