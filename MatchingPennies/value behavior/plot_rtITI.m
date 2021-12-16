function plot_rtITI(ITI_array, rt_array)

% plot the response time against previous ITI

% concatenate them first
ITI_vector = [];
rt_vector = [];
for ii = 1:length(ITI_array)
    ITI_vector = [ITI_vector, ITI_array{ii}];
    rt_vector = [rt_vector;rt_array{ii}(2:end)];
end

startPoint = floor(min(ITI_vector));
endPoint = ceil(max(ITI_vector));
windows = 0.5;

% calculate the bin
ITI_interval = startPoint:windows:endPoint-windows;
for kk = 1:length(ITI_interval)
    if kk < length(ITI_interval)
        rt_bin{kk} = rt_vector(ITI_vector > ITI_interval(kk) & ITI_vector < ITI_interval(kk+1));
    else
        rt_bin{kk} = rt_vector(ITI_vector > ITI_interval(kk));
    end
end

% get the mean and standard deviation
meanRT_bin = cellfun(@nanmean,rt_bin);
stdRT_bin = cellfun(@nanstd,rt_bin);

figure;
plot(ITI_interval, meanRT_bin,'black')
hold on;
%errorbar(ITI_interval, meanRT_bin, stdRT_bin);
errorshade(ITI_interval(1:22), meanRT_bin(1:22)+stdRT_bin(1:22), meanRT_bin(1:22)-stdRT_bin(1:22),[0.7 0.7 0.7]);
hold on;
plot(ITI_interval, meanRT_bin,'black')
xlabel('Intertrial interval(s)');
ylabel('Response time(s)');

