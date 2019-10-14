function benchRandom

addpath('../matlabHelpers/')
% create figure to show data
figure(1);
clf;
hold on;

legendStr={};

dataFile = loadHDF5('/Users/gonciarz/allRandomGraphs_new4.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'random'];


figure(1);
hold on

% size(data.randomTime)
% data.randomTime = diff(data.randomTime);
% data.random=data.random(2:end)
% data.exact=data.exact(2:end)
% data.exactTime=data.exactTime(2:end)
% data.gr=data.gr(2:end)
% data.grTime=data.grTime(2:end)
% size(data.randomTime)

s = data.random ./ data.exact;
[min(s), mean(s), max(s)] %print range of solutions
h1 = histogram(s);
hold on;
d = data.gr ./ data.exact;
[min(d), mean(d), max(d)]
h2 = histogram(d);

data

[fs, idx] = sort(data.exact);
et=data.exactTime(idx);
rt=data.randomTime(idx);

et(et==0)=0.00001;

size(fs);
[vals, startIdx] = unique(fs);

times = [];
dev=[];
for i = 1:size(startIdx, 1)
    si = startIdx(i);
    
    if i == size(startIdx, 1)
        ei = size(et, 1);
    else
        ei = startIdx(i + 1) - 1;
    end

    m = mean(rt(si:ei)./et(si:ei));
    d = nanstd(rt(si:ei)./et(si:ei));
    [m, d];
    times = [times, m];
    dev=[dev, d];
end

figure(3)
clf;
plot(vals, times')
hold on
p = polyfit(vals(1:end), times(1:end)',6);
p
f1 = polyval(p, vals(1:end)');
plot(vals(1:end), f1', 'r')


% errorbar(vals, times', dev')

% figure(1)
% errorRandom = data.random ./ data.exact;
% hist(errorRandom)
% figure(2)
% errorGR = data.gr ./ data.exact;
% hist(errorGR)
% 
% plotOne(data.vertices, [data.randomTime]);
% 
% 
% l = legend(legendStr);
% l.FontSize = 10;
% title('Time benchmark');
% xlabel('#v');
% ylabel('time(s)');
% 
% set(gca,'YScale','log')
% 
% function plotOne(dataRef, dataIn)
%     xVals = unique(dataRef)'
% 
%     output result with size (numberOfDataPoints, numberOfDataSets)
%     average = zeros(size(xVals, 2), size(dataIn, 2));
% 
%     Go through each data column separately
%     for inputData = 1:size(dataIn, 2)
%         currData=dataIn(:, inputData);
% 
%         Calculate mean per each x value
%         averagec = zeros(length(xVals), 1);
%         idx = 1;
%         for i = xVals
%             averagec(idx) = mean(currData(dataRef == i));
%             idx = idx + 1;
%         end
%         average(:,inputData) = averagec;
%     end
% 
%     for i = 1:size(average, 2)
%         plot(xVals, average(:,i), '-*', 'LineWidth', 3);
%         f=fit(xVals',average(:,i), 'poly3')
%         plot(f, xVals, average(:,i));
%     end
% end

end