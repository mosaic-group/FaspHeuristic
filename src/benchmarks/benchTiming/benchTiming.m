function benchTiming

addpath('../matlabHelpers/')
% create figure to show data
figure(1);
clf;
hold on;

legendStr={};

dataFile = loadHDF5('data/TimingConstWeightVarFaspConstVE_v_30_e_120_f_1-20_s_20_lin__r_200.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'v=30 e=120 gr', 'v=30 e=120 delta', 'v=30 e=120 random'];
plotOne(data.exact, [data.grTime, data.deltaTime, data.randomTime]);
dataFile = loadHDF5('data/TimingConstWeightVarFaspConstVE_v_40_e_120_f_1-20_s_20_lin__r_200.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'v=40 e=120 gr', 'v=40 e=120 delta', 'v=40 e=120 random'];
plotOne(data.exact, [data.grTime, data.deltaTime, data.randomTime]);
dataFile = loadHDF5('data/TimingConstWeightVarFaspConstVE_v_60_e_120_f_1-20_s_20_lin__r_200.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'v=60 e=120 gr', 'v=60 e=120 delta', 'v=60 e=120 random'];
plotOne(data.exact, [data.grTime, data.deltaTime, data.randomTime]);


l = legend(legendStr);
l.FontSize = 10;
title('Time benchmark');
xlabel('FASP size');
ylabel('time(s)');


function plotOne(dataRef, dataIn)
    xVals = unique(dataRef)';

    % output result with size (numberOfDataPoints, numberOfDataSets)
    average = zeros(size(xVals, 2), size(dataIn, 2));

    % Go through each data column separately
    for inputData = 1:size(dataIn, 2)
        currData=dataIn(:, inputData);

        % Calculate mean per each x value
        averagec = zeros(length(xVals), 1);
        idx = 1;
        for i = xVals
            averagec(idx) = mean(currData(dataRef == i));
            idx = idx + 1;
        end
        average(:,inputData) = averagec;
    end

    for i = 1:size(average, 2)
        plot(xVals, average(:,i), '-*', 'LineWidth', 4);
%         f=fit(xVals',average(:,i), 'poly3')
%         plot(f, xVals, average(:,i));
    end
end

end