function benchHeuristics

addpath('../matlabHelpers/')
% create figure to show data
figure(1);
clf;
hold on;

legendStr={};

% dataFile = loadHDF5('data/HeuristicsConstWeightVarFaspConstVE_v_15_e_30_f_1-15_s_15_log__r_1000.h5');
% data = dataFile.Analysis_data; legendStr=[legendStr, 'v=15 e=30 gr', 'v=15 e=30 delta', 'v=15 e=30 random'];
% plotOne(data.exact, [data.gr, data.delta, data.random]);
% dataFile = loadHDF5('data/HeuristicsConstWeightVarFaspConstVE_v_15_e_55_f_1-15_s_15_log__r_1000.h5');
% data = dataFile.Analysis_data; legendStr=[legendStr, 'v=15 e=55 gr', 'v=15 e=55 delta', 'v=15 e=55 random'];
% plotOne(data.exact, [data.gr, data.delta, data.random]);
% dataFile = loadHDF5('data/HeuristicsConstWeightVarFaspConstVE_v_15_e_80_f_1-15_s_15_log__r_1000.h5');
% data = dataFile.Analysis_data; legendStr=[legendStr, 'v=15 e=80 gr', 'v=15 e=80 delta', 'v=15 e=80 random'];
% plotOne(data.exact, [data.gr, data.delta, data.random]);
% dataFile = loadHDF5('data/HeuristicsConstWeightVarFaspConstVE_v_15_e_105_f_1-15_s_15_log__r_1000.h5');
% data = dataFile.Analysis_data; legendStr=[legendStr, 'v=15 e=105 gr', 'v=15 e=105 delta', 'v=15 e=105 random'];
% plotOne(data.exact, [data.gr, data.delta, data.random]);


% 
% dataFile = loadHDF5('data/HeuristicsConstWeightVarFaspConstVE_v_30_e_60_f_1-30_s_20_log__r_1000.h5');
% data = dataFile.Analysis_data; legendStr=[legendStr, 'v=30 e=60 gr', 'v=30 e=60 delta', 'v=30 e=60 random'];
% plotOne(data.exact, [data.gr, data.delta, data.random]);
% dataFile = loadHDF5('data/HeuristicsConstWeightVarFaspConstVE_v_30_e_120_f_1-30_s_20_log__r_1000.h5');
% data = dataFile.Analysis_data; legendStr=[legendStr, 'v=30 e=120 gr', 'v=30 e=120 delta', 'v=30 e=120 random'];
% plotOne(data.exact, [data.gr, data.delta, data.random]);
% dataFile = loadHDF5('data/HeuristicsConstWeightVarFaspConstVE_v_30_e_180_f_1-30_s_20_log__r_1000.h5');
% data = dataFile.Analysis_data; legendStr=[legendStr, 'v=30 e=180 gr', 'v=30 e=180 delta', 'v=30 e=180 random'];
% plotOne(data.exact, [data.gr, data.delta, data.random]);
% dataFile = loadHDF5('data/HeuristicsConstWeightVarFaspConstVE_v_30_e_300_f_1-30_s_20_log__r_1000.h5');
% data = dataFile.Analysis_data; legendStr=[legendStr, 'v=30 e=120 gr', 'v=30 e=120 delta', 'v=30 e=120 random'];
% plotOne(data.exact, [data.gr, data.delta, data.random]);

dataFile = loadHDF5('data/HeuristicsConstWeightVarFaspConstVE_v_50_e_100_f_1-30_s_20_log__r_1000.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'v=50 e=100 gr', 'v=50 e=100 delta', 'v=50 e=100 random'];
plotOne(data.exact, [data.gr, data.delta, data.random]);
dataFile = loadHDF5('data/HeuristicsConstWeightVarFaspConstVE_v_50_e_200_f_1-30_s_20_log__r_1000.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'v=50 e=200 gr', 'v=50 e=200 delta', 'v=50 e=200 random'];
plotOne(data.exact, [data.gr, data.delta, data.random]);

l = legend(legendStr);
l.FontSize = 30;
title('Super Algorithm efficiency');
xlabel('FASP size');
ylabel('efficiency');


function plotOne(dataRef, dataIn)
    xVals = unique(dataRef)';

    % output result with size (numberOfDataPoints, numberOfDataSets)
    average = zeros(size(xVals, 2), size(dataIn, 2));

    % Go through each data column separately
    for inputData = 1:size(dataIn, 2)
        currData=dataIn(:, inputData);

        % Calculate how good we ded comparing to exact solution
        quality = currData ./ dataRef;

        % Calculate mean per each x value
        averagec = zeros(length(xVals), 1);
        idx = 1;
        for i = xVals
            averagec(idx) = mean(quality(dataRef == i));
            idx = idx + 1;
        end
        average(:,inputData) = averagec;
    end

    for i = 1:size(average, 2)
        plot(xVals, average(:,i), '-*', 'LineWidth', 8);
    end
end

end