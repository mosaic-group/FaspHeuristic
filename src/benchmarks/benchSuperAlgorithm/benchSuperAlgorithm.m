function benchSuperAlgorithm

addpath('../matlabHelpers/')
% create figure to show data
figure(1);
clf;
hold on;

legendStr={};

dataFile = loadHDF5('data/SuperAlgorithmConstWeightVarFaspConstVE_v_30_e_60_f_1-30_s_20_r_1000.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'v=30 e=60'];
plotOne(data.exact, [data.sa]);
dataFile = loadHDF5('data/SuperAlgorithmConstWeightVarFaspConstVE_v_30_e_120_f_1-30_s_20_r_1000.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'v=30 e=120'];
plotOne(data.exact, [data.sa]);
dataFile = loadHDF5('data/SuperAlgorithmConstWeightVarFaspConstVE_v_30_e_180_f_1-30_s_20_r_1000.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'v=30 e=180'];
plotOne(data.exact, [data.sa]);
dataFile = loadHDF5('data/SuperAlgorithmConstWeightVarFaspConstVE_v_30_e_240_f_1-30_s_20_r_1000.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'v=30 e=240'];
plotOne(data.exact, [data.sa]);
dataFile = loadHDF5('data/SuperAlgorithmConstWeightVarFaspConstVE_v_30_e_300_f_1-30_s_20_r_1000.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'v=30 e=300'];
plotOne(data.exact, [data.sa]);


dataFile = loadHDF5('data/SuperAlgorithmConstWeightVarFaspConstVE_v_15_e_30_f_1-15_s_15_r_1000.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'v=15 e=30'];
plotOne(data.exact, [data.sa]);
dataFile = loadHDF5('data/SuperAlgorithmConstWeightVarFaspConstVE_v_15_e_55_f_1-15_s_15_r_1000.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'v=15 e=55'];
plotOne(data.exact, [data.sa]);
dataFile = loadHDF5('data/SuperAlgorithmConstWeightVarFaspConstVE_v_15_e_80_f_1-15_s_15_r_1000.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'v=15 e=80'];
plotOne(data.exact, [data.sa]);
dataFile = loadHDF5('data/SuperAlgorithmConstWeightVarFaspConstVE_v_15_e_105_f_1-15_s_15_r_1000.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'v=15 e=105'];
plotOne(data.exact, [data.sa]);
% 
% 
dataFile = loadHDF5('data/SuperAlgorithmConstWeightVarFaspConstVE_v_50_e_100_f_1-30_s_20_r_1000.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'v=50 e=100'];
plotOne(data.exact, [data.sa]);
dataFile = loadHDF5('data/SuperAlgorithmConstWeightVarFaspConstVE_v_50_e_200_f_1-30_s_20_r_1000.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'v=50 e=200'];
plotOne(data.exact, [data.sa]);
dataFile = loadHDF5('data/SuperAlgorithmConstWeightVarFaspConstVE_v_50_e_300_f_1-30_s_20_r_1000.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'v=50 e=300'];
plotOne(data.exact, [data.sa]);
dataFile = loadHDF5('data/SuperAlgorithmConstWeightVarFaspConstVE_v_50_e_400_f_1-30_s_20_r_1000.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'v=50 e=400'];
plotOne(data.exact, [data.sa]);
dataFile = loadHDF5('data/SuperAlgorithmConstWeightVarFaspConstVE_v_50_e_500_f_1-30_s_20_r_1000.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'v=50 e=500'];
plotOne(data.exact, [data.sa]);




l = legend(legendStr);
l.FontSize = 20;
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
        plot(xVals, average(:,i), '-*');
    end
end

end