function benchTimingConstFaspVarVE

addpath('../matlabHelpers/')
% create figure to show data
figure(1);
clf;
hold on;

legendStr={};

dataFile = loadHDF5('data/TimingConstWeightAndFaspVarVE_v_20-200_d_2_f_5_s_19lin_r_200.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'e=40-400 f=5 random'];
plotOne(data.vertices, [data.randomTime]);

dataFile = loadHDF5('data/TimingConstWeightAndFaspVarVE_v_20-200_d_2_f_10_s_19lin_r_200.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'e=40-400 f=10 random'];
plotOne(data.vertices, [data.randomTime]);

dataFile = loadHDF5('data/TimingConstWeightAndFaspVarVE_v_20-200_d_2_f_15_s_19lin_r_200.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'e=40-400 f=15 random'];
plotOne(data.vertices, [data.randomTime]);



% dataFile = loadHDF5('data/TimingConstWeightAndFaspVarVE_v_20-200_d_3_f_5_s_19lin_r_200.h5');
% data = dataFile.Analysis_data; legendStr=[legendStr, 'e=60-600 f=5 random'];
% plotOne(data.vertices, [data.randomTime]);
% 
% dataFile = loadHDF5('data/TimingConstWeightAndFaspVarVE_v_20-200_d_3_f_10_s_19lin_r_200.h5');
% data = dataFile.Analysis_data; legendStr=[legendStr, 'e=60-600 f=10 random'];
% plotOne(data.vertices, [data.randomTime]);
% 
% 
% 
% dataFile = loadHDF5('data/TimingConstWeightAndFaspVarVE_v_20-200_d_4_f_5_s_19lin_r_200.h5');
% data = dataFile.Analysis_data; legendStr=[legendStr, 'e=80-800 f=5 random'];
% plotOne(data.vertices, [data.randomTime]);
% 
% 
% dataFile = loadHDF5('data/TimingConstWeightAndFaspVarVE_v_20-200_d_5_f_5_s_19lin_r_200.h5');
% data = dataFile.Analysis_data; legendStr=[legendStr, 'e=100-1000 f=5 random'];
% plotOne(data.vertices, [data.randomTime]);
% 
% dataFile = loadHDF5('data/TimingConstWeightAndFaspVarVE_v_20-200_d_6_f_5_s_19lin_r_200.h5');
% data = dataFile.Analysis_data; legendStr=[legendStr, 'e=120-1200 f=5 random'];
% plotOne(data.vertices, [data.randomTime]);

l = legend(legendStr);
l.FontSize = 10;
title('Time benchmark');
xlabel('#v');
ylabel('time(s)');

% set(gca,'YScale','log')

function plotOne(dataRef, dataIn)
    xVals = unique(dataRef)'

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
        plot(xVals, average(:,i), '-*', 'LineWidth', 3);
%         f=fit(xVals',average(:,i), 'poly3')
%         plot(f, xVals, average(:,i));
    end
end

end