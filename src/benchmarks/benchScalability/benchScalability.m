function benchScalability

addpath('../matlabHelpers/')

figure(1)
clf;
d = [];
d = addData(loadHDF5('./data/DensityVsFaspSizeTime/NewTimingConstWeightVarFaspConstVE_v_300_e_375_f_1-41_s_21_lin__r_100.h5'), d);
d = addData(loadHDF5('./data/DensityVsFaspSizeTime/NewTimingConstWeightVarFaspConstVE_v_300_e_450_f_1-41_s_21_lin__r_100.h5'), d);
d = addData(loadHDF5('./data/DensityVsFaspSizeTime/NewTimingConstWeightVarFaspConstVE_v_300_e_525_f_1-41_s_21_lin__r_100.h5'), d);
d = addData(loadHDF5('./data/DensityVsFaspSizeTime/NewTimingConstWeightVarFaspConstVE_v_300_e_600_f_1-41_s_21_lin__r_100.h5'), d);
d = addData(loadHDF5('./data/DensityVsFaspSizeTime/NewTimingConstWeightVarFaspConstVE_v_300_e_675_f_1-41_s_21_lin__r_100.h5'), d);
d = addData(loadHDF5('./data/DensityVsFaspSizeTime/NewTimingConstWeightVarFaspConstVE_v_300_e_750_f_1-41_s_21_lin__r_100.h5'), d);
d = addData(loadHDF5('./data/DensityVsFaspSizeTime/NewTimingConstWeightVarFaspConstVE_v_300_e_825_f_1-41_s_21_lin__r_100.h5'), d);
d = addData(loadHDF5('./data/DensityVsFaspSizeTime/NewTimingConstWeightVarFaspConstVE_v_300_e_900_f_1-41_s_21_lin__r_100.h5'), d);
d = addData(loadHDF5('./data/DensityVsFaspSizeTime/NewTimingConstWeightVarFaspConstVE_v_300_e_975_f_1-41_s_21_lin__r_100.h5'), d);
d = addData(loadHDF5('./data/DensityVsFaspSizeTime/NewTimingConstWeightVarFaspConstVE_v_300_e_1050_f_1-41_s_21_lin__r_100.h5'), d);
d = addData(loadHDF5('./data/DensityVsFaspSizeTime/NewTimingConstWeightVarFaspConstVE_v_300_e_1125_f_1-41_s_21_lin__r_100.h5'), d);
d = addData(loadHDF5('./data/DensityVsFaspSizeTime/NewTimingConstWeightVarFaspConstVE_v_300_e_1200_f_1-41_s_21_lin__r_100.h5'), d);
d = addData(loadHDF5('./data/DensityVsFaspSizeTime/NewTimingConstWeightVarFaspConstVE_v_300_e_1275_f_1-41_s_21_lin__r_100.h5'), d);
d = addData(loadHDF5('./data/DensityVsFaspSizeTime/NewTimingConstWeightVarFaspConstVE_v_300_e_1350_f_1-41_s_21_lin__r_100.h5'), d);
d = addData(loadHDF5('./data/DensityVsFaspSizeTime/NewTimingConstWeightVarFaspConstVE_v_300_e_1425_f_1-41_s_21_lin__r_100.h5'), d);
[d, xvals] = addData(loadHDF5('./data/DensityVsFaspSizeTime/NewTimingConstWeightVarFaspConstVE_v_300_e_1500_f_1-41_s_21_lin__r_100.h5'), d);

[X, Y] = meshgrid(xvals, 1.25:0.25:5)
s = mesh(X,Y,d);
s.FaceColor = 'flat';
colorbar
xlabel('|FASP|')
ylabel('density = |E| / |V|')
zlabel('time')
title('|V| = 300')
xlim([1 41])
ylim([1.25 5])

figure(2)
clf;
d = [];
d = addData2(loadHDF5('./data/DensityVsVerticesTime/TimingConstWeightVarFaspConstVE_v_1000-4000_d_1_25_f_40_s_7lin_r_100.h5'), d);
d = addData2(loadHDF5('./data/DensityVsVerticesTime/TimingConstWeightVarFaspConstVE_v_1000-4000_d_1_50_f_40_s_7lin_r_100.h5'), d);
d = addData2(loadHDF5('./data/DensityVsVerticesTime/TimingConstWeightVarFaspConstVE_v_1000-4000_d_1_75_f_40_s_7lin_r_100.h5'), d);
d = addData2(loadHDF5('./data/DensityVsVerticesTime/TimingConstWeightVarFaspConstVE_v_1000-4000_d_2_00_f_40_s_7lin_r_100.h5'), d);
d = addData2(loadHDF5('./data/DensityVsVerticesTime/TimingConstWeightVarFaspConstVE_v_1000-4000_d_2_25_f_40_s_7lin_r_100.h5'), d);
d = addData2(loadHDF5('./data/DensityVsVerticesTime/TimingConstWeightVarFaspConstVE_v_1000-4000_d_2_50_f_40_s_7lin_r_100.h5'), d);
d = addData2(loadHDF5('./data/DensityVsVerticesTime/TimingConstWeightVarFaspConstVE_v_1000-4000_d_2_75_f_40_s_7lin_r_100.h5'), d);
d = addData2(loadHDF5('./data/DensityVsVerticesTime/TimingConstWeightVarFaspConstVE_v_1000-4000_d_3_00_f_40_s_7lin_r_100.h5'), d);
d = addData2(loadHDF5('./data/DensityVsVerticesTime/TimingConstWeightVarFaspConstVE_v_1000-4000_d_3_25_f_40_s_7lin_r_100.h5'), d);
d = addData2(loadHDF5('./data/DensityVsVerticesTime/TimingConstWeightVarFaspConstVE_v_1000-4000_d_3_50_f_40_s_7lin_r_100.h5'), d);
d = addData2(loadHDF5('./data/DensityVsVerticesTime/TimingConstWeightVarFaspConstVE_v_1000-4000_d_3_75_f_40_s_7lin_r_100.h5'), d);
[d, xvals] = addData2(loadHDF5('./data/DensityVsVerticesTime/TimingConstWeightVarFaspConstVE_v_1000-4000_d_4_00_f_40_s_7lin_r_100.h5'), d);

[X, Y] = meshgrid(xvals, 1.25:0.25:4)
s = mesh(X,Y,d);
s.FaceColor = 'flat';
colorbar
xlabel('|V|')
ylabel('density = |E| / |V|')
zlabel('time')
title('|FASP| = 40')
ylim([1.25 4])

    function [dd, uniqeValues] = addData2(dataFile, dd)
        data = dataFile.Analysis_data;
        [data, uniqeValues, startIdx] = sortValues(data.vertices, data);
        dd = [dd; calcMeanOnUniqueValues(data.random, startIdx)];
    end


    function [dd, uniqeValues] = addData(dataFile, dd)
        data = dataFile.Analysis_data;
        [data, uniqeValues, startIdx] = sortValues(data.exact, data);
        dd = [dd; calcMeanOnUniqueValues(data.random, startIdx)];
    end


    function [sortedData, uniqueVals, startIndices]=sortValues(key, data)
        keyCopy = key;
        [sortedKey, indices] = sort(keyCopy);
        data.exact=data.exact(indices);
        data.exactTime=data.exactTime(indices);
        data.random=data.random(indices);
        data.randomTime=data.randomTime(indices);
        [uniqueVals, startIndices] = unique(sortedKey);
        sortedData = data;
    end

    function result = calcMeanOnUniqueValues(data, startIdxOfUniqueData)
        result = zeros(1, size(startIdxOfUniqueData, 1));
        for i = 1:size(startIdxOfUniqueData, 1)
            si = startIdxOfUniqueData(i);

            if i == size(startIdxOfUniqueData, 1)
                ei = size(data, 1);
            else
                ei = startIdxOfUniqueData(i + 1) - 1;
            end
            
            m = mean(data(si:ei));
            result(i) = m;
        end
    end
end
