function benchmark

addpath('../matlabHelpers/')

legendStr={};
dataFile = loadHDF5('/Users/gonciarz/faspResults/random.h5');
data = dataFile.Analysis_data; legendStr=[legendStr, 'random'];
data



% ----------------- histogram of results ----------------------------------
figure(1);
clf;
hold on;

s = data.random ./ data.exact;
[min(s), mean(s), max(s)] %print range of solutions
h1 = histogram(s);

d = data.gr ./ data.exact;
[min(d), mean(d), max(d)]
h2 = histogram(d);


% -------------------------------------------------------------------------

[data, uniqueVals, startIdx] = sortValues(data.exactTime, data);
et=data.exactTime;
rt=data.randomTime;
gt=data.grTime;

et(et==0)=0.00001;
rt(rt==0)=0.00001;
gt(gt==0)=0.00001;

es=data.exact;
rs=data.random;
gs=data.gr;


figure(4)
clf;
% showData(uniqueVals, calcMeanOnUniqueValues(rs./es, startIdx)', 6);
showData(uniqueVals, calcMeanOnUniqueValues(rt./et, startIdx)', 4);
% showData(uniqueVals, calcMeanOnUniqueValues(es, startIdx)', 6);


    function showData(xVals, yVals, fitPolyDegree)
        plot(xVals, yVals);
        hold on;
        fitFunc = polyfit(xVals, yVals, fitPolyDegree);
        fitFuncVals = polyval(fitFunc, xVals');
        plot(xVals, fitFuncVals', 'r')
    end

    function [sortedData, uniqueVals, startIndices]=sortValues(key, data)
        keyCopy = key;
        [sortedKey, indices] = sort(keyCopy);
        data.exact=data.exact(indices);
        data.exactTime=data.exactTime(indices);
        data.random=data.random(indices);
        data.randomTime=data.randomTime(indices);
        data.gr=data.gr(indices);
        data.grTime=data.grTime(indices);
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