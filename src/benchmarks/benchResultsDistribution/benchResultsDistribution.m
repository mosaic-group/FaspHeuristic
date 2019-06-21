function benchResultsDistribution

addpath('../matlabHelpers/')
% create figure to show data
figure(1);
clf;
hold on;


dataFile = loadHDF5('../benchHeuristics/data/HeuristicsConstWeightVarFaspConstVE_v_50_e_200_f_1-30_s_20_log__r_1000.h5');

data = dataFile.Analysis_data;


figure(1);
hold on

s = data.random ./ data.exact;
[min(s), mean(s), max(s)] %print range of solutions
h1 = histogram(s);
hold on;
d = data.gr ./ data.exact;
[min(d), mean(d), max(d)]
h2 = histogram(d);

% set(gca,'YScale','log')

h1.Normalization = 'probability';
h2.Normalization = 'probability';
h1.BinWidth = 0.1;
h2.BinWidth = 0.1;
h1.FaceAlpha = 0.3;
h2.FaceAlpha = 0.1;

% plotHistOverlay(h1, 'b')
plotHistOverlay(h2, 'r')

    function plotHistOverlay(histogram, color)
        binRange = histogram.BinEdges(2:end);
        hcx = histogram.Values;
        plot(binRange - histogram.BinWidth/2, hcx, color, 'LineWidth', 1);
    end


end