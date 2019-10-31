function xl
addpath('../matlabHelpers/')

gs=[];
gsd=[];
rs=[];
rsd=[];
v=[];
e=[];

for i=100:100:1000
        file = strcat('./dataDensity5percent/NewTimingConstWeightVarFaspConstVE_v_', num2str(i),'_e_', num2str(i*(i-1)*0.05),'_f_20-20_s_1_lin__r_10.h5');
        dataFile = loadHDF5(file);
        data = dataFile.Analysis_data;
        gs = [gs; mean(data.gr)];
        gsd = [gsd; std(data.gr)];
        rs = [rs; mean(data.random)];
        rsd = [rsd; std(data.random)];
        v = [v; data.vertices(1)];
        e = [e; data.edges(1)];
end
    
    gs'
    rs'
    v'
    e'
    figure(1);
    clf;
    hold on;
    errorbar(e, gs, gsd);
    errorbar(e, rs, rsd);
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')
end
