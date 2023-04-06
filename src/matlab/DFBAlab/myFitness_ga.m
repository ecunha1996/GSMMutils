function fitness = myFitness_ga(x0)
% Simulate the kinetic model with the given parameters
matrix = readExpMatrix('C:\Users\Bisbii\PythonProjects\ExpGSMM\data\experimental\Matriz- DCCR Dunaliella salina_dfba.xlsx');
disp(x0);
INFO.ro1 = x0(1);
INFO.ro0 = x0(2);
INFO.wPmin = x0(3);
INFO.wPopt = x0(4);
INFO.a0 = x0(5);
INFO.a1 = x0(6);
INFO.a2 = x0(7);
INFO.a3 = x0(8);
INFO.l = x0(9);
INFO.smoothing_factor = x0(10);
INFO.vhpo4max = x0(11);
INFO.ExA = x0(12);
INFO.vcarmax = x0(13);
ssq = 0;
for z=35:35
    [resT, Yt] = main(matrix, z, INFO);
    res = [resT Yt];
    %res =  [sort(4*rand(100,1), 'ascend'), rand(100,1), sort(linspace(0.1, 0.9, 100)', 'ascend'), rand(100,15)];
    asTable = array2table(res);
    asTable.Properties.VariableNames(1:18) = {'Time (d)', 'Volume (L)', 'DW', 'P (mmol)', 'N (mmol)', 'Active Biomass (g/L)', 'Starch Concentration (g/L)', 'Carotene Concentration', 'TAG concentration',...
        'glycerol concentration', 'nitrogen quota', 'chl quota', 'starch quota', 'glycerol quota', 'Caro', 'TAG quota', 'P quota', 'Penalty'};
    asTable = unique(asTable);
    time = matrix{z}.("Time (d)");
    if ~isfloat(time) && ~strcmp(class(time), 'double')
        if all(cellfun(@ischar, time))
            time = str2double(matrix{z}.("Time (d)"));
        end
    end
    time_interp = interp1(asTable.('Time (d)')', 1:numel(asTable.('Time (d)')'), time, 'nearest', 'extrap');
    to_evaluate = ["DW", "Caro"];
    for i=1:length(to_evaluate)
        v = to_evaluate(i);
        exp_value = matrix{z}(:,to_evaluate(i));
        sim_value = asTable(time_interp, to_evaluate(i));
        relative_error = abs(sim_value - exp_value) ./ exp_value;
        sum_square = sum(relative_error.^2);
        sum_square = sum_square{1,1};
        ssq = ssq + sum_square;
    end
end
disp(ssq);
fitness = ssq;
end