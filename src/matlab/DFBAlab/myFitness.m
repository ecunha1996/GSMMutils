function fitness = myfitness(x, x0)
% Simulate the kinetic model with the given parameters
matrix = readExpMatrix('C:\Users\Bisbii\PythonProjects\ExpGSMM\data\experimental\Matriz- DCCR Dunaliella salina_dfba.xlsx');
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
for z=34:34
    [resT, Yt] = main(matrix, z, INFO);
    res = [resT Yt];
    % time_temp = sort(4*rand(100,1), 'ascend');
    % values = sort(linspace(0.1, 0.9, 100)', 'ascend');
    % res =  [time_temp, rand(100,1), values, rand(100,15)];
    asTable = array2table(res);
    asTable.Properties.VariableNames(1:18) = {'Time (d)', 'Volume (L)', 'Biomass (g/L)', 'P (mmol)', 'N (mmol)', 'Active Biomass (g/L)', 'Starch Concentration (g/L)', 'Carotene Concentration', 'TAG concentration',...
        'glycerol concentration', 'nitrogen quota', 'chl quota', 'starch quota', 'glycerol quota', 'carotene quota', 'TAG quota', 'P quota', 'Penalty'};

    time = matrix{z}.("Time (d)");
    if ~isfloat(time) && ~strcmp(class(time), 'double')
        if all(cellfun(@ischar, time))
            time = str2double(matrix{z}.("Time (d)"));
        end
    end
    time_interp = interp1(asTable.('Time (d)')', 1:numel(asTable.('Time (d)')'), time, 'nearest', 'extrap');
    exp_biomass = matrix{z}(:,{'DW'}).DW;
    sim_biomass = asTable(time_interp, 'Biomass (g/L)');

    relative_error = abs(sim_biomass - exp_biomass) ./ exp_biomass;

    ssq = sum(relative_error.^2);

    fitness = -ssq;
end