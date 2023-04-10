function ssq = evaluateTrial(matrix, z, INFO)
    [resT, Yt] = main(matrix, z, INFO);
    res = [resT Yt];
    ssq = 0;
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
    to_evaluate = {'DW'; 'Caro'};
    for i=1:length(to_evaluate)
        exp_value = rmmissing(matrix{z}(:,{'Time (d)', to_evaluate{i}}));
        time_value = str2double(exp_value.("Time (d)"));
        exp_value = exp_value(:,to_evaluate{i});
        time_interp = interp1(asTable.("Time (d)"), 1:numel(asTable.("Time (d)")), time_value, 'nearest', 'extrap');
        sim_value = asTable(time_interp, to_evaluate(i));
        relative_error = abs(sim_value - exp_value) ./ exp_value;
        sum_square = sum(relative_error.^2);
        sum_square = sum_square{1,1};
        ssq = ssq + sum_square;
    end
end