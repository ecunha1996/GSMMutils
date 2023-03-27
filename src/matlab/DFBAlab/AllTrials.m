

matrix = readExpMatrix('C:\Users\Bisbii\PythonProjects\ExpAlgae\data\experimental\Matriz- DCCR Dunaliella salina_dfba.xlsx');



for z=1:length(matrix)-2
    [resT, Yt] = main(matrix, z);
    res = [resT Yt];
    asTable = array2table(res);
    asTable.Properties.VariableNames(1:6) = {'Time (d)', 'Volume (L)', 'Biomass (g/L)', 'P (mmol)', 'N (mmol)', 'Penalty'};
    filename = strcat(strcat('Trial_',int2str(z)), '.csv');
    writetable(asTable, filename);

    
    %% Plotting
    trial = readtable(filename, 'VariableNamingRule', 'preserve');
    [~,sheet_name]=xlsfinfo('C:\Users\Bisbii\PythonProjects\ExpAlgae\data\experimental\Matriz- DCCR Dunaliella salina_dfba.xlsx');
    T = trial.("Time (d)");
    time = str2double(matrix{z}.("Time (d)"));
    f1 = figure(1);
    hold on
    yyaxis left
    a = plot(T,trial.('P (mmol)'));
    xlabel('Time (d)')
    ylabel('Concentration  (mmol/L)');
    yyaxis right
    b = plot(T,trial.('N (mmol)'));
    ylabel('Concentration (mmol/L)');
    hold off;
    legend([a; b] , "HPO4", "NO3");
    outname = strcat(strcat(sheet_name{z}, '_Concentrations'), '.fig');
    savefig(outname);
    close(f1);
    f2 = figure(2);
    plot(T,trial.("Biomass (g/L)"));
    xlabel('Time (d)');
    ylabel('Biomass (g/L)');
    hold on
    scatter(time,matrix{z}(:,{'DW'}).DW,'filled')
    hold off
    legend('in silico', 'Experimental');
    outname = strcat(strcat(sheet_name{z}, '_Biomass'), '.fig');
    savefig(outname);
    message = strcat(sheet_name{z}, ' is over!');
    disp(message);
    close(f2);
end



