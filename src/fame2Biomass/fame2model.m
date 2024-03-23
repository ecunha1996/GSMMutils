function return_val = fame2model(lipidClass, lipid_name, b)

disp(lipid_name);

matrix_filename = strcat('model_plutheri/', strcat(lipid_name, '_metmat.mat'));
%res_filename = strcat('model_ngaditana_lipids/', strcat(lipid_name, '_b_vector.mat'));

metmat = double(load(matrix_filename).metmat);
%res_vector = double(load(res_filename).b_vector);

%b = res_vector';
b = b';

[nExp, nMets] = size(metmat);

mdl2 = [];
mdl2.A = [
    sparse(metmat), speye(nExp), -speye(nExp);
    ];
mdl2.rhs =  b; 


mdl2.lb = zeros(nMets+nExp+nExp,1);
mdl2.ub = [ones(nMets, 1); 1000*ones(nExp*2,1)];

% ADD UPPER LIMITS FOR KNOWN LIPID CONFIGURATIONS

const = "C:\Users\Bisbii\PythonProjects\gsmmutils\data\fame2biomass\model_plutheri\constraints.json";
jsonStr = fileread(const);
data = jsondecode(jsonStr);
fieldNames = fieldnames(data);

for i = 1:numel(fieldNames)
    fieldName = fieldNames{i};
    fieldValue = data.(fieldName);
    mdl2.lb(find(strcmp(lipidClass, fieldName))) = fieldValue(1);
    mdl2.ub(find(strcmp(lipidClass, fieldName))) = fieldValue(2);

    if isempty(mdl2.lb(find(strcmp(lipidClass, fieldName)))) && startsWith(fieldName, lipidClass)
        disp("##########");
        disp('holy shit');
        disp(fieldName);
        disp("##########");
    end
end


q = [1./b; 1./b];
mdl2.obj = [zeros(nMets,1); q];

error = 0.2;
mdl2.sense = '=';
% Constrain the errors based on known experimental error/noise
% x must sum to 1
mdl2.lb( [nMets+1, nMets+nExp+1]) = 0;
mdl2.ub( [nMets+1, nMets+nExp+1]) = 0;
%
relErrStoich = 0.1; %allow 5% deviation from sum of x
relErrFA = b(2:end);
relErrFA = (relErrFA.*error);% adds the error from the error argument
relErr = [relErrStoich;relErrFA];

mdl2.ub( nMets + 1:end) = [relErr;relErr];


% Solve
sol = gurobi(mdl2);


stoichs = sol.x(1:nMets);

ep = sol.x( nMets + (1:nExp) ); %upper error
en = sol.x( nMets + nExp + (1:nExp)); %lower error
errs = ep - en; 

LipidStoich = sol.x(1:length(lipidClass));
metIndex = [];
for i=1:length(LipidStoich)
    in = LipidStoich(i);
    if in ~= 0
        metIndex = [metIndex;i];
    end
end
StoichOut = LipidStoich(metIndex);
MetsOut = lipidClass(metIndex)';
TableOut = horzcat(MetsOut, num2cell(StoichOut));
outfilename = strcat('C:\Users\Bisbii\PythonProjects\gsmmutils\data\fame2biomass\model_plutheri\', lipid_name);
outfilename = strcat(outfilename, '_opt_results.tsv');
writecell(TableOut, outfilename, 'FileType', 'text','Delimiter','tab');
return_val = 'return';
end