function return_val = fame2model(lipidClass, lipid_name)

disp(lipid_name);

matrix_filename = strcat(lipid_name, '_metmat.mat');
res_filename = strcat(lipid_name, '_metmat_b.mat');

metmat = double(load(matrix_filename).metmat);
res_vector = double(load(res_filename).metmat_b);

b = 2*res_vector';

disp(b);

[nExp, nMets] = size(metmat);

mdl2 = [];
mdl2.A = [
    sparse(metmat), speye(nExp), -speye(nExp);
    ];
mdl2.rhs =  b; 


mdl2.lb = zeros(nMets+nExp+nExp,1);
mdl2.ub = [ones(nMets, 1); 1000*ones(nExp*2,1)];

% ADD UPPER LIMITS FOR KNOWN LIPID CONFIGURATIONS

mdl2.lb(find(strcmp(lipidClass, 'PG__16_0__14_0'))) = 0.004;
mdl2.ub(find(strcmp(lipidClass, 'PG__16_0__14_0'))) = 0.482;

mdl2.lb(find(strcmp(lipidClass, 'PG__16_0__18_2'))) = 0.039;
mdl2.ub(find(strcmp(lipidClass, 'PG__16_0__18_2'))) = 0.482;
 
mdl2.lb(find(strcmp(lipidClass, 'PG__18_1__16_1'))) = 0.031;
mdl2.ub(find(strcmp(lipidClass, 'PG__18_1__16_1'))) = 0.428;

mdl2.lb(find(strcmp(lipidClass, 'PG__18_2__16_1'))) = 0.039;
mdl2.ub(find(strcmp(lipidClass, 'PG__18_2__16_1'))) = 0.428;

mdl2.ub(find(strcmp(lipidClass, 'PG__16_0__18_3'))) = 0.482;

mdl2.lb(find(strcmp(lipidClass, 'PG__18_3__16_1'))) = 0.006;
mdl2.ub(find(strcmp(lipidClass, 'PG__18_3__16_1'))) = 0.576;


q = [1./b; 1./b];
mdl2.obj = [zeros(nMets,1); q];

error = 0.1;
mdl2.sense = '=';
% Constrain the errors based on known experimental error/noise
% x must sum to 1
mdl2.lb( [nMets+1, nMets+nExp+1]) = 0;
mdl2.ub( [nMets+1, nMets+nExp+1]) = 0;
%
relErrStoich = 0.05; %allow 5% deviation from sum of x
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
writecell(TableOut, 'C:\Users\Bisbii\OneDrive - Universidade do Minho\Algae\Models\Dsalina\lipids\pg_opt_results.tsv', 'FileType', 'text','Delimiter','tab');
return_val = 'return';
end