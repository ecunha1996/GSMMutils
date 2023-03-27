function return_val = lipids_pc(lipidClass)

metmat = double(load('pc_metmat.mat').metmat);

% lipidClass = {'PC__16_1__18_3v2', 'PC__14_0__18_3', 'PC__18_1__18_3', 'PC__16_2__16_1', 'PC__18_3__18_1', 'PC__18_1__18_0', 'PC__16_0__16_1', 'PC__18_3v2__18_3', 'PC__18_2__16_4', 'PC__16_2__16_2', 'PC__16_2__16_0', 'PC__14_0__14_0', 'PC__16_1__18_2', 'PC__16_1__14_0', 'PC__18_3v2__14_0', 'PC__18_0__14_0', 'PC__18_3__16_0', 'PC__18_0__16_2', 'PC__16_2__18_2', 'PC__18_0__16_0', 'PC__14_0__18_1', 'PC__16_0__16_0', 'PC__18_3__16_1', 'PC__18_3__18_0', 'PC__18_1__18_1', 'PC__16_2__18_1', 'PC__18_3v2__18_0', 'PC__14_0__16_0', 'PC__14_0__18_2', 'PC__18_1__18_2', 'PC__18_3v2__18_1', 'PC__18_0__18_1', 'PC__18_2__18_3v2', 'PC__16_1__16_1', 'PC__16_0__18_3v2', 'PC__18_3v2__16_0', 'PC__18_3__14_0', 'PC__16_0__18_1', 'PC__16_2__18_3v2', 'PC__16_0__18_3', 'PC__16_1__16_0', 'PC__16_1__18_1', 'PC__18_0__18_3', 'PC__16_1__16_2', 'PC__16_2__14_0', 'PC__18_2__18_1', 'PC__18_1__14_0', 'PC__18_1__16_2', 'PC__16_1__18_0', 'PC__14_0__18_3v2', 'PC__18_2__16_0', 'PC__18_3__18_3', 'PC__18_0__16_1', 'PC__14_0__16_2', 'PC__18_3v2__16_2', 'PC__18_1__16_1', 'PC__18_1__16_0', 'PC__18_3__18_2', 'PC__18_2__14_0', 'PC__16_1__18_3', 'PC__18_0__18_0', 'PC__16_0__18_0', 'PC__18_3v2__18_3v2', 'PC__18_1__16_4', 'PC__18_2__18_0', 'PC__18_3v2__16_1', 'PC__18_3__16_4', 'PC__18_2__18_2', 'PC__18_2__18_3', 'PC__18_2__16_1', 'PC__14_0__18_0', 'PC__16_0__14_0', 'PC__18_0__18_3v2', 'PC__18_2__16_2', 'PC__16_0__16_2', 'PC__16_0__18_2', 'PC__18_1__18_3v2', 'PC__18_3__18_3v2', 'PC__16_2__18_3', 'PC__18_0__18_2', 'PC__16_2__18_0', 'PC__18_3v2__18_2', 'PC__14_0__16_1'};
% 

% metmat = [
% 
% 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
% 0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,1,1,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,2,0,1,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0;
% 1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,1;
% 0,0,2,1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0;
% 0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,1,0,0,1,1,1,0,2,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0;
% 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0;
% 0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,2,0,0,1,1,0,1,0,1,0,1,0,1,0,1,1,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1;
% 0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,1,0,1,0,0,1,1,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0;
% 0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,1,1,0,0,1,0,0;
% 0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,1,0,0,0,1,0,0,1,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,1,0;
% 1,1,0,0,0,0,0,0,2,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0;
% ];

b = 2*[0.5, 0.027, 0.32, 0.024, 0.013, 0.093, 0.082, 0.138, 0.112, 0.084, 0.107]';

[nExp, nMets] = size(metmat);

mdl2 = [];
mdl2.A = [
    sparse(metmat), speye(nExp), -speye(nExp);
    ];
mdl2.rhs =  b; 


mdl2.lb = zeros(nMets+nExp+nExp,1);
mdl2.ub = [ones(nMets, 1); 1000*ones(nExp*2,1)];

% ADD UPPER LIMITS FOR KNOWN LIPID CONFIGURATIONS

mdl2.lb(find(strcmp(lipidClass, 'PC__16_0__14_0'))) = 0.0045;
mdl2.ub(find(strcmp(lipidClass, 'PC__16_0__14_0'))) = 0.0290;

mdl2.lb(find(strcmp(lipidClass, 'PC__16_0__18_1'))) = 0.0047;
mdl2.ub(find(strcmp(lipidClass, 'PC__16_0__18_1'))) = 0.0189;

mdl2.lb(find(strcmp(lipidClass, 'PC__16_0__18_2'))) = 0.1058;
mdl2.ub(find(strcmp(lipidClass, 'PC__16_0__18_2'))) = 0.4233;

mdl2.lb(find(strcmp(lipidClass, 'PC__16_0__18_3'))) = 0.0296;
mdl2.ub(find(strcmp(lipidClass, 'PC__16_0__18_3'))) = 0.1185;

mdl2.lb(find(strcmp(lipidClass, 'PC__18_3__16_1'))) = 0.0037;
mdl2.ub(find(strcmp(lipidClass, 'PC__18_3__16_1'))) = 0.0149;

mdl2.lb(find(strcmp(lipidClass, 'PC__18_1__18_1'))) = 0.0004;
mdl2.ub(find(strcmp(lipidClass, 'PC__18_1__18_1'))) = 0.0017;

mdl2.lb(find(strcmp(lipidClass, 'PC__18_0__18_2'))) = 0.0004;
mdl2.ub(find(strcmp(lipidClass, 'PC__18_0__18_2'))) = 0.0017;

mdl2.lb(find(strcmp(lipidClass, 'PC__18_1__18_2'))) = 0.0054;
mdl2.ub(find(strcmp(lipidClass, 'PC__18_1__18_2'))) = 0.0216;

mdl2.lb(find(strcmp(lipidClass, 'PC__18_0__18_3'))) = 0.0054;
mdl2.ub(find(strcmp(lipidClass, 'PC__18_0__18_3'))) = 0.0216;

mdl2.lb(find(strcmp(lipidClass, 'PC__18_2__18_2'))) = 0.0275;
mdl2.ub(find(strcmp(lipidClass, 'PC__18_2__18_2'))) = 0.1101;

mdl2.lb(find(strcmp(lipidClass, 'PC__18_1__18_3'))) = 0.0056;
mdl2.ub(find(strcmp(lipidClass, 'PC__18_1__18_3'))) = 0.0224;

mdl2.lb(find(strcmp(lipidClass, 'PC__18_1__18_3v2'))) = 0.0056;
mdl2.ub(find(strcmp(lipidClass, 'PC__18_1__18_3v2'))) = 0.0224;

mdl2.lb(find(strcmp(lipidClass, 'PC__18_2__18_3v2'))) = 0.0133;
mdl2.ub(find(strcmp(lipidClass, 'PC__18_2__18_3v2'))) = 0.0534;

mdl2.lb(find(strcmp(lipidClass, 'PC__18_2__18_3'))) = 0.0023;
mdl2.ub(find(strcmp(lipidClass, 'PC__18_2__18_3'))) = 0.0092;

mdl2.lb(find(strcmp(lipidClass, 'PC__18_3__18_3'))) = 0.0020;
mdl2.ub(find(strcmp(lipidClass, 'PC__18_3__18_3'))) = 0.0081;



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
writecell(TableOut, 'C:\Users\Bisbii\OneDrive - Universidade do Minho\Algae\Models\Dsalina\lipids\pc_opt_results.tsv', 'FileType', 'text','Delimiter','tab');
return_val = 'return';
end