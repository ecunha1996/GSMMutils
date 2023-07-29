function return_val = lipids_pe(lipidClass)

metmat = double(load('pe_metmat.mat').metmat);



% lipidClass = {'PE__18_1__16_0', 'PE__16_1__18_3v2', 'PE__16_2__16_0', 'PE__16_1__18_0', 'PE__18_3__18_1', 'PE__16_1__16_2', 'PE__18_3__18_3', 'PE__14_0__16_0', 'PE__18_0__16_2', 'PE__14_0__18_0', 'PE__16_0__14_0', 'PE__16_0__16_2', 'PE__18_3__18_0', 'PE__18_0__18_3v2', 'PE__18_1__18_3', 'PE__16_0__18_3', 'PE__16_1__16_0', 'PE__18_3v2__16_1', 'PE__16_1__16_1', 'PE__18_1__18_0', 'PE__18_2__18_0', 'PE__18_1__14_0', 'PE__18_1__18_3v2', 'PE__18_2__18_3v2', 'PE__16_2__16_1', 'PE__14_0__16_2', 'PE__16_2__16_2', 'PE__16_0__18_3v2', 'PE__18_0__18_3', 'PE__18_3v2__18_3', 'PE__16_0__16_1', 'PE__18_3__14_0', 'PE__16_0__18_2', 'PE__16_1__18_3', 'PE__18_3__16_0', 'PE__18_3__16_1', 'PE__14_0__18_1', 'PE__18_0__18_1', 'PE__18_1__18_1', 'PE__18_3__16_4', 'PE__18_3v2__14_0', 'PE__18_1__18_2', 'PE__18_2__16_4', 'PE__16_2__18_0', 'PE__16_1__14_0', 'PE__16_1__18_1', 'PE__18_2__18_2', 'PE__14_0__14_0', 'PE__16_2__18_2', 'PE__18_3v2__18_1', 'PE__18_0__18_0', 'PE__16_0__16_0', 'PE__14_0__18_3', 'PE__18_2__18_3', 'PE__18_2__14_0', 'PE__18_3__18_2', 'PE__18_0__16_0', 'PE__18_3__18_3v2', 'PE__14_0__16_1', 'PE__18_2__18_1', 'PE__18_3v2__16_2', 'PE__18_0__14_0', 'PE__18_1__16_1', 'PE__18_0__16_1', 'PE__14_0__18_3v2', 'PE__16_2__18_3', 'PE__16_1__18_2', 'PE__18_3v2__16_0', 'PE__18_3v2__18_3v2', 'PE__16_0__18_0', 'PE__18_1__16_2', 'PE__18_1__16_4', 'PE__18_3v2__18_2', 'PE__18_2__16_0', 'PE__16_2__14_0', 'PE__18_0__18_2', 'PE__18_3v2__18_0', 'PE__16_0__18_1', 'PE__16_2__18_3v2', 'PE__18_2__16_1', 'PE__14_0__18_2'}
% 
% 
% 
% metmat = [
%             1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
%             0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,2,0,0,0,0,1,0,1,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1;
%             1,0,1,0,0,0,0,1,0,0,1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,1,0,0,0;
%             0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1,1,2,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0;
%             0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0;
%             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0;
%             0,0,0,1,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,2,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0;
%             1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,2,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0;
%             0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,2,0,1,0,0,0,0,1,1,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,1,0,0,0,1,1;
%             0,0,0,0,1,0,2,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,1,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
%             0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,1,2,0,0,0,1,0,0,0,1,0,1,0,0;
%             ];

b = 2*[0.5, 0.0690, 0.2552, 0.0414, 0.0517, 0.1165, 0.0710, 0.1788, 0.1051, 0.0324, 0.0789]';

[nExp, nMets] = size(metmat);

mdl2 = [];
mdl2.A = [
    sparse(metmat), speye(nExp), -speye(nExp);
    ];
mdl2.rhs = [ b; ];


mdl2.lb = zeros(nMets+nExp+nExp,1);
mdl2.ub = [ones(nMets, 1); 1000*ones(nExp*2,1)];

% ADD UPPER LIMITS FOR KNOWN LIPID CONFIGURATIONS

mdl2.lb(find(strcmp(lipidClass, 'PE__16_0__14_0'))) = 0.0047;
mdl2.ub(find(strcmp(lipidClass, 'PE__16_0__14_0'))) = 0.0846;

mdl2.lb(find(strcmp(lipidClass, 'PE__16_0__18_1'))) = 0.0106;
mdl2.ub(find(strcmp(lipidClass, 'PE__16_0__18_1'))) = 0.0423;


mdl2.lb(find(strcmp(lipidClass, 'PE__16_0__18_2'))) = 0.0174;
mdl2.ub(find(strcmp(lipidClass, 'PE__16_0__18_2'))) = 0.0696;


mdl2.lb(find(strcmp(lipidClass, 'PE__16_0__18_3'))) = 0.0017;
mdl2.ub(find(strcmp(lipidClass, 'PE__16_0__18_3'))) = 0.0067;

mdl2.lb(find(strcmp(lipidClass, 'PE__18_1__18_1'))) = 0.0073;
mdl2.ub(find(strcmp(lipidClass, 'PE__18_1__18_1'))) = 0.0294;


mdl2.lb(find(strcmp(lipidClass, 'PE__18_0__18_2'))) = 0.0073;
mdl2.ub(find(strcmp(lipidClass, 'PE__18_0__18_2'))) = 0.0294;


mdl2.lb(find(strcmp(lipidClass, 'PE__18_1__18_2'))) = 0.0209;
mdl2.ub(find(strcmp(lipidClass, 'PE__18_1__18_2'))) = 0.0835;


mdl2.lb(find(strcmp(lipidClass, 'PE__18_0__18_3'))) = 0.0209;
mdl2.ub(find(strcmp(lipidClass, 'PE__18_0__18_3'))) = 0.0835;


mdl2.lb(find(strcmp(lipidClass, 'PE__18_2__18_2'))) = 0.0190;
mdl2.ub(find(strcmp(lipidClass, 'PE__18_2__18_2'))) = 0.0758;

mdl2.lb(find(strcmp(lipidClass, 'PE__18_1__18_3'))) = 0.0190;
mdl2.ub(find(strcmp(lipidClass, 'PE__18_1__18_3'))) = 0.0758;

mdl2.lb(find(strcmp(lipidClass, 'PE__18_1__18_3v2'))) = 0.0190;
mdl2.ub(find(strcmp(lipidClass, 'PE__18_1__18_3v2'))) = 0.0758;

mdl2.lb(find(strcmp(lipidClass, 'PE__18_2__18_3v2'))) = 0.0177;
mdl2.ub(find(strcmp(lipidClass, 'PE__18_2__18_3v2'))) = 0.0707;

mdl2.lb(find(strcmp(lipidClass, 'PE__18_3__18_3'))) = 0.0024;
mdl2.ub(find(strcmp(lipidClass, 'PE__18_3__18_3'))) = 0.0097;



q = [1./b; 1./b];
mdl2.obj = [zeros(nMets,1); q];

error = 0.1;
mdl2.sense = '=';
% Constrain the errors based on known experimental error/noise
% x must sum to 1
mdl2.lb( [nMets+1, nMets+nExp+1]) = 0;
mdl2.ub( [nMets+1, nMets+nExp+1]) = 0;
%
relErrStoich = 0.05; %allow 1% deviation from sum of x
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
writecell(TableOut, 'C:\Users\Bisbii\OneDrive - Universidade do Minho\Algae\Models\Dsalina\lipids\pe_opt_results.tsv', 'FileType', 'text','Delimiter','tab');
return_val = 'return';
end