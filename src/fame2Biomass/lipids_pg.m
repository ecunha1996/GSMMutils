function return_val = lipids_pg(lipidClass, b)


metmat = double(load('pg_metmat.mat').metmat);


% lipidClass = {'PG__18_0__18_3v2', 'PG__16_0__18_3v2', 'PG__18_3v2__18_1', 'PG__18_2__14_0', 'PG__18_3v2__18_3v2', 'PG__18_0__18_3', 'PG__14_0__14_0', 'PG__16_0__18_0', 'PG__18_3__16_0', 'PG__18_1__18_0', 'PG__18_0__14_0', 'PG__18_3__18_0', 'PG__16_0__18_3', 'PG__16_1__16_0', 'PG__18_1__14_0', 'PG__16_1__18_3', 'PG__18_2__16_0', 'PG__18_2__18_2', 'PG__16_0__18_2', 'PG__14_0__16_1', 'PG__16_1__18_0', 'PG__18_3__16_1', 'PG__18_0__18_0', 'PG__16_1__18_1', 'PG__18_0__18_2', 'PG__14_0__18_1', 'PG__18_1__18_3v2', 'PG__18_0__16_1', 'PG__16_0__18_1', 'PG__18_2__18_3', 'PG__14_0__18_3', 'PG__18_2__16_4', 'PG__18_1__18_3', 'PG__18_2__18_0', 'PG__14_0__18_2', 'PG__16_1__18_3v2', 'PG__18_3__18_1', 'PG__18_1__16_0', 'PG__18_1__16_4', 'PG__14_0__18_3v2', 'PG__18_1__18_1', 'PG__18_3v2__16_1', 'PG__18_3v2__18_0', 'PG__18_3__14_0', 'PG__18_3__16_4', 'PG__16_1__16_1', 'PG__16_0__16_1', 'PG__18_3v2__18_2', 'PG__18_0__16_0', 'PG__16_1__14_0', 'PG__18_2__18_3v2', 'PG__16_0__16_0', 'PG__18_2__18_1', 'PG__14_0__16_0', 'PG__18_1__18_2', 'PG__18_3v2__16_0', 'PG__18_3__18_2', 'PG__18_0__18_1', 'PG__16_0__14_0', 'PG__18_3v2__18_3', 'PG__18_3__18_3', 'PG__18_2__16_1', 'PG__18_3v2__14_0', 'PG__18_1__16_1', 'PG__14_0__18_0'};



% 
% metmat = [
%    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
% 0,0,2,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,1,0,0,0,1,1,0,0,0;
% 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,2,1,0,0,1,0,1,0,1,1,0,1,0,1,1,0,0,1,0,0,0,1,0,0,0,0,0;
% 0,1,0,0,2,1,1,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0;
% 0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;
% 1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,1,1,0,0,0,0,2,0,0,1,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,1,0;
% 0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,2,0,0,0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0;
% 0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,1,0,0,1,0,0,1,0,0,2,0,0,0,0,0,0,0,1,0,1;
% 0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,1,1,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,2,0,0,0,1,0,0,0;
% 1,0,0,1,0,0,0,1,0,1,1,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1;
%     ];

%b = 2*[0.5, 0.0501, 0.2317, 0.2825, 0.0745, 0.0305, 0.0805, 0.1091, 0.0281, 0.1129]';

[nExp, nMets] = size(metmat);

mdl2 = [];
mdl2.A = [
    sparse(metmat), speye(nExp), -speye(nExp);
    ];
mdl2.rhs =  b; 


mdl2.lb = zeros(nMets+nExp+nExp,1);
mdl2.ub = [ones(nMets, 1); 1000*ones(nExp*2,1)];

% ADD UPPER LIMITS FOR KNOWN LIPID CONFIGURATIONS

mdl2.lb(find(strcmp(lipidClass, 'PG__16_0__14_0'))) = 0.0032;
mdl2.ub(find(strcmp(lipidClass, 'PG__16_0__14_0'))) = 0.0213;

mdl2.lb(find(strcmp(lipidClass, 'PG__18_1__16_1'))) = 0.0161;
mdl2.ub(find(strcmp(lipidClass, 'PG__18_1__16_1'))) = 0.0322;

mdl2.lb(find(strcmp(lipidClass, 'PG__16_0__18_2'))) = 0.0161;
mdl2.ub(find(strcmp(lipidClass, 'PG__16_0__18_2'))) = 0.0322;


mdl2.lb(find(strcmp(lipidClass, 'PG__18_2__16_1'))) = 0.1367;
mdl2.ub(find(strcmp(lipidClass, 'PG__18_2__16_1'))) = 0.2734;

mdl2.lb(find(strcmp(lipidClass, 'PG__16_0__18_3'))) = 0.0070;
mdl2.ub(find(strcmp(lipidClass, 'PG__16_0__18_3'))) = 0.0141;

mdl2.lb(find(strcmp(lipidClass, 'PG__18_3__16_1'))) = 0.2288;
mdl2.ub(find(strcmp(lipidClass, 'PG__18_3__16_1'))) = 0.4575;


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