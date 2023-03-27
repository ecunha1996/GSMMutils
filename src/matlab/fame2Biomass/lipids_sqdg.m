function return_val = lipids_sqdg(lipidClass)

metmat = double(load('sqdg_metmat.mat').metmat);


% lipidClass = {'SQDG__16_2__18_3', 'SQDG__18_3__16_0', 'SQDG__16_0__16_2', 'SQDG__18_0__18_0', 'SQDG__18_1__18_2', 'SQDG__18_2__16_1', 'SQDG__16_0__16_0', 'SQDG__16_0__14_0', 'SQDG__16_0__18_0', 'SQDG__14_0__16_2', 'SQDG__14_0__16_0', 'SQDG__16_1__18_2', 'SQDG__16_1__16_2', 'SQDG__16_2__18_1', 'SQDG__18_3__18_3', 'SQDG__18_0__18_1', 'SQDG__18_2__16_4', 'SQDG__18_1__16_0', 'SQDG__14_0__14_0', 'SQDG__18_1__18_0', 'SQDG__16_1__16_0', 'SQDG__16_0__16_1', 'SQDG__18_2__18_2', 'SQDG__16_1__18_0', 'SQDG__18_2__18_0', 'SQDG__16_2__16_2', 'SQDG__16_1__18_3', 'SQDG__16_2__18_2', 'SQDG__18_3__16_4', 'SQDG__16_2__18_0', 'SQDG__14_0__16_1', 'SQDG__18_2__16_0', 'SQDG__18_1__18_1', 'SQDG__18_1__16_1', 'SQDG__18_0__18_2', 'SQDG__16_1__16_1', 'SQDG__16_1__18_1', 'SQDG__18_1__16_4', 'SQDG__18_2__18_3'};
% 
% 
% 
% 
% metmat = [
%         1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
%         0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0;
%         0,1,1,0,0,0,2,1,1,0,1,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0;
%         0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,1,0,1,0,0,1,0,0,0,1,0,0,1,0,2,1,0,0;
%         1,0,1,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,2,0,1,0,1,0,0,0,0,0,0,0,0,0;
%         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0;
%         0,0,0,2,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0;
%         0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,2,1,0,0,1,1,0;
%         0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,2,0,1,0,0,1,0,0,0,1,0,0,1,0,0,0,1;
%         1,1,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1;
%             ];

b = 2*[0.5, 0.0245, 0.5236, 0.0300, 0.0189, 0.0515, 0.0798, 0.0774, 0.0723, 0.1221]';

[nExp, nMets] = size(metmat);

mdl2 = [];
mdl2.A = [
    sparse(metmat), speye(nExp), -speye(nExp);
    ];
mdl2.rhs = [ b; ];


mdl2.lb = zeros(nMets+nExp+nExp,1);
mdl2.ub = [ones(nMets, 1); 1000*ones(nExp*2,1)];

% ADD UPPER LIMITS FOR KNOWN LIPID CONFIGURATIONS
% % 
mdl2.lb(find(strcmp(lipidClass, 'SQDG__16_0__16_0'))) = 0.2910;
mdl2.ub(find(strcmp(lipidClass, 'SQDG__16_0__16_0'))) = 0.5819;

% 
mdl2.lb(find(strcmp(lipidClass, 'DGDG__18_1__16_0'))) = 0.0354;
mdl2.ub(find(strcmp(lipidClass, 'DGDG__18_1__16_0'))) = 0.0708;

% % 
mdl2.lb(find(strcmp(lipidClass, 'SQDG__18_3__16_0'))) = 0.0472;
mdl2.ub(find(strcmp(lipidClass, 'SQDG__18_3__16_0'))) = 0.0944;
% 
mdl2.lb(find(strcmp(lipidClass, 'SQDG__18_2__16_0'))) = 0.0241;
mdl2.ub(find(strcmp(lipidClass, 'SQDG__18_2__16_0'))) = 0.0483;


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

writecell(TableOut, 'C:\Users\Bisbii\OneDrive - Universidade do Minho\Algae\Models\Dsalina\lipids\sqdg_opt_results.tsv', 'FileType', 'text','Delimiter','tab');
return_val = 'return';
end