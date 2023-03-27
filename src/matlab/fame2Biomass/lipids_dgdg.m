function return_val = lipids_dgdg(lipidClass)

metmat = double(load('dgdg_metmat.mat').metmat);


% lipidClass = {'DGDG__18_0__18_0', 'DGDG__18_3__18_0', 'DGDG__16_0__16_2', 'DGDG__16_0__18_2', 'DGDG__16_1__16_1', 'DGDG__16_2__18_3', 'DGDG__18_3__16_0', 'DGDG__18_2__18_1', 'DGDG__16_1__18_0', 'DGDG__18_1__18_3', 'DGDG__18_2__18_2', 'DGDG__18_3__16_1', 'DGDG__16_1__16_0', 'DGDG__18_1__16_0', 'DGDG__18_1__14_0', 'DGDG__16_1__16_2', 'DGDG__18_0__18_3', 'DGDG__16_2__18_0', 'DGDG__18_3__18_2', 'DGDG__18_2__16_1', 'DGDG__18_2__16_4', 'DGDG__18_0__16_0', 'DGDG__18_0__14_0', 'DGDG__18_2__18_3', 'DGDG__18_3__16_4', 'DGDG__16_0__18_0', 'DGDG__18_2__14_0', 'DGDG__18_0__18_2', 'DGDG__18_1__16_1', 'DGDG__18_1__18_1', 'DGDG__18_2__16_0', 'DGDG__16_1__18_3', 'DGDG__16_2__18_1', 'DGDG__16_1__18_2', 'DGDG__16_2__16_2', 'DGDG__18_1__16_4', 'DGDG__18_1__18_0', 'DGDG__16_1__14_0', 'DGDG__16_0__16_1', 'DGDG__18_3__18_1', 'DGDG__18_0__16_1', 'DGDG__16_0__18_1', 'DGDG__18_0__18_1', 'DGDG__16_0__18_3', 'DGDG__16_2__18_2', 'DGDG__18_2__18_0', 'DGDG__18_3__18_3', 'DGDG__16_0__14_0', 'DGDG__16_0__16_0', 'DGDG__18_1__18_2', 'DGDG__18_3__14_0', 'DGDG__16_1__18_1'};
%  
% metmat = [
%         1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
%         0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0;
%         0,0,1,1,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,1,2,0,0,0;
%         0,0,0,0,2,0,0,0,1,0,0,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,1;
%         0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,2,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0;
%         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
%         2,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,1,1,0,0,1,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0;
%         0,0,0,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,0,0,1,0,0,1,1,0,0,1,0,1,1,0,0,0,0,0,0,1,0,1;
%         0,0,0,1,0,0,0,1,0,0,2,0,0,0,0,0,0,0,1,1,1,0,0,1,0,0,1,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0;
%         0,1,0,0,0,1,1,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,2,0,0,0,1,0;
%             ];

b = 2*[0.5, 0.0742, 0.3392, 0.0281, 0.1136, 0.0455, 0.0505, 0.0860, 0.0970, 0.1660]';

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
mdl2.lb(find(strcmp(lipidClass, 'DGDG__18_3__16_4'))) = 0.0111;
mdl2.ub(find(strcmp(lipidClass, 'DGDG__18_3__16_4'))) = 0.0222;
% 

mdl2.lb(find(strcmp(lipidClass, 'DGDG__18_3__16_2'))) = 0.0724;
mdl2.ub(find(strcmp(lipidClass, 'DGDG__18_3__16_2'))) = 0.1448;


mdl2.lb(find(strcmp(lipidClass, 'DGDG__18_2__16_2'))) = 0.0634;
mdl2.ub(find(strcmp(lipidClass, 'DGDG__18_2__16_2'))) = 0.1267;


mdl2.lb(find(strcmp(lipidClass, 'DGDG__18_3__16_0'))) = 0.1576;
mdl2.ub(find(strcmp(lipidClass, 'DGDG__18_3__16_0'))) = 0.3152;
% 

mdl2.lb(find(strcmp(lipidClass, 'DGDG__18_2__16_0'))) = 0.0549;
mdl2.ub(find(strcmp(lipidClass, 'DGDG__18_2__16_0'))) = 0.1097;


mdl2.lb(find(strcmp(lipidClass, 'DGDG__18_1__16_0'))) = 0.0103;
mdl2.ub(find(strcmp(lipidClass, 'DGDG__18_1__16_0'))) = 0.0207;


mdl2.lb(find(strcmp(lipidClass, 'DGDG__18_3__18_3'))) = 0.0318;
mdl2.ub(find(strcmp(lipidClass, 'DGDG__18_3__18_3'))) = 0.0636;


mdl2.lb(find(strcmp(lipidClass, 'DGDG__18_2__16_1'))) = 0.0181;
mdl2.ub(find(strcmp(lipidClass, 'DGDG__18_2__16_1'))) = 0.0362;

mdl2.lb(find(strcmp(lipidClass, 'DGDG__18_1__16_2'))) = 0.0181;
mdl2.ub(find(strcmp(lipidClass, 'DGDG__18_1__16_2'))) = 0.0362;

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

writecell(TableOut, 'C:\Users\Bisbii\OneDrive - Universidade do Minho\Algae\Models\Dsalina\lipids\dgdg_opt_results.tsv', 'FileType', 'text','Delimiter','tab');
return_val = 'return';
end