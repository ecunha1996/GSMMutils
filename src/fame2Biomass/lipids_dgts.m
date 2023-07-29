function return_val = lipids_dgts(lipidClass)

metmat = double(load('dgts_metmat.mat').metmat);

% lipidClass = {'DGTS__16_0__16_2', 'DGTS__14_0__18_3', 'DGTS__16_1__18_1', 'DGTS__18_3__18_3v2', 'DGTS__14_0__18_0', 'DGTS__14_0__18_1', 'DGTS__18_0__18_3', 'DGTS__16_2__16_2', 'DGTS__18_2__16_4', 'DGTS__18_3__16_4', 'DGTS__16_0__18_3', 'DGTS__14_0__16_2', 'DGTS__14_0__18_2', 'DGTS__18_1__18_3v2', 'DGTS__14_0__16_0', 'DGTS__16_2__18_0', 'DGTS__18_2__18_2', 'DGTS__16_2__18_1', 'DGTS__16_1__16_1', 'DGTS__16_0__14_0', 'DGTS__18_2__18_3', 'DGTS__16_0__18_0', 'DGTS__18_0__18_0', 'DGTS__16_1__18_0', 'DGTS__18_1__18_2', 'DGTS__16_1__16_2', 'DGTS__14_0__14_0', 'DGTS__16_0__18_1', 'DGTS__16_0__16_1', 'DGTS__18_2__18_3v2', 'DGTS__18_0__18_2', 'DGTS__18_1__16_4', 'DGTS__14_0__16_1', 'DGTS__18_3__18_3', 'DGTS__18_3v2__16_4', 'DGTS__16_2__18_3', 'DGTS__16_0__18_3v2', 'DGTS__18_1__18_3', 'DGTS__18_1__18_1', 'DGTS__16_2__18_2', 'DGTS__16_1__18_2', 'DGTS__16_1__18_3', 'DGTS__18_0__18_1', 'DGTS__16_0__16_0', 'DGTS__16_0__18_2'};
% 
% metmat = [
%     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1;
%     0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,2,0,0,1,0,0,1,0,0,1,0;
%     1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,2,1,0,0,0,1,1,0,0,0,0,0,0;
%     1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0;
%     0,0,0,0,0,0,0,0,2,1,0,0,0,0,0,1,0,0,0,1,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0;
%     0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0;
%     0,0,0,0,0,0,0,0,0,0,0,1,0,2,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1;
%     0,1,1,0,0,2,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1;
%     0,1,0,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,2,0,0;
%     0,0,0,1,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,2,0,0,0,0,0,0,1,0,0,0,0,1,0;
%     0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
%             ];

b = 2*[0.5, 0.0253, 0.3300, 0.0036, 0.0764, 0.1734, 0.0457, 0.1018, 0.1358, 0.0389, 0.0691]';

[nExp, nMets] = size(metmat);

mdl2 = [];
mdl2.A = [
    sparse(metmat), speye(nExp), -speye(nExp);
    ];
mdl2.rhs =  b; 


mdl2.lb = zeros(nMets+nExp+nExp,1);
mdl2.ub = [ones(nMets, 1); 1000*ones(nExp*2,1)];

% ADD UPPER LIMITS FOR KNOWN LIPID CONFIGURATIONS
% 
mdl2.lb(find(strcmp(lipidClass, 'DGTS__16_0__14_0'))) = 0.0010;
mdl2.ub(find(strcmp(lipidClass, 'DGTS__16_0__14_0'))) = 0.0042;

mdl2.lb(find(strcmp(lipidClass, 'DGTS__18_2__18_3'))) = 0.0071;
mdl2.ub(find(strcmp(lipidClass, 'DGTS__18_2__18_3'))) = 0.0286;

mdl2.lb(find(strcmp(lipidClass, 'DGTS__18_3__18_3'))) = 0.0036;
mdl2.ub(find(strcmp(lipidClass, 'DGTS__18_3__18_3'))) = 0.146;

mdl2.lb(find(strcmp(lipidClass, 'DGTS__18_0__18_3'))) = 0.0019;
mdl2.ub(find(strcmp(lipidClass, 'DGTS__18_0__18_3'))) = 0.0075;

mdl2.lb(find(strcmp(lipidClass, 'DGTS__18_1__18_2'))) = 0.0100;
mdl2.ub(find(strcmp(lipidClass, 'DGTS__18_1__18_2'))) = 0.0398;
% 
mdl2.lb(find(strcmp(lipidClass, 'DGTS__16_0__18_3v2'))) = 0.0215;
mdl2.ub(find(strcmp(lipidClass, 'DGTS__16_0__18_3v2'))) = 0.0859;
% 
mdl2.lb(find(strcmp(lipidClass, 'DGTS__16_0__18_3'))) = 0.0128;
mdl2.ub(find(strcmp(lipidClass, 'DGTS__16_0__18_3'))) = 0.0511;
% 
mdl2.lb(find(strcmp(lipidClass, 'DGTS__16_0__18_2'))) = 0.0463;
mdl2.ub(find(strcmp(lipidClass, 'DGTS__16_0__18_2'))) = 0.1851;

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

writecell(TableOut, 'C:\Users\Bisbii\OneDrive - Universidade do Minho\Algae\Models\Dsalina\lipids\dgts_opt_results.tsv', 'FileType', 'text','Delimiter','tab');
return_val = 'return';
end