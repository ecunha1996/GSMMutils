
% Define the objective function


% Set the genetic algorithm parameters
options = optimoptions('ga', 'PopulationSize',2,'Generations',1,'MutationFcn',@mutationadaptfeasible, 'UseParallel',true);
%options = optimset('MaxIter',10, 'Display', 'iter');
% Define the search space
% INFO.ro1 = x0(1);
% INFO.ro0 = x0(2);
% INFO.wPmin = x0(3);
% INFO.wPopt = x0(4);
% INFO.a0 = x0(5);
% INFO.a1 = x0(6);
% INFO.a2 = x0(7);
% INFO.a3 = x0(8);
% INFO.l = x0(9);
% INFO.smoothing_factor = x0(10);
%INFO.vhpo4max

x0 = [45, 0, 0.12, 0.17, 6.5 * 10^-2, 1 * 10^-5, 50, 40, 2, 4, 0.05];
lb = [0 0 0 0 0 0 0 0 0 0 0];
ub = [200 1 1 1 1 1 200 200 10 10 0.3];


% Run the genetic algorithm
%[x,fval] = fminsearch(@myFitness, x0, options);
[x,fval] = ga(@myFitness_ga, length(x0),[],[],[],[],lb,ub,[],options);

% Display the optimized parameters
disp(x);