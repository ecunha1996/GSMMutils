
% Define the objective function


% Set the genetic algorithm parameters
options = optimoptions('ga', 'PopulationSize',100,'Generations',50,'MutationFcn',@mutationadaptfeasible);

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

x0 = [45, 0, 0.12, 0.17, 6.5 * 10^-2, 1 * 10^-5, 50, 40, 2, 4];
lb = [0 0 0 0 0 0 0 0 0 0];
ub = [100 1 1 1 1 1 100 100 10 10];


% Run the genetic algorithm
[x,fval] = ga(@(x) myFitness(x, x0),length(x0),[],[],[],[],lb,ub,[],options);

% Display the optimized parameters
disp(x);