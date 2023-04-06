
% Define the objective function


% Set the genetic algorithm parameters

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
% INFO.vhpo4max

x0 = [45, 0, 0.10, 0.17, 6.5 * 10^-2, 1 * 10^-5, 40, 50, 2, 4, 0.05, 420/1000*24, 18*10^-3*24];
lb = [0 -0.1 0.05 0 0 0 0 0 0 0 0 0 0];
ub = [200 1 1 1 1 1 200 200 10 10 0.3 300 1];

%options = optimoptions('ga', 'PopulationSize',1,'Generations',1,'MutationFcn',@mutationadaptfeasible, 'UseParallel',true); %
%options.InitialPopulationMatrix = x0;
options = optimset('MaxIter',10, 'Display', 'iter');
%options = optimoptions('fmincon', 'MaxIter',2,'Display', 'iter','UseParallel',true);


% Run the genetic algorithm
[x,fval] = fminsearch(@myFitness_ga, x0, options);
%[opt_params, obj_val] = fmincon(@myFitness_ga, x0, [], [], [], [], lb, ub, [], options);
%[x,fval] = ga(@myFitness_ga, length(x0),[],[],[],[],lb,ub,[],options);

% Display the optimized parameters
disp(x);