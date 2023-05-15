
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

x0 = [100, 0,    0.105,    0.17,    0.0650,    0.0000,   40.0000,   50.0000,    2.0000,    4.0000, 10.0800, 0.4320, 1e-5];
lb = [50 -0.1 0.05 0 0 0 0 0 0 0 0 0 0 0];
ub = [200 1 1 1 1 1 200 200 10 10 0.3 300 1 0.001];


%
%options = optimoptions('ga', 'PopulationSize',1,'Generations',1,'MutationFcn',@mutationadaptfeasible, 'UseParallel',true); %
%options.InitialPopulationMatrix = x0;
options = optimset('MaxIter',2, 'Display', 'iter');
%options = optimoptions('fmincon', 'MaxIter',2,'Display', 'iter','UseParallel',true);


% Run the genetic algorithm
[x,fval] = fminsearch(@myFitness_ga, x0, options);
%[opt_params, obj_val] = fmincon(@myFitness_ga, x0, [], [], [], [], lb, ub, [], options);
%[x,fval] = ga(@myFitness_ga, length(x0),[],[],[],[],lb,ub,[],options);

% Display the optimized parameters
disp(x);