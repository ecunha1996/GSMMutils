function fitness = myFitness_ga(x0)
% Simulate the kinetic model with the given parameters
matrix = readExpMatrix('C:\Users\Bisbii\PythonProjects\ExpGSMM\data\experimental\Matriz- DCCR Dunaliella salina_dfba.xlsx');
disp(x0);
INFO.ro1 = x0(1);
INFO.ro0 = x0(2);
INFO.wPmin = x0(3);
INFO.wPopt = x0(4);
INFO.a0 = x0(5);
INFO.a1 = x0(6);
INFO.a2 = x0(7);
INFO.a3 = x0(8);
INFO.l = x0(9);
INFO.smoothing_factor = x0(10);
INFO.vhpo4max = 0.001;
INFO.ExA = x0(11);
INFO.vcarmax = x0(12);
INFO.c0 = x0(13);
ssq = 0;
parpool(4);
parfor z=1:3
    try
        res = evaluateTrial(matrix, z, INFO);
    catch
        disp('Problem infeasible');
        res = 1e10;
    end
    ssq = ssq + res;
end
delete(gcp);
disp('SSQ:')
disp(ssq);
fitness = ssq;
end