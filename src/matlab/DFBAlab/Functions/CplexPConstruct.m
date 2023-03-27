%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DFBAlab: Dynamic Flux Balance Analysis laboratory                       %
% Process Systems Engineering Laboratory, Cambridge, MA, USA              %
% July 2014                                                               %
% Written by Jose A. Gomez                                                %
% Revised by Kai Höffner                                                  %
%                                                                         % 
% This code can only be used for academic purposes. When using this code  %
% please cite:                                                            %
%                                                                         %
% Gomez, J.A., Höffner, K. and Barton, P. I. (2014).                      %
% DFBAlab: A fast and reliable MATLAB code for Dynamic Flux Balance       %
% Analysis. BMC Bioinformatics, 15:409                                    % 
%                                                                         %
% COPYRIGHT (C) 2014 MASSACHUSETTS INSTITUTE OF TECHNOLOGY                %
%                                                                         %
% Read the LICENSE.txt file for more details.                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CPLEXobjs = CplexPConstruct(model,nmodel,C,ncost,sense,tol,tolPh1,solver,disp,error)
global CPLEXobjs
for k=1:nmodel

    A = model{k}.A;
    b = model{k}.b;
    lb = model{k}.lb;
    ub = model{k}.ub;
    
    CPLEXobjs{k,1} = Cplex();
    CPLEXobjs{k,1}.Model.sense='minimize';
    CPLEXobjs{k,1}.addCols(C{k}(:,1),[],lb,ub);
    CPLEXobjs{k,1}.addRows(b,A,b);
    CPLEXobjs{k,1}.DisplayFunc =[];
    CPLEXobjs{k,1}.Param.feasopt.tolerance.Cur =tol;
    CPLEXobjs{k,1}.Param.simplex.tolerances.feasibility.Cur = tol;
    CPLEXobjs{k,1}.Param.simplex.tolerances.optimality.Cur = tol;
    CPLEXobjs{k,1}.Param.barrier.convergetol.Cur = tol;
    
    for i=2:ncost(k)
    if(sense{k}(i)==-1)
        C{k}(:,i) = -C{k}(:,i);
    end
    CPLEXobjs{k,i} = Cplex();
    CPLEXobjs{k,i}.Model.sense='minimize';
    CPLEXobjs{k,i}.addCols(C{k}(:,i),[],lb,ub);
    CPLEXobjs{k,i}.addRows([b; zeros(i-1,1)],[A;(C{k}(:,1:i-1))'],[b; zeros(i-1,1)]);
    CPLEXobjs{k,i}.DisplayFunc =[];
    CPLEXobjs{k,i}.Param.feasopt.tolerance.Cur =tol;
    CPLEXobjs{k,i}.Param.simplex.tolerances.feasibility.Cur = tol;
    CPLEXobjs{k,i}.Param.simplex.tolerances.optimality.Cur = tol;
    CPLEXobjs{k,i}.Param.barrier.convergetol.Cur = tol;
    end

    [flux] = CplexLexicographicPSolve(k,b,ncost(k),tolPh1,error,0);
    if disp == 1
        for i=2:ncost(k)
           flux(i) = sense{k}(i)*flux(i); 
        end
        display(flux)
    end

    for i=1:ncost
       CPLEXobjs{k,i}.Param.lpmethod.Cur=solver;  
    end
end
return
