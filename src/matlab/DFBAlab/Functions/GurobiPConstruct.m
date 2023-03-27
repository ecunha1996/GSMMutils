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

function CPLEXobjs = GurobiPConstruct(model,nmodel,C,ncost,sense,tol,tolPh1,solver,disp,error)
global CPLEXobjs
for k=1:nmodel

    A = model{k}.A;
    b = model{k}.b;
    lb = model{k}.lb;
    ub = model{k}.ub;
    
    % Pack LP parameters
    params.FeasibilityTol = tol;
    params.OptimalityTol = tol;
    params.outputflag = 0;
    params.method = -1;
    
    CPLEXobjs{k,1}.params = params;
    CPLEXobjs{k,1}.obj = C{k}(:,1);
    CPLEXobjs{k,1}.lb = lb;
    CPLEXobjs{k,1}.ub = ub;
    CPLEXobjs{k,1}.A = sparse(A);
    CPLEXobjs{k,1}.rhs = b;
    CPLEXobjs{k,1}.sense = ['='];

    for i=2:ncost(k)
    if(sense{k}(i)==-1)
        C{k}(:,i) = -C{k}(:,i);
    end
    CPLEXobjs{k,i}.obj = C{k}(:,i);
    CPLEXobjs{k,i}.lb = lb;
    CPLEXobjs{k,i}.ub = ub;
    CPLEXobjs{k,i}.A = sparse([A;(C{k}(:,1:i-1))']);
    CPLEXobjs{k,i}.rhs = [b; zeros(i-1,1)];
    CPLEXobjs{k,i}.sense = ['='];
    end

    [flux] = GurobiLexicographicPSolve(k,b,ncost(k),tolPh1,error,0);
    if disp == 1
        for i=2:ncost(k)
           flux(i) = sense{k}(i)*flux(i); 
        end
        display(flux)
    end
    CPLEXobjs{k,1}.params.method = solver;
end
return
