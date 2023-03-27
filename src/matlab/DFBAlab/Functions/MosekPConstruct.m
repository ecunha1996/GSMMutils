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

function CPLEXobjs = MosekPConstruct(model,nmodel,C,ncost,sense,tol,tolPh1,solver,disp,error)
global CPLEXobjs
for k=1:nmodel

    A = model{k}.A;
    b = model{k}.b;
    lb = model{k}.lb;
    ub = model{k}.ub;
    
    % Pack LP parameters
    param.MSK_IPAR_LOG = 0; % No display
    param.MSK_DPAR_INTPNT_TOL_PFEAS = tol; % Primal feasibility   
    param.MSK_DPAR_INTPNT_TOL_DFEAS = tol; % Dual feasibility
    param.MSK_DPAR_BASIS_TOL_X = tol; % Primal feasibility
    param.MSK_DPAR_BASIS_TOL_S = tol; % Dual feasibility
    param.MSK_IPAR_OPTIMIZER = 0; 
    
    CPLEXobjs{k,1}.param = param;
    CPLEXobjs{k,1}.c = C{k}(:,1);
    CPLEXobjs{k,1}.blx = lb;
    CPLEXobjs{k,1}.bux = ub;
    CPLEXobjs{k,1}.a = sparse(A);
    CPLEXobjs{k,1}.blc = b;
    CPLEXobjs{k,1}.buc = b;

    for i=2:ncost(k)
    if(sense{k}(i)==-1)
        C{k}(:,i) = -C{k}(:,i);
    end
    CPLEXobjs{k,i}.c = C{k}(:,i);
    CPLEXobjs{k,i}.blx = lb;
    CPLEXobjs{k,i}.bux = ub;
    CPLEXobjs{k,i}.a = sparse([A;(C{k}(:,1:i-1))']);
    CPLEXobjs{k,i}.blc = [b; zeros(i-1,1)];
    CPLEXobjs{k,i}.buc = [b; zeros(i-1,1)];
    end

    [flux] = MosekLexicographicPSolve(k,b,ncost(k),tolPh1,error,0);
    if disp == 1
        for i=2:ncost(k)
           flux(i) = sense{k}(i)*flux(i); 
        end
        display(flux)
    end
    CPLEXobjs{k,1}.param.MSK_IPAR_OPTIMIZER = solver; 
end
return
