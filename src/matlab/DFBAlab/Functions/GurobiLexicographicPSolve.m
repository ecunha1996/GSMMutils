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

function [ flux ] = GurobiLexicographicPSolve(k, b, ncost, tolPh1, error, t )
global CPLEXobjs
flux = zeros(ncost,1); % This array will hold the fluxes information.
CPLEXobjs{k,1}.rhs = b;
for i=1:ncost
    result = gurobi(CPLEXobjs{k,i},CPLEXobjs{k,1}.params);
    if strcmp(result.status, 'OPTIMAL')
        CPLEXobjs{k,i}.vbasis = result.vbasis;
        CPLEXobjs{k,i}.cbasis = result.cbasis;
    else
        solver = CPLEXobjs{k,1}.params.method;
            if solver == 1
               CPLEXobjs{k,1}.params.method=0;
            else
               CPLEXobjs{k,1}.params.method=1;
            end
       result = gurobi(CPLEXobjs{k,i},CPLEXobjs{k,1}.params);
       if strcmp(result.status, 'OPTIMAL')
         CPLEXobjs{k,i}.vbasis = result.vbasis;
         CPLEXobjs{k,i}.cbasis = result.cbasis;
       else
          fprintf('Error in model %i, level %i. \n',k,i);
          fprintf('Simulation time = %d. \n',t);
          if error == 1           
             error('LPsolve:status', 'Solution of LP not optimal'); 
          end
       end
       CPLEXobjs{k,1}.params.method=solver;
    end
    flux(i)=result.objval;
    if i==1 && abs(flux(i))<tolPh1
        flux(i)=0;
    end
    if i<ncost
        CPLEXobjs{k,i+1}.rhs = [b; flux(1:i)];
    end
    
end

return

