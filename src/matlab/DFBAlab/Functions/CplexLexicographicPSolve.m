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
% Analysis. BMC Bioinformatics, 15:409                                    % %
%                                                                         %
% COPYRIGHT (C) 2014 MASSACHUSETTS INSTITUTE OF TECHNOLOGY                %
%                                                                         %
% Read the LICENSE.txt file for more details.                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ flux ] = CplexLexicographicPSolve(k, b, ncost, tolPh1, error, t )
global CPLEXobjs
flux = zeros(ncost,1); % This array will hold the fluxes information.
CPLEXobjs{k,1}.Model.rhs = b;
CPLEXobjs{k,1}.Model.lhs = b;
for i=1:ncost
    CPLEXobjs{k,i}.solve();
        if CPLEXobjs{k,i}.Solution.status ~= 1
            solver = CPLEXobjs{k,i}.Param.lpmethod.Cur;
            if solver == 2
               CPLEXobjs{k,i}.Param.lpmethod.Cur=1;
            else
               CPLEXobjs{k,i}.Param.lpmethod.Cur=2;
            end
            CPLEXobjs{k,i}.solve();
               if CPLEXobjs{k,i}.Solution.status ~= 1
                    fprintf('Error in model %i, level %i. \n',k,i);
                    fprintf('Simulation time = %d. \n',t);
                    fprintf('LP solution status: %i. \n',CPLEXobjs{k,i}.Solution.status);
                    fprintf('Status string: ');
                    display(CPLEXobjs{k,i}.Solution.statusstring);
                    if error == 1           
                        error('LPsolve:status', 'Solution of LP not optimal'); 
                    end
               end
            CPLEXobjs{k,i}.Param.lpmethod.Cur=solver;
        end
    flux(i)=CPLEXobjs{k,i}.Solution.objval;
    if i==1 && abs(flux(i))<tolPh1
        flux(i)=0;
    end

    if i<ncost
        CPLEXobjs{k,i+1}.Model.rhs = [b; flux(1:i)];
        CPLEXobjs{k,i+1}.Model.lhs = [b; flux(1:i)];
    end
    
end

return

