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

function [ flux ] = MosekLexicographicPSolve(k, b, ncost, tolPh1, error, t )
global CPLEXobjs
flux = zeros(ncost,1); % This array will hold the fluxes information.
CPLEXobjs{k,1}.blc = b;
CPLEXobjs{k,1}.buc = b;
for i=1:ncost
    [r,res] = mosekopt('minimize statuskeys(1) echo(0)',CPLEXobjs{k,i},CPLEXobjs{k,1}.param); 
    CPLEXobjs{k,i}.sol.bas = res.sol.bas;
    if res.sol.bas.prosta ~= 1
        solver = CPLEXobjs{k,1}.param.MSK_IPAR_OPTIMIZER;
            if solver == 4
               CPLEXobjs{k,1}.param.MSK_IPAR_OPTIMIZER=3;
            else
               CPLEXobjs{k,1}.param.MSK_IPAR_OPTIMIZER=4;
            end
       [r,res] = mosekopt('minimize statuskeys(1) echo(0)',CPLEXobjs{k,i},CPLEXobjs{k,1}.param); 
        CPLEXobjs{k,i}.sol.bas = res.sol.bas;
       if res.sol.bas.prosta ~= 1
          fprintf('Error in model %i, level %i. \n',k,i);
          fprintf('Simulation time = %d. \n',t);
          if error == 1           
             error('LPsolve:status', 'Solution of LP not optimal'); 
          end
       end
       CPLEXobjs{k,1}.param.MSK_IPAR_OPTIMIZER=solver;
    end
    flux(i)=res.sol.bas.pobjval;
    if i==1 && abs(flux(i))<tolPh1
        flux(i)=0;
    end

    if i<ncost
        CPLEXobjs{k,i+1}.blc = [b; flux(1:i)];
        CPLEXobjs{k,i+1}.buc = [b; flux(1:i)];
    end
end

return

