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

function [ rflux, penalty ] = ModelObjsSolve( INFO, lbx, ubx )
global CPLEXobjs
ncost = INFO.ncost;
lexID = INFO.lexID;
nmodel = INFO.nmodel;
tolPh1 = INFO.tolPh1;
Method = INFO.Method;
error = INFO.error;
bmodel = INFO.b;
sense = INFO.sense;
LlexID = INFO.LlexID;
lbct = INFO.lbct;
LPsolver = INFO.LPsolver;
t = INFO.t;

flux = zeros(nmodel,max(ncost));

if Method == 1
    for i = 1:nmodel
        b = bmodel{i};
        lb = lbx(i,1:lexID(i));
        ub = ubx(i,1:lexID(i));
        lb(find(lb==-Inf))=[];
        ub(find(ub==Inf))=[];
        b(1:length(lb)) = lb;
        b(length(lb)+lbct(i)+1:length(lb)+lbct(i)+length(ub)) = ub;
        if LPsolver == 0 % CPLEX
            [flux(i,1:ncost(i))] = CplexLexicographicPSolve(i,b,ncost(i),tolPh1,error,t);
        elseif LPsolver == 1 % Gurobi
            [flux(i,1:ncost(i))] = GurobiLexicographicPSolve(i,b,ncost(i),tolPh1,error,t);
        elseif LPsolver == 2 % Mosek
            [flux(i,1:ncost(i))] = MosekLexicographicPSolve(i,b,ncost(i),tolPh1,error,t);
        else
            error('Some error occurred: solver not matched correctly');
        end
            
         for j=2:ncost(i)
            flux(i,j) = sense{i}(j)*flux(i,j); 
        end
    end
else
    for i = 1:nmodel
        b = bmodel{i};
        lb = lbx(i,1:lexID(i));
        ub = ubx(i,1:lexID(i));
        lb(find(lb==-Inf))=[];
        ub(find(ub==Inf))=[];
        b(1:length(lb)) = lb;
        b(length(lb)+lbct(i)+1:length(lb)+lbct(i)+length(ub)) = ub;
        if LPsolver == 0 % CPLEX
            [flux(i,1:ncost(i))] = CplexLexicographicSolve(i,b,ncost(i),tolPh1,error,t);
        elseif LPsolver == 1 % Gurobi
            [flux(i,1:ncost(i))] = GurobiLexicographicSolve(i,b,ncost(i),tolPh1,error,t);
        elseif LPsolver == 2 % Mosek
            [flux(i,1:ncost(i))] = MosekLexicographicSolve(i,b,ncost(i),tolPh1,error,t);
        else
            error('Some error occurred: solver not matched correctly');
        end
        flux(i,:) = -1*flux(i,:);
        for j=2:ncost(i)
            flux(i,j) = sense{i}(j)*flux(i,j); 
        end
    end
end
   penalty = flux(:,1);
   rflux = flux(:,2:size(flux,2));
end

