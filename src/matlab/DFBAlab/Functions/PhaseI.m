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

function [model,C,pair] = PhaseI(model,nmodel,exID,ncost,optionslacks)

for k=1:nmodel
    for i = 1:length(exID{k})
        model{k}.lb(exID{k}(i)) = -Inf;
        model{k}.ub(exID{k}(i)) = Inf;
    end

    [model{k},pair{k}]=LPstandardform(model{k},exID{k});
    % Phase I LP construction

    A = model{k}.A;
    b = model{k}.b;
    cost = model{k}.c;

    [m,n] = size(A);
    if optionslacks(k) ==1 % Full set of slacks
        for i = 1:m
           A(i,n+2*i-1) = 1;
           A(i,n+2*i) = -1;
           A(m+1,n+2*i-1) = 1;
           A(m+1,n+2*i) = 1;
           cost(n+2*i-1) = 0;
           cost(n+2*i) = 0;
        end
           A(m+1,n+2*m+1) = -1;
           cost(n+2*m+1) = 0;
           b(m+1) = 0;  

    else % Reduced set of slacks. 
        for i = 1:length(exID{k})
           A(i,n+i) = 1;
           A(m+1,n+i) = 1;
           A(length(exID{k})+i,n+length(exID{k})+i) = -1;
           A(m+1,n+length(exID{k})+i) = 1;
           cost(n+i) = 0;
           cost(n+length(exID{k})+i) = 0;
        end
        ct = 2*length(exID{k})+1;
% These slacks are not needed: Jose A. Gomez (Aug. 25, 2014)
        for i = 2*length(exID{k})+1:m
           if b(i)<0
                A(i,n+ct) = -1;
                A(m+1,n+ct) = 1;
                cost(n+ct) = 0;
                ct = ct + 1;
           elseif b(i)>0
                A(i,n+ct) = 1;
                A(m+1,n+ct) = 1;
                cost(n+ct) = 0;
                ct = ct + 1;
           end
        end
           A(m+1,n+ct) = -1;
           cost(n+ct) = 0;
           b(m+1) = 0;
    end

    lb = zeros(size(A,2),1);
    ub = Inf*ones(size(A,2),1);
    C{k} = zeros(size(A,2),ncost(k));
    % Cost matrix for phase I LP
    C{k}(size(A,2),1) = 1;
    model{k}.A = A;
    model{k}.b = b;
    model{k}.c = cost;
    model{k}.lb = lb;
    model{k}.ub = ub;
end
return