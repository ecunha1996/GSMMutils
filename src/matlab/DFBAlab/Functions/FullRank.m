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

function [model] = FullRank(model,nmodel,DB)

for i=1:nmodel
    Ain=model{i}.S;
    bin=model{i}.b;
    cin=model{i}.c;
    lbin=model{i}.lb;
    ubin=model{i}.ub;

    lbin(find(lbin==-DB(i)))=-Inf;
    ubin(find(ubin==DB(i)))=Inf;

    % Reformulating Ax = b for an equivalent system using QR factorization.
    % AE = QR, E = E', Q = Q'
    % QRE'x = b
    % RE'x = Q'b
    % Since R is an upper triangular matrix, only the first r rows of R (where
    % r is the rank of R) will contain the relevant information.
    [Q,R,E]=qr(Ain);
    rnk=rank(full(R)); % this could be replaced by an approximation
    Ain=R(1:rnk,:)*E';
    bin=Q'*bin;
    bin=bin(1:rnk);

    model{i}.A = Ain;
    model{i}.b = bin;
    model{i}.c = cin;
    model{i}.lb = lbin;
    model{i}.ub = ubin;
end
return

