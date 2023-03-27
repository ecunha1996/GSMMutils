%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DFBAlab: Dynamic Flux Balance Analysis laboratory                       %
% Process Systems Engineering Laboratory, Cambridge, MA, USA              %
% July 2014                                                               %
% Written by Kai Höffner                                                  %
% Modified by Jose A. Gomez                                               %
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

function [model,pair] = LPstandardform(model,exID)
%LPSTANDARDFORM Summary of this function goes here
%   Detailed explanation goes here


Ain = model.A;
bin = model.b;
lbin = model.lb;
ubin = model.ub;

[m,n] = size(Ain);
A1=zeros(m,2*n);
A2=zeros(1,2*n);
b1=zeros(m,1);
b2=0;

rhslb=zeros(m,1);
rhsub=zeros(m,1);
for i = 1:n
   pair(i) = 0; % Only those fluxes that get redefined will get a pair 
                % value equal to the pair variable.
end
k=0;
l = 0;
% There are several different possibilities
for i=1:n
    A1(:,i) = Ain(:,i);    
    if (lbin(i)==-Inf)
        if (ubin(i)==Inf) 
            % variable is unconstrained
            l = l+1;
            A1(:,n+l) = -Ain(:,i);
            pair(i) = n+l;
        else
            % change sign of the variable and shift upper bound to zero            
            A1(:,i) = -Ain(:,i);        
            rhsub = rhsub - Ain(:,i)*ubin(i);
        end
    else % lb~=-Inf
        % shift finite lower bound to zero
        if (ubin(i)==Inf)
            % shift lower bound to zero
            rhslb = rhslb - Ain(:,i)*lbin(i);
% We include in this case the situation where both lower bound and
% upper bound are equal. This avoids introducing linearly dependent rows.
        elseif (ubin(i)>=lbin(i))
           rhslb = rhslb - Ain(:,i)*lbin(i);
           ubin(i) = ubin(i) - lbin(i);
           % add upper bound constraint
           k = k + 1;
           l = l + 1;
           A2(k,i) = 1;
           A2(k,n+l) = 1;
           b2(k) = ubin(i); 
        end
    end
end
if A2 == 0
    B = A1;
else
    B = [A1 ; A2];
end
A = B(:,1:n+l);
b = [ bin + rhslb + rhsub ; b2'];

nvar = length(exID);
[m,n] = size(A);
A = [zeros(2*nvar,n); A];
b = [zeros(2*nvar,1);b];
lb = zeros(size(A,2),1);
ub = Inf*ones(size(A,2),1);
for i = 1:nvar
   A(i,exID(i)) = 1;
   A(i,pair(exID(i))) = -1;
   A(i,n+i) = -1;
   b(i) = lbin(exID(i));
   A(nvar+i,exID(i)) = 1;
   A(nvar+i,pair(exID(i))) = -1;
   A(nvar+i,n+nvar+i) = 1;
   b(nvar+i) = ubin(exID(i));
end

model.A = A;
model.b = b;
model.lb = lb;
model.ub = ub;
end

