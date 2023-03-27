%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DFBAlab: Dynamic Flux Balance Analysis laboratory                       %
% Process Systems Engineering Laboratory, Cambridge, MA, USA              %
% November 2014                                                           %
% Written by Jose A. Gomez                                                %
%                                                                         % 
% This code can only be used for academic purposes. When using this code  %
% please cite:                                                            %
%                                                                         %
% Gomez, J.A., Höffner, K. and Barton, P. I. (2014).                      %
% DFBAlab: A fast and reliable MATLAB code for Dynamic Flux Balance       %
% Analysis. BMC Bioinformatics, 15:409                                    %                                               %
%                                                                         %
% COPYRIGHT (C) 2014 MASSACHUSETTS INSTITUTE OF TECHNOLOGY                %
%                                                                         %
% Read the LICENSE.txt file for more details.                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [INFO] = bupdate(tint,Y0,INFO)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[lbx,ubx]=RHS(tint,Y0,INFO);
   
lexID = INFO.lexID;
nmodel = INFO.nmodel;
bmodel = INFO.b;
lbct = INFO.lbct;
indlb = INFO.indlb;
indub = INFO.indub;

for i=1:nmodel
    b = bmodel{i};
    lb = lbx(i,1:lexID(i));
    ub = ubx(i,1:lexID(i));
    lb(indlb{i}) = [];
    ub(indub{i})=[];
    b(1:length(lb)) = lb;
    b(length(lb)+lbct(i)+1:length(lb)+lbct(i)+length(ub)) = ub;
    INFO.b{i} = b;
end

end

