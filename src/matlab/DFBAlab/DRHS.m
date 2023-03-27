%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DFBAlab: Dynamic Flux Balance Analysis laboratory                       %
% Process Systems Engineering Laboratory, Cambridge, MA, USA              %
% July 2014                                                               %
% Written by Jose A. Gomez and Kai Höffner                                %
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

function dy = DRHS(t, y, INFO)
% Initial conditions
% Y1 = Volume (L)
% Y2 = Biomass (gDW/L)
% Y3 = CO2
% Y4 = Acetate
% Y5 = O2
% Y6 = Penalty

% Assign values from states
nmodel = INFO.nmodel;

% Feed rates
Fin = 0.006;
Fout = 0.006;



q = y(10) / 6.35;
n = 1 - (q / (q + 3.5));  %0.049


%% Update bounds and solve for fluxes
INFO.t = t;
[flux, penalty] = solveModel(t, y, INFO);
%%

%% Dynamics

dy = zeros(16,1);    % a column vector
dy(1) = Fin-Fout;    % Volume

for i=1:nmodel
    growth_rate_active = flux(i,1)*(1-n);
    dy(3) = flux(i,8)*y(2); % P mM
    dy(4) = flux(i,9)*y(2); % N mM
    dy(5) = growth_rate_active*y(2); % Active Biomass    
    dy(10) = -flux(i,9) - growth_rate_active*y(10); %nitrogen quota
    dy(11) = (flux(i, 3) * 893.49/1000 + flux(i, 4)*907.49/1000) - growth_rate_active*y(11); % chlorophyll quota
    dy(12) = flux(i,2)*48660.195/1000  - growth_rate_active*y(12);
    dy(13) = flux(i,7)*92.09/1000 - growth_rate_active*y(13);  % glycerol quota
    dy(14) = flux(i,5)*536.87/1000 - growth_rate_active*y(14); % carotene quota
    dy(15) = flux(i,6)*904.78/1000 - growth_rate_active*y(15); % tag quota
    dy(6) = dy(12)*y(2); %Starch /1000 (mmol/gdWd -> g/gDWd) *1000 (obj. wts) -> g/L
    dy(7) = dy(14)*y(2); % carotene concentration (g/L)
    dy(8) = dy(15)*y(2); % TAG concentration
    dy(9) = dy(13)*y(2); % glycerol concentration
    dy(2) = dy(5) + dy(6) + dy(7) + dy(8) + dy(9);
    dy(16) = dy(16) + penalty(i);
end
%dy(8) = dy(8) + Sfeed(1)*Fin + MT(1);
end