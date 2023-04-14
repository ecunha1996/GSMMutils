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


q = y(10) / 6.78;
n = 1 - (q / (q + 0.049));  %0.049


%% Update bounds and solve for fluxes
INFO.t = t;
[flux, penalty] = solveModel(t, y, INFO);
%%

%% Dynamics

dy = zeros(17,1);    % a column vector
dy(1) = Fin-Fout;    % Volume

for i=1:nmodel
    growth_rate_active = flux(i,1);
    chl_production  =  (flux(i, 3) * 893.49/1000 + flux(i, 4)*907.49/1000);
    starch_production = flux(i,2)*48660.195/1000;
    glycerol_production = flux(i,7)*92.09/1000;
    tag_production = flux(i,6)*904.78/1000;
    caro_production = flux(i,5)*536.87/1000;
    total_growth_rate = growth_rate_active + chl_production + starch_production + glycerol_production + tag_production + caro_production;

    dy(3) = flux(i,8)*y(2); % P mM
    dy(4) = flux(i,9)*y(2); % N mM
    dy(5) = growth_rate_active*y(2); % Active Biomass
    dy(10) = -flux(i,9) - growth_rate_active*y(10); %nitrogen quota 
    dy(11) = chl_production - total_growth_rate*y(11); % chlorophyll quota
    dy(12) = starch_production  - total_growth_rate*y(12);
    dy(13) = glycerol_production - total_growth_rate*y(13);  % glycerol quota
    dy(14) = caro_production - total_growth_rate*y(14); % carotene quota
    dy(15) = tag_production - total_growth_rate*y(15); % tag quota
    dy(16) = -flux(i,8) - total_growth_rate*y(16); % P quota     
    dy(6) = starch_production*y(2); %Starch /1000 (mmol/gdWd -> g/gDWd) *1000 (obj. wts) -> g/L
    dy(7) = caro_production*y(2); % carotene concentration (g/L)
    dy(8) = tag_production*y(2); % TAG concentration
    dy(9) = glycerol_production*y(2); % glycerol concentration
    dy(2) = dy(5) + dy(6) + dy(7) + dy(8) + dy(9);
    dy(17) = dy(17) + penalty(i);
end
%dy(8) = dy(8) + Sfeed(1)*Fin + MT(1);
end