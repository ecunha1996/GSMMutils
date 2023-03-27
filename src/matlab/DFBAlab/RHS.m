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

function [lb,ub] = RHS( ~,y,INFO )
nmodel = INFO.nmodel;

% Initial conditions
% Y1 = Volume (L)
% Y2 = Biomass (gDW/L)
% Y3 = CO2
% Y4 = O2
% Y5 = Penalty

% This subroutine updates the upper and lower bounds for the fluxes in the
% exID arrays in main. The output should be two matrices, lb and ub. The lb matrix
% contains the lower bounds for exID{i} in the ith row in the same order as
% exID. The same is true for the upper bounds in the ub matrix.
% Infinity can be used for unconstrained variables, however, it should be 
% fixed for all time. 

lb = zeros(nmodel, 2);
ub = zeros(nmodel, 2);
    for j = 1:nmodel
        %Biomass auto
        lb(j, 1) = 0;
        ub(j, 1) = Inf;
    
        % Light
        L = 0.2; % meters depth of pond
        biomass = y(2);
             
        %Ke = Ke1 + Ke2*(biomass);
        Ke = 11.5*biomass*y(11);
        Io = 178;
        lb(j,2) = 0;
        E = Io*(1-exp(-L*Ke))/(Ke*L);
        ro1 = 10;%3.2;
        ro0 = 0;
        ro = (ro1*y(11) + ro0);
        light_uptake = ro / (biomass*L) * E * 5.7;

        lb(j,2) = 0;
        ub(j,2) = light_uptake;

        % P
        if y(1+nmodel+1) < 0
            lb(j,3) = 0;
        else
            lb(j,3) = -0.04*y(3)/(0.0015+y(3));
        end
        %lb(j,3) = -1000;
        ub(j,3) = Inf;

        % N
        if y(1+nmodel+1) < 0
            lb(j,4) = 0;
        else
            lb(j,4) = -4.07*y(4)/(0.011+y(4))*(1-(y(10)/6.78)); % y(10)
        end
        %lb(j,4) = -1000;
        ub(j,4) = Inf;


        % Starch accumulation
        q = y(10) / 6.35;
        n = 1 - (q / (q + 0.049));
        Tmax = 1.74;
        xstarch = y(12);
        T = 1/(1-xstarch);
        z = (T-1) / (Tmax-1);
        Ks = 0.034;
        vmax = 0.66/48660.195 * 1000;
        lb(j,5) = -vmax * y(6) / (y(6) + Ks)*z;
        ub(j,5) = 0.001 * (1-z);
    

        % chla
        ymax = 0.37; % max chl content: 0.0118 ; min N quota : 2.29 -> 0.0118 / (2.29*14/1000)
        KE = 12.5 * 5.7;
        yE = ymax * (KE/(light_uptake + KE));
        nitrogen_mass_quota = y(10) * 14.007 / 1000;
        sum_chl = yE - (y(11)/nitrogen_mass_quota);
        
        lb(j, 6) = sum_chl * 1.73/2.73;
        ub(j, 6) = Inf;

        lb(j, 7) = sum_chl/2.73;
        ub(j, 7) = Inf;

        %  Carotene
        l = 2;
        Exn = light_uptake^l;
        ExnA = (420/1000*24)^l;
        vcarmax = 18 * 10^-3 * 24;      %0.8 * 10^-3 * 24;
        vcargen = vcarmax * (Exn / (Exn + ExnA));
        a0 = 6.5 * 10^-2;
        a1 = 1 * 10^-5;
        vcar = vcargen * phi(a1 * light_uptake + a0 - y(10));
        lb(j, 8) = vcar;
        ub(j, 8) = Inf;

        % TAG
        rtag = 49/y(5)/904.78;
        lb(j, 9) = rtag * n;
        ub(j, 9) = Inf;


        % glycerol
        
        nacl = INFO.nacl; %g/L
        y0 = 0.966; % g glycerol / g chl
        m = 0.15; %  g glycerol / g chl * gnacl/L
        %nacl_production = (y0 + nacl*m) * y(11); %g gly / gDW
        
        %lb(j, 10) = nacl_production / 92.09*1000 * n; % mmol gly/ gDW
        wgly_max = 0.17; %https://doi.org/10.1016/j.biortech.2008.02.042
        max_production = (1e-5*nacl^2 + 0.002*nacl + 0.112) / y(2) * (1- y(13)/wgly_max);

        lb(j, 10) = max_production;
        ub(j, 10) = Inf;

        % CO2 consumption
        r_co2_max = 12.28; %
        Km = 0.873;
        %lb(j, 11) = -r_co2_max * (y(8) / (y(8)+ Km));
        lb(j, 11) = -r_co2_max * (1-z);
        %lb(j, 11) = -100;
        ub(j, 11) = Inf;

        % Nitrate intracellular metabolization
        vno3max = 0.19*24;
        wnmin = 2.29;
        lb(j, 12) = -vno3max * (1- wnmin/y(10));
        ub(j, 12) = inf;
    end
end

