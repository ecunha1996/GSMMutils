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

function [CPLEXobjs,INFO] = ModelSetup(model,Y0,INFO)
global CPLEXobjs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subroutine to create full row rank model.
% Inputs
% model: Biological model.
% DB: Default bound on model (DB = Infinity)
% Outputs:
% model: Full row rank model with Infinity in place of DB. 
% Unpack parameters
nmodel = INFO.nmodel;
DB = INFO.DB;
Cost = INFO.C;
exID = INFO.exID;
optionslacks = INFO.optionslacks;
solver = INFO.solver;
LPsolver = INFO.LPsolver;
error = INFO.error;
Method = INFO.Method;
tol = INFO.tol;
tolPh1 = INFO.tolPh1;
disp = INFO.disp;

[model] = FullRank(model, nmodel, DB);

lexID = zeros(1,nmodel);
ncost = zeros(1,nmodel);
for i=1:nmodel
   ncost(i) = length(Cost{i})+1;
   lexID(i) = length(exID{i});
   a = zeros(length(model{i}.rxns),1);
   for j=1:lexID(i)
      a(exID{i}(j)) = 1; 
   end
   % This part is just adding elements of the cost vector that are not in
   % exID to the exID array such that LPstandardform does not redefine
   % them.
   LlexID = lexID; % LlexID will be the real length of exID now.
   for j=1:ncost(i)-1
      for k=1:length(Cost{i}(j).rxns)
        if a(Cost{i}(j).rxns(k))==0
           exID{i} = [exID{i}  Cost{i}(j).rxns(k)];
           LlexID(i) = LlexID(i) + 1;
        end
      end
   end
   % Store upper and lower bounds of the model.
   lb{i} = model{i}.lb;
   ub{i} = model{i}.ub;
end

% INFO structure: add more parameters
INFO.ncost = ncost;
INFO.lexID = lexID;
INFO.LlexID = LlexID;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subroutine to construct Phase I problem.
% A = Standard form stoichiometry matrix (equivalent to S).
% b = Standard form right-hand side.
% cost = Standard form cost vector (negative of the one in the model).
% C = Cost matrix. Each column contains a cost vector for each level of the
% lexicographic optimization. Only the first column has been set up.
% lb = lower bound on varibles. Vector of zeros because problem is in
% standard form.
% ub = upper bound on variables. Infinity because it is in standard form.
% pair = Vector containing pairings for variables that can be positive or
% negative. This vector is useful for cost vectors involving reactions on
% exID.
[model,C,pair] = PhaseI(model,nmodel,exID,ncost,optionslacks);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Restore upper and lower bounds to default values for reactions that were
% not originally in exID. Count non infinities added.
for i=1:nmodel
    lbct(i) = 0;
    ubct(i) = 0;
   for j=lexID(i)+1:length(exID{i})
        model{i}.b(j) = lb{i}(exID{i}(j));
        if model{i}.b(j) ~= -Inf
           lbct(i) = lbct(i) + 1; 
        end
        model{i}.b(length(exID{i})+j) = ub{i}(exID{i}(j));
        if model{i}.b(length(exID{i})+j) ~= Inf
           ubct(i) = ubct(i) + 1; 
        end
   end
end 
% Pack counters
INFO.lbct = lbct;
INFO.ubct = ubct;
% Cost vectors for problem in standard form. 
% The first column of the cost vector contains the Phase I objective 
% (minimizing infeasibilities).
% Structure:
% C{model} is a matrix. Each column contains a cost vector. First column
% contains the phase I const vector and next are the following objectives
% in order. The pair value are for those variables that get split due to 
% standard form transformations. 
for i=1:nmodel
    sense{i}(1) = 1; % Phase I minimizes
   for j=1:ncost(i)-1
       sense{i}(j+1) = Cost{i}(j).sense;
       for k=1:length(Cost{i}(j).rxns)
          C{i}(Cost{i}(j).rxns(k),j+1) = Cost{i}(j).wts(k);
          if pair{i}(Cost{i}(j).rxns(k))~=0
                C{i}(pair{i}(Cost{i}(j).rxns(k)),j+1) = -Cost{i}(j).wts(k);
          end
       end
   end
end

% Modify fixed bounds of exID. LPstandardform assigns the default bounds
% on the model. You can use the vex vector to update bounds. 

[lbx,ubx]=RHS(0,Y0,INFO);
for j = 1:nmodel
    for i=1:lexID(j)
       model{j}.b(i) = lbx(j,i);
       model{j}.b(length(exID{j})+i) = ubx(j,i);
    end
end

for i=1:nmodel
   clear('ind');
   ind = find(model{i}.b==Inf);
   model{i}.A(ind,:) = [];
   model{i}.b(ind) = [];
   clear('ind');
   ind = find(model{i}.b==-Inf);
   model{i}.A(ind,:) = [];
   model{i}.b(ind) = [];
   b{i} = model{i}.b;
end

% Solver dependent section
if LPsolver == 0
    if solver == 3 % Barrier method
        solver = 4;
    elseif solver == 0 || solver == 1 || solver == 2
    else
        fprintf('Solver code not supported. Continuing with CPLEX automatic solver')
        solver = 0;
    end
    if Method==1 % Solve the primal problem
        CPLEXobjs = CplexPConstruct(model,nmodel,C,ncost,sense,tol,tolPh1,solver,disp,error);
    else % Solve the dual problem
        CPLEXobjs = CplexConstruct(model,nmodel,C,ncost,sense,tol,tolPh1,solver,disp,error);
    end
elseif LPsolver == 1 % Gurobi
    solver = solver - 1;
    if solver == -1 || solver == 0 || solver == 1 || solver == 2
    else
       fprintf('Solver code not supported. Continuing with Gurobi automatic solver')
       solver = -1;
    end     
    if Method==1 % Solve the primal problem
        CPLEXobjs = GurobiPConstruct(model,nmodel,C,ncost,sense,tol,tolPh1,solver,disp,error);
    else % Solve the dual problem
        CPLEXobjs = GurobiConstruct(model,nmodel,C,ncost,sense,tol,tolPh1,solver,disp,error);
    end
elseif LPsolver == 2 % Mosek
    if solver == 0
    elseif solver == 1
        solver = 3;
    elseif solver == 2
        solver = 4;
    elseif solver == 3
        solver = 1;
    else
       fprintf('Solver code not supported. Continuing with Mosek automatic solver')
       solver = 0;
    end    
    if Method==1 % Solve the primal problem
        CPLEXobjs = MosekPConstruct(model,nmodel,C,ncost,sense,tol,tolPh1,solver,disp,error);
    else % Solve the dual problem
        CPLEXobjs = MosekConstruct(model,nmodel,C,ncost,sense,tol,tolPh1,solver,disp,error);
    end
else
    error('Solver code not supported. Pick a solver: CPLEX = 0, Gurobi = 1, Mosek = 2');
end

% Finally, b cell and sense cell.
INFO.sense = sense;
INFO.b = b;
INFO.pair = pair;
end

