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
% Analysis. BMC Bioinformatics, 15:409                                    % 
%                                                                         %
% COPYRIGHT (C) 2014 MASSACHUSETTS INSTITUTE OF TECHNOLOGY                %
%                                                                         %
% Read the LICENSE.txt file for more details.                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [INFO] = LexicographicOpt(model,INFO)
% This function implements Algorithm 2 in Harwood, Höffner, and Barton:
% "Solution of ordinary differential equations with a linear program 
% embedded: the right-hand side case".
if INFO.flagbasis == 0 
    for k = 1:INFO.nmodel
        A = model{k}.A;
        B{k} = A;
        C = model{k}.C;
        CB{k} = C;
        obj = zeros(INFO.ncost(k),1);
        CPLEXlex = Cplex();
        CPLEXlex.Model.sense='minimize';
        CPLEXlex.addCols(model{k}.C(:,1),[],model{k}.lb,model{k}.ub);
        CPLEXlex.addRows(INFO.b{k},A,INFO.b{k});
        CPLEXlex.DisplayFunc =[];
        CPLEXlex.Param.simplex.tolerances.feasibility.Cur = INFO.tol;
        CPLEXlex.Param.simplex.tolerances.optimality.Cur = INFO.tol;
        CPLEXlex.Param.lpmethod.Cur = 1;
        lexid = length(INFO.exID{k});
        rows = size(A,1);
        cols = size(A,2);
        basis = zeros(cols,1);
        null = [];
        j = cols;
            for i=rows:-1:2*lexid+1
                flag = 0;
                while flag ==0
                    if A(i,j) ~= 0
                        basis(j) = 1;
                        j = j-1;
                        flag = 1;
                    else
                        j = j-1;
                    end
                end
            end 
            for i =1:2*lexid
               ind = find(A(i,:)); 
               if (A(i,ind(end))*INFO.b{k}(i))<0
                  basis(ind(end-1))=1;
               else
                  basis(ind(end)) = 1;
               end
            end

        CPLEXlex.Start.basis.colstat = basis;
        CPLEXlex.Start.basis.rowstat = zeros(rows,1);
        CPLEXlex.solve();
        if CPLEXlex.Solution.status ~= 1
           fprintf('Error in model %i, level 1. \n',k);
           fprintf('LP solution status: %i. \n',CPLEXlex.Solution.status);
           fprintf('Status string: ');
           display(CPLEXlex.Solution.statusstring);
           error('Simulation stopped.');
        end
        CPLEXlex.Param.lpmethod.Cur = 1;
        basis = CPLEXlex.Solution.basis.colstat;
        reducedcost = CPLEXlex.Solution.reducedcost;
        obj(1) = CPLEXlex.Solution.objval;

    % Three cases: unique solution, linearly dependent cost vector, or need
    % to continue solving LPs;
        totind = [];
        for i=2:INFO.ncost(k)
            ind = find(reducedcost>INFO.tol);
            l = length(ind);
            totlength = length(totind);
            if l == 0 % Cost vector linearly dependent
                CPLEXlex.Model.obj = C(:,i);
            elseif l == cols - rows - totlength % Unique solution
                break;
            else % Need to add constraint and eliminate null variables
                [~,j] = max(reducedcost);
                ind(find(ind==j))=[];
                totind = [totind ind'];
                A(:,ind) = zeros(size(A,1),l-1);
                C(ind,:) = zeros(l-1,INFO.ncost(k));
                null(end+1) = j;
                basis(j) = 1;
                CPLEXlex.Model.obj = C(:,i);
                CPLEXlex.Model.A = A;
                A(rows+1,:) = C(:,i-1)';
                rows = rows + 1;
                CPLEXlex.addRows(obj(i-1),C(:,i-1)',obj(i-1));
                CPLEXlex.Start.basis.rowstat = zeros(size(A,1),1);
                CPLEXlex.Start.basis.colstat = basis;
            end

            CPLEXlex.solve();
            if CPLEXlex.Solution.status ~= 1
                fprintf('Error in model %i, level %i. \n',k,i);
                fprintf('LP solution status: %i. \n',CPLEXlex.Solution.status);
                fprintf('Status string: ');
                display(CPLEXlex.Solution.statusstring);
                error('Simulation stopped.');
            end
            basis = CPLEXlex.Solution.basis.colstat;
            reducedcost = CPLEXlex.Solution.reducedcost;
            obj(i) = CPLEXlex.Solution.objval;
        end
    % Store the basis matrix
        for i=1:length(null)
           basis(null(i)) = 0; 
        end
        bs{k} = basis;
        ind = find(basis==0);
        B{k}(:,ind) = [];
        CB{k}(ind,:) = [];
        [L{k},U{k},P{k},Q{k}] = lu(B{k});
        % Store vectors to multiply times right hand-side and obtain fluxes.
        for i=1:INFO.ncost(k)
            zzz = CB{k}(:,i)'*Q{k}/U{k}/L{k};
            plambda{k}(:,i) = (zzz*P{k})';
        end
        clear CPLEXlex
    end

    % Store important information
    INFO.basis = bs;
    INFO.B = B;
    INFO.CB = CB;
    INFO.L = L;
    INFO.U = U;
    INFO.P = P;
    INFO.Q = Q;
    INFO.plambda = plambda;
else
    k = INFO.flagbasis;
    A = model{k}.A;
    INFO.B{k} = A;
    C = model{k}.C;
    costvec = size(C,2);
    INFO.CB{k} = C;
    obj = zeros(INFO.ncost(k),1);
    CPLEXlex = Cplex();
    CPLEXlex.Model.sense='minimize';
    CPLEXlex.addCols(model{k}.C(:,1),[],model{k}.lb,model{k}.ub);
    CPLEXlex.addRows(INFO.b{k},A,INFO.b{k});
    CPLEXlex.DisplayFunc =[];
    CPLEXlex.Param.simplex.tolerances.feasibility.Cur = INFO.tol;
    CPLEXlex.Param.simplex.tolerances.optimality.Cur = INFO.tol;
    
    rows = size(A,1);
    cols = size(A,2);
    A = [A;zeros(costvec,cols)];
    null = zeros(costvec,1);
    basis = INFO.basis{k};
    CPLEXlex.Param.lpmethod.Cur = 2; % Dual feasible basis

    CPLEXlex.Start.basis.rowstat = zeros(rows,1);
    CPLEXlex.Start.basis.colstat = basis;
    CPLEXlex.solve();
    if CPLEXlex.Solution.status ~= 1
        INFO.flagbasis = 0; 
        INFO = LexicographicOpt(model,INFO);
        return;         
    end
    CPLEXlex.Param.lpmethod.Cur = 1;
    basis = CPLEXlex.Solution.basis.colstat;
    reducedcost = CPLEXlex.Solution.reducedcost;
    obj(1) = CPLEXlex.Solution.objval;

% Three cases: unique solution, linearly dependent cost vector, or need
% to continue solving LPs;
    totind = [];
    for i=2:INFO.ncost(k)
        ind = find(reducedcost>INFO.tol);
        l = length(ind);
        totlength = length(totind);
        if l == 0 % Cost vector linearly dependent
            CPLEXlex.Model.obj = C(:,i);
        elseif l == cols - rows - totlength % Unique solution
            break;
        else % Need to add constraint and eliminate null variables
            [~,j] = max(reducedcost);
            ind(find(ind==j))=[];
            totind = [totind ind'];
            A(:,ind) = zeros(size(A,1),l-1);
            C(ind,:) = zeros(l-1,INFO.ncost(k));
            null(null(1)+2) = j;
            null(1) = null(1) + 1;
            basis(j) = 1;
            CPLEXlex.Model.obj = C(:,i);
            CPLEXlex.Model.A = A(1:rows,:);
            A(rows+1,:) = C(:,i-1)';
            rows = rows + 1;
            costvec = costvec - 1;
            CPLEXlex.addRows(obj(i-1),C(:,i-1)',obj(i-1));
            CPLEXlex.Start.basis.rowstat = zeros(rows,1);   
            CPLEXlex.Start.basis.colstat = basis;
        end
        CPLEXlex.solve();
        if CPLEXlex.Solution.status ~= 1
            INFO.flagbasis = 0; 
            INFO = LexicographicOpt(model,INFO);
            return;         
        end
        basis = CPLEXlex.Solution.basis.colstat;
        reducedcost = CPLEXlex.Solution.reducedcost;
        obj(i) = CPLEXlex.Solution.objval;
    end
% Store the basis matrix
    for i=2:null(1)+1
       basis(null(i)) = 0; 
    end
    INFO.basis{k} = basis;
    ind = find(basis==0);
    INFO.B{k}(:,ind) = [];
    INFO.CB{k}(ind,:) = [];
    [INFO.L{k},INFO.U{k},INFO.P{k},INFO.Q{k}] = lu(INFO.B{k});
    % Store vectors to multiply times right hand-side and obtain fluxes.
    for i=1:INFO.ncost(k)
         zzz = INFO.CB{k}(:,i)'*INFO.Q{k}/INFO.U{k}/INFO.L{k};
         INFO.plambda{k}(:,i) = (zzz*INFO.P{k})';
    end
    clear CPLEXlex 
end

