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

function [INFO] = LexicographicOptG(model,INFO)
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
        
        INFO.params.FeasibilityTol = INFO.tol;
        INFO.params.OptimalityTol = INFO.tol;
        INFO.params.outputflag = 0;
        INFO.params.method = 0;
        
        INFO.CPLEXlex{k}.obj = full(C(:,1));
        INFO.CPLEXlex{k}.lb = model{k}.lb;
        INFO.CPLEXlex{k}.ub = model{k}.ub;
        INFO.CPLEXlex{k}.A = A;
        INFO.CPLEXlex{k}.rhs = INFO.b{k};
        INFO.CPLEXlex{k}.sense = ['='];
        
        CPLEXlex = INFO.CPLEXlex{k};
        
        lexid = length(INFO.exID{k});
        rows = size(A,1);
        cols = size(A,2);
        basis = -1*ones(cols,1);
        null = [];
        j = cols;
            for i=rows:-1:2*lexid+1
                flag = 0;
                while flag ==0
                    if A(i,j) ~= 0
                        basis(j) = 0;
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
                  basis(ind(end-1))=0;
               else
                  basis(ind(end)) = 0;
               end
            end

        CPLEXlex.vbasis = basis;
        CPLEXlex.cbasis = -1*ones(rows,1);
        INFO.CPLEXlex{k}.cbasis = CPLEXlex.cbasis;
        result = gurobi(CPLEXlex,INFO.params);
        basis = result.vbasis;
        reducedcost = result.rc;
        obj(1) = result.objval;

    % Three cases: unique solution, linearly dependent cost vector, or need
    % to continue solving LPs;
         totind = [];
         for i=2:INFO.ncost(k)
            ind = find(reducedcost>INFO.tol);
            l = length(ind);
            totlength = length(totind);
            if l == 0 % Cost vector linearly dependent
               CPLEXlex.obj = full(C(:,i));
            elseif l == cols - rows - totlength % Unique solution
                break;
            else % Need to add constraint and eliminate null variables
                [~,j] = max(reducedcost);
                ind(find(ind==j))=[];
                totind = [totind ind'];
                A(:,ind) = zeros(size(A,1),l-1);
                C(ind,:) = zeros(l-1,INFO.ncost(k));
                null(end+1) = j;
                basis(j) = 0;
                CPLEXlex.obj = full(C(:,i));
                A(end+1,:) = C(:,i-1)';
                CPLEXlex.A = A;
                CPLEXlex.rhs(end+1) = obj(i-1);
                CPLEXlex.cbasis = -1*ones(size(A,1),1);
                CPLEXlex.vbasis = basis;
            end
            result = gurobi(CPLEXlex,INFO.params);
            basis = result.vbasis;
            reducedcost = result.rc;
            obj(i) = result.objval;
        end
    % Store the basis matrix
        for i=1:length(null)
           basis(null(i)) = -1; 
        end
        INFO.basis{k} = basis;
        ind = find(basis==-1);
        B{k}(:,ind) = [];
        CB{k}(ind,:) = [];
        [L{k},U{k},P{k},Q{k}] = lu(B{k});
        % Store vectors to multiply times right hand-side and obtain fluxes.
        for i=1:INFO.ncost(k)
            zzz = CB{k}(:,i)'*Q{k}/U{k}/L{k};
            plambda{k}(:,i) = (zzz*P{k})';
        end
    end

    % Store important information
    INFO.B = B;
    INFO.CB = CB;
    INFO.L = L;
    INFO.U = U;
    INFO.P = P;
    INFO.Q = Q;
    INFO.plambda = plambda;
else
    k = INFO.flagbasis;
    CPLEXlex = INFO.CPLEXlex{k};
    CPLEXlex.rhs = INFO.b{k};
    A = model{k}.A;
    INFO.B{k} = A;
    C = model{k}.C;
    costvec = size(C,2);
    INFO.CB{k} = C;
    obj = zeros(INFO.ncost(k),1);
    
    rows = size(A,1);
    cols = size(A,2);
    A = [A;zeros(costvec,cols)];
    null = zeros(costvec,1);
    basis = INFO.basis{k};
    INFO.params.method = 1; % Dual feasible basis

    CPLEXlex.vbasis = basis;
    result = gurobi(CPLEXlex,INFO.params);
    TF = strcmp(result.status,'OPTIMAL');
    if TF == 0
        INFO.flagbasis = 0; 
        INFO = LexicographicOptG(model,INFO);
        return;            
    end
    basis = result.vbasis;
    reducedcost = result.rc;
    obj(1) = result.objval;
    INFO.params.method = 0; % Back to primal simplex

% Three cases: unique solution, linearly dependent cost vector, or need
% to continue solving LPs;
     totind = [];
     for i=2:INFO.ncost(k)
        ind = find(reducedcost>INFO.tol);
        l = length(ind);
        totlength = length(totind);
        if l == 0 % Cost vector linearly dependent
           CPLEXlex.obj = full(C(:,i));
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
            basis(j) = 0;
            CPLEXlex.obj = full(C(:,i));
            A(rows+1,:) = C(:,i-1)';
            rows = rows + 1;
            CPLEXlex.A = A(1:rows,:);
            CPLEXlex.rhs(end+1) = obj(i-1);
            costvec = costvec - 1;
            CPLEXlex.cbasis = -1*ones(rows,1);
            CPLEXlex.vbasis = basis;
        end
        result = gurobi(CPLEXlex,INFO.params);
        TF = strcmp(result.status,'OPTIMAL');
        if TF == 0
            INFO.flagbasis = 0; 
            INFO = LexicographicOptG(model,INFO);
            return;            
        end
        basis = result.vbasis;
        reducedcost = result.rc;
        obj(i) = result.objval;
    end
% Store the basis matrix
    for i=2:null(1)+1
       basis(null(i)) = -1; 
    end
    INFO.basis{k} = basis;
    ind = find(basis==-1);
    INFO.B{k}(:,ind) = [];
    INFO.CB{k}(ind,:) = [];
    [INFO.L{k},INFO.U{k},INFO.P{k},INFO.Q{k}] = lu(INFO.B{k});
    % Store vectors to multiply times right hand-side and obtain fluxes.
    for i=1:INFO.ncost(k)
         zzz = INFO.CB{k}(:,i)'*INFO.Q{k}/INFO.U{k}/INFO.L{k};
         INFO.plambda{k}(:,i) = (zzz*INFO.P{k})';
    end
end

