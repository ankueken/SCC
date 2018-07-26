function [x,Tol] = fractional_LP(N,Tol_Lower,Tol_Upper,Min_t,Tol_t,bs,idx1,idx2,obj)
%% fractional linear programming to find ratio between two reactions r_idx1/r_idx2
% Input:  N: stoichiometric matrix (mxn)
%         Tol_Lower: lower bounds on variables (nx1)
%         Tol_Upper: upper bounds on  variables (nx1)
%         Min_t: lower bound on t (1x1)
%         Tol_t: upper bound on t (1X1)
%         bs: minimum flux below (-1) or above 1 (1) (nx1)
%         idx1: index reaction 1
%         idx2: index reaction 2
%         obj: minimize objective (1), maximize objective (-1)

% Output: x: flux distribution
%         Tol: solver status
        
y_i = zeros(1,size(N,2));
y_i(idx1) = 1;

N_opt = [N zeros(size(N,1),1);  % steady state
    y_i bs(find(y_i==1))  % y_j + bt = 1, b=0
    ];

X_LB = [ones(size(N,2),1).*Tol_Lower; Min_t];
X_UB = [ones(size(N,2),1).*Tol_Upper; Tol_t];
b_LB = [zeros(size(N,1),1);
    1];
b_UB = [zeros(size(N,1),1);
    1];
c_opt = zeros(size(N,2)+1,1);
c_opt(idx2) = 1;

[y, ~, ~,~, ~, ~, ~, Tol] = cplex(obj*c_opt, N_opt, X_LB, X_UB, b_LB, b_UB);
if y(end)==0
    x = y;
else
    x = (1/y(end)).*y(1:size(N,2));
end
end
                        