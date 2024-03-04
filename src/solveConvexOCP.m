function [d_x,d_u] = solveConvexOCP(xCand,uCand,Aeq,beq,Aineq,bineq,opts)
%SOLVECONVEXOCP - this function sets up the cost function and calls the
% optimal control problem solver Mosek
%
% Syntax:
%       [d_x,d_u] = solveConvexOCP(xCand,uCand,Aeq,beq,Aineq,bineq,opts)
%
% Input Arguments:
%       - xCand:            current candidate for the state vector
%       - uCand:            current candidate for the input vector
%       - Aeq:              constraint matrix - (linear) equality constraint 
%                           (sparse real matrix of dimension [nConEq,idxVars.nVars], 
%                           where nConEq denotes the number of equality constraints)
%       - beq:              right-hand side - (linear) equality constraints
%                           (real vector of dimension [nConEq,1])
%       - Aineq:            constraint matrix - (linear) inequality constraint
%                           (sparse real matrix of dimension 
%                           [nConIneq,idxVars.nVars] where nConIneq denotes the 
%                           number of inequality constraints)
%       - bineq:            right-hand side - (linear) inequality constraints
%                           (real vector of dimension [nConIneq,1])
%                         
%       - opts:             structure containing the algorithm settings
%
% Output Arguments:
%       - d_x:              optimal state trajectory
%       - d_u:              optimal inputs
%
%Other Varibles Used:-
%       - y: Number of varibles except states and inputs
%       - H_size: Size of H to be filled with Q and R
%       - H__fill: Part of H to be filled with Q and R
%       - H: Quadratic cost matrix of the cost function
%       - h: Linear term of the cost function
%       - z: Stacked vector of optimized variables obtained from Mosek
% ------------------------------------------------------------------------

y = opts.idxVars.nVars - (opts.idxVars.n_x + opts.idxVars.n_u);

% Only define part of H that needs to be filled with Q and R
H_size = (opts.nx + opts.nu) * opts.N;
H_fill = zeros(H_size,H_size);
    
H_fill(1:opts.nu, 1:opts.nu) = opts.R;
H_fill((H_size - opts.nx + 1):H_size, (H_size - opts.nx + 1):H_size) = opts.Q;

% Defining h vector
h = zeros(1,opts.idxVars.nVars);
h(opts.idxVars.penalty_slack_terminalCon) = opts.mu;
h(opts.idxVars.trustRegionRad) = opts.lambda;

for kk = 1:opts.N

    % Calculating h according to its indices of the optimization vector
    h(opts.idxVars.x(:,kk)) = -2*xCand(:,kk)'*opts.Q; 
    h(opts.idxVars.u(:,kk)) = -2*uCand(:,kk)'*opts.R;

    h(opts.idxVars.penalty_vrtCtrl_pos(:,kk)) = opts.mu;
    h(opts.idxVars.penalty_vrtCtrl_neg(:,kk)) = opts.mu;
    h(opts.idxVars.penalty_slack_stageCon(:,kk)) = opts.mu;

    % Filling rest of the H matrix
    if kk < opts.N
        start_index = opts.nu + (kk-1)*(opts.nx + opts.nu) + 1;
        H_fill(start_index:start_index+opts.nx-1, start_index:start_index+opts.nx-1) = opts.Q;
        start_index = start_index + opts.nx;
        H_fill(start_index:start_index+opts.nu-1, start_index:start_index+opts.nu-1) = opts.R;

    end
end

% Fill H with zeros to match dimensionality
H = blkdiag(H_fill,zeros(y,y));

% Get optimized variables through Mosek
z = wrapperMosekSOCP(H,h,Aeq,beq,Aineq,bineq,opts.lb,opts.ub,opts.idxVars);

% Segregate the states and inputs
d_x = z(opts.idxVars.x);
d_u = z(opts.idxVars.u);

end

