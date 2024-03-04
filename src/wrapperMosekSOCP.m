function z = wrapperMosekSOCP(H,h,Aeq,beq,Aineq,bineq,lb,ub,idxVars)
%WRAPPERMOSEKSOCP - this function reformulates the conic quadratic
% approximation of the non-convex OCP in the standard form of a SOCP as
% specified in https://docs.mosek.com/10.0/toolbox/tutorial-cqo-shared.html
% and interfaces Mosek.
%
% Syntax:
%       opts = initRealTimeMPC(benchmark,opts,params)
%
% Input Arguments:
%   - H:            quadratic cost matrix
%                   (sparse real matrix of dimension 
%                   [idxVars.nVars,idxVars.nVars])
%   - h:            coefficient vector of the linear term in the cost 
%                   function 
%                   ((dense) real vector of dimension [idxVars.nVars,1])
%   - Aeq:          constraint matrix - (linear) equality constraint 
%                   (sparse real matrix of dimension [nConEq,idxVars.nVars], 
%                   where nConEq denotes the number of equality constraints)
%   - beq:          right-hand side - (linear) equality constraints
%                   (real vector of dimension [nConEq,1])
%   - Aineq:        constraint matrix - (linear) inequality constraint
%                   (sparse real matrix of dimension 
%                   [nConIneq,idxVars.nVars] where nConIneq denotes the 
%                   number of inequality constraints)
%   - bineq:        right-hand side - (linear) inequality constraints
%                   (real vector of dimension [nConIneq,1])
%   - lb:           lower bound on optimization variables
%                   (real vector of dimension [idxVars.nVars,1])
%   - ub:           upper bound on optimization variables
%                   (real vector of dimension [idxVars.nVars,1])
%   - idxVars:      structure containting the indices of the optimization
%                   variables in the corresponding vector
%
% Output Arguments:
%   - z:    (interior-point) solution of the SOCP
%           (real vector of dimension [idxVars.nVars,1])
%
% ------------------------------------------------------------------------


% coefficient vector of the cost function -> include 
h(idxVars.trustRegionRad) = 0.5;  %% todo: replace by corresponding weight
prob.c = h;

% setup constraint matrix
nCon = size(Aeq,1) + size(Aineq,1);
prob.a = [Aeq; Aineq;sparse(1,idxVars.trustRegionRad,-2,1,idxVars.nVars)];
prob.blc = [beq; -inf(size(Aineq,1),1);-inf]; 
prob.buc = [beq; bineq;0];

% bounds on optimization variables
prob.blx = lb;
prob.bux = ub;

% quadratic cost (conversion to sparse since you provide dense matrices as
% an input)
Htril = tril(sparse(H));
[row,col,val] = find(Htril);
prob.qosubi = row;
prob.qosubj = col;
prob.qoval = 2*val;

% quadratic constraint - soft trust region penalty
prob.qcsubi = [idxVars.x(:)',idxVars.u(:)'];  % row indices
prob.qcsubj = [idxVars.x(:)',idxVars.u(:)'];  % column indices
prob.qcval = ones(1,idxVars.n_x+idxVars.n_u); % entries -> all 1
prob.qcsubk = (nCon+1)*ones(1,idxVars.n_x+idxVars.n_u); 
% number of the constraint the quadratic function is associated with (see
% Mosek Opt. Toolbox for Matlab 10.0.46 Optimization Tutorials Sec. 6.9)

% solve SOCP (for now: modeled as a (convex) QCQP)
[~,res] = mosekopt('minimize',prob);
z = res.sol.itr.xx;

end

