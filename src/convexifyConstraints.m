function [Aeq,beq,Aineq,bineq] = convexifyConstraints(x0,x_cand,u_cand,opts)
%CONVEXIFYCONSTRAINTS - evalutes the system dynamics and constraint functions 
% and the respective jacobians at the current iterate with the endgoal of
% convexifying the respective constraints
% 
%
% Syntax:
%       [Aeq,beq,Aineq,bineq] = convexifyConstraints(x0,x_cand,u_cand,opts)
%
% Input Arguments:
%       - x0:               initial state
%       - x_cand:           current candidate for the state vector
%       - u_cand:           current candidate for the input vector
%       - opts:             a structure containing the algorithm settings 
%                           (see function getOptimizationVariables for more info)
%
% Output Arguments:
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


% fill the following matrices Aeq,Aineq and vectors beq,bineq with
% linearized constraints (see report on how exactly this is done)
Aeq = [];
beq = [];
Aineq = [];
bineq = [];



% get easier access to system dynamics and constraint functions
jacobianDynamics = opts.constraints.jac_dynamics;
dynamics = opts.constraints.dynamics;

stageCon = opts.constraints.stageCon;
jac_stageCon = opts.constraints.jac_stageCon;


%% take care of case x_0

% Aeq and beq:
J = full(jacobianDynamics(x0,u_cand(:,1)));

B_k = J(:,opts.nx+(1:opts.nu));

% initialize first row with zeros
matrix = zeros(opts.nx,opts.idxVars.nVars);

% fill first row (same principle for remaining rows of Aeq,beq,Aineq,bineq)
matrix(:,opts.idxVars.x(:,1)) = -eye(opts.nx);
matrix(:,opts.idxVars.u(:,1)) = B_k;
matrix(:,opts.idxVars.vrtCtrl_neg(:,1)) = -eye(opts.nx);
matrix(:,opts.idxVars.vrtCtrl_pos(:,1)) = eye(opts.nx);

vector = full(dynamics(x0,u_cand(:,1),x_cand(:,1)));

beq = [beq;-vector];

Aeq = [Aeq;matrix];

% Aineq and bineq:
% initialize with zeros
matrix = zeros(opts.constraints.num_stageConstraints,opts.idxVars.nVars);
vector = zeros(opts.constraints.num_stageConstraints,1);
    
% fill rows
for ii = 1:opts.constraints.num_stageConstraints
    K = full(stageCon{ii}(x0,u_cand(:,1)));
    J = full(jac_stageCon{ii}(x0,u_cand(:,1)));

    C_k = J(:,1:opts.nx);
    D_k = J(:,opts.nx+(1:opts.nu));

    % Aineq
    matrix(ii,opts.idxVars.x(:,1)) = C_k;
    matrix(ii,opts.idxVars.u(:,1)) = D_k;
    matrix(ii,opts.idxVars.slack_stageCon(:,1)) = -eye(1,opts.constraints.num_stageConstraints);
    % bineq
    vector(ii,1) = -K;
end

Aineq = [Aineq;matrix];
bineq = [bineq;vector];


%% setup Aeq and beq
for k = 2:opts.N

    J = full(jacobianDynamics(x_cand(:,k-1),u_cand(:,k)));
    
    A_k = J(:,1:opts.nx);
    B_k = J(:,opts.nx+(1:opts.nu));

    % initialize with zeros
    matrix = zeros(opts.nx,opts.idxVars.nVars);
  
    % fill row for Aeq
    matrix(:,opts.idxVars.x(:,k-1)) = A_k;
    matrix(:,opts.idxVars.x(:,k)) = -eye(opts.nx);
    matrix(:,opts.idxVars.u(:,k)) = B_k;
    matrix(:,opts.idxVars.vrtCtrl_neg(:,k)) = -eye(opts.nx);
    matrix(:,opts.idxVars.vrtCtrl_pos(:,k)) = eye(opts.nx);

    Aeq = [Aeq;matrix];

    % fill row for beq    
    vector = full(dynamics(x_cand(:,k-1),u_cand(:,k),x_cand(:,k)));

    beq = [beq;-vector];

end

%% setup Aineq and bineq:
for k=2:opts.N-1
    
    % initialize with zeros
    matrix = zeros(opts.constraints.num_stageConstraints,opts.idxVars.nVars);
    vector = zeros(opts.constraints.num_stageConstraints,1);
    
    % fill rows
    for ii = 1:opts.constraints.num_stageConstraints
        K = full(stageCon{ii}(x_cand(:,k-1),u_cand(:,k)));
        J = full(jac_stageCon{ii}(x_cand(:,k-1),u_cand(:,k)));
        
        C_k = J(:,1:opts.nx);
        D_k = J(:,opts.nx+(1:opts.nu));
        
        % fill row of Aineq
        matrix(ii,opts.idxVars.x(:,k)) = C_k;
        matrix(ii,opts.idxVars.u(:,k)) = D_k;
        matrix(ii,opts.idxVars.slack_stageCon(:,k)) = -eye(1,opts.constraints.num_stageConstraints);
        % fill row of bineq
        vector(ii,1) = -K;
    end

    Aineq = [Aineq;matrix];

    bineq = [bineq;-vector];

end

%% fill last row of Aineq,bineq according terminal constraint
matrix = zeros(opts.constraints.num_terminalConstraints,opts.idxVars.nVars);
vector = zeros(opts.constraints.num_terminalConstraints,1);

terminalCon = opts.constraints.terminalCon;
jac_terminalCon = opts.constraints.jac_terminalCon;

for ii = 1:opts.constraints.num_terminalConstraints

    K = full(terminalCon{ii}(x_cand(:,opts.N)));
    J = full(jac_terminalCon{ii}(x_cand(:,opts.N)));

    E_k = J(:,1:opts.nx);

    matrix(ii,opts.idxVars.x(:,opts.N)) = E_k;
    matrix(ii,opts.idxVars.slack_terminalCon(:,ii)) = -1;
    vector(ii,1) = -K;
end

Aineq = [Aineq;matrix];

bineq = [bineq;vector];


%% add reformulation of exact penalty terms to Aineq
for ii=1:opts.N

    matrix = zeros(opts.nx,opts.idxVars.nVars);

    matrix(:,opts.idxVars.vrtCtrl_pos(:,ii)) = eye(opts.nx);
    matrix(:,opts.idxVars.vrtCtrl_neg(:,ii)) = eye(opts.nx);
    matrix(:,opts.idxVars.penalty_vrtCtrl_pos(ii)) = -ones(opts.nx,1);
    matrix(:,opts.idxVars.penalty_vrtCtrl_neg(ii)) = -ones(opts.nx,1);

    Aineq = [Aineq;matrix];

end

for ii=1:opts.N-1
    matrix = zeros(opts.constraints.num_stageConstraints,opts.idxVars.nVars);

    matrix(:,opts.idxVars.slack_stageCon(:,ii)) = eye(opts.constraints.num_stageConstraints);
    
    matrix(:,opts.idxVars.penalty_slack_stageCon(ii)) = -ones(opts.constraints.num_stageConstraints,1);

    Aineq = [Aineq;matrix];

end

matrix = zeros(opts.constraints.num_terminalConstraints,opts.idxVars.nVars);

matrix(:,opts.idxVars.slack_terminalCon) = eye(opts.constraints.num_terminalConstraints);
matrix(:,opts.idxVars.penalty_slack_terminalCon) = -ones(opts.constraints.num_terminalConstraints,1);

Aineq = [Aineq;matrix];

% add 0s to bineq to uphold dimension between Aineq and bineq

a = size(Aineq);
b = size(bineq);

z = a(1) - b(1);
for i=1:z
    bineq = [bineq;0];
end

end
