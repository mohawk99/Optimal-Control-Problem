function idxVars = getOptimizationVariables(opts)
%GETOPTIMIZATIONVARIABLES computes the indices of the states,
% control inputs, virtual controls and slack variables in the vector of 
% optimization variabels
%
% Input Arguments:
%       - opts: structure with the fields
%           .nx:           dimension of the state space
%           .nu:           dimension of the input space
%           .N:            number of time steps
%           .constraints: structure with the fields
%               .num_stageConstraints:      number of constraints enforced
%                                           at 0 <= k < N
%               .num_terminalConstraints:   number of terminal constraints
%
%
% Output Arguments:
%       - idxVars: structure containing the following fields
%           .x:                             indices of the states (the 
%                                           initial state is not modeled as
%                                           an optimization variable)
%                                           (matrix of dimension 
%                                           [opts.nx,opts.N])
%           .n_x:                           number of entries of .x
%                                           variables modeling the states
%           .u:                             indices of the control inputs
%                                           (matrix of dimension 
%                                           [opts.nu,opts.N])
%           .n_u:                           number of entries of .u
%           .vrtCtrl_pos:                   indices of the positive virtual
%                                           control inputs 
%                                           (matrix of dimension
%                                           [opts.nx,opts.N])
%           .vrtCtrl_neg:                   indices of the negative virtual
%                                           control inputs 
%                                           (matrix of dimension
%                                           [opts.nx,opts.N])        
%           .slack_stageCon:                indices of the slack variables
%                                           of constraints enforced at 
%                                           0 <= k < N
%                                           (matrix of dimension
%                                           [opts.constraints.num_stageConstraints,opts.N])
%           .slack_terminalCon:             indices of the slack variables
%                                           of the terminal constraints 
%                                           (vector of length 
%                                           opts.constraints.num_terminalConstraints)
%           .penalty_vrtCtrl_pos:           indices of the corresponding
%                                           penalties
%                                           (vector of length opts.N)
%           .penalty_vrtCtrl_neg:           indices of the corresponding
%                                           penalties
%                                           (vector of length opts.N)
%           .penalty_slack_stageCon:        indices of the corresponding
%                                           penalties
%                                           (vector of length opts.N)
%           .penalty_slack_terminalCon:     index of the corresponding
%                                           penalties
%                                           (vector of length 1)
%           .trustRegionRad:                index of the soft trust region
%                                           penalty (vector of length 1)
%           .nVars:                         total number of optimization 
%                                           variables
%
% ------------------------------------------------------------------------


% vector of (nominal) variables: [u_{0},x_{1},....x_{N-1},u_{N-1},x_{N}]

% 0 <= k <= N-1
tmp = reshape(1:opts.N*(opts.nx+opts.nu),[],opts.N);
idxVars.u = tmp(1:opts.nu,:);
idxVars.x = tmp(opts.nu+(1:opts.nx),:);
idxVars.n_u = length(idxVars.u(:));
idxVars.n_x = length(idxVars.x(:));

countVars = opts.N*(opts.nx+opts.nu);

% virtual control inputs
tmp = reshape(countVars + (1:opts.nx*opts.N),[],opts.N);
idxVars.vrtCtrl_pos = tmp;
tmp = tmp + opts.nx*opts.N;
idxVars.vrtCtrl_neg = tmp;

countVars = countVars + 2*opts.nx*opts.N;

% slack variables
% -> stage constraints
idxVars.slack_stageCon = ...
    reshape(countVars + (1:opts.constraints.num_stageConstraints*opts.N),[],opts.N);
countVars = countVars + opts.constraints.num_stageConstraints*opts.N;
% -> terminal constraints
idxVars.slack_terminalCon = ...
    countVars + (1:opts.constraints.num_terminalConstraints);
countVars = countVars + opts.constraints.num_terminalConstraints;

% penalties
idxVars.penalty_vrtCtrl_pos = countVars + (1:opts.N);
countVars = countVars + opts.N;
idxVars.penalty_vrtCtrl_neg = countVars + (1:opts.N);
countVars = countVars + opts.N;
idxVars.penalty_slack_stageCon = countVars + (1:opts.N);
countVars = countVars + opts.N;
idxVars.penalty_slack_terminalCon = countVars + 1;
countVars = countVars + 1;

% soft penalties
idxVars.trustRegionRad = countVars + 1;
countVars = countVars + 1;

idxVars.nVars = countVars;

end
