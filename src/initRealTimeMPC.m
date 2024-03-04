function opts = initRealTimeMPC(benchmark,opts,params)
%INITREALTIMEMPC - initializes the optimization vector, values needed for 
% the cost function in the optimal control problem, as well as the system 
% dynamics, stage and terminal constraint functions and boundary conditions
%
% Syntax:
%       opts = initRealTimeMPC(benchmark,opts,params)
%
% Input Arguments:
%       - benchmark:        name of the considered benchmark model (see
%                           "aroc/benchmarks/dynamics/...")
%       - opts:             a structure containing the algorithm settings 
%                           (see function getOptimizationVariables for more info)
%       - params:           a structure containing the benchmark parameters
%
% Output Arguments:
%       - opts:             a structure containing the algorithm settings 
%                           (see function getOptimizationVariables for more info)


% setup values needed for the cost function in the optimal control problem
opts.mu = 100;
opts.lambda = 1;

opts.toleranceConvergence = 1e-6;
opts.toleranceFeasability = 10e-6;


% system dynamcis --------------------------------------------------------

str = ['funHandle = @(x,u,w)',benchmark,'(x,u,w);'];
eval(str);
% get number of states, inputs, and disturbances
[count,out] = inputArgsLength(funHandle,3);
opts.nx = out; % dimension of state space
opts.nu = count(2); % dimension of input space
opts.nw = count(3); % dimension of disturbance space

% nominal dynamics
w_ref = zeros(opts.nw,1);
system = @(x,u) unicycle(x,u,w_ref);

% optimal control --------------------------------------------------------
import casadi.*

% convert to symbolic and back to a function handle (without this step, 
% some operations would not be defined for casadi)
xx = sym('x',[opts.nx,1],'real');
uu = sym('u',[opts.nu,1],'real');
ff = system(xx,uu);
dyn = matlabFunction(ff,'Vars',{xx,uu});

% integrate dynamics \dot(x) = f(x,p,....) using casadi
x_ = MX.sym('x_',opts.nx,1);
u_ = MX.sym('u_',opts.nu,1);
odeOpts.t0 = 0;
odeOpts.tf = opts.dt;
ode = struct('x',x_,'p',u_,'ode',dyn(x_,u_));
Dx = integrator('Dx','rk',ode,odeOpts);

% setup dynamics constraint: x_{k+1} - f(x_{k},u_{k},0) = 0
x_k = MX.sym('x_k',opts.nx,1); % statet at time step k
u_k = MX.sym('u_k',opts.nu,1); % input at time step k
x_int = getfield(Dx('x0',x_k,'p',u_k),'xf'); % integrate system until time k+1
x_next = MX.sym('x_next',opts.nx,1); % state at time step k+1
% -> constraint function
opts.constraints.dynamics = Function('con_dynamics',{x_k,u_k,x_next}, ...
    {x_next-x_int});
opts.constraints.jac_dynamics = Function('jac_dynamics',{x_k,u_k}, ...
    {jacobian(x_int,[x_k;u_k])});

% setup stage and terminal constraints functions
opts.constraints.num_stageConstraints = length(opts.stageConstraints);
for ii = 1:opts.constraints.num_stageConstraints
    stageCon = opts.stageConstraints{ii}(xx,uu);
    stageCon = matlabFunction(stageCon,'Vars',{xx,uu});
    opts.constraints.stageCon{ii} = Function('stageCon',{x_k,u_k}, ...
        {stageCon(x_k,u_k)});
    opts.constraints.jac_stageCon{ii} = Function('jac_stageCon',{x_k,u_k}, ...
        {jacobian(stageCon(x_k,u_k),[x_k;u_k])});
end

opts.constraints.num_terminalConstraints = length(opts.terminalConstraints);
for ii = 1:opts.constraints.num_terminalConstraints
    terminalCon = opts.terminalConstraints{ii}(xx);
    terminalCon = matlabFunction(terminalCon,'Vars',{xx});
    opts.constraints.terminalCon{ii} = Function('terminalCon',{x_k}, ...
        {terminalCon(x_k)});
    opts.constraints.jac_terminalCon{ii} = Function('jac_terminalCon',{x_k}, ...
        {jacobian(terminalCon(x_k),x_k)});
end

% setup vector of optimization variables, i.e. compute indices of 
% optimization variables
opts.idxVars = getOptimizationVariables(opts);


% setup input and state boundary constraints and
% extract lower bounds on optimization variables
opts.ub = inf(opts.idxVars.nVars,1);
opts.lb = zeros(opts.idxVars.nVars,1);

opts.constraints.lbu = params.U.inf;
opts.constraints.ubu = params.U.sup;

for k=1:opts.idxVars.n_u
    if mod(k,2) == 1
        opts.ub(opts.idxVars.u) = 4;
        opts.lb(opts.idxVars.u) = -4;
    else
        opts.ub(opts.idxVars.u) = 2.5;
        opts.lb(opts.idxVars.u) = -2.5;
    end
end
if isfield(params,'X')
    opts.constraints.lbx = params.X.inf;
    opts.constraints.ubx = params.X.sup;

    opts.ub(opts.idxVars.x) = inf;
    opts.lb(opts.idxVars.x) = -inf;
else
    opts.constraints.lbx = -inf(opts.nx,1);
    opts.constraints.ubx = inf(opts.nx,1);

    opts.ub(opts.idxVars.x) = inf;
    opts.lb(opts.idxVars.x) = -inf;
end
end