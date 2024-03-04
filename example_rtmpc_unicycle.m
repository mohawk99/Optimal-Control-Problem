%% optimal control of a unicycle
% example based on: https://arxiv.org/pdf/2209.03535v1.pdf

%% Run this script in order to execute the code

% Currently the code does not provide the expected solution of the optimal
% control problem (second plot) compared to desired plot as seen in the report.


clear all

% define benchmark, parameters, and settings
benchmark = 'unicycle';

% parameters
params.X = interval(-inf(3,1),inf(3,1));
width = [4;2.5];
params.U = interval(-width,width);
width = [1; 1];
params.W = interval(-width,width);
params.x0 = [0;0;0];
params.xf = [5;5;0];

% algorihtm settings
opts.N = 30;
opts.dt = 0.1;
opts.tSim = opts.N*opts.dt;  % only "one-shot-planning"
% cost matrices (see the cost function in Sec. IV) 
opts.Q = zeros(3);
opts.R = eye(2);
% reference trajectories (see the cost function in Sec. IV) 
opts.xref = zeros(3,opts.N+1);  % since Q=0, any arbitrary choice is admissible
opts.uref = zeros(2,opts.N); 
% max. number of iterations
opts.maxIter = 20;

% obstacle avoidance (stage constraints)
E1 = diag([0.75^2,1.5^2]);
c1 = [1;2];
opts.obstacles{1} = ellipsoid(E1,c1); % just for plotting 
% -> use plot(opts.obstacles{1},[xx,yy]) where xx denotes the dimension of
% the state space to be shown on the abscissa and yy denotes the dimension
% of the state space to be shown on the ordinate
opts.stageConstraints{1} = @(x,u) 1 - ([x(1);x(2)]-c1)'*inv(E1)*([x(1);x(2)]-c1); % constraint function
E2 = diag([0.75^2,1.5^2]);
c2 = [4;3];
opts.obstacles{2} = ellipsoid(E2,c2);
opts.stageConstraints{2} = @(x,u) 1- ([x(1);x(2)]-c2)'*inv(E2)*([x(1);x(2)]-c2);
% terminal constraints - box enclosure of the bounday condition
% ellipsoid(xf,Qf) defined in Sec. IV
opts.terminalConstraints{1} = @(x) 4.5-x(1);
opts.terminalConstraints{2} = @(x) x(1) - 5.5;
opts.terminalConstraints{3} = @(x) 4.5 - x(2);
opts.terminalConstraints{4} = @(x) x(2) - 5.5;
opts.terminalConstraints{5} = @(x) -(20*pi)/180 - x(3);
opts.terminalConstraints{6} = @(x) x(3) - (20*pi)/180;

% real-time MPC
res = realTimeMPC(benchmark,opts,params);
