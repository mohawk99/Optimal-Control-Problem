function sol = realTimeMPC(benchmark,opts,params)
%REALTIMEMPC - this function initializes the problem, calls the straight
%line interpolation function and starts the optimal control problem at hand
%
% Syntax:
%       sol = realTimeMPC(benchmark,opts,params)
%       
%
% Input Arguments
%  
%       -benchmark:    name of the considered benchmark model (see
%                      "aroc/benchmarks/dynamics/...")
%       -params:             a structure containing the benchmark parameters
%
%           -.x0:           initial state
%           -.xf:           goal state (optional)
%           -.U:            set of admissible control inputs (class:
%                           interval)
%           -.W:            set of disturbances (class: zonotope)
%
%       -opts:              a structure containing the algorithm settings
%
%           -.tSim:         final time for the closed-loop simulation
%                           (>= dt*N)
%           -.dt:           sampling time [{0.1} / positive real number]
%           -.N:            prediction horizon
%                           [{10} / positive integer]
%           -.Q:            state weighting matrix for the cost function of
%                           optimal control problem (reference trajectory)
%           -.R:            input weighting matrix for the cost function of
%                           optimal control problem (reference trajectory)
%           -.xref:         reference trajectory (dim: [nx,ceil(tSim/dt)+1])
%           -.uref:         reference trajectory (dim: [nx,ceil(tSim/dt)])    
%           -.maxIter:      maximum number of iterations for the optimal
%                           control problem (PTR iterations)
%                           [{10} / positive integer]
%           -.stageConstraints:     constraints to be enforced along the
%                                   prediction horizon (cell array of
%                                   function handles)
%           -.terminalConstraints   constraints to be enforced at the last 
%                                   step of the prediction horizon
%                                   (cell array of function handles)
%
% Output Arguments:
%
%       -sol:       results object storing the computed state trajectory
%
% ------------------------------------------------------------------------


% initialization
opts = initRealTimeMPC(benchmark,opts,params);

% linear interpolation between x0 and xf
initialGuess = initialGuessRTMPC(params.x0,params.xf,opts);

% start the optimal control problem
sol = mainOptimalControl(params.x0,initialGuess,opts);


% plots linear interpolation
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
figure; hold on
plot(opts.obstacles{1},[1,2],'color',[1,0,0])
plot(opts.obstacles{2},[1,2],'color',[1,0,0])
plot(initialGuess.x(1,:),initialGuess.x(2,:),'kx-')
xlabel('$x_1$ in $m$')
ylabel('$x_2$ in $m$')

% plot iterates of solution
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
figure; hold on
plot(opts.obstacles{1},[1,2],'color',[1,0,0])
plot(opts.obstacles{2},[1,2],'color',[1,0,0])
plot(sol.x(1,:),sol.x(2,:),'kx-')
xlabel('$x_1$ in $m$')
ylabel('$x_2$ in $m$')

end
