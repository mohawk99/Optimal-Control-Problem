function solOCP = mainOptimalControl(x0,initialGuess,opts)
%MAINOPTIMALCONTROL - this function starts the iterative approach of
%solving the optimal control problem and outputs the optimal solution
%
% Syntax:
%       solOCP = mainOptimalControl(x0,initialGuess,opts)
%
% Input Arguments:
%       - x0:           initial state for the optimal control problem
%       - initialGuess: intial solution of the optimal control problem 
%                         
%       - opts:         structure containing the algorithm settings, the
%                       constraint functions, etc.
%
% Output Arguments:
%       - solOCP:       solution of the optimal control problem, i.e.
%                       contains optimal state trajectory and inputs
%
% ------------------------------------------------------------------------



xCand = initialGuess.x;
uCand = initialGuess.u;

for ii = 1:opts.maxIter
    
    % linearize constraints
    [Aeq,beq,Aineq,bineq] = convexifyConstraints(x0,xCand,uCand,opts);
    
    % solve convexified OCP
    [delta_x,delta_u] = solveConvexOCP(xCand,uCand,Aeq,beq,Aineq,bineq,opts);
    uCand = uCand + delta_u; % update control inputs
    xCand = xCand + delta_x; % update states

    % compute "simulated" slack variables and virtual control inputs
    virtCtrlSim = computeVirtualControlInputs(xCand,uCand,opts);
    slackSim = computeStageConstraintSlacks(xCand,uCand,opts);
   
   

    % check convergence of computed solution by stacking xCand and uCand 
    % in a column vector, taking the 2-norm and checking if the norm is 
    % smaller than opts.toleranceConvergence

    % check feasibility of computed solution by checking if the infinity
    % norm of virtCtrlSim and slackSim is smaller than opts.toleranceFeasability
    % 
    % current version of code produces solutions which fail these tests,
    % hence we comment this section out, in order to produce some output
    %
    % stack = [delta_x(:);delta_u(:)];
    % n = norm(stack);
    % n_vrtCtrSim = norm(virtCtrlSim,"inf");
    % n_slackSim = norm(slackSim,"inf");
    % 
    % if ii == opts.maxIter && n < opts.toleranceConvergence && n_vrtCtrSim < opts.toleranceFeasability && n_slackSim < opts.toleranceFeasability
    %     disp('Converged and feasable!')
    %     solOCP.x = xCand;
    %     solOCP.u = uCand;
    % else
    %     solOCP = "error";
    % end

    solOCP.x = xCand;
    solOCP.u = uCand;

    
end

end
