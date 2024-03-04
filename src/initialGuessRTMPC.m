function initialGuess = initialGuessRTMPC(x0,xf,opts)
% INITIALGUESRTMPC - provides the inital solution guess between points x0
% and xf based on the straight line interpolation for the optimization problem
%
% Syntax:
%       initialGuess = initialGuessRTMPC(x0,xf,opts)
%
% Input Arguments:
%       - x0:               initial state for the optimal control problem, i.e,
%                           the (measured/estimated) system state at time k
%       - xf:               final state
%       - opts:             a structure containing the algorithm settings 
%                           (see function getOptimizationVariables for more info)
%
% Output Arguments:
%       - initialGuess:     straight line interpolation between x0 and xf


    N = opts.N;

    % If there is no final state given, use the reference trajectory
    if (isempty(xf) || all(xf==0))
        x_final = opts.xref(:,opts.nx);   
    else
        x_final = xf;   
    end

    dims = size(x0);

    x_inter = zeros(dims(1), N);
    
    % Calculating the straight line trajectory
    for k=1:N
        const = (N - k)/(N - 1);
        x_inter(:, k) = const*x0 + (1-const)*x_final;
    end

    initialGuess.x = x_inter;
    initialGuess.u = zeros(opts.nu,N);

end