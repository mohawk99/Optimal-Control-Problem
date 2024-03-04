function virtCtrlSim = computeVirtualControlInputs(xCand,uCand,opts)
%COMPUTEVIRTUALCONTROLINPUTS - this function computes the virutal control
%inputs
%
% Syntax:
%       virtCtrlSim = computeVirtualControlInputs(xCand,uCand,opts)
%
% Input Arguments:
%       - xCand:            current candidate for the state vector
%       - uCand:            current candidate for the input vector
%                         
%       - opts:             structure containing the algorithm settings
%
% Output Arguments:
%       - virtCtrlSim:      computed virtual control inputs
%
% ------------------------------------------------------------------------

virtCtrlSim = zeros(opts.nx,opts.N);

% compute virtual control inputs
for j=2:opts.N
    
    virtCtrlSim(:,j-1) = abs(full(opts.constraints.dynamics(xCand(:,j-1),uCand(:,j),xCand(:,j))));

end

end

