function slackSim = computeStageConstraintSlacks(xCand,uCand,opts)
%COMPUTESTAGECONSTRAINTSLACKS - this function computes the stage constraint
%slacks
%
% Syntax:
%       slackSim = computeStageConstraintSlacks(xCand,uCand,opts)
%
% Input Arguments:
%       - xCand:            current candidate for the state vector
%       - uCand:            current candidate for the input vector
%                         
%       - opts:             structure containing the algorithm settings
%
% Output Arguments:
%       - slackSim:         computed stage constraint slacks
%
% ------------------------------------------------------------------------

slackSim = zeros(opts.constraints.num_stageConstraints);

% compute stage constraint slacks
for j=2:opts.N-1
    for ii = 1:opts.constraints.num_stageConstraints
        slackSim(ii,j-1) = max(0,full(opts.constraints.stageCon{ii}(xCand(:,j-1),uCand(:,j))));
    end
    
end


end

