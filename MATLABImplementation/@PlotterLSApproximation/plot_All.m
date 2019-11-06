function plot_All(obj)

if (~isempty(obj.handles_approximations)); delete(obj.handles_approximations(logical(obj.handles_approximations))); end
obj.handles_approximations = zeros(length(obj.m_LeastSquaresApproximations), 1);
% obj.show_SamplePoints;
obj.show_Approximation(1:length(obj.m_LeastSquaresApproximations));

end