function show_Approximation(obj, approximationIndices)

figure(obj.handle_figure); hold all;

for approximationIndex = approximationIndices
    lineWidth = 1;
%    if approximationIndex == 2; lineWidth = 3; end;
    if (length(obj.handles_approximations)>approximationIndex); if (obj.handles_approximations(approximationIndex))>0; delete(obj.handles_approximations(approximationIndex)); end; end
    obj.handles_approximations(approximationIndex) = plot(obj.m_LeastSquaresApproximations(approximationIndex).m_evaluationPoints, obj.m_LeastSquaresApproximations(approximationIndex).m_evaluations, 'LineWidth', lineWidth);
end

end