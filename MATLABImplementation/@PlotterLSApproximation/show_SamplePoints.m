function show_SamplePoints(obj)

    figure(obj.handle_figure); hold all;
    if (~isempty(obj.handle_samplePoints)); delete(obj.handle_samplePoints); end
    obj.handle_samplePoints = plot(obj.m_LeastSquaresApproximations(1).m_samplePoints, obj.m_LeastSquaresApproximations(1).m_sampleValues, 'o');

end