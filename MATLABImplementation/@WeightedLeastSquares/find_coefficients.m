function find_coefficients(obj)
    
[obj.m_coefficients, obj.m_conditionNumber] = findCoefficients_WeightedLeastSquares(obj.m_samplePoints, obj.m_sampleValues, @(x)obj.m_polynomial.evaluate(x), obj.m_solver, obj.m_condNumEvaluator, obj.m_weightPerStationaryPoint);

end