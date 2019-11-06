function evaluationMethod(obj)

obj.m_evaluations = evaluatePoints_WeightedLeastSquares_globalLCS(obj.m_evaluationPoints, obj.m_stationaryPoints, @(x)obj.m_polynomial.evaluate(x), @(d)obj.m_weightingFunction.evaluate(d), obj.m_coefficients);

end