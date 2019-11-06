function derivativesEvaluationMethod(obj)

obj.m_derivatives = findDerivatives_MovingLeastSquaresLCS(obj.m_samplePoints, obj.m_stationaryPoints, obj.m_sampleValues, @(x)obj.m_polynomial.evaluate(x), obj.m_weightPerStationaryPoint);

end