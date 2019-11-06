function evaluationMethod(obj)

obj.m_evaluations = evaluatePoints_ClassicLeastSquares(obj.m_evaluationPoints, @(x)obj.m_polynomial.evaluate(x), obj.m_coefficients);

end