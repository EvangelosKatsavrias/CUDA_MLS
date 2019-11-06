function values = evaluatePoints_MovingLeastSquaresLCS(coefficients_local)

numOfEvaluationPoints = size(coefficients_local, 1);
values = zeros(numOfEvaluationPoints, 1);

for evaluationPointIndex = 1:numOfEvaluationPoints
    values(evaluationPointIndex) = coefficients_local(evaluationPointIndex, 1);
end

end