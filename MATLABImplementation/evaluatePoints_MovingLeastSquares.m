function values = evaluatePoints_MovingLeastSquares(evaluationPoints, polynomial, coefficients_local)

numOfEvaluationPoints = size(evaluationPoints, 1);
values = zeros(numOfEvaluationPoints, 1);

for evaluationPointIndex = 1:numOfEvaluationPoints
    values(evaluationPointIndex) = coefficients_local(evaluationPointIndex, :)*polynomial(evaluationPoints(evaluationPointIndex, :));
end

end