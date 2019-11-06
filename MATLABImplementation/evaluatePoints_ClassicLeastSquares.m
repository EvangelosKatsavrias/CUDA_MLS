function evaluations = evaluatePoints_ClassicLeastSquares(evaluationPoints, polynomial, coefficients)

    numOfEvaluationPoints = size(evaluationPoints, 1);
    evaluations = zeros(numOfEvaluationPoints, 1);

    for evaluationPointIndex = 1:numOfEvaluationPoints
        evaluations(evaluationPointIndex) = polynomial(evaluationPoints(evaluationPointIndex, :))'*coefficients;
    end


end