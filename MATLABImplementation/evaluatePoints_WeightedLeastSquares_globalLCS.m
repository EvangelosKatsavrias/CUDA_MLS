function values = evaluatePoints_WeightedLeastSquares_globalLCS(evaluationPoints, stationaryPoints, polynomial, weightingFunction, coefficients_local)

numOfEvaluationPoints = size(evaluationPoints, 1);
numOfStationaryPoints = size(stationaryPoints, 1);
values = zeros(numOfEvaluationPoints, 1);

for evaluationPointIndex = 1:numOfEvaluationPoints

    weights         = zeros(numOfStationaryPoints, 1);
    weights_sum     = 0;
    coefficients    = 0;
    
    for stationaryPointIndex = 1:numOfStationaryPoints
        weights(stationaryPointIndex) = weightingFunction(norm(evaluationPoints(evaluationPointIndex,:)-stationaryPoints(stationaryPointIndex,:), 2));
        weights_sum = weights_sum + weights(stationaryPointIndex);
    end
    
    phi = weights/weights_sum;
    
    for stationaryPointIndex = 1:numOfStationaryPoints
        coefficients = coefficients + phi(stationaryPointIndex)*coefficients_local(stationaryPointIndex, :)*polynomial(evaluationPoints(evaluationPointIndex, :)-stationaryPoints(stationaryPointIndex, :));
    end
    
    values(evaluationPointIndex) = coefficients;
    
end

end