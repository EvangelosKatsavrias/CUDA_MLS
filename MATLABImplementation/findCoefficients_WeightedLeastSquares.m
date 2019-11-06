function [coefficients, conditionNumber] = findCoefficients_WeightedLeastSquares(samplePoints, sampleValues, polynomial, solver, condNumEvaluator, weightPerStationaryPoint)

numOfStationaryPoints   = size(weightPerStationaryPoint, 1);
numOfSamplePoints       = size(samplePoints, 1);
numOfMonomials          = length(polynomial(samplePoints(1,:)));
coefficients            = zeros(numOfStationaryPoints, numOfMonomials);
conditionNumber         = zeros(numOfStationaryPoints, 1);

for stationaryPointIndex = 1:numOfStationaryPoints

    Sum_bbT = zeros(numOfMonomials, numOfMonomials);
    Sum_bf  = zeros(numOfMonomials, 1);
    for samplePointIndex = 1:numOfSamplePoints
        b       = polynomial(samplePoints(samplePointIndex,:));
        theta_b = weightPerStationaryPoint(stationaryPointIndex, samplePointIndex)*b;
        Sum_bbT = Sum_bbT   + theta_b*b';
        Sum_bf  = Sum_bf    + theta_b*sampleValues(samplePointIndex);
    end

    coefficients(stationaryPointIndex, :) = solver(Sum_bbT,Sum_bf);
    conditionNumber(stationaryPointIndex) = condNumEvaluator(Sum_bbT);

end

end