function weightPerStationaryPoint = findWeights_WeightedLeastSquares(samplePoints, stationaryPoints, weightingFunction)

numOfSamplePoints = size(samplePoints, 1);
numOfStationaryPoints = size(stationaryPoints, 1);
weightPerStationaryPoint = zeros(numOfStationaryPoints, numOfSamplePoints);

for stationaryPointIndex = 1:numOfStationaryPoints
    
    for samplePointIndex = 1:numOfSamplePoints
        weightPerStationaryPoint(stationaryPointIndex, samplePointIndex) = weightingFunction(norm(samplePoints(samplePointIndex, :) -stationaryPoints(stationaryPointIndex, :), 2));
    end
    
end

end