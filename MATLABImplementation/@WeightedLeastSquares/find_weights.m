function find_weights(obj)
    
    obj.m_weightPerStationaryPoint = findWeights_WeightedLeastSquares(obj.m_samplePoints, obj.m_stationaryPoints, @(d)obj.m_weightingFunction.evaluate(d));
    
end