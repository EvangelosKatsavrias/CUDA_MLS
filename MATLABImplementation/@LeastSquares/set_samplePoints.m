function set_samplePoints(obj, newSamplePoints)
    obj.m_samplePoints = newSamplePoints;
    obj.m_numOfSamplePoints = size(obj.m_samplePoints, 1);
end