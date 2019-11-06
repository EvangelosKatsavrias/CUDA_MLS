function set_sampleValues(obj, newSampleValues)
    if (obj.m_numOfSamplePoints == size(newSampleValues, 1))
        obj.m_sampleValues = newSampleValues;
    else error('The number of sample values does not correspond to the number of sample points.');
    end
end