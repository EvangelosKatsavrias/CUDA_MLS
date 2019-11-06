function copy_properties(obj, sourceObject)
    obj.set_samplePoints(sourceObject.m_samplePoints);
    obj.set_sampleValues(sourceObject.m_sampleValues);
    obj.set_evaluationPoints(sourceObject.m_evaluationPoints);
    obj.set_polynomialBasis(sourceObject.m_polynomial);
    obj.m_solver = sourceObject.m_solver;
    obj.m_condNumEvaluator = sourceObject.m_condNumEvaluator;
end