function set_evaluationPoints(obj, newEvaluationPoints)
    obj.m_evaluationPoints = newEvaluationPoints;
    obj.m_numOfEvaluationPoints = size(obj.m_evaluationPoints, 1);
end