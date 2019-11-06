classdef MovingLeastSquares_LCS < WeightedLeastSquares_LCS

    methods
        function obj = MovingLeastSquares_LCS(varargin)
            obj = obj@WeightedLeastSquares_LCS(varargin{:});
            obj.m_stationaryPoints = obj.m_evaluationPoints;
        end
    end 

    methods (Access = protected)
        evaluationMethod(obj);
        derivativesEvaluationMethod(obj);
    end
end