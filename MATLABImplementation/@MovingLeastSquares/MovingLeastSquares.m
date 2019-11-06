classdef MovingLeastSquares < WeightedLeastSquares
    
    methods
        function obj = MovingLeastSquares(varargin)
            obj = obj@WeightedLeastSquares(varargin{:});
            obj.m_stationaryPoints = obj.m_evaluationPoints;
        end
    end

	methods (Access = protected)
        evaluationMethod(obj);
%         derivativesEvaluationMethod(obj);
    end
end