classdef WeightedLeastSquares_LCS < WeightedLeastSquares

	methods
        function obj = WeightedLeastSquares_LCS(varargin)
            obj = obj@WeightedLeastSquares(varargin{:});
        end
    end
    
    methods (Access = protected)
        find_coefficients(obj);
        evaluationMethod(obj);
%         derivativesEvaluationMethod(obj);
    end
    
end