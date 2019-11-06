classdef ClassicLeastSquares < LeastSquares

    methods
        function obj = ClassicLeastSquares(varargin)
            obj = obj@LeastSquares(varargin{:});
        end
    end
    
	methods (Access = protected)
        checkEvaluationInputData(obj);
        find_coefficients(obj);
        evaluationMethod(obj);
        derivativesEvaluationMethod(obj);
    end
    
end