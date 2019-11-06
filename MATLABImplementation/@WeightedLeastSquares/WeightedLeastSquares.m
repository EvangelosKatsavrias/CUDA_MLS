classdef WeightedLeastSquares < LeastSquares

    properties
        m_weightingFunction = WeightingFunction.Wendland
        m_stationaryPoints = linspace(0, 30, 101)'
        m_numOfStationaryPoints = 101
        m_weightPerStationaryPoint
    end
    
    methods
        function obj = WeightedLeastSquares(varargin)
            obj = obj@LeastSquares(varargin{:});
            for arginIndex = 1:nargin
                if isa(varargin{arginIndex}, 'WeightedLeastSquares'); obj.copy_weightProperties(varargin{arginIndex}); end
            end
            for arginIndex = 1:nargin
                if isa(varargin{arginIndex}, 'WeightingFunction'); obj.set_weightingFunction(varargin{arginIndex}); end
                if strcmp(varargin{arginIndex}, 'stationaryPoints'); obj.set_stationaryPoints(varargin{arginIndex+1}); end
            end
        end
        
        function set_weightingFunction(obj, newWeightingFunction)
            obj.m_weightingFunction = newWeightingFunction;
        end
        
        function set_stationaryPoints(obj, newStationaryPoints)
            obj.m_stationaryPoints = newStationaryPoints;
            obj.m_numOfEvaluationPoints = size(obj.m_stationaryPoints, 1);
        end
        
        function copy_weightProperties(obj, sourceObject)
            obj.m_weightingFunction = sourceObject.m_weightingFunction;
            obj.set_stationaryPoints(sourceObject.m_stationaryPoints);
        end
    end
    
    methods (Access = protected)
        checkEvaluationInputData(obj);
        find_coefficients(obj);
        find_weights(obj);
        evaluationMethod(obj);
        derivativesEvaluationMethod(obj);
    end

end