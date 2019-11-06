classdef LeastSquares < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    % min|| c_i*N_i(x_k) - f(x_k) ||,  
    % f  : sample values of the uknown distribution, at specific sample points x_k
    % c_i: unknown coefficients (control parameters of the uknown distribution)
    % N_i: basis functions (components of the unknown distribution), evaluated at specific sample points x_k
    
    properties
        m_samplePoints = [0 1 3 6 7 9 15 20 21 23 24 28 30]'
        m_sampleValues = [1 2 0 3 4 8 5 6 9 10 5 2 1]'
        m_evaluationPoints = linspace(0, 30, 101)'
        m_polynomial = Polynomial
        m_numOfSamplePoints = 13
        m_numOfEvaluationPoints = 101
        m_solver = @(A, b)linsolve(A, b)
        m_condNumEvaluator = @(A)cond(A)
        m_coefficients
        m_evaluations
        m_conditionNumber
        m_derivatives
    end

    methods
        function obj = LeastSquares(varargin); obj = constructorProcesses(obj, varargin{:}); end
        set_samplePoints(obj, newSamplePoints);
        set_sampleValues(obj, newSampleValues);
        set_evaluationPoints(obj, newEvaluationPoints);
        set_polynomialBasis(obj, newPolynomialBasis);
        copy_properties(obj, sourceObject);
        evaluatePoints(obj);
    end
    
    methods (Abstract, Access = protected)
        checkEvaluationInputData(obj);
        find_coefficients(obj);
        evaluationMethod(obj);
        derivativesEvaluationMethod(obj);
    end
    
end


function obj = constructorProcesses(obj, varargin)
    for arginIndex = 1:nargin-1
        if isa(varargin{arginIndex}, 'LeastSquares'); obj.copy_properties(varargin{arginIndex}); end
    end
    for arginIndex = 1:nargin-1
        if strcmp(varargin{arginIndex}, 'samplePoints'); obj.set_samplePoints(varargin{arginIndex+1}); end
        if strcmp(varargin{arginIndex}, 'sampleValues'); obj.set_sampleValues(varargin{arginIndex+1}); end
        if strcmp(varargin{arginIndex}, 'evaluationPoints'); obj.set_evaluationPoints(varargin{arginIndex+1}); end
        if isa(varargin{arginIndex}, 'Polynomial'); obj.set_polynomialBasis(varargin{arginIndex}); end
    end
end