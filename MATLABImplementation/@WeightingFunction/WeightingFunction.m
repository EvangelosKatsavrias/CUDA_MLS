classdef WeightingFunction < handle

    properties (Access = private)
        m_Span = 10
        m_eval_Funct
        m_Exponent = 4
    end

    properties (Hidden, Access = private)
        functionHandle
    end

    enumeration
        Wendland(@(h, d, exponent)WendlandDistribution(h, d, exponent)),
        Gauss(@(h, d, exponent)gaussianDistribution(h, d, exponent)),
        Radial(@(h, d, exponent)radialDistribution(h, d, exponent)),
        Constant(@(h, d, exponent)ConstantDistribution(h, d, exponent))
    end

    methods
        function obj = WeightingFunction(varargin)
            if nargin > 1; obj.m_Span = varargin{2}; end
            obj.functionHandle = varargin{1};
            obj.m_eval_Funct = @(d)varargin{1}(obj.m_Span, d);
        end

        function obj = set_Span(obj, newSpanValue)
            obj.m_Span = newSpanValue;
            obj.m_eval_Funct = @(d)obj.functionHandle(newSpanValue, d, obj.m_Exponent);
        end
        
        function obj = set_Exponent(obj, newValue)
            obj.m_Exponent = newValue;
            obj.m_eval_Funct = @(d)obj.functionHandle(obj.m_Span, d, obj.m_Exponent);
        end

        function spanValue = get_Span(obj); spanValue = obj.m_Span; end

        function evaluation = evaluate(obj, point); evaluation = obj.m_eval_Funct(point); end

    end

end