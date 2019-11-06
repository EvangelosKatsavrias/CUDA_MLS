classdef Polynomial
    
    properties
        m_dimensions = 1
        m_degree = 1
    end
    
    methods
        
        function obj = Polynomial(varargin)
            if nargin > 0; obj.m_dimensions = varargin{1}; end
            if nargin > 1; obj.m_degree = varargin{2}; end
        end
        
        function values = evaluate(obj, x)
            values = obj.univariateMonomials(obj.m_degree, x);
        end
        
        function values = evaluate2(obj, x)
            for index = 2:obj.m_dimensions
    %            values = evaluate2(obj, x);
                values = univariateMonomials(obj.m_degree, x(1));
            end
        end
        
    end
    
    methods (Static)
        function values = univariateMonomials(degree, x)
            values = zeros(degree+1, 1);
            values(1) = 1;
            for index = 2:degree+1
                values(index) = values(index-1)*x;
            end
        end
    end
    
end