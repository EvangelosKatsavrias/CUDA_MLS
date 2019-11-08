%%  CUDA_MLS Framework
%
%   Copyright 2017-2018 Evangelos D. Katsavrias, Luxembourg
%
%   This file is part of the CUDA_MLS Framework.
%
%   CUDA_MLS Framework is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License version 3 as published by
%   the Free Software Foundation.
%
%   CUDA_MLS Framework is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with CUDA_MLS Framework.  If not, see <https://www.gnu.org/licenses/>.
%
%   Contact Info:
%   Evangelos D. Katsavrias
%   email/skype: vageng@gmail.com
% -----------------------------------------------------------------------

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