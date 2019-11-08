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