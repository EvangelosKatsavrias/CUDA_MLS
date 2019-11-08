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

function [coefficients, conditionNumber] = findCoefficients_WeightedLeastSquares(samplePoints, sampleValues, polynomial, solver, condNumEvaluator, weightPerStationaryPoint)

numOfStationaryPoints   = size(weightPerStationaryPoint, 1);
numOfSamplePoints       = size(samplePoints, 1);
numOfMonomials          = length(polynomial(samplePoints(1,:)));
coefficients            = zeros(numOfStationaryPoints, numOfMonomials);
conditionNumber         = zeros(numOfStationaryPoints, 1);

for stationaryPointIndex = 1:numOfStationaryPoints

    Sum_bbT = zeros(numOfMonomials, numOfMonomials);
    Sum_bf  = zeros(numOfMonomials, 1);
    for samplePointIndex = 1:numOfSamplePoints
        b       = polynomial(samplePoints(samplePointIndex,:));
        theta_b = weightPerStationaryPoint(stationaryPointIndex, samplePointIndex)*b;
        Sum_bbT = Sum_bbT   + theta_b*b';
        Sum_bf  = Sum_bf    + theta_b*sampleValues(samplePointIndex);
    end

    coefficients(stationaryPointIndex, :) = solver(Sum_bbT,Sum_bf);
    conditionNumber(stationaryPointIndex) = condNumEvaluator(Sum_bbT);

end

end