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

function values = evaluatePoints_WeightedLeastSquares_global(evaluationPoints, stationaryPoints, polynomial, weightingFunction, coefficients_local)

numOfEvaluationPoints = size(evaluationPoints, 1);
numOfStationaryPoints = size(stationaryPoints, 1);
values = zeros(numOfEvaluationPoints, 1);

for evaluationPointIndex = 1:numOfEvaluationPoints

    weights         = zeros(numOfStationaryPoints, 1);
    weights_sum     = 0;
    coefficients    = 0;
    
    for stationaryPointIndex = 1:numOfStationaryPoints
        weights(stationaryPointIndex) = weightingFunction(norm(evaluationPoints(evaluationPointIndex,:)-stationaryPoints(stationaryPointIndex,:), 2));
        weights_sum = weights_sum + weights(stationaryPointIndex);
    end
    
    phi = weights/weights_sum;
    
    for stationaryPointIndex = 1:numOfStationaryPoints
        coefficients = coefficients + phi(stationaryPointIndex)*coefficients_local(stationaryPointIndex, :);
    end
    
    values(evaluationPointIndex) = polynomial(evaluationPoints(evaluationPointIndex, :))'*coefficients';
    
end

end