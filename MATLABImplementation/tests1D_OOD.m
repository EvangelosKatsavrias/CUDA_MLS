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

%%  Data
A = dlmread('./data3');
samplePoints        = A(:,1);
sampleValues        = A(:,2);
evaluationPoints    = linspace(samplePoints(1),samplePoints(end), 101)';
evaluationPoints = [evaluationPoints; samplePoints];
evaluationPoints = sort(evaluationPoints);
evaluationPoints = unique(evaluationPoints);
localPnt = evaluationPoints(42);
%evaluationPoints = evaluationPoints(evaluationPoints <= localPnt);
polynomialBasis     = Polynomial(1, 4);

weightingFunction_1 = WeightingFunction.Radial; weightingFunction_1.set_Span(0.2);
stationaryPoints = localPnt;

weightingFunction_1.set_Exponent(4);

%%  Test objects
classicLSApprox     = ClassicLeastSquares('samplePoints', samplePoints, 'sampleValues', sampleValues, 'evaluationPoints', evaluationPoints, polynomialBasis);
weightedLSApprox    = WeightedLeastSquares(classicLSApprox, 'stationaryPoints', stationaryPoints, weightingFunction_1);
% weightedLSApproxLCS = WeightedLeastSquares_LCS(weightedLSApprox);
% movingLSApprox      = MovingLeastSquares(weightedLSApprox);
movingLSApproxLCS   = MovingLeastSquares_LCS(weightedLSApprox);
% 
classicLSApprox.evaluatePoints; 
weightedLSApprox.evaluatePoints;
% weightedLSApproxLCS.evaluatePoints;
% movingLSApprox.evaluatePoints;
movingLSApproxLCS.evaluatePoints;
% 
%%  Plotter
testPlotter = PlotterLSApproximation;
testPlotter.m_LeastSquaresApproximations = ...
   [classicLSApprox weightedLSApprox movingLSApproxLCS]; 
%    weightedLSApprox; weightedLSApproxLCS; 
%    movingLSApprox; movingLSApproxLCS];

testPlotter.show_SamplePoints;
testPlotter.plot_All;
ylim([-0.2 1.2]);
% legend({'Data points' 'Classic LS' 'Weighted LS' 'Weighted LS LCS' 'Moving LS' 'Moving LS LCS'});
% legend({'Data points' 'Classic LS' 'Local LS'});
%%  Interpolated surface plot
%V = interp1(samplePoints(:,1), sampleValues, evaluationPoints, 'pchip');
%plot(evaluationPoints, V);
%legend({'Data points' 'Piecewise cubic interpolation'});
