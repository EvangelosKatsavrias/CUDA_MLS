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
A = dlmread('../giannakoglou/data');
samplePoints = [0 1 3 6 7 9 15 20 21 23 24 28 30]';
sampleValues = [1 2 0 3 4 8 5 6 9 10 5 2 1]';
samplePoints = A(:,1); 
sampleValues = A(:,2);

evaluationPoints = linspace(samplePoints(1),samplePoints(end),101)';
stationaryPoints = evaluationPoints;
numOfStationaryPoints = size(stationaryPoints, 1);
weightingFunction_span = 0.1;
numOfSamplePoints = size(samplePoints, 1);


figure(1); plot(samplePoints, sampleValues, 'o'); hold all; axis equal;


%%  Interpolated surface plot
V = interp1(samplePoints(:,1),sampleValues,evaluationPoints,'pchip');
figure(1);
plot(evaluationPoints, V);
legend({'Data points' 'Piecewise cubic interpolation'});


%%  Curve plot - Least squares
evaluations_LS = zeros(size(evaluationPoints,1),1);
[coefficients_LS, conditionNumber_LS] = findCoefficients_ClassicLeastSquares(samplePoints, sampleValues, @(x)polynomialD1M2(x), @(A,b)linsolve(A,b), @(A)cond(A));
for coordX = 1:size(evaluationPoints,1)
        evaluations_LS(coordX) = polynomialD1M2(evaluationPoints(coordX))'*coefficients_LS;
end

figure(1);
plot(evaluationPoints, evaluations_LS);
legend({'Data points' 'Piecewise cubic interpolation' 'Simple least squares, 2nd order'});


%%  Curve plot - Weighted least squares - global Gaussian weighting functions
weightPerStationaryPoint = zeros(numOfStationaryPoints, numOfSamplePoints);
for stationaryPointIndex=1:numOfStationaryPoints
    for samplePointIndex = 1:numOfSamplePoints
        weightPerStationaryPoint(stationaryPointIndex, samplePointIndex) = gaussianDistribution(weightingFunction_span, norm(samplePoints(samplePointIndex, :)-stationaryPoints(stationaryPointIndex, :), 2));
    end
end

coefficients_local = findCoefficients_WeightedLeastSquares(samplePoints, sampleValues, @(x)polynomialD1M2(x), @(A,b)linsolve(A,b), @(A)cond(A), weightPerStationaryPoint);
evaluations_WLSGlobal_Gaussian = evaluatePoints_WeightedLeastSquares_global(evaluationPoints, stationaryPoints, @(x)polynomialD1M2(x), @(d)gaussianDistribution(weightingFunction_span, d), coefficients_local);

figure(1);
plot(evaluationPoints, evaluations_WLSGlobal_Gaussian);
legend({'Data points' 'Piecewise cubic interpolation' 'Simple least squares, 2nd order' 'Global weighted least squares - Gauss function'});


%%  Surface plot - Weighted least squares - global Wendland weighting functions
weightPerStationaryPoint = zeros(numOfStationaryPoints, numOfSamplePoints);
for stationaryPointIndex=1:numOfStationaryPoints
    for samplePointIndex = 1:numOfSamplePoints
        weightPerStationaryPoint(stationaryPointIndex, samplePointIndex) = WendlandDistribution(weightingFunction_span, norm(samplePoints(samplePointIndex, :)-stationaryPoints(stationaryPointIndex, :), 2));
    end
end

coefficients_local = findCoefficients_WeightedLeastSquares(samplePoints, sampleValues, @(x)polynomialD1M2(x), @(A,b)linsolve(A,b), @(A)cond(A), weightPerStationaryPoint);
evaluations_WLSGlobal_Wendland = evaluatePoints_WeightedLeastSquares_global(evaluationPoints, stationaryPoints, @(x)polynomialD1M2(x), @(d)WendlandDistribution(weightingFunction_span, d), coefficients_local);

figure(1);
plot(evaluationPoints, evaluations_WLSGlobal_Wendland);
legend({'Data points' 'Piecewise cubic interpolation' 'Simple least squares, 2nd order' 'Global weighted least squares - Gauss function' 'Global weighted least squares - Wendland function'});


%%  Surface plot - Weighted least squares - global radial weighting functions
weightPerStationaryPoint = zeros(numOfStationaryPoints, numOfSamplePoints);
for stationaryPointIndex=1:numOfStationaryPoints
    for samplePointIndex = 1:numOfSamplePoints
        weightPerStationaryPoint(stationaryPointIndex, samplePointIndex) = radialDistribution(weightingFunction_span, norm(samplePoints(samplePointIndex, :)-stationaryPoints(stationaryPointIndex, :), 2));
    end
end

coefficients_local = findCoefficients_WeightedLeastSquares(samplePoints, sampleValues, @(x)polynomialD1M2(x), @(A,b)linsolve(A,b), @(A)cond(A), weightPerStationaryPoint);
evaluations_WLSGlobal_radial = evaluatePoints_WeightedLeastSquares_global(evaluationPoints, stationaryPoints, @(x)polynomialD1M2(x), @(d)radialDistribution(weightingFunction_span, d), coefficients_local);

figure(1);
plot(evaluationPoints, evaluations_WLSGlobal_radial);
legend({'Data points' 'Piecewise cubic interpolation' 'Simple least squares, 2nd order' 'Global weighted least squares - Gauss function' 'Global weighted least squares - Wendland function' 'Global weighted least squares - Radial function'});


%%  Curve plot - Weighted least squares - global Gaussian weighting functions - LCS
weightPerStationaryPoint = zeros(numOfStationaryPoints, numOfSamplePoints);
for stationaryPointIndex=1:numOfStationaryPoints
    for samplePointIndex = 1:numOfSamplePoints
        weightPerStationaryPoint(stationaryPointIndex, samplePointIndex) = gaussianDistribution(weightingFunction_span, norm(samplePoints(samplePointIndex, :)-stationaryPoints(stationaryPointIndex, :), 2));
    end
end

coefficients_local = findCoefficients_WeightedLeastSquaresLCS(samplePoints, stationaryPoints, sampleValues, @(x)polynomialD1M2(x), @(A,b)linsolve(A,b), @(A)cond(A), weightPerStationaryPoint);
evaluations_WLSGlobal_Gaussian = evaluatePoints_WeightedLeastSquares_globalLCS(evaluationPoints, stationaryPoints, @(x)polynomialD1M2(x), @(d)gaussianDistribution(weightingFunction_span, d), coefficients_local);

figure(1);
plot(evaluationPoints, evaluations_WLSGlobal_Gaussian);
legend({'Data points' 'Piecewise cubic interpolation' 'Simple least squares, 2nd order' 'Global weighted least squares - Gauss function' 'Global weighted least squares - Wendland function' 'Global weighted least squares - Radial function' 'Global weighted least squares in LCS - Gauss function'});


%%  Surface plot - Weighted least squares - global Wendland weighting functions - LCS
weightPerStationaryPoint = zeros(numOfStationaryPoints, numOfSamplePoints);
for stationaryPointIndex=1:numOfStationaryPoints
    for samplePointIndex = 1:numOfSamplePoints
        weightPerStationaryPoint(stationaryPointIndex, samplePointIndex) = WendlandDistribution(weightingFunction_span, norm(samplePoints(samplePointIndex, :)-stationaryPoints(stationaryPointIndex, :), 2));
    end
end

coefficients_local = findCoefficients_WeightedLeastSquaresLCS(samplePoints, stationaryPoints, sampleValues, @(x)polynomialD1M2(x), @(A,b)linsolve(A,b), @(A)cond(A), weightPerStationaryPoint);
evaluations_WLSGlobal_Wendland = evaluatePoints_WeightedLeastSquares_globalLCS(evaluationPoints, stationaryPoints, @(x)polynomialD1M2(x), @(d)WendlandDistribution(weightingFunction_span, d), coefficients_local);

figure(1);
plot(evaluationPoints, evaluations_WLSGlobal_Wendland);
legend({'Data points' 'Piecewise cubic interpolation' 'Simple least squares, 2nd order' 'Global weighted least squares - Gauss function' 'Global weighted least squares - Wendland function' 'Global weighted least squares - Radial function' 'Global weighted least squares in LCS - Gauss function' 'Global weighted least squares in LCS - Wendland function'});


%%  Surface plot - Weighted least squares - global radial weighting functions - LCS
weightPerStationaryPoint = zeros(numOfStationaryPoints, numOfSamplePoints);
for stationaryPointIndex=1:numOfStationaryPoints
    for samplePointIndex = 1:numOfSamplePoints
        weightPerStationaryPoint(stationaryPointIndex, samplePointIndex) = radialDistribution(weightingFunction_span, norm(samplePoints(samplePointIndex, :)-stationaryPoints(stationaryPointIndex, :), 2));
    end
end

coefficients_local = findCoefficients_WeightedLeastSquaresLCS(samplePoints, stationaryPoints, sampleValues, @(x)polynomialD1M2(x), @(A,b)linsolve(A,b), @(A)cond(A), weightPerStationaryPoint);
evaluations_WLSGlobal_radial = evaluatePoints_WeightedLeastSquares_globalLCS(evaluationPoints, stationaryPoints, @(x)polynomialD1M2(x), @(d)radialDistribution(weightingFunction_span, d), coefficients_local);

figure(1);
plot(evaluationPoints, evaluations_WLSGlobal_radial);
legend({'Data points' 'Piecewise cubic interpolation' 'Simple least squares, 2nd order' 'Global weighted least squares - Gauss function' 'Global weighted least squares - Wendland function' 'Global weighted least squares - Radial function' 'Global weighted least squares in LCS - Gauss function' 'Global weighted least squares in LCS - Wendland function' 'Global weighted least squares in LCS - Radial function'});


%%  Curve plot - Moving least squares - Gaussian weighting functions
weightPerStationaryPoint = zeros(numOfStationaryPoints, numOfSamplePoints);
for stationaryPointIndex=1:numOfStationaryPoints
    for samplePointIndex = 1:numOfSamplePoints
        weightPerStationaryPoint(stationaryPointIndex, samplePointIndex) = gaussianDistribution(weightingFunction_span, norm(samplePoints(samplePointIndex, :)-stationaryPoints(stationaryPointIndex, :), 2));
    end
end

coefficients_local = findCoefficients_WeightedLeastSquares(samplePoints, sampleValues, @(x)polynomialD1M2(x), @(A,b)linsolve(A,b), @(A)cond(A), weightPerStationaryPoint);
evaluations_MLS_Gaussian = evaluatePoints_MovingLeastSquares(evaluationPoints, @(x)polynomialD1M2(x), coefficients_local);

figure(1);
plot(evaluationPoints, evaluations_MLS_Gaussian);
legend({'Data points' 'Piecewise cubic interpolation' 'Simple least squares, 2nd order' 'Global weighted least squares - Gauss function' 'Global weighted least squares - Wendland function' 'Global weighted least squares - Radial function' 'Global weighted least squares in LCS - Gauss function' 'Global weighted least squares in LCS - Wendland function' 'Global weighted least squares in LCS - Radial function' 'Moving least squares - Gauss function'});


%%  Curve plot - Moving least squares - Wendland weighting functions
weightPerStationaryPoint = zeros(numOfStationaryPoints, numOfSamplePoints);
for stationaryPointIndex=1:numOfStationaryPoints
    for samplePointIndex = 1:numOfSamplePoints
        weightPerStationaryPoint(stationaryPointIndex, samplePointIndex) = WendlandDistribution(weightingFunction_span, norm(samplePoints(samplePointIndex, :)-stationaryPoints(stationaryPointIndex, :), 2));
    end
end

coefficients_local = findCoefficients_WeightedLeastSquares(samplePoints, sampleValues, @(x)polynomialD1M2(x), @(A,b)linsolve(A,b), @(A)cond(A), weightPerStationaryPoint);
evaluations_MLS_Wendland = evaluatePoints_MovingLeastSquares(evaluationPoints, @(x)polynomialD1M2(x), coefficients_local);

figure(1);
plot(evaluationPoints, evaluations_MLS_Wendland);
legend({'Data points' 'Piecewise cubic interpolation' 'Simple least squares, 2nd order' 'Global weighted least squares - Gauss function' 'Global weighted least squares - Wendland function' 'Global weighted least squares - Radial function' 'Global weighted least squares in LCS - Gauss function' 'Global weighted least squares in LCS - Wendland function' 'Global weighted least squares in LCS - Radial function' 'Moving least squares - Gauss function' 'Moving least squares - Wendland function'});



%%  Curve plot - Moving least squares - Radial weighting functions
weightPerStationaryPoint = zeros(numOfStationaryPoints, numOfSamplePoints);
for stationaryPointIndex=1:numOfStationaryPoints
    for samplePointIndex = 1:numOfSamplePoints
        weightPerStationaryPoint(stationaryPointIndex, samplePointIndex) = radialDistribution(weightingFunction_span, norm(samplePoints(samplePointIndex, :)-stationaryPoints(stationaryPointIndex, :), 2));
    end
end

coefficients_local = findCoefficients_WeightedLeastSquares(samplePoints, sampleValues, @(x)polynomialD1M2(x), @(A,b)linsolve(A,b), @(A)cond(A), weightPerStationaryPoint);
evaluations_MLS_Radial = evaluatePoints_MovingLeastSquares(evaluationPoints, @(x)polynomialD1M2(x), coefficients_local);

figure(1);
plot(evaluationPoints, evaluations_MLS_Radial);
legend({'Data points' 'Piecewise cubic interpolation' 'Simple least squares, 2nd order' 'Global weighted least squares - Gauss function' 'Global weighted least squares - Wendland function' 'Global weighted least squares - Radial function' 'Global weighted least squares in LCS - Gauss function' 'Global weighted least squares in LCS - Wendland function' 'Global weighted least squares in LCS - Radial function' 'Moving least squares - Gauss function' 'Moving least squares - Wendland function' 'Moving least squares - Radial function'});


%%  Curve plot - Moving least squares in LCS - Gaussian weighting functions
weightPerStationaryPoint = zeros(numOfStationaryPoints, numOfSamplePoints);
for stationaryPointIndex=1:numOfStationaryPoints
    for samplePointIndex = 1:numOfSamplePoints
        weightPerStationaryPoint(stationaryPointIndex, samplePointIndex) = gaussianDistribution(weightingFunction_span, norm(samplePoints(samplePointIndex, :)-stationaryPoints(stationaryPointIndex, :), 2));
    end
end

coefficients_local = findCoefficients_WeightedLeastSquaresLCS(samplePoints, stationaryPoints, sampleValues, @(x)polynomialD1M2(x), @(A,b)linsolve(A,b), @(A)cond(A), weightPerStationaryPoint);

figure(1);
plot(evaluationPoints, coefficients_local(:, 1));
legend({'Data points' 'Piecewise cubic interpolation' 'Simple least squares, 2nd order' 'Global weighted least squares - Gauss function' 'Global weighted least squares - Wendland function' 'Global weighted least squares - Radial function' 'Global weighted least squares in LCS - Gauss function' 'Global weighted least squares in LCS - Wendland function' 'Global weighted least squares in LCS - Radial function' 'Moving least squares - Gauss function' 'Moving least squares - Wendland function' 'Moving least squares - Radial function' 'Moving least squares in LCS - Gauss function'});


%%  Curve plot - Moving least squares in LCS - Wendland weighting functions
weightPerStationaryPoint = zeros(numOfStationaryPoints, numOfSamplePoints);
for stationaryPointIndex=1:numOfStationaryPoints
    for samplePointIndex = 1:numOfSamplePoints
        weightPerStationaryPoint(stationaryPointIndex, samplePointIndex) = WendlandDistribution(weightingFunction_span, norm(samplePoints(samplePointIndex, :)-stationaryPoints(stationaryPointIndex, :), 2));
    end
end

coefficients_local_wendland = findCoefficients_WeightedLeastSquaresLCS(samplePoints, stationaryPoints, sampleValues, @(x)polynomialD1M2(x), @(A,b)linsolve(A,b), @(A)cond(A), weightPerStationaryPoint);

figure(1);
plot(evaluationPoints, coefficients_local_wendland(:, 1));
legend({'Data points' 'Piecewise cubic interpolation' 'Simple least squares, 2nd order' 'Global weighted least squares - Gauss function' 'Global weighted least squares - Wendland function' 'Global weighted least squares - Radial function' 'Global weighted least squares in LCS - Gauss function' 'Global weighted least squares in LCS - Wendland function' 'Global weighted least squares in LCS - Radial function' 'Moving least squares - Gauss function' 'Moving least squares - Wendland function' 'Moving least squares - Radial function' 'Moving least squares in LCS - Gauss function' 'Moving least squares in LCS - Wendland function'});


%%  Curve plot - Moving least squares in LCS - Radial weighting functions
weightPerStationaryPoint = zeros(numOfStationaryPoints, numOfSamplePoints);
for stationaryPointIndex=1:numOfStationaryPoints
    for samplePointIndex = 1:numOfSamplePoints
        weightPerStationaryPoint(stationaryPointIndex, samplePointIndex) = radialDistribution(weightingFunction_span, norm(samplePoints(samplePointIndex, :)-stationaryPoints(stationaryPointIndex, :), 2));
    end
end

coefficients_local = findCoefficients_WeightedLeastSquaresLCS(samplePoints, stationaryPoints, sampleValues, @(x)polynomialD1M2(x), @(A,b)linsolve(A,b), @(A)cond(A), weightPerStationaryPoint);

figure(1);
plot(evaluationPoints, coefficients_local(:, 1));
legend({'Data points' 'Piecewise cubic interpolation' 'Simple least squares, 2nd order' 'Global weighted least squares - Gauss function' 'Global weighted least squares - Wendland function' 'Global weighted least squares - Radial function' 'Global weighted least squares in LCS - Gauss function' 'Global weighted least squares in LCS - Wendland function' 'Global weighted least squares in LCS - Radial function' 'Moving least squares - Gauss function' 'Moving least squares - Wendland function' 'Moving least squares - Radial function' 'Moving least squares in LCS - Gauss function' 'Moving least squares in LCS - Wendland function' 'Moving least squares in LCS - Radial function'});

% 
% B = dlmread('../giannakoglou/rmlsl');
% plot(B(:,1), B(:,2));
% legend({'Data points' 'Piecewise cubic interpolation' 'Simple least squares, 2nd order' 'Global weighted least squares - Gauss function' 'Global weighted least squares - Wendland function' 'Global weighted least squares - Radial function' 'Global weighted least squares in LCS - Gauss function' 'Global weighted least squares in LCS - Wendland function' 'Global weighted least squares in LCS - Radial function' 'Moving least squares - Gauss function' 'Moving least squares - Wendland function' 'Moving least squares - Radial function' 'Moving least squares in LCS - Gauss function' 'Moving least squares in LCS - Wendland function' 'Moving least squares in LCS - Radial function' 'Giannakoglou'});
% 
