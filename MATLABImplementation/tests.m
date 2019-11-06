%%  Data
points = ...
[0 0
 5 0
 5 5
 0 5
 1 1
 4 4
 4 0
 0 4
 2 2
 3 3
 2 3
 3 2
];
fieldValues = [0 0 0 0 1 1 1 1 5 5 5 5]';
[X,Y] = meshgrid(0:0.2:5,0:0.2:5);


%%  Solutions
[coefficients_LS, conditionNumber_LS]                       = findCoefficients_LeastSquares(points, fieldValues, @(x)polynomialD2M2(x));
[coefficients_WLS_Gaussian, conditionNumber_WLS_Gaussian]   = findCoefficients_WeightedLeastSquares(points, fieldValues, @(x)polynomialD2M2(x), @(x_center, x_point)gaussianDistribution(5, x_center, x_point));
[coefficients_WLS_Wendland, conditionNumber_WLS_Wendland]   = findCoefficients_WeightedLeastSquares(points, fieldValues, @(x)polynomialD2M2(x), @(x_center, x_point)WendlandDistribution(5, x_center, x_point));
[coefficients_WLS_radial, conditionNumber_WLS_radial]       = findCoefficients_WeightedLeastSquares(points, fieldValues, @(x)polynomialD2M2(x), @(x_center, x_point)radialDistribution(5, x_center, x_point));

coefficients_WLS_global_Gaussian    = zeros(length(coefficients_LS),size(X,1),size(X,2));
coefficients_WLS_global_Wendland    = zeros(length(coefficients_LS),size(X,1),size(X,2));
coefficients_WLS_global_radial      = zeros(length(coefficients_LS),size(X,1),size(X,2));
for coordX = 1:size(X,1)
    for coordY = 1:size(X,2)
        coefficients_WLS_global_Gaussian(:,coordX,coordY)     = findCoefficients_WeightedLeastSquares_global(points, fieldValues, @(x)polynomialD2M2(x), @(x_center, x_point)gaussianDistribution(5, x_center, x_point), 2);
        coefficients_WLS_global_Wendland(:,coordX,coordY)     = findCoefficients_WeightedLeastSquares_global(points, fieldValues, @(x)polynomialD2M2(x), @(x_center, x_point)WendlandDistribution(5, x_center, x_point), 2);
        coefficients_WLS_global_radial(:,coordX,coordY)       = findCoefficients_WeightedLeastSquares_global(points, fieldValues, @(x)polynomialD2M2(x), @(x_center, x_point)radialDistribution(5, x_center, x_point), 2);
    end
end


%%  Interpolated surface plot
figure(1); axis equal
V=griddata(points(:,1),points(:,2),fieldValues,X,Y,'cubic');
surf(X,Y,V);


%%  Surface plot - Least squares
V1=zeros(size(X,1),size(X,2));
for coordX = 1:size(X,1)
    for coordY = 1:size(X,2)
        V1(coordX,coordY) = polynomialD2M2([X(coordX,coordY); Y(coordX,coordY)])'*coefficients_LS;
    end
end

figure(2); axis equal
surf(X,Y,V1);


%%  Surface plot - Weighted least squares - global Gaussian weighting functions
V2=zeros(size(X,1),size(X,2));
for coordX = 1:size(X,1)
    for coordY = 1:size(X,2)
        V2(coordX,coordY) = polynomialD2M2([X(coordX,coordY); Y(coordX,coordY)])'*coefficients_WLS_global_Gaussian(:,coordX,coordY);
    end
end

figure(3); axis equal
surf(X,Y,V2);


%%  Surface plot - Weighted least squares - global Wendland weighting functions
V3=zeros(size(X,1),size(X,2));
for coordX = 1:size(X,1)
    for coordY = 1:size(X,2)
        V3(coordX,coordY) = polynomialD2M2([X(coordX,coordY); Y(coordX,coordY)])'*coefficients_WLS_global_Wendland(:,coordX,coordY);
    end
end

figure(4); axis equal
surf(X,Y,V3);


%%  Surface plot - Weighted least squares - global radial weighting functions
V4=zeros(size(X,1),size(X,2));
for coordX = 1:size(X,1)
    for coordY = 1:size(X,2)
        V4(coordX,coordY) = polynomialD2M2([X(coordX,coordY); Y(coordX,coordY)])'*coefficients_WLS_global_radial(:,coordX,coordY);
    end
end

figure(5); axis equal
surf(X,Y,V4);


%%  Surface plot - Weighted least squares - local Wendland weighting functions
V5=zeros(size(X,1),size(X,2));
for coordX = 1:size(X,1)
    for coordY = 1:size(X,2)
        V5(coordX,coordY) = polynomialD2M2([X(coordX,coordY); Y(coordX,coordY)])'*coefficients_WLS_Wendland(:,coordX,coordY);
    end
end

figure(5); axis equal
surf(X,Y,V5);