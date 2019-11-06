function [coefficients, conditionNumber] = findCoefficients_ClassicLeastSquares(samplePoints, sampleValues, polynomial, solver, condNumEvaluator)

numOfSamplePoints   = size(samplePoints, 1);
numOfMonomials      = length(polynomial(samplePoints(1,:)));
Sum_bbT             = zeros(numOfMonomials, numOfMonomials);
Sum_bf              = zeros(numOfMonomials, 1);


for samplePointIndex = 1:numOfSamplePoints
    b       = polynomial(samplePoints(samplePointIndex,:));
    Sum_bbT = Sum_bbT   + b*b';
    Sum_bf  = Sum_bf    + b*sampleValues(samplePointIndex);
end


coefficients    = solver(Sum_bbT,Sum_bf);
conditionNumber = condNumEvaluator(Sum_bbT);

end