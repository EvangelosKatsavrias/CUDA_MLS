function dudx = findDerivatives_MovingLeastSquaresLCS(samplePoints, stationaryPoints, sampleValues, polynomial, weightPerStationaryPoint)

numOfStationaryPoints   = size(weightPerStationaryPoint, 1);
numOfSamplePoints       = size(samplePoints, 1);
numOfMonomials          = length(polynomial(samplePoints(1,:)));
dudx                    = zeros(numOfStationaryPoints, 1);

% eval_dw = @(x) -4/(x^5+1e-12);

eval_dw = @(x) -4/(x^5+1e-12);


% for stationaryPointIndex = 1:numOfStationaryPoints
% 
%     A  = zeros(numOfMonomials, numOfMonomials);
%     dA = zeros(numOfMonomials, numOfMonomials);
%     
%     for samplePointIndex = 1:numOfSamplePoints
%         b       = polynomial(samplePoints(samplePointIndex, :)-stationaryPoints(stationaryPointIndex, :));
%         w       = weightPerStationaryPoint(stationaryPointIndex, samplePointIndex);
%         theta_b = w*b; 
%         db      = dbdx( samplePoints(samplePointIndex, :) -stationaryPoints(stationaryPointIndex, :), numOfMonomials );
%         dw      = eval_dw( samplePoints(samplePointIndex, :) -stationaryPoints(stationaryPointIndex, :) );
% 
%         A       = A  + theta_b*b';
%         dA      = dA + db*theta_b' +theta_b*db' +theta_b*b'*dw;
%     end
% 
%     inv_A = inv(A);
%     dinv_A = -inv_A*dA*inv_A;
%     
%     for samplePointIndex = 1:numOfSamplePoints
%         b   = polynomial(samplePoints(samplePointIndex, :)-stationaryPoints(stationaryPointIndex, :));
%         w   = weightPerStationaryPoint(stationaryPointIndex, samplePointIndex);
%         db  = dbdx( samplePoints(samplePointIndex, :) -stationaryPoints(stationaryPointIndex, :), numOfMonomials );
%         dw  = eval_dw( samplePoints(samplePointIndex, :) -stationaryPoints(stationaryPointIndex, :) );
% 
%         dudx(stationaryPointIndex) = dudx(stationaryPointIndex) + ( db'*inv_A*b*w +b'*(dinv_A*b*w +inv_A*( b*dw +db*w ) ) )*sampleValues(samplePointIndex);
%     end
%     
% end



for stationaryPointIndex = 1:numOfStationaryPoints
    
    sum_di = 0;
    
    for samplePointIndex = 1:numOfSamplePoints

        d_i = samplePoints(samplePointIndex, :)-stationaryPoints(stationaryPointIndex, :);

        if d_i == 0
            dudx(stationaryPointIndex) = 0;
            continue;
        end
        
        dx = d_i;
        
        sum_di = sum_di +d_i^-2;
        
        for samplePointIndex_j = 1:numOfSamplePoints
            
            if samplePointIndex_j == samplePointIndex; continue; end;
            
            d_j = samplePoints(samplePointIndex_j, :)-stationaryPoints(stationaryPointIndex, :);
        
            dz = sampleValues(samplePointIndex) -sampleValues(samplePointIndex_j);

            dudx(stationaryPointIndex) = dudx(stationaryPointIndex) + d_i^-4*d_j^-2*dx*dz*sampleValues(samplePointIndex);
            
        end
        
    end
    
    dudx(stationaryPointIndex) = dudx(stationaryPointIndex)/sum_di;
    
end


end



function db = dbdx(x, n)

db = zeros(n, 1);
if (n > 0); db(2) = 1; end

for monomial = 3:n
    db(monomial) = (monomial-1)*db(monomial-1)*x;
end

end