degVariants = [ 0 ]; %[0 1 2 3 4]
expVariants = [ 4 ] ;%[ 1.5 2 2.5 3 3.5 4 4.5 ];
spanVariants = [ 0.04 ];%[ 0.05 0.1 0.5 1 2 ];
n_deg = length(degVariants);
n_exp = length(expVariants);
n_span = length(spanVariants);
error = zeros(n_deg, n_exp, n_span);
for iterDeg=1:n_deg
    for iterExp=1:n_exp
        for iterSpan=1:n_span
%%  Data
A = dlmread('./data');
samplePoints        = A(:,1);
sampleValues        = A(:,2);
evaluationPoints    = linspace(samplePoints(1),samplePoints(end), 101)';
polynomialBasis     = Polynomial(1, degVariants(iterDeg));

weightingFunction_1 = WeightingFunction.Radial; weightingFunction_1.set_Span(spanVariants(iterSpan));
stationaryPoints = evaluationPoints;

weightingFunction_1.set_Exponent(expVariants(iterExp));
%%  Test objects2
movingLSApproxLCS   = MovingLeastSquares_LCS('samplePoints', samplePoints, 'sampleValues', sampleValues, 'evaluationPoints', evaluationPoints, polynomialBasis, 'stationaryPoints', stationaryPoints, weightingFunction_1);
movingLSApproxLCS.evaluatePoints;
% movingLSApproxLCS.evaluateDerivatives;
%plot(evaluationPoints, weightedLSApproxLCS.m_conditionNumber);
testPlotter = PlotterLSApproximation;
testPlotter.m_LeastSquaresApproximations = [ movingLSApproxLCS]; 
if ( iterDeg==1 && iterExp==1 && iterSpan==1 ); testPlotter.show_SamplePoints; end
testPlotter.plot_All;

for i = 1:size(samplePoints, 1)
    coord = samplePoints(i);
    pos = []; tol = 1e-6;
    while (isempty(pos))
        pos = find((coord-tol < evaluationPoints & evaluationPoints < coord+tol));
        tol = tol*2;
    end
    val = movingLSApproxLCS.m_evaluations(pos(1));
    error(iterDeg, iterExp, iterSpan) = error(iterDeg, iterExp, iterSpan) +((val-sampleValues(i))/sampleValues(i))^2;
end
error(iterDeg, iterExp, iterSpan) = sqrt(error(iterDeg, iterExp, iterSpan));
% error
% disp('The interpolation error is: '); disp(error);
        end
    end
end

% refined parametric for r<2, 1.5 < exp <4.5, 
%%
% hold all; n_surfvar = n_exp; n_var1 = n_span; n_var2 = n_deg;
% var1 = spanVariants; var2 = degVariants;
% for surfvar = 1:n_surfvar
% surf( var1, var2, squeeze(error(:,surfvar,:)), ...
%         surfvar/n_surfvar*ones(n_var1,n_var2)','FaceAlpha',surfvar/n_surfvar,'EdgeColor',surfvar/n_surfvar/5*[surfvar/n_surfvar/5 1 1]);
% end

%     xlabel('Support radius r');
%     ylabel('Degree n');
%     zlabel('Approximation relative error');
%     titl='Parametric surface plot';
%     title(titl);
%     legend('r=.02','r=.05','r=.1', 'r=.5', 'r=1', 'r=2')
%     set(gca,'XTick', var1)
%     set(gca,'XTickLabel',genxticklabl)
%     set(gca,'YTick', var2)
%     set(gca,'YTickLabel',{'1st' '2nd' '3rd' '4th'})
%     grid on
%     view(110,30)
%     set(gcf,'renderer','opengl')
%      set(gcf,'renderer','painters')
%     print(figure(1),'-tiff','text.tiff','-opengl')
%     hgexport(figure(1),'text.tiff',hgexport('factorystyle'),'format','tiff')
% saveas(figure(1), 'parametric_final.png')
%     view(120,10)
%     saveas(figure(1),'hp_parambwsurf2.tiff')

% surf( var1, var2, squeeze(error(:,1,:)), squeeze(error(:,1,:)), ...
%         'FaceAlpha', 1,'EdgeColor', 1/5*[1/5 1 1]);


%%  Test objects
% classicLSApprox     = ClassicLeastSquares('samplePoints', samplePoints, 'sampleValues', sampleValues, 'evaluationPoints', evaluationPoints, polynomialBasis);
% weightedLSApprox    = WeightedLeastSquares(classicLSApprox, 'stationaryPoints', stationaryPoints, weightingFunction_1);
% weightedLSApproxLCS = WeightedLeastSquares_LCS(weightedLSApprox);
% movingLSApprox      = MovingLeastSquares(weightedLSApprox);
% movingLSApproxLCS   = MovingLeastSquares_LCS(weightedLSApprox);
% 
% classicLSApprox.evaluatePoints; weightedLSApprox.evaluatePoints;
% weightedLSApproxLCS.evaluatePoints; movingLSApprox.evaluatePoints;
% movingLSApproxLCS.evaluatePoints;
% 
% %%  Plotter
% testPlotter = PlotterLSApproximation;
% testPlotter.m_LeastSquaresApproximations = ...
%    [classicLSApprox]; 
%    weightedLSApprox; weightedLSApproxLCS; 
%    movingLSApprox; movingLSApproxLCS];
% testPlotter.plot_All;
% 
% legend({'Data points' 'Classic LS' 'Weighted LS' 'Weighted LS LCS' 'Moving LS' 'Moving LS LCS'});
%%  Interpolated surface plot
%V = interp1(samplePoints(:,1), sampleValues, evaluationPoints, 'pchip');
%plot(evaluationPoints, V);
%legend({'Data points' 'Piecewise cubic interpolation'});
