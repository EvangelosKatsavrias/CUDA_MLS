%% 1D
function du = derivativesOfMLS(x, x_u, u)
wx_eval = @(x) -4/(x.^5+1e-12);
w_eval = @(x) 1/(x.^4+1e-12);

du = zeros( length(x), 1 );

for evalPoint = 1:length(x)

    p1 = 0;
    p2 = 0;
    sum_w = 0;
    sum_wx = 0;
    
    for samplePoint = 1:length(u)
        d = x_u(samplePoint) -x(evalPoint);
        wx = wx_eval( d );
        w = w_eval( d );
        b   = polynomial( d );
        p1 = p1 + wx*u(samplePoint);
        p2 = p2 + w*u(samplePoint);
        sum_w = sum_w + w;
        sum_wx = sum_wx +wx;
    end

    du(evalPoint) = p1/sum_w -sum_wx/sum_w^2;

end

end
