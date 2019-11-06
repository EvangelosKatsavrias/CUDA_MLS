function theta = gaussianDistribution(h, d, varargin)

theta = exp(-(d/h).^2);

theta(d>h) = 0;

end