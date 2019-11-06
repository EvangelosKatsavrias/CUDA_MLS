function theta = ConstantDistribution(h, d, varargin)

theta = 1;
theta(d>h) = 0;

end