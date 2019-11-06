function theta = radialDistribution(h, d, varargin)

if nargin > 2; exponent = varargin{1};
else exponent = 4;
end


theta = 1./( (d/h).^exponent +1e-12 );

% theta = 1 /( (d/h)^4 +1e-5 );

% theta = exp(-(d/h^2)^2 ) /( sqrt(pi*h) );
% theta = exp(-(d/h^2)^2 );

% theta = 1/((d/h)^2 + 1e-12);

% theta = exp(-(d/h)^4 );

%theta(d>h) = 0;

end