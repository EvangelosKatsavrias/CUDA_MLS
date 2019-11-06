function theta = WendlandDistribution(h, d, varargin)

if nargin > 2; exponent = varargin{1};
else exponent = 4;
end

theta = (1-(d/h)).^4.*(4.*d/h+1)./( (d/h).^exponent +1e-12 );
% theta = (1-(d/h)).^6.*(35*(d/h).^2+18*d/h+3);

theta(d>h) = 0;

end