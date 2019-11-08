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