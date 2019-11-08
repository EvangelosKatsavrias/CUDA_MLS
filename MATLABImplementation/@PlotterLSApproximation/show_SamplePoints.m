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

function show_SamplePoints(obj)

    figure(obj.handle_figure); hold all;
    if (~isempty(obj.handle_samplePoints)); delete(obj.handle_samplePoints); end
    obj.handle_samplePoints = plot(obj.m_LeastSquaresApproximations(1).m_samplePoints, obj.m_LeastSquaresApproximations(1).m_sampleValues, 'o');

end