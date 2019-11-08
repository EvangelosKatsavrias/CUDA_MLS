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

function show_Approximation(obj, approximationIndices)

figure(obj.handle_figure); hold all;

for approximationIndex = approximationIndices
    lineWidth = 1;
%    if approximationIndex == 2; lineWidth = 3; end;
    if (length(obj.handles_approximations)>approximationIndex); if (obj.handles_approximations(approximationIndex))>0; delete(obj.handles_approximations(approximationIndex)); end; end
    obj.handles_approximations(approximationIndex) = plot(obj.m_LeastSquaresApproximations(approximationIndex).m_evaluationPoints, obj.m_LeastSquaresApproximations(approximationIndex).m_evaluations, 'LineWidth', lineWidth);
end

end