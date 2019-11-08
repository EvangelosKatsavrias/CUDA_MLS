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


%datpath = '../WLS_cpp_cuda/exe/';
datpath = '../WLS_cpp_cuda/exe/doublePrecision/';

fileID = fopen([datpath 'minQualities_30.dat']);
min_30 = textscan(fileID, '%f');
min_30 = min_30{1};
fclose(fileID);

fileID = fopen([datpath 'maxQualities_30.dat']);
max_30 = textscan(fileID, '%f');
max_30 = max_30{1};
fclose(fileID);

fileID = fopen([datpath 'meanQualities_30.dat']);
mean_30 = textscan(fileID, '%f');
mean_30 = mean_30{1};
fclose(fileID);

fileID = fopen([datpath 'stdDevQualities_30.dat']);
stddev_30 = textscan(fileID, '%f');
stddev_30 = stddev_30{1};
fclose(fileID);

fileID = fopen([datpath 'badElements_30.dat']);
badEl_30 = textscan(fileID, '%f');
badEl_30 = badEl_30{1};
fclose(fileID);


fileID = fopen([datpath 'minQualities_45.dat']);
min_45 = textscan(fileID, '%f');
min_45 = min_45{1};
fclose(fileID);

fileID = fopen([datpath 'maxQualities_45.dat']);
max_45 = textscan(fileID, '%f');
max_45 = max_45{1};
fclose(fileID);

fileID = fopen([datpath 'meanQualities_45.dat']);
mean_45 = textscan(fileID, '%f');
mean_45 = mean_45{1};
fclose(fileID);

fileID = fopen([datpath 'stdDevQualities_45.dat']);
stddev_45 = textscan(fileID, '%f');
stddev_45 = stddev_45{1};
fclose(fileID);

fileID = fopen([datpath 'badElements_45.dat']);
badEl_45 = textscan(fileID, '%f');
badEl_45 = badEl_45{1};
fclose(fileID);


fileID = fopen([datpath 'minQualities_60.dat']);
min_60 = textscan(fileID, '%f');
min_60 = min_60{1};
fclose(fileID);

fileID = fopen([datpath 'maxQualities_60.dat']);
max_60 = textscan(fileID, '%f');
max_60 = max_60{1};
fclose(fileID);

fileID = fopen([datpath 'meanQualities_60.dat']);
mean_60 = textscan(fileID, '%f');
mean_60 = mean_60{1};
fclose(fileID);

fileID = fopen([datpath 'stdDevQualities_60.dat']);
stddev_60 = textscan(fileID, '%f');
stddev_60 = stddev_60{1};
fclose(fileID);

fileID = fopen([datpath 'badElements_60.dat']);
badEl_60 = textscan(fileID, '%f');
badEl_60 = badEl_60{1};
fclose(fileID);



fileID = fopen([datpath 'minQualities_80.dat']);
min_80 = textscan(fileID, '%f');
min_80 = min_80{1};
fclose(fileID);

fileID = fopen([datpath 'maxQualities_80.dat']);
max_80 = textscan(fileID, '%f');
max_80 = max_80{1};
fclose(fileID);

fileID = fopen([datpath 'meanQualities_80.dat']);
mean_80 = textscan(fileID, '%f');
mean_80 = mean_80{1};
fclose(fileID);

fileID = fopen([datpath 'stdDevQualities_80.dat']);
stddev_80 = textscan(fileID, '%f');
stddev_80 = stddev_80{1};
fclose(fileID);

fileID = fopen([datpath 'badElements_80.dat']);
badEl_80 = textscan(fileID, '%f');
badEl_80 = badEl_80{1};
fclose(fileID);


%%

degVariants = [ 0 1 2 3 4 5];
expVariants = [ 0 1 2 3 4 5 6 ];
spanVariants = [ 6 8 10 12 15 18 20 25 30 40];
n_deg = length(degVariants);
n_exp = length(expVariants);
n_span = length(spanVariants);

error = reshape(max_80, [n_span, n_exp, n_deg] );
error = error(1:end,:,:);

hold all; n_surfvar = n_deg; n_var1 = n_exp; n_var2 = n_span;
var1 = expVariants; var2 = spanVariants;
figure('rend','painters', 'pos', [10 10 900 600]);
for surfvar = 1:n_surfvar
    surf(gca,var1, var2, squeeze(error(:,:,surfvar)));
    title(['Mean value of the mesh quality,' ' for polynomial degree ' num2str(surfvar-1)]);
    xlabel('Interpolation parameter a');
    ylabel('Support size r');
    zlabel('Mean value of mesh quality');
    set(gca,'XTick', var1);
    set(gca,'YTick', var2);
    view(0,90);
    caxis([0 1]);
    %set(gcf,'renderer','painters');
%     set(gca,'ZScale', 'log');
    colorbar
%     grid on
    saveas(gca, ['stddev_60_' 'deg' num2str(surfvar-1) '.png']);
    figure('rend','painters', 'pos', [10 10 900 600]);
end
%surfvar/n_surfvar*ones(n_var1,n_var2)'
% 
%     xlabel('Interpolation parameter a');
%     ylabel('Support size r');
%     zlabel('No. of degenerated elements');
%     titl='Parametric surface plot, polynomial degrees 0-4';
%     title(titl);
%     legend('n=0', 'n=1', 'n=2', 'n=3', 'n=4')
%     set(gca,'XTick', var1)
% %    set(gca,'XTickLabel',genxticklabl)
%     set(gca,'YTick', var2)
% %    set(gca,'YTickLabel',{'1st' '2nd' '3rd' '4th'})
%     grid on
%     view(110,30)
% %      set(gcf,'renderer','opengl')
%      set(gcf,'renderer','painters')
%    % print(figure(1),'-tiff','text.tiff','-opengl')
%    % hgexport(figure(1),'text.tiff',hgexport('factorystyle'),'format','tiff')
% %saveas(figure(1), 'parametric_final.png')
%     view(120,10)
%  %   saveas(figure(1),'par1.tiff')

 
%  figure
%  
%  plot(spanVariants, error(:,1,1))
%  xlabel('Support size r')
%  ylabel('No. of degenerated elements')
%  title('No. of degenerated elements for 0-degree polynomial')
%  grid on
 