% A script to produce figures for WAXS
close all;
clear;
global wavelength pixelSize sDist beamX beamZ
addpath(genpath('C:\Documents and Settings\Owner\My Documents\work\MATLAB_UserFunctions\Functions\Downloaded'));
addpath('C:\Documents and Settings\Owner\My Documents\work\MATLAB_UserFunctions\Functions\waxs_kiyo');
addpath(genpath('C:\Documents and Settings\Owner\My Documents\work\chess11'));

addpath(genpath('~/WinC/MATLAB_UserFunctions/Functions/Downloaded'));
addpath('~/WinC/MATLAB_UserFunctions/Functions/waxs_kiyo');
addpath('~/WinC/MATLAB_UserFunctions/Functions');
addpath(genpath('~/WinE/chess11'));

wavelength = 1.176;
pixelSize = 0.07113;
sDist = 169.8;
beamX = 32;
beamZ = 110; % same as X_cen = 1024 - 110 + 1 = 915;

%% ripple phase data
a = slurp('ripple_060_cz.tif', 'c');
b = slurp('ripple_061_cz.tif', 'c');
tmp = a - b;
tmp = rotateAround(tmp, 915, 32, -0.97, 'bicubic');
tmp = flipud(tmp);
waxs = transform_ccd2q(tmp, [0.0 2.01], [0 2.01], 0.0022, 0.0022, 0.2, beamX, beamZ);
waxs.Int = int64(waxs.Int);
figure
qshow(waxs, [0 500]);
saveTightFigure(gcf, 'ripple_waxs_low-res.pdf');
%axis([1 1.8 0 1])

figure
[qr, Int] = qrplot_q(waxs, [0.11 0.13]);
dlmwrite('ripple_waxs_low-res_weak.dat', [qr Int]);

figure
[qr, Int] = qrplot_q(waxs, [0.19 0.21]);
dlmwrite('ripple_waxs_low-res_strong.dat', [qr Int]);

%% ripple phase data
% a = slurp('ripple_060_cz.tif', 'c');
% a = flipud(a);
% fig1 = figure;
% imagesc(a, [0 1000])
% colormap(gray)
% set(gca, 'TickDir', 'out', 'FontSize', 16, 'FontName', 'Times New Roman')
% axis image
% axis xy
% saveTightFigure(fig1, 'ripple_060.pdf')
% 
% b = slurp('ripple_061_cz.tif', 'c');
% b = flipud(b);
% fig2 = figure;
% imagesc(b, [0 1000])
% colormap(gray)
% set(gca, 'TickDir', 'out', 'FontSize', 16, 'FontName', 'Times New Roman')
% axis image
% axis xy
% saveTightFigure(fig2, 'ripple_061.pdf')
% 
% tmp = a - b;
% fig3 = figure;
% imagesc(tmp, [0 300])
% colormap(gray)
% set(gca, 'TickDir', 'out', 'FontSize', 16, 'FontName', 'Times New Roman')
% axis image
% axis xy
% saveTightFigure(fig3, 'ripple_060_ccd.pdf')
% 
% a = slurp('ripple_060_cz.tif', 'c');
% b = slurp('ripple_061_cz.tif', 'c');
% tmp = a - b;
% tmp = rotateAround(tmp, 915, 32, -0.97, 'bicubic');
% waxs1 = transform_ccd2q(tmp, [0 2.2], [0 2], 0.0022, 0.0022, 0.2);
% waxs1.Int = int64(waxs1.Int);
% fig4 = figure;
% qshow(waxs1, [0 300]);
% saveTightFigure(fig4, 'ripple_060_q.pdf');
% %axis([1 1.8 0 1])
% %saveas(fig4, 'ripple_060.pdf');
% %saveTightFigure(fig4, 'ripple_060.pdf')

%% fluid phase data
% a = slurp('fluid_036_cz.tif', 'c');
% a = flipud(a);
% figure
% imagesc(a, [0 3000])
% colormap(gray)
% set(gca, 'TickDir', 'out', 'FontSize', 16, 'FontName', 'Times New Roman')
% xlabel('p_x', 'interpreter', 'tex');
% ylabel('p_z', 'interpreter', 'tex');
% axis image
% axis xy
% saveTightFigure(gcf, 'fluid_036.pdf')
% 
% b = slurp('fluid_039_cz.tif', 'c');
% b = flipud(b);
% figure
% imagesc(b, [0 3000])
% colormap(gray)
% set(gca, 'TickDir', 'out', 'FontSize', 16, 'FontName', 'Times New Roman')
% xlabel('p_x', 'interpreter', 'tex');
% ylabel('p_z', 'interpreter', 'tex');
% axis image
% axis xy
% saveTightFigure(gcf, 'fluid_039.pdf')
% 
% tmp = a - b;
% figure
% imagesc(tmp, [0 1000])
% colormap(gray)
% set(gca, 'TickDir', 'out', 'FontSize', 16, 'FontName', 'Times New Roman')
% xlabel('p_x', 'interpreter', 'tex');
% ylabel('p_z', 'interpreter', 'tex');
% axis image
% axis xy
% saveTightFigure(gcf, 'fluid_036_ccd.pdf')
% 
% tmp = flipud(tmp);
% tmp = rotateAround(tmp, 915, 32, -0.97, 'bicubic');
% waxs2 = transform_ccd2q(tmp, [0 2.2], [0 2], 0.0022, 0.0022, 0.5);
% waxs2.Int = int64(waxs2.Int);
% figure
% qshow(waxs2, [0 1000]);
% saveTightFigure(gcf, 'fluid_036_q.pdf')
% 
% k = 0;
% for i = 0.09:0.1:0.6
%   figure 
%   [qr, Int] = qrplot_q(waxs1, [i i+0.02]);
%   avg = i + 0.01;
%   filename = strcat('qr_ripple_', num2str(k), '.dat');
%   dlmwrite(filename, [qr Int])
%   k = k + 1;
% end
% 
% k = 0;
% for i = 0.09:0.1:0.6
%   figure 
%   [qr, Int] = qrplot_q(waxs2, [i i+0.02]);
%   avg = i + 0.01;
%   filename = strcat('qr_fluid_', num2str(k), '.dat');
%   dlmwrite(filename, [qr Int])
%   k = k + 1;
% end


% k = 0;
% for i = 0.09:0.1:0.6
%   figure 
%   [qr, Int] = qrplot_q(waxs1, [i i+0.02], 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k');
%   hold on
%   qrplot_q(waxs2, [i i+0.02], 'LineStyle', 'none', 'Marker', 'o', 'Color', 'r');
%   hold off
%   set(gca,'XMinorTick','on','YMinorTick','on')
%   xlabel(strcat('q_r (', char(197), '^{-1})'), 'interpreter', 'tex', ...
%          'FontName', 'Times New Roman', 'FontSize', 18);
%   axis([1 1.8 0 1300]);
%   avg = i + 0.01;
%   str = strcat(num2str(avg), {' '}, char(197), '^{-1}');
%   legend(str, 'FontName', 'Times New Roman', 'FontSize', 18, 'Location', 'NorthWest');
%   saveTightFigure(gcf, strcat('qr_fluid&ripple', num2str(k), '.pdf'))
%   %saveas(gcf, strcat('qr_fluid&ripple', num2str(k), '.pdf'))
%   k = k + 1;
% end

% figure;
% tmp = waxs1;
% tmp.Int = waxs1.Int-waxs2.Int;
% qrplot_q(tmp, [0.19 0.21], 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k');
% xlabel(strcat('q_r (', char(197), '^{-1})'), 'interpreter', 'tex', ...
%        'FontName', 'Times New Roman', 'FontSize', 18);
% axis([1 1.8 0 1300]);
% avg = 0.20;
% str = strcat(num2str(avg), {' '}, char(197), '^{-1}');
% legend(str, 'FontName', 'Times New Roman', 'FontSize', 18, 'Location', 'NorthWest');
% saveTightFigure(gcf, strcat('qr_ripple-fluid.pdf'))