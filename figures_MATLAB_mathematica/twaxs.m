% A script to produce figures for WAXS
close all;
clear;
global wavelength pixelSize sDist beamX beamZ
addpath(genpath('C:\Documents and Settings\Owner\My Documents\work\MATLAB_UserFunctions\Functions\Downloaded'));
addpath('C:\Documents and Settings\Owner\My Documents\work\MATLAB_UserFunctions\Functions\waxs_kiyo');
addpath(genpath('C:\Documents and Settings\Owner\My Documents\work\chess11\tripple'));

addpath(genpath('~/WinC/MATLAB_UserFunctions/Functions/Downloaded'));
addpath('~/WinC/MATLAB_UserFunctions/Functions/waxs_kiyo');
addpath('~/WinC/MATLAB_UserFunctions/Functions');
addpath(genpath('~/WinE/chess11'));

wavelength = 1.176;
pixelSize = 0.07113;
%sDist = 158.6;
sDist = 161.8;

%% Ripple phase 
a = slurp('waxs_003_cz.tif', 'c');
b = slurp('bkgd_002_cz.tif', 'c');
a = a - 110;
b = b - 100;
b = b * 6.9;
tmp = a - b;
tmp = flipud(tmp);
beamX = 29;
beamZ = 510;
%tmp = rotateAround(tmp, 929, -245, -0.57, 'bicubic');
twaxs1 = transform_ccd2q(tmp, [0.6 1.8], [-0.6 1], 0.0024, 0.0024, -45, beamX, beamZ);
twaxs1.Int = int64(twaxs1.Int);
figure
qshow(twaxs1, [0 2000]);
axis([1.3 1.7 -0.6 0.6])
%saveas(fig1, 'waxs_003.pdf');

figure
[qz, Int] = qzplot_q(twaxs1, [1.465 1.481]);
dlmwrite('twaxs_ripple_qz_weak.dat', [qz Int]);

figure
[qz, Int] = qzplot_q(twaxs1, [1.483 1.5]);
dlmwrite('twaxs_ripple_qz_strong.dat', [qz Int]);

figure
[qz, Int] = qzplot_q(twaxs1, [1.465 1.5]);
dlmwrite('twaxs_ripple_qz_both.dat', [qz Int]);

figure
[qr, Int] = qrplot_q(twaxs1, [0.11 0.13]);
dlmwrite('twaxs_ripple_qr_weak.dat', [qr Int]);

figure
[qr, Int] = qrplot_q(twaxs1, [0.19 0.21]);
dlmwrite('twaxs_ripple_qr_strong.dat', [qr Int]);

% k = 0;
% for i = -0.05:0.02:0.21
%   figure;
%   qrplot_q(waxs1, [i i+0.02], 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k');
%   xlabel(strcat('q_r (', char(197), '^{-1})'), 'interpreter', 'tex', ...
%          'FontName', 'Times New Roman', 'FontSize', 18);
% 
%   axis([1.38 1.58 0 2000]);
%   avg = i + 0.01;
%   str = strcat(num2str(avg), {' '}, char(197), '^{-1}');
%   legend(str, 'FontName', 'Times New Roman', 'FontSize', 18);
%   %saveTightFigure(gcf, strcat('twaxs_qrplot', num2str(k), '.pdf'))
%   k = k + 1;
% end
% % do qz=0 again
% figure;
% qrplot_q(waxs1, [-0.01 0.01], 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k');
% xlabel(strcat('q_r (', char(197), '^{-1})'), 'interpreter', 'tex', ...
%        'FontName', 'Times New Roman', 'FontSize', 18);
% axis([1.38 1.58 0 2000]);
% str = strcat('0', {' '}, char(197), '^{-1}');
% legend(str, 'FontName', 'Times New Roman', 'FontSize', 18);
% saveTightFigure(gcf, strcat('twaxs_qrplot2.pdf'))
  

%% Gel phase
a = slurp('gel_020_cz.tif', 'c');
b = slurp('bkgd_005_cz.tif', 'c');
a = flipud(a);
b = flipud(b);
tmp = a - b;
beamX = 30;
beamZ = 431;
tmp = rotateAround(tmp, beamZ, beamX, -0.57, 'bicubic');
figure
imagesc(tmp, [0 100]);
colormap(gray)
axis image
axis xy
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
xlabel('p_x', 'interpreter', 'tex');
ylabel('p_z', 'interpreter', 'tex');
set(gca,'tickdir','out');
saveTightFigure(gcf, 'twaxs_gel_ccd.pdf');
twaxs2 = transform_ccd2q(tmp, [0 2.4], [-1 1], 0.0024, 0.0024, -45, beamX, beamZ);
twaxs2.Int = int64(twaxs2.Int);
figure
qshow(twaxs2, [0 100]);
%axis([1.3 1.7 -0.6 0.6])
saveTightFigure(gcf, 'twaxs_gel_q.pdf');

figure
[qz, Int] = qzplot_q(twaxs2, [1.47 1.51]);
dlmwrite('twaxs_gel_qz_20.dat', [qz Int]);

figure
[qr, Int] = qrplot_q(twaxs2, [-0.05 0.05]);
dlmwrite('twaxs_gel_qr_20.dat', [qr Int]);