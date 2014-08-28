% A script to produce figures for WAXS
close all;
clear;
global wavelength pixelSize sDist beamX beamZ
addpath(genpath('C:\Documents and Settings\Owner\My Documents\work\MATLAB_UserFunctions\Functions\Downloaded'));
addpath('C:\Documents and Settings\Owner\My Documents\work\MATLAB_UserFunctions\Functions\waxs_kiyo');
addpath(genpath('C:\Documents and Settings\Owner\My Documents\work\chess13a'));

addpath(genpath('~/WinC/MATLAB_UserFunctions/Functions/Downloaded'));
addpath('~/WinC/MATLAB_UserFunctions/Functions/waxs_kiyo');
addpath('~/WinC/MATLAB_UserFunctions/Functions');
addpath(genpath('~/WinE/chess13a'));

wavelength = 1.175;
pixelSize = 0.07113;
sDist = 220.6;
beamX = -241;
beamZ = 96; % same as X_cen = 1024 - 96 + 1 = 929;

%% ripple pahse
% D = 59.2
a = slurp('dmpc1_046_cz.tif', 'c');
b = slurp('dmpc1_047_cz.tif', 'c');
tmp = a - b;
tmp = rotateAround(tmp, 929, -245, -0.57, 'bicubic');
%tmp = rotateAround(tmp, beamZ, beamX, -0.57, 'bicubic');
tmp = 0.7 * tmp;
tmp = flipud(tmp);
waxs1 = transform_ccd2q(tmp, [0.4 2.1], [0 1.6], 0.0017, 0.0017, 0.2, beamX, beamZ);
waxs1.Int = int64(waxs1.Int);
fig1 = figure;
qshow(waxs1, [0 50]);
%saveTightFigure(fig1, 'dmpc1_046.pdf');

figure
[q, I] = sector_q(waxs1, [1 2], [5 15]);
dlmwrite('sector1.dat', [q I]);
figure
[q, I] = sector_q(waxs1, [1 2], [15 25]);
dlmwrite('sector2.dat', [q I]);
figure
[q, I] = sector_q(waxs1, [1 2], [25 35]);
dlmwrite('sector3.dat', [q I]);
figure
[q, I] = sector_q(waxs1, [1 2], [35 45]);
dlmwrite('sector4.dat', [q I]);
figure
[q, I] = sector_q(waxs1, [1 2], [45 55]);
dlmwrite('sector5.dat', [q I]);
figure
[q, I] = sector_q(waxs1, [1 2], [55 65]);
dlmwrite('sector6.dat', [q I]);


% D = 60.8
a = slurp('dmpc1_052_cz.tif', 'c');
b = slurp('dmpc1_053_cz.tif', 'c');
tmp = a - b;
a = slurp('dmpc1_060_cz.tif', 'c');
b = slurp('dmpc1_061_cz.tif', 'c');
tmp2 = a - b;
tmp = (tmp + tmp2) / 2;
tmp = rotateAround(tmp, 929, -241, -0.57, 'bicubic');
tmp = flipud(tmp);
waxs2 = transform_ccd2q(tmp, [0.4 2.1], [0 1.6], 0.0017, 0.0017, 0.2, beamX, beamZ);
waxs2.Int = int64(waxs2.Int);
fig2 = figure;
qshow(waxs2, [0 50]);
%saveas(fig2, 'dmpc1_052-060.pdf');
fig3 = figure;
qshow(waxs2, [0 100]);
axis([1.3 1.6 0 0.8])
%saveas(fig3, 'dmpc1_enlarge.pdf');

% difference betweene two D's
waxs3 = waxs1;
waxs3.Int = waxs1.Int - waxs2.Int;
figure
qshow(waxs3, [-5, 5]);
%saveTightFigure(gcf, '046_vs_052.pdf');

% k = 0;
% for i = 0.03:0.02:0.31
%   figure;
%   [qr, Int] = qrplot_q(waxs2, [i i+0.02], 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k');
%   xlabel(strcat('q_r (', char(197), '^{-1})'), 'interpreter', 'tex', ...
%          'FontName', 'Times New Roman', 'FontSize', 18);
% 
%   axis([1.38 1.58 0 110]);
%   avg = i + 0.01;
%   str = strcat(num2str(avg), {' '}, char(197), '^{-1}');
%   legend(str, 'FontName', 'Times New Roman', 'FontSize', 18);
%   filename = strcat('qrplot', num2str(k));
%   %saveTightFigure(gcf, strcat(filename, '.pdf'))
%   k = k + 1;
%   dlmwrite(strcat(filename, '.dat'), [qr Int])
% end

% gel phase data
a = slurp('dmpc1_107_cz.tif', 'c');
b = slurp('dmpc1_108_cz.tif', 'c');
tmp = a - b;
tmp = rotateAround(tmp, 929, -245, -0.57, 'bicubic');
tmp = flipud(tmp);
waxs4 = transform_ccd2q(tmp, [0.4 2.1], [0 1.6], 0.0017, 0.0017, 0.2, beamX, beamZ);
waxs4.Int = int64(waxs4.Int);
fig4 = figure;
qshow(waxs4, [0 50]);
axis([1.2 1.7 0 1])
% saveas(fig4, 'dmpc1_107.pdf');
 
% k = 0;
% for i = 0.05:0.02:0.07
%   figure;
%   qrplot_q(waxs4, [i i+0.02], 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k');
%   xlabel(strcat('q_r (', char(197), '^{-1})'), 'interpreter', 'tex', ...
%          'FontName', 'Times New Roman', 'FontSize', 18);
% 
%   axis([1.38 1.58 0 80]);
%   avg = i + 0.01;
%   str = strcat(num2str(avg), {' '}, char(197), '^{-1}');
%   legend(str, 'FontName', 'Times New Roman', 'FontSize', 18);
%   saveTightFigure(gcf, strcat('gel_qrplot', num2str(k), '.pdf'))
%   k = k + 1;
% end
