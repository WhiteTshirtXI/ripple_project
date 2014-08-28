% A script to produce a figure for LAXS
close all;
clear;

addpath('C:\Documents and Settings\Owner\My Documents\work\MATLAB_UserFunctions\Functions\waxs_kiyo');
addpath(genpath('C:\Documents and Settings\Owner\My Documents\work\MATLAB_UserFunctions\Functions\Downloaded'));
addpath('C:\Documents and Settings\Owner\My Documents\work\chess13a\dmpclr');

addpath('~/WinC/MATLAB_UserFunctions/Functions/waxs_kiyo');
addpath('~/WinC/MATLAB_UserFunctions/Functions');
addpath(genpath('~/WinC/MATLAB_UserFunctions/Functions/Downloaded'));
addpath('~/WinE/chess13a/dmpclr');


a = slurp('ripple_083_cz.tif', 'c');
pedestal = 100;
a = a - pedestal;
a = imrotate(a, -0.46, 'bicubic');
fig1 = figure;
%tview(a, [0 100]);
imagesc(a, [0 100]);
colormap(gray)
%colorbar
axis image
axis([322 722 0 1024])
qr = 0:0.1:1;
qz = 0:0.1:1;
wavelength = 1.175;
beamX = 522;
beamZ = 948;
sDist = 359.3;
pixelSize = 0.07113;
add_qr_qz_labels(qr, qz, wavelength, beamX, beamZ, sDist, pixelSize);
color_log([0 1000])
%saveTightFigure(fig1, 'ripple083.pdf');
saveas(fig1, 'ripple083.pdf')

b = slurp('ripple_085_cz.tif', 'c');
pedestal = 112;
b = b - pedestal;
b = imrotate(b, -0.46, 'bicubic');
fig2 = figure;
imagesc(b, [0 300]);
colormap(gray)
%colorbar
axis image
axis([222 822 0 1024])
qr = 0:0.1:1;
qz = 0:0.1:1;
wavelength = 1.175;
beamX = 522;
beamZ = 948;
sDist = 359.3;
pixelSize = 0.07113;
add_qr_qz_labels(qr, qz, wavelength, beamX, beamZ, sDist, pixelSize);
color_log([0 2000])
%saveTightFigure(fig2, 'ripple085.pdf');
saveas(fig2, 'ripple085.pdf')
