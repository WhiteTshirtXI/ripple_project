% A script to produce figures
close all;
clear;

addpath('C:\Documents and Settings\Owner\My Documents\work\MATLAB_UserFunctions\Functions\waxs_kiyo');
addpath(genpath('C:\Documents and Settings\Owner\My Documents\work\MATLAB_UserFunctions\Functions\Downloaded'));
addpath('C:\Documents and Settings\Owner\My Documents\work\chess13a\lastday');

addpath('~/WinC/MATLAB_UserFunctions/Functions/waxs_kiyo');
addpath('~/WinC/MATLAB_UserFunctions/Functions');
addpath(genpath('~/WinC/MATLAB_UserFunctions/Functions/Downloaded'));
addpath('~/WinE/chess13a/lastday');

a = slurp('olddopc_044_c.tif', 'c');
pedestal = 103;
a = a - pedestal;
a = flipud(a);
fig1 = figure;
imagesc(a, [0 100]);
colormap(gray)
axis image
axis xy
axis([422 622 1 400]);
px = 100:100:1000;
pz = 100:100:1000;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
set(gca, 'xtick', px);
set(gca, 'ytick', pz);
set(gca, 'xticklabel', px);
set(gca, 'yticklabel', pz);
set(gca, 'TickDir', 'out');
xlabel('p_x', 'interpreter', 'tex', 'FontName', 'Times New Roman');
ylabel('p_z', 'interpreter', 'tex', 'FontName', 'Times New Roman');
saveTightFigure(fig1, 'olddopc044.pdf');

b = slurp('olddopc_045_c.tif', 'c');
pedestal = 102;
b = b - pedestal;
b = flipud(b);
fig2 = figure;
imagesc(b, [0 100]);
colormap(gray)
axis image
axis xy
axis([422 622 1 400]);
px = 100:100:1000;
pz = 100:100:1000;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
set(gca, 'xtick', px);
set(gca, 'ytick', pz);
set(gca, 'xticklabel', px);
set(gca, 'yticklabel', pz);
set(gca, 'TickDir', 'out');
xlabel('p_x', 'interpreter', 'tex', 'FontName', 'Times New Roman');
ylabel('p_z', 'interpreter', 'tex', 'FontName', 'Times New Roman');
saveTightFigure(fig2, 'olddopc045.pdf');

figure;
qzplot(a, [529 533], 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', ...
       'MarkerFaceColor', 'k', 'MarkerSize', 4);
hold on
qzplot(b, [529 533], 'LineStyle', 'none', 'Marker', 'o', 'Color', 'r', ...
       'MarkerFaceColor', 'r', 'MarkerSize', 4);
hold off
axis([0 400 -100 4000])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('p_z', 'interpreter', 'tex', 'FontName', 'Times New Roman');
saveTightFigure(gcf, 'attenuator1.pdf')

b = 6.9 * b;
figure;
qzplot(a, [529 533], 'LineStyle', 'none', 'Marker', 'o', 'Color', 'k', ...
       'MarkerFaceColor', 'k', 'MarkerSize', 4);
hold on
qzplot(b, [529 533], 'LineStyle', 'none', 'Marker', 'o', 'Color', 'r', ...
       'MarkerFaceColor', 'r', 'MarkerSize', 4);
hold off
axis([0 400 -100 4000])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('p_z', 'interpreter', 'tex', 'FontName', 'Times New Roman');
saveTightFigure(gcf, 'attenuator2.pdf')
