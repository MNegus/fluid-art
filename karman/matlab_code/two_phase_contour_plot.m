%%
addpath("src");

clear;


% Parent directory
parent_dir = "/home/michael/Documents/Generative art/karman_flow/two_phase_data";

TIMESTEP = 10;

%%
% Reads in data
read_mat = dlmread(sprintf("%s/f_output_%d.txt", parent_dir, TIMESTEP), ' ', 1);

% Loads in vectors
x_read = read_mat(:, 1);
y_read = read_mat(:, 2);    
upper_read = read_mat(:, 3);
lower_read = read_mat(:, 4);

%%
% % Finds the dimensions
dimen = sqrt(length(x_read))
x = reshape(x_read, dimen, dimen);
y = reshape(y_read, dimen, dimen);
upper = reshape(upper_read, dimen, dimen);
lower = reshape(lower_read, dimen, dimen);

%% (Optional) Filters out in x direction
% max_x = 3.5;
% max_idx = 2809^2;
% x_filter = x_read(1 : max_idx);
% y_filter = y_read(1 : max_idx);
% upper_filter = upper_read(1 : max_idx);
% lower_filter = lower_read(1 : max_idx);
% 
% % Finds the dimensions
% length(x_filter)
% length(y_filter)
% length(upper_filter)
% length(lower_filter)
% dimen = sqrt(length(x_filter))
% x = reshape(x_filter, dimen, []);
% y = reshape(y_filter, dimen, dimen);
% upper = reshape(upper_filter, dimen, dimen);
% lower = reshape(lower_filter, dimen, dimen);
% 

%% Plots the contours
close all;
fig = figure(1);
hold on;

% Plots contours
contour(x, y, upper, [2], 'color', 'black');
contour(x, y, lower, [2], 'color', 'black');

% surf(x, y, upper, 'edgecolor', 'none')

ang=0:0.01:2*pi; 
r = 0.125/2;
xp=r*cos(ang);
yp=r*sin(ang);
plot(xp, yp, 'color', 'black');


%%
% pbaspect([1 1 1]);


% xlim([-0.5, 3.5]);
% ylim([-1, 1]);
% pbaspect([1 (2 / 4) 1]);

hAxes = gca;
hAxes.XRuler.Axle.LineStyle = 'none'; 
hAxes.YRuler.Axle.LineStyle = 'none'; 
set(gca,'xtick',[]);
set(gca,'ytick',[])

% set(gca, 'visible', 'off');
% % set(gca,'LooseInset',get(gca,'TightInset'));
% % saveas(gcf, "contours.svg");


plot2svg('two_phase_contours.svg');



