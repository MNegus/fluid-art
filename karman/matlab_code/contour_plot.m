%%
addpath("src");

clear;
close all;

% Parent directory
parent_dir = "/home/michael/Documents/Generative art/karman_flow/data";


fig = figure(1);
hold on;


% %%
for interp = [0.01 : 0.01 : 0.99]
    % Reads in text file
    read_mat = dlmread(sprintf("%s/f_output_%g.txt", parent_dir, interp), ' ', 1);

    % Loads in vectors
    x_read = read_mat(:, 1);
    y_read = read_mat(:, 2);    
    z_read = read_mat(:, 3);

    % Finds the dimensions
    dimen = sqrt(length(x_read));
    x = reshape(x_read, dimen, dimen);
    y = reshape(y_read, dimen, dimen);
    z = reshape(z_read, dimen, dimen);

    % Plots contours
    contour(x, y, z, [1], 'color', 'black');
    
    % Plots reversed contours
%     contour(x, -y, z, [1], 'color', 'black');
    
end

%%
ang=0:0.01:2*pi; 
r = 0.125;
xp=r*cos(ang);
yp=r*sin(ang);
plot(xp, yp, 'color', 'black');


%%
pbaspect([1 1 1]);

hAxes = gca;
hAxes.XRuler.Axle.LineStyle = 'none'; 
hAxes.YRuler.Axle.LineStyle = 'none'; 
set(gca,'xtick',[]);
set(gca,'ytick',[])

% set(gca, 'visible', 'off');
% set(gca,'LooseInset',get(gca,'TightInset'));
% saveas(gcf, "contours.svg");
plot2svg('contours_alt.svg');



