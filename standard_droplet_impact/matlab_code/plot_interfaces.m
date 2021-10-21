addpath("src");




%% Data directories and parameters
data_dir = "/home/michael/repos/fluid-art/standard_droplet_impact/level_11";
max_timestep = 1500;
interval = 75;
plate_tol = 1e-2;

%% Loops and plots
close all;
figure(1);
% hold on;
for k = 0 : 37 : 1500
   % Name of interface file
    filename ...
        = sprintf('%s/interface_%d.txt', data_dir, k);
    
    % Reads in the start and end points of the line segments, with x along
    % the horizontal axis and y along the vertical
    [start_points, end_points] = read_interface_points(filename, false);
    
    % Finds unique values of y in all the points
    all_points = [start_points; end_points];
    [~, uniq_idxs, ~] = uniquetol(all_points(:, 1), 1e-4);
    uniq_points = all_points(uniq_idxs, :);
    
    % Removes all points that are below the plate_tol
    uniq_points = uniq_points(uniq_points(:, 1) > plate_tol, :);
    
    % Sorts the resulting vector in decreasing y order
    [~, sorted_idxs] = sort(uniq_points(:, 1), 'descend');
    sorted_points = uniq_points(sorted_idxs, :);
    
    % Saves xs and ys
    ys = sorted_points(:, 1);
    xs = sorted_points(:, 2);
    
    previous_x = xs(1);
    q = 2;
    while (q <= length(xs))
        if abs(xs(q) - previous_x) > 0.1
            xs(q) = [];
            ys(q) = [];
        else
            previous_x = xs(q);
        end
        q = q + 1;
    end
    
    % Full xs and ys
    full_xs = [-flip(xs); 0; xs];
    full_ys = [flip(ys); max(ys); ys];
    
    % Plots
    plot(full_xs, full_ys, 'color', 'black');
    
    xlim([-4, 4]);
    ylim([0, 4]);
    pbaspect([2 1 1]);
    
    drawnow;
    pause(0.1);
    
    % Outputs as text file
    full_data = [full_xs'; full_ys'];
    fileID = fopen(sprintf("interface_data/interface_%d.txt", k),'w');
    formatSpec = '%.4f, %.4f\n';
    fprintf(fileID, formatSpec, full_data);
    fclose(fileID);
    
    
end




%% Saves plot
lim = max(max(xs), max(ys));
xlim([-lim, lim]);
ylim([0, lim]);
pbaspect([2 1 1]);

hAxes = gca;
hAxes.XRuler.Axle.LineStyle = 'none'; 
hAxes.YRuler.Axle.LineStyle = 'none'; 
set(gca,'xtick',[]);
set(gca,'ytick',[])
set(gca, 'color', 'none');

plot2svg('droplet_interfaces.svg');