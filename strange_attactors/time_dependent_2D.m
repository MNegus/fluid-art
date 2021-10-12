%% Parameters
found_initial_solution = false;
NMAX = 1000;

% range = [-1, 0, 1];
% range = -0.6 : 0.05 : 0.6;
range = -0.6 : 0.1 : 0.6;
no_params = 18;

t = -2;
tmax = 2;
dt = 1e-2;

%% Finds initial configuration
while found_initial_solution == false
    
    % Randomise parameters
    as = zeros(no_params, 1);
    for q = 1 : no_params
        pos = randi(length(range));
        as(q) = range(pos);
    end
    as
    
    % Performs the first timestep 
    [ORIG_X_POINTS, ORIG_Y_POINTS, stopped] = initialise(as, NMAX, t);
    
    % Plots if found solution
    if (stopped == false)
        found_initial_solution = true;
    end
end



%% Plots
t = -2;
X_POINTS = ORIG_X_POINTS;
Y_POINTS = ORIG_Y_POINTS;

% Finds ranges
xmin = min(X_POINTS);
xmax = max(X_POINTS);
ymin = min(X_POINTS);
ymax = max(X_POINTS);

close all;
close(figure(1));
figure(1);
% hold on;
g = arrayfun(@(x) animatedline(), 1 : NMAX)

all_x_points = [];
all_y_points = [];

hWaitbar = waitbar(0, 'Iteration 1', 'Name', 'Solving problem','CreateCancelBtn','delete(gcbf)');
while t < tmax
    % Update time
    t = t + dt
    
    % Find new points
    xdiff = 0;
    for n = 1 : NMAX
       [XNEW, YNEW] = update(ORIG_X_POINTS(n), ORIG_Y_POINTS(n), t, as);
       xdiff = max(xdiff, abs(XNEW - X_POINTS(n)));
       X_POINTS(n) = XNEW;
       Y_POINTS(n) = YNEW;
    end
    
    %% Scatter method
    % Filter points
    filtered_x_points = X_POINTS(abs(X_POINTS) < 1e3);
    filtered_y_points = Y_POINTS(abs(X_POINTS) < 1e3);
    
    % Plotted points
    all_x_points = [all_x_points; filtered_x_points];
    all_y_points = [all_y_points; filtered_y_points];
    
%     scatter(all_x_points, all_y_points, 0.1);
    
    %% Animated lines method
    for n = 1 : 50
        x = X_POINTS(n);
        y = Y_POINTS(n);
        if (abs(x) >= 1e3 && abs(y) >= 1e3)
            clearpoints(g(n)); 
        elseif (t > -2)
            addpoints(g(n), x, y);
        end
    end
    
    drawnow;
    if ~ishandle(hWaitbar)
        % Stop the if cancel button was pressed
        disp('Stopped by user');
        break;
    end
    
%     xlim([min(0.5 * xmin, 2 * xmin), max(0.5 * xmax, 2 * xmax)]);
%     ylim([min(0.5 * ymin, 2 * ymin), max(0.5 * ymax, 2 * ymax)]);
    pause(0.001);
    
    
    
    
    xdiff
    
end

%% Function declarations
function [X_POINTS, Y_POINTS, stopped] = initialise(as, NMAX, t)

    % Initialise points
    X_POINTS = zeros(NMAX, 1);
    Y_POINTS = zeros(NMAX, 1);
    
    X = t;
    Y = t;

    % Initialise Variables
    LSUM = 0;
    NL = 0;
    
    XE = X + 0.000001; YE = Y;

    stopped = false;

    for N = 1 : NMAX
        N
        X_POINTS(N) = X;
        Y_POINTS(N) = Y;
        
        % Update X and Y
        [XNEW, YNEW] = update(X, Y, t, as);

        % Stops if unbounded points
        if (abs(XNEW) + abs(YNEW) > 1e6) 
            disp("Unbounded");
            stopped = true;
            break;
        end
        
        % Stops if fixed point
        if (abs(XNEW - X) + abs(YNEW - Y) < 1e-6) 
            disp("Fixed point");
            stopped = true;
            break;
        end

        % Updates X and Y
        X = XNEW;
        Y = YNEW;
        
        %% Calculate Lyapunov exponent
        XSAVE = XNEW; YSAVE = YNEW; 
        X = XE; Y = YE;
        
        % Find new values of XNEW and YNEW
        [XNEW, YNEW] = update(X, Y, t, as);
        
        % Update DLX and DLY
        DLX = XNEW - XSAVE;
        DLY = YNEW - YSAVE;
        DL2 = DLX^2 + DLY^2;
        
        if (DL2 > 0)
            DF = 1e12 * DL2;
            RS = 1 / sqrt(DF);
            XE = XSAVE + RS * (XNEW - XSAVE);
            YE = YSAVE + RS * (YNEW - YSAVE);
            X = XSAVE;
            Y = YSAVE;
            
            LSUM = LSUM + log2(DF);
            NL = NL + 1;
            L = LSUM / NL;
            
            % Stops if limit cycle
            if ((N > 100) && (L < 0.005))
               disp("Limit cycle");
               stopped = true;
               break;
            end
        end
        

    end
end

function [XNEW, YNEW] = update(X, Y, t, as)
    % Determine new values of X and Y
    XNEW = as(1) * X^2 + as(2) * Y^2 + as(3) * t^2 + as(4) * X * Y ...
        + as(5) * X * t + as(6) * Y * t + as(7) * X + as(8) * Y + as(9) * t;
    YNEW = as(10) * X^2 + as(11) * Y^2 + as(12) * t^2 + as(13) * X * Y ...
        + as(14) * X * t + as(15) * Y * t + as(16) * X + as(17) * Y + as(18) * t;
end