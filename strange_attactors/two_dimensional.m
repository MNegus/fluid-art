%% Parameters
found_solution = false;
NMAX = 1e4;
X_VALS = zeros(NMAX, 1);
Y_VALS = zeros(NMAX, 1);

range = -1.2 : 0.1 : 1.2;

while found_solution == false
    
    as = zeros(13, 1);
    for q = 1 : 13
        pos = randi(length(range));
        as(q) = range(pos);
    end
    
    LSUM = 0;
    NL = 0;

    X = 0.5;
    Y = 0.5;
    
    XE = X + 0.000001; YE = Y;

    stopped = false;

    for N = 1 : NMAX
        N
        X_VALS(N) = X;
        Y_VALS(N) = Y;
        
        % Update X and Y
        [XNEW, YNEW] = update(X, Y, as);

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
%         N = N - 1;
        
        % Find new values of XNEW and YNEW
        [XNEW, YNEW] = update(X, Y, as);
        
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
        end
        
        % Stops if limit cycle
        if ((N > 100) && (L < 0.005))
           disp("Limit cycle");
           stopped = true;
           break;
        end

    end
    
    % Plots if found solution
    if (stopped == false)
        found_solution = true;
    end
end
%%
close(figure(1));
figure(1);
hold on;
scatter(X_VALS, Y_VALS, 1);



function [XNEW, YNEW] = update(X, Y, as)
    % Determine new values of X and Y
    XNEW = as(1) + as(2) * X + as(3) * X^2 + as(4) * X * Y + as(5) * Y + as(6) * Y^2;
    YNEW = as(7) + as(8) * X + as(9) * X^2 + as(10) * X * Y + as(11) * Y + as(12) * Y^2;
end