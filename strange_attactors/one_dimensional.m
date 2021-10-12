
%% Parameters
found_solution = false;
NMAX = 11000;
X_VALS = zeros(NMAX, 1);

while found_solution == false
    as = 9 * rand(3, 1) - 4.5;
    
    LSUM = 0;
    NL = 0;

    X = 0;
    

    stopped = false;

    for N = 1 : NMAX
        N
        X_VALS(N) = X;
        XNEW = as(1) + as(2) * X + as(3) * X^2;

        % Return if determined to be unbounded
        if (abs(XNEW) > 1e6) 
            disp("Unbounded");
            stopped = true;
            break;
        end

        % Calculate Lyapunov exponent
        DF = abs(as(2) + 2 * as(3) * X);
        if (DF > 0)
           LSUM = LSUM + log2(DF);
           NL = NL + 1;
           L = LSUM / NL;
        end

        % Stops if fixed point
        if (abs(XNEW - X) < 1e-6) 
            disp("Fixed point");
            stopped = true;
            break;
        end

        % Stops if limit cycle
        if ((N > 100) && (L < 0.005))
           disp("Limit cycle");
           stopped = true;
           break;
        end

        X = XNEW

    end
    
    % Plots if found solution
    if (stopped == false)
        found_solution = true;
    end
end
%%
close(figure(1));
figure(1);
interval = 200;

gap = 5;
scatter(X_VALS(1 : end - gap), X_VALS(1 + gap: end) );
