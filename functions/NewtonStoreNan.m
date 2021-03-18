function [x, bestx, tolreached] =   NewtonStore(yTarg, x, tolerance, trunc, maxiters, maxX)  
% Multivariate Newtonian solver.  Function and Jacobian are evaluated using
% "MomsAndJacobianFlexi".
if nargin < 6
    maxX= 25;
end

if nargin < 5
    maxiters = 250;
end

if nargin < 4
    trunc = 40;
end

if nargin < 3
    tolerance = 1e-3;
end
    % Initialise: ensure correct dimensions, and create placeholder for f
    yTarg = yTarg(:);     
    [bestx, x] = deal(x(:));       
    y0 = MomsAndJacobianFlexi(x,trunc);
    f = y0(:)-yTarg;
    [bestNormF, normF] = deal(norm(f));
%     tic
    count = 0;  %%%%%
     while normF > tolerance & count <= maxiters
        [y,J] = MomsAndJacobianFlexi(x,trunc);
        f = y(:)-yTarg;
        
        normF = norm(f);
        if normF < bestNormF
            bestNormF = normF;
            bestx = x;          
        end
        if rcond(J) > 1e-14 % i.e. J is invertible
            x = x - J\f;
            x(1) = abs(x(1)); 
            x(2) = abs(x(2));
            if x(1) > maxX 
                x(1) = NaN;
            end
            if x(2) > maxX
                x(2) = NaN;
            end
        else 
            x = NaN*x; 
            % This violates the "while" loop condition, so iteration stops
            % and NaN returned as result.
        end
        count = count + 1;

     end
     tolreached = normF <= tolerance;
end
     