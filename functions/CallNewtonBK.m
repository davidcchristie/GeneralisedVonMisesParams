    function x = CallNewton(TargMoms, x0, tolerance, trunc)
    
    if nargin<4
        trunc = 50;
    end
    
    FinalTarget = TargMoms;
    needtosolve = true;
    x0Orig = x0;

    y0 = MomsAndJacobianFlexi(x0, trunc);  
    yTarg = FinalTarget;   
    isFinalTarget = true;
    movingPointCount = 0; 
    % If Newton's method fails, use an intermediate target point.  
    % movingPointCount ensures this is only attempted a finite number of
    % times


    while (needtosolve && movingPointCount<200)
        x = Newton(yTarg,x0, tolerance);
        needtosolve = any(isnan(x)) | ~isFinalTarget;
        if isnan(x) % i.e. Newton's method failed
            yTarg = (y0+yTarg)/2; % Bring target closer to known values
            isFinalTarget = false;
            movingPointCount = movingPointCount + 1;
        elseif ~isFinalTarget % Newton succeeded for an intermediate target
                x0 = x; % Use new known values to define new starting point
                y0 = MomsAndJacobianFlexi(x0); 
                
                yTarg = FinalTarget; % Next time, attempt to solve for final target
                isFinalTarget = true;
        end

    end
    
    if any(isnan(x))
        x = x0Orig;
        disp(['Orig x0 returned after ', num2str(movingPointCount), ' iterations']);
    end
        
    end