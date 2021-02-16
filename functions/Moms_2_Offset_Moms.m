function t0m1m2n2 = Moms_2_Offset_Moms(a1b1a2b2)    
%     Converts a1 b1 a2 b2 to centred moments theta0, m1, m2, n2
    if size(a1b1a2b2,2) == 1
        a1b1a2b2 = a1b1a2b2';
    end

    a1 = a1b1a2b2(:, 1);
    b1 = a1b1a2b2(:, 2);
    a2 = a1b1a2b2(:, 3);
    b2 = a1b1a2b2(:, 4);
    
    t0 = wrapToPi( atan2(b1, a1));
    m1 = sqrt(a1.^2 + b1.^2);
    
    
    m2 = a2.*cos(2*t0) + b2.*sin(2*t0);
    n2 = -a2.*sin(2*t0) + b2.*cos(2*t0);
    
    t0m1m2n2 = [t0, m1, m2, n2];
end