function a1b1a2b2 = Offset_Moms_2_Moms(t0m1m2n2)    
%     Converts a1 b1 a2 b2 to centred moments theta0, m1, m2, n2
    if size(t0m1m2n2, 2) == 1
        t0m1m2n2 = t0m1m2n2';
    end

    if size(t0m1m2n2,2)>=4
            t0 = t0m1m2n2(:,1);
        else
            t0 = zeros(size(t0m1m2n2,1), 1);
    end

    m1 = t0m1m2n2(:,end-2);
    m2 = t0m1m2n2(:,end-1);
    n2 = t0m1m2n2(:,end);
    
    a1 = m1.*cos(t0);
    b1 = m1.*sin(t0);
    a2 = m2.*cos(2*t0) - n2.*sin(2*t0);
    b2 = m2.* sin(2*t0) + n2.*cos(2*t0);

    
    a1b1a2b2 = [a1, b1, a2, b2];
end