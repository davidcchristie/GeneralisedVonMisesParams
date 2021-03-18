function fourierOK = KrogstadTest(a1b1a2b2)

% Are the Fourier coefficients suitable?  The inequalities below are
% necessary and sufficient for the supplied coefficients to be Fourier
% coefficients of a non-negative distribution (which can be modelled by the
% Generalised von Mises), a consequence of Bochner's theorem.
% If not, then the iteration procedure for GvM parameters will fail as no
% suitable GvM parameters will reproduce a distribution with the required
% moment set. 

% Refs:
% Barstow, Stephen F., et al. Measuring and analysing the directional 
% spectrum of ocean waves. COST Office, 2005.

% Kuik, A. J., G. Ph Van Vledder, and L. H. Holthuijsen. 
% "A method for the routine analysis of pitch-and-roll buoy wave data." 
% Journal of physical oceanography 18.7 (1988): 1020-1034.
    a1 = a1b1a2b2(:,1);
    b1 = a1b1a2b2(:,2);
    a2 = a1b1a2b2(:,3);
    b2 = a1b1a2b2(:,4);
    
    c1 = a1 + 1i * b1;
    c2 = a2 + 1i * b2;
%     c1 = a1b1a2b2(:,1) + 1i*a1b1a2b2(:,2);
%     c2 = a1b1a2b2(:,3) + 1i*a1b1a2b2(:,4);
    t0m1m2n2 = Moms_2_Offset_Moms(a1b1a2b2);
    m1 = t0m1m2n2(:,2);
    m2 = t0m1m2n2(:,3);
    n2 = t0m1m2n2(:,4);
    
    ineq1 = abs(c1) > 0;
    ineq2 = abs(c1.^2-c2) <= 1- abs(c1).^2;
    

    
%     ineq3 = 1-a1b1a2b2(:,3).^2 - a1b1a2b2(:,4).^2 - 2*m1.^2.*(1-m1.^2) >= 0;        
    ineq3 = 1-m2.^2 - n2.^2 - 2*m1.^2.*(1-m2) >= 0; 
    % Mardia (1972) 3.3.4 (Looks like there's a typo in Kuik)
    ineq4 = abs(c1) <= 1;
    ineq5 = abs(c2) <= 1;
    
    
    fourierOK = ineq1 & ineq2 & ineq3 & ineq4 & ineq5;

end