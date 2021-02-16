function [A_kap1_mu1_kap2_mu2, Errs] = moms2params(a1_b1_a2_b2, tolerance, trunc)
    
%     Given the first two pairs of trig moments, returns the GvM
%     parameters: A (divide by this to normalise), kappa1 (concentration
%     parameter 1), mu1 (location parameter 1), kappa2, mu2.  Can also
%     specify tolerance (or use 0.01 default value)
%     

    addpath([pwd, filesep, 'functions']);

momsOK = KrogstadTest(a1_b1_a2_b2);


if any(~momsOK) 
    warning('Some supplied moments are not moments of a non-negative distribution: will return NaNs');
    a1_b1_a2_b2(~momsOK,:) = NaN;
    
end


if nargin < 3 
    trunc = 40;
    if nargin<2
        tolerance = 0.01;
    end
end


    
    
    t0_m1_m2_n2 = Moms_2_Offset_Moms(a1_b1_a2_b2);

    
    
    [A_kap1_mu1_kap2_mu2, Errs] = Offset_Moms_2_GvM_Parameters(t0_m1_m2_n2, tolerance, trunc);
    


end