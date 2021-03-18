function  [F0_K1mu1_K2mu2, Errs, IterCounts] = Offset_Moms_2_GvM_Parameters(TargSpecMoms, tolerance, trunc)


% Takes moment list in the form [...(PSD) (MWD(rads)), m1, m2, n2]
    
if nargin < 3 
    trunc = 40;
    if nargin<2
        tolerance = 0.001;
    end
end


global SIGrid;  
% This is the 3D interpolation grid set used to find initial guess for GvM
% parameters.  Delaring as global means that you only have to generate it
% once, even if function is called multiple times (eg if you are
% calculating multiple spectra)

% If it hasn't already been generated, do so...
if isempty(SIGrid)
    tic
    disp('Generating interpolation grid for first approximation of GvM parameters...');
    SIGrid = GetInterpGrid;
    toc
end

    TargCentredMoms = TargSpecMoms(:,(end-2):end);

    
    if size(TargSpecMoms, 2) >= 4
        theta0 = wrapToPi(TargSpecMoms(:,end-3));
    else
        theta0 = 0;
    end
      
    
    
    
    OKParams = ~isnan(prod(TargCentredMoms, 2));

    [GvMParamsRel, InitialK1K2Psi,GvMLocsD0, Errs] = deal(NaN(size(TargCentredMoms)));
    IterCounts = zeros(size(TargCentredMoms,1), 2);
    

        
%  Get Starting Parameter Set
    for paramindx = 1:3
        InitialK1K2Psi(OKParams,paramindx) = SIGrid{paramindx}(TargCentredMoms(OKParams,1), TargCentredMoms(OKParams,2), TargCentredMoms(OKParams,3));
    end
    
    for fn = 1:size(TargCentredMoms,1)
        
       if OKParams(fn)
           [GvMParamsRel(fn,:), IterCounts(fn,:)] = CallNewton(TargCentredMoms(fn,:), InitialK1K2Psi(fn,:), tolerance, trunc);
           GvMLocsD0(fn,:) = MomsAndJacobianFlexi(GvMParamsRel(fn,:), trunc, true); % returns location params
           if nargout >= 2
               Errs(fn,:) = MomsAndJacobianFlexi(GvMParamsRel(fn,:))-TargCentredMoms(fn,:);
           end
       end
    end
    
    F0_K1mu1_K2mu2 = NaN(size(GvMParamsRel,1), 5);    
    F0_K1mu1_K2mu2(:,1) = GvMLocsD0(:,3); % F0
    F0_K1mu1_K2mu2(:,[2 4]) = GvMParamsRel(:,1:2); % kappa1 kappa2
    F0_K1mu1_K2mu2(:,3) = GvMLocsD0(:,1) + theta0; % mu 1
    F0_K1mu1_K2mu2(:,5) = GvMLocsD0(:,2) + theta0; % mu 2
    

end

   