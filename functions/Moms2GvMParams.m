function  [K1mu1_K2mu2, Errs] = Moms2GvMParams(TargSpecMoms, tolerance)


% Takes moment list in the form [...(PSD) (MWD(rads)), m1, m2, n2]
if nargin < 2
    tolerance = 0.01;
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
    if size(TargSpecMoms, 2) >= 5
        PSD = TargSpecMoms(:,end-4);
    else
        PSD = 1;
    end
    
    if size(TargSpecMoms, 2) >= 4
        MWD = wrapToPi(TargSpecMoms(:,end-3));
    else
        MWD = 0;
    end
      
    
    
    
    OKParams = ~isnan(prod(TargCentredMoms, 2));

    [GvMParamsRel, InitialK1K2Psi,GvMLocsD0, Errs] = deal(NaN(size(TargCentredMoms)));

    

        
%  Get Starting Parameter Set
    for paramindx = 1:3
        InitialK1K2Psi(OKParams,paramindx) = SIGrid{paramindx}(TargCentredMoms(OKParams,1), TargCentredMoms(OKParams,2), TargCentredMoms(OKParams,3));
    end
    
    for fn = 1:size(TargCentredMoms,1)
        
       if OKParams(fn)
           GvMParamsRel(fn,:) = CallNewton(TargCentredMoms(fn,:), InitialK1K2Psi(fn,:), tolerance);
           GvMLocsD0(fn,:) = MomsAndJacobianFlexi(GvMParamsRel(fn,:), 50, true); % returns location params
           if nargout == 2
               Errs(fn,:) = MomsAndJacobianFlexi(GvMParamsRel(fn,:))-TargCentredMoms(fn,:);
           end
       end
    end
    
    K1mu1_K2mu2 = NaN(size(GvMParamsRel,1), 4);
    K1mu1_K2mu2(:,[1 3]) = GvMParamsRel(:,1:2); % kappa1 kappa2
    K1mu1_K2mu2(:,2) = GvMLocsD0(:,1) + MWD; % mu 1
    K1mu1_K2mu2(:,4) = GvMLocsD0(:,2) + MWD; % mu 2
end

   