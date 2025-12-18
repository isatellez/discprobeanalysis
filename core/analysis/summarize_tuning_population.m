function pop = summarize_tuning_population(allUnitTuning, config)
% Collect preferred directions and widths across units

    nUnits = numel(allUnitTuning);
    phi0 = nan(nUnits,1);
    fwhm = nan(nUnits,1);
    R2   = nan(nUnits,1);

    for ii = 1:nUnits
        ut = allUnitTuning{ii};
        if isempty(ut) || isempty(ut.fit)
            continue;
        end
        phi0(ii) = ut.fit.phi0_deg;
        fwhm(ii) = ut.fit.fwhm_deg;
        R2(ii)   = ut.fit.R2;
    end

    good = R2 >= 0.5;   % tweak threshold if necessary

    pop.phi0_deg = phi0;
    pop.fwhm_deg = fwhm;
    pop.R2       = R2;
    pop.good     = good;
end
