function mask = get_legacy_trial_mask(U, config)

mask = true(U.nTrials,1);

if ~isfield(config,'legacy') || ~isfield(config.legacy,'useLegacySats')
    return;
end
if ~config.legacy.useLegacySats
    return;
end

if ~isfield(U,'trials') || ~isfield(U.trials,'satID')
    return;
end

sats = config.legacy.sats(:).';
tol  = config.legacy.satTol;

satID = U.trials.satID;
maskS = false(size(satID));

for s = sats
    maskS = maskS | abs(satID - s) < tol;
end

mask = mask & maskS;

end
