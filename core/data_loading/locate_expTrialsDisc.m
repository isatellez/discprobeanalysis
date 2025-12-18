function expFile = locate_expTrialsDisc(PATHS, MONKEY, date6)
% find the ExpTrialsDisc .mat file for a date

    %location should be /DiscProbe/250513/Jacomo_250513_ExpTrialsDisc.mat
    cand1 = fullfile(PATHS.baseDiscProbeLocal, date6, ...
                     sprintf('%s_%s_ExpTrialsDisc.mat', MONKEY, date6));
    if isfile(cand1)
        expFile = cand1;
        return;
    end

    % but you can also search anywhere under baseDiscProbeLocal
    pat = sprintf('*%s*ExpTrialsDisc.mat', date6);
    D = dir(fullfile(PATHS.baseDiscProbeLocal, '**', pat));

    if ~isempty(D)
        [~, ix] = max([D.datenum]);
        expFile = fullfile(D(ix).folder, D(ix).name);
        fprintf('[info] Using ExpTrialsDisc at: %s\n', expFile);
        return;
    end

    error('Could not locate an ExpTrialsDisc .mat for date %s (base=%s)', ...
          date6, PATHS.baseDiscProbeLocal);
end
