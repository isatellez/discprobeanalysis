function dates = normalize_dates(D)
% Normalize various date input formats into a cellstr of 'YYMMDD' char rows.

% single char, like '250513'
if ischar(D)
    dates = {char(D)};
    return;
end

% string or string array
if isstring(D)
    if isscalar(D)
        dates = {char(D)};
    else
        dates = cellstr(D);
    end
    return;
end

% cell array of anything string-like
if iscell(D)
    % convert each element to char
    dates = cellfun(@char, D, 'UniformOutput', false);
    return;
end

% numeric or other weird types – try to be forgiving
if isnumeric(D)
    % assume rows or elements are datenum-like or YYMMDD-ish, convert via num2str
    D = num2cell(D);
    dates = cellfun(@(x) char(string(x)), D, 'UniformOutput', false);
    return;
end

% fallback
dates = {char(string(D))};

end
