function DATES = normalize_dates(inputDATES, config)
% Converts 'all', strings, or cells into a standard cell array of date strings.
    if nargin < 1 || isempty(inputDATES)
        inputDATES = 'all';
    end
    
    if (ischar(inputDATES) || isstring(inputDATES)) && strcmpi(strtrim(inputDATES), 'all')
        % Look in the base data directory
        d = dir(config.paths.base);
        names = {d([d.isdir]).name};
        % Find folders containing 6 digits (the date)
        validIdx = cellfun(@(s) ~isempty(regexp(s, '\d{6}', 'once')), names) & ...
                   ~ismember(names, {'.','..'});
        DATES = sort(names(validIdx));
    elseif ischar(inputDATES) || isstring(inputDATES)
        DATES = {char(inputDATES)};
    else
        DATES = cellstr(inputDATES);
    end
end