function info = parse_session_id(id)
% parse session identifier into components
% can be:
%   '250513'
%   'Jacomo_250513'
%   'Jacomo_250513_anything" in case there is two sessions
%
% Returns struct with fields:
%   sessionID  - full string as given
%   dateStr    - 6-digit YYMMDD if found, otherwise original string
%   monkey     - leading chunk before the 6-digit date, if any
%   suffix     - anything after the date, if any

    s = string(id);

    info = struct( ...
        'sessionID', s, ...
        'dateStr',   "", ...
        'monkey',    "", ...
        'suffix',    "" );

    %usual pattern is Monkey_250513_optionalStuff
    tok = regexp(char(s), ...
        '^(?<monkey>[^_]+)_(?<date>\d{6})(?<suffix>.*)$', ...
        'names');

    if ~isempty(tok)
        info.monkey  = string(tok.monkey);
        info.dateStr = string(tok.date);
        info.suffix  = string(tok.suffix);
        return;
    end

    % otherwise try to find any 6-digit date inside
    tok2 = regexp(char(s), '(?<date>\d{6})', 'names', 'once');
    if ~isempty(tok2)
        info.dateStr = string(tok2.date);
    else
        % worst case treat whole string as "dateStr"
        info.dateStr = s;
    end
end
