function actualName = resolve_session_folder(basePath, inputStr)
    actualName = '';
    if exist(fullfile(basePath, inputStr), 'dir')
        actualName = inputStr;
        return;
    end
    if regexp(inputStr, '^\d{6}$')
        d = dir(fullfile(basePath, ['*' inputStr '*']));
        isDir = [d.isdir] & ~ismember({d.name}, {'.','..'});
        matches = d(isDir);
        if ~isempty(matches), actualName = matches(1).name; end
    end
end