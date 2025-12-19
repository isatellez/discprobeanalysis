function info = find_sta_file(sessionDir, unitType, phyID)
    info = struct('fullpath','','name','');
    staDir = fullfile(sessionDir, 'STAs');
    if ~exist(staDir,'dir'), return; end
    
    phyTag = sprintf('phy%05d', phyID);
    umTag = iff(strcmpi(unitType, 'SU'), 'SU', 'MU');
    
    d = dir(fullfile(staDir, sprintf('*%s_%s*_STA.mat', phyTag, umTag)));
    if isempty(d)
        d = dir(fullfile(staDir, sprintf('*%s*_STA.mat', phyTag)));
    end
    
    if ~isempty(d)
        info.fullpath = fullfile(staDir, d(1).name);
        info.name = d(1).name;
    end
end
function val = iff(cond,a,b), if cond, val=a; else, val=b; end, end