classdef spatial_transforms
    % Utilities for mapping hues/colors into angle / DKL space.

    methods (Static)

        function theta_deg = hue_to_angle(hues, COL, config) %like index to dkl for example
        if nargin < 2 || isempty(COL)
            COL = struct();
        end
        if nargin < 3
            config = struct();
        end
    
        hues = double(hues);
        theta_deg = nan(size(hues));
    
        mode = 'index';
        if isfield(config, 'space') && isfield(config.space, 'mode')
            mode = lower(config.space.mode);
        end
    
        useDKL = strcmp(mode, 'dkl') && isfield(COL, 'dkl') && ~isempty(COL.dkl);
    
        if useDKL
            % pick saturation plane
            if isfield(config, 'space') && isfield(config.space, 'satIndexForAngle')
                satIdx = config.space.satIndexForAngle;
            elseif isfield(COL, 'nSat')
                satIdx = COL.nSat;
            else
                sz = size(COL.dkl);
                satIdx = sz(1);
            end
    
                       sz = size(COL.dkl);
            h = round(hues);
            valid = ~isnan(h) & h >= 1 & h <= sz(2);

            theta_deg(~valid) = NaN;
            if any(valid)
                % DKL is [LD RG YV] from calib script / dkl.csv
                D = squeeze(COL.dkl(satIdx, :, :));   % [nHue x 3]

                hv = h(valid);
                v  = D(hv, 2:3);   % columns 2–3: [RG YV]
                vx = v(:,1);       % RG
                vy = v(:,2);       % YV

                theta = atan2(vy, vx);               % angle in RG–YV plane
                theta_deg(valid) = mod(rad2deg(theta), 360);
            end
        else
            % index mode as before
            if isfield(COL, 'nHue')
                nHueAll = COL.nHue;
            elseif isfield(config, 'space') && isfield(config.space, 'nHue')
                nHueAll = config.space.nHue;
            else
                nHueAll = max(hues(~isnan(hues)));
            end
    
            zeroHue = 1;
            if isfield(config, 'space') && isfield(config.space, 'zeroHue')
                zeroHue = config.space.zeroHue;
            end
    
            theta_deg = (hues - zeroHue) .* (360 / nHueAll);
            theta_deg = mod(theta_deg, 360);
        end
    end


        function [L, M, S] = index_to_dkl(idx, COL)
            % Map "master row" indices into DKL coordinates (LD, RG, YV).
        
            L = nan(size(idx));
            M = nan(size(idx));
            S = nan(size(idx));
        
            if ~isfield(COL, 'dklRows') || isempty(COL.dklRows) || isempty(idx)
                return;
            end
        
            D  = COL.dklRows;      % [nProbes x 3]
            nD = size(D,1);
        
            idx   = double(idx);
            valid = ~isnan(idx) & idx >= 1 & idx <= nD;
            if ~any(valid)
                return;
            end
        
            useIdx    = idx(valid);
            L(valid)  = D(useIdx, 1);
            M(valid)  = D(useIdx, 2);
            S(valid)  = D(useIdx, 3);
        end



        function [theta_deg, r] = lm_to_angle(L, M)
            % Convert L/M coordinates into polar angle and radius.

            L = double(L);
            M = double(M);

            theta = atan2(M, L);
            theta_deg = mod(rad2deg(theta), 360);
            r = hypot(L, M);
        end

    end
end
