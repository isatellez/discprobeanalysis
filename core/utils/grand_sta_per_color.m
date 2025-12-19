function [w_chan, chanNames] = grand_sta_per_color(staPath)
    S = load(staPath);
    if ~isfield(S,'STA'), error('STA variable not found in %s', staPath); end
    STA = S.STA;
    if numel(STA)>1, STA = STA(1); end
    sta_flat = STA.sta_flat;
    chanNames = cellstr(STA.color_channels);
    nColor = numel(chanNames);
    [~, nCols] = size(sta_flat);
    
    if nCols == nColor
        w_chan = mean(sta_flat,1);
    elseif mod(nCols,nColor)==0
        nPix = nCols/nColor;
        sta_resh = reshape(sta_flat,[],nPix,nColor);
        w_chan = squeeze(sum(sum(sta_resh,1),2));
        w_chan = w_chan(:)';
    else
        w_chan = mean(sta_flat,1);
    end
    
    names_upper = upper(strtrim(chanNames(:)));
    if (nColor==3) && all(ismember(names_upper,{'R','G','B'}))
        [tL,tM,tS] = get_gun_LMS_weights();
        wL_val = 0; wM_val = 0; wS_val = 0;
        for i=1:3
            nm = names_upper{i}; val = w_chan(i);
            switch nm
                case 'R', wL_val=wL_val+val*tL(1); wM_val=wM_val+val*tM(1); wS_val=wS_val+val*tS(1);
                case 'G', wL_val=wL_val+val*tL(2); wM_val=wM_val+val*tM(2); wS_val=wS_val+val*tS(2);
                case 'B', wL_val=wL_val+val*tL(3); wM_val=wM_val+val*tM(3); wS_val=wS_val+val*tS(3);
            end
        end
        w_chan = [wL_val, wM_val, wS_val];
        chanNames = {'L','M','S'};
    end
end

function [tL,tM,tS] = get_gun_LMS_weights()
    rx=0.6280; ry=0.3310; gx=0.3059; gy=0.5826; bx=0.1639; by=0.0617;
    Wr=18.6469/2; Wg=75.8449/2; Wb=10.5313/2;
    R=cie2lms(rx,ry); G=cie2lms(gx,gy); B=cie2lms(bx,by);
    M = [R;G;B]' * diag([Wr Wg Wb]);
    tL=M(1,:)'; tM=M(2,:)'; tS=M(3,:)';
end

function v = cie2lms(a,b)
    x=a; y=b; z=1-x-y;
    M = [.15514 .54316 -.03286; -.15514 .45684 .03286; 0 0 .01608];
    lms = (M*[x;y;z])';
    den = lms(1)+lms(2);
    v = [lms(1)/den, lms(2)/den, lms(3)/den];
end