function COL = add_lms_axes(COL)
% Build L-M / S / Lum axes from COL.lms (columns = [L M S])

if ~isfield(COL,'lms')
    error('COL.lms missing. Make sure lms.csv is loaded into COL.lms.');
end

L = COL.lms(:,1);
M = COL.lms(:,2);
S = COL.lms(:,3);

COL.axis_LM  = L - M;
COL.axis_S   = S;
COL.axis_Lum = L + M;
end
