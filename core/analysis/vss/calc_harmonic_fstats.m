function [F, p_val, df_num, df_den] = calc_harmonic_fstats(varExp, nHues, nParams)
if nargin < 3 || isempty(nParams)
    nParams = 2;
end

R2 = varExp;

if isnan(R2)
    F = NaN;
    p_val = NaN;
    df_num = NaN;
    df_den = NaN;
    return;
end

df_num = nParams;
df_den = nHues - nParams - 1;

if df_den <= 0
    F = NaN;
    p_val = NaN;
    df_num = NaN;
    df_den = NaN;
    return;
end

if R2 >= 1
    F = Inf;
    p_val = 0;
    return;
end

F = (R2 / df_num) / ((1 - R2) / df_den);
p_val = 1 - fcdf(F, df_num, df_den);
end
