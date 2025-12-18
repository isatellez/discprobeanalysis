% calibration/generate_dkl_lms_tables.m
function generate_dkl_lms_tables()

S = load('/Users/tellezi2/Documents/DiscProbesCode/masterRGBUpdated.mat');
fullRGB = double(S.masterRGB);

if max(fullRGB(:)) > 1
    RGB = fullRGB/255;
else
    RGB = fullRGB;
end

rx=0.6280; ry=0.3310; gx=0.3059; gy=0.5826; bx=0.1639; by=0.0617;
Wr=18.6469/2; Wg=75.8449/2; Wb=10.5313/2;

R = cie2lms(rx,ry); G = cie2lms(gx,gy); B = cie2lms(bx,by);
Lr=R(1); Lg=G(1); Lb=B(1);
Mr=R(2); Mg=G(2); Mb=B(2);
Sr=R(3); Sg=G(3); Sb=B(3);

M_rgb2lms = [Lr Lg Lb; Mr Mg Mb; Sr Sg Sb] * diag([Wr Wg Wb]);
LMS = RGB * M_rgb2lms.';

[ldrgyvMAT, lmrgyvMAT, CIS] = calibMethod2_fb_01_10_2024();
Minv = inv(ldrgyvMAT);

DKL = 2*(RGB - 0.5) * Minv.';

dklname = fullfile(pwd,'dkl.csv');
lmsname = fullfile(pwd,'lms.csv');
rgbname = fullfile(pwd,'rgb.csv');

writematrix(DKL, dklname);
writematrix(LMS, lmsname);
writematrix(RGB, rgbname);

save('calib_outputs.mat','RGB','LMS','DKL','M_rgb2lms','ldrgyvMAT','lmrgyvMAT','CIS');

end

function v = cie2lms(a,b)
x=a; y=b; z=1-x-y;
M=[ .15514 .54316 -.03286; -.15514 .45684 .03286; 0 0 .01608];
lms = (M*[x;y;z]).';
den = lms(1)+lms(2);
v = [lms(1)/den, lms(2)/den, lms(3)/den];
end
