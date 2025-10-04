function [vv, hh, x] = NMM3D_LUT_NRCS_40degree_interp(rms0,eps0,ratio0)
% interpolate the NMM3D look up table of 40 degree incidence angle to find
% backscatter of bare soil surface. Interpolation is performed in dB scale.
% In the LUT (look up table) file, each row is a special case, and contains
% eight columns, listing the incidence angle in degree (1), ratio (2), 
% real part of permittivity (3), imag part of permittivity (4), 
% normalized rms height with respect to wavelength in top medium (5), 
% vv (6), hh (7), and crosspol (8) backscatter, respectively.
% 
% inputs:
%   rms0, normalized rms height over wavelength
%   eps0, real part of perm2/perm1, 
%       perm1 is permittivity of the top medium, perm2 is permittivity of
%       the bottom mediium.
%   ratio0, correlation length/rms height
% 
% outputs:
%   vv, hh and cross-pol backscatter in dB scale.
% 

LUT = load('NMM3D_LUT_NRCS_40degree.dat');
LUT = LUT(:,2:end);
ratio = [4 7 10 15];
if ratio0 < ratio(1) || ratio0 > ratio(end)
    error('ratio out of range');
end
r = find(ratio == ratio0, 1);
if ~isempty(r)
    [vv,hh,x] = rblock_interp(LUT,ratio0, eps0, rms0);
else
    rs = find(ratio < ratio0, 1);
    [vvs,hhs,xs] = rblock_interp(LUT,ratio(rs), eps0, rms0);
    rb = find(ratio > ratio0, 1);   
    [vvb,hhb,xb] = rblock_interp(LUT,ratio(rb), eps0, rms0);
    vv = r_interp(ratio(rs),ratio(rb),vvs,vvb,ratio0);
    hh = r_interp(ratio(rs),ratio(rb),hhs,hhb,ratio0);
    x = r_interp(ratio(rs),ratio(rb),xs,xb,ratio0);
end
    
end

function [vv,hh,x] = rblock_interp(LUT,ratio, eps0, rms0)
% construct 2D block for a particular ratio

idx = find(LUT(:,1) == ratio);
N = length(idx);
blk = LUT(idx,2:end);
idx = find(blk(:,1) == blk(1,1));
nrms = length(idx);
rms = blk(idx,3);
if rms0 < rms(1) || rms0 > rms(end)
    error('RMS out of range');
end

neps = N/nrms;
eps = blk(1:nrms:N,1);
if eps0 < eps(1) || eps0 > eps(end)
    error('permittivity out of range');
end

vv_blk = reshape(blk(:,4),[nrms neps]);
hh_blk = reshape(blk(:,5),[nrms neps]);
x_blk = reshape(blk(:,6),[nrms neps]);

vv = interp2(eps,rms,vv_blk,eps0,rms0);
hh = interp2(eps,rms,hh_blk,eps0,rms0);
x = interp2(eps,rms,x_blk,eps0,rms0);
end

function val = r_interp(rs,rb,vs,vb,r0)
% linear interpolation over ratio, need to consider -inf.
% 
if isfinite(vs) && isfinite(vb)
    val = interp1([rs rb], [vs vb], r0);  
else
    val = -inf;
end
end
