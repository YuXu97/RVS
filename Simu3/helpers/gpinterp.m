function [hz, extr] = gpinterp(z,zp,hp)
% GPINTERP  Helper function for Gaussian Process evaluation.
%   HZ = GPINTERP(Z,ZP,HP) evaluates the function h(.) specified at the
%   grid {ZP, HP} at the points Z using interpolation.

if(min(z) < min(zp) || max(z) > max(zp))
%     warning('Extrapolating');
    extr = 1;
else
    extr = 0;
end

hz = interp1(zp(:),hp(:),z(:),'linear','extrap')';