function [I] = multiply_folder(vec,mut,gamma,mu,step)

% first multiply gamma and vec with length ndeg2 
ndeg = length(mut);
ndeg2 = length(mu);
I2 = gamma*vec;
% then interpolate back to mut, by step
I2d = reshape(I2,step,ndeg2);
I = spline(mu',I2d,mut');
I = reshape(I,step*ndeg,1);

end