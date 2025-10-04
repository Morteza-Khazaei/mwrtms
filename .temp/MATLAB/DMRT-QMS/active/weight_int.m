function [I] = weight_int(vec,mut,gamma,mu,aa,step)

% first interpolate vec from mut to mu, by step
ndeg = length(mut);
ndeg2 = length(mu);
vec2d = reshape(vec,step,ndeg);
vec2 = spline(mut',vec2d,mu');
% weight
for i = 1:ndeg2
    vec2(:,i) = vec2(:,i)*aa(i);
end
vec2 = reshape(vec2,step*ndeg2,1);
I2 = gamma*vec2;
% interpolate back to mut
I2d = reshape(I2,step,ndeg2);
I = spline(mu',I2d,mut');
I = reshape(I,step*ndeg,1)/2;

end