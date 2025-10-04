function vec = weight(vec0,aa,s)
% weight of vector, 
% assert(): length(aa)*s = length(vec0)
vec = vec0;
for i = 1:length(aa)
    vec((i - 1)*s + 1:i*s) = vec((i - 1)*s + 1:i*s)*aa(i);
end
end