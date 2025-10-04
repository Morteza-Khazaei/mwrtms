function y = App(x,n)
% append or trucate to make sure y has the same length as x.

nx = length(x);
y = zeros(1,n);
if nx > n
    y = x(1:n);
elseif nx < n
    y(1:nx) = x;
    y(nx + 1:n) = x(end);
else
    y = x;
end