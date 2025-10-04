function I1U = CalcIU(coem,vecup,vecum,eigvm,ke,d,ctinc,part1u,part2u,tup)
% reconstruct the specific intensity at the top boundary in the air

% vecup and vecum are eigen vectors in the top layer
% eigvm are eigenvalues in the top layer
% ke is the top layer extinction coefficient
% ctinc is the top lay incident angle cosine
% d is the top layer depth
% part1u,part2u are the particular solution at the top layer
% tup is tup10.

[dimm,npol] = size(coem);
dimh = dimm/2;

I1U = zeros(dimh,npol);

for ipol = 1:npol
    % homogeneous solution
    for jcol=1:dimh
        I1UP = coem(jcol,ipol) * vecup(:,jcol);
        I1U(:,ipol) = I1U(:,ipol) + I1UP;
        arg = eigvm(jcol) * d;
        if (abs(real(arg))<20)
            I1UM = coem(jcol + dimh,ipol) * vecum(:,jcol);
            I1U(:,ipol) = I1U(:,ipol) + I1UM * exp(arg);
        end
    end
    % particular solution
    I1U(:,ipol) = I1U(:,ipol) + part2u(:,ipol);
    arg2 = ke * d / ctinc; 
    if (abs(real(arg2))<20)
        I1U(:,ipol) = I1U(:,ipol) + part1u(:,ipol) * exp(-arg2);
    end
    I1U(:,ipol) = tup * I1U(:,ipol);
end

    

end
