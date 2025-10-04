function [Am0, Am] = CalcAm(mu,aa,P0,PC,PS)

[Mmax,ntai,ntas,~,~] = size(PC);
Am0 = zeros(4*ntas,4*ntai);
Am = zeros(8*ntas,8*ntai,Mmax);

for idegi=1:ntai
    idi_ = block(idegi,4);
    idi = block(idegi,8);
    for idegs=1:ntas 
        ids_ = block(idegs,4);
        ids = block(idegs,8);  
        coe = pi*aa(idegi)/mu(idegs);        
        Am0(ids_,idi_)=2*coe*P0(idegs,idegi,:,:);                     
        for m = 1:Mmax
            Am(ids(1:4),idi(1:4),m)= coe*PC(m,idegs,idegi,:,:);
            Am(ids(1:4),idi(5:8),m)=-coe*PS(m,idegs,idegi,:,:);
            Am(ids(5:8),idi(1:4),m)= coe*PS(m,idegs,idegi,:,:);
            Am(ids(5:8),idi(5:8),m)= coe*PC(m,idegs,idegi,:,:);
        end
    end 
end

end

function id = block(x,s)
id = (x-1)*s + 1:x*s;
end
