clear;

n = 5;


c = rand; lam = rand; rho = 2*rand-1;

Knew = zeros(n,n,n,n);

for t1 = 1:n
    for t2 = 1:n
        for s1 = 1:n
            for s2 = 1:n
                Knew(t1,t2,s1,s2) = lam^(t1+t2+s1+s2)*rho^(abs(t1+t2-s1-s2));
            end
        end
    end
end

Knew = c*Knew;

K1 = Knew(1,2,:,:);

K2 = Knew(2,1,:,:);


Kcd = zeros(n,n,n,n);

for t1 = 1:n
    for t2 = 1:n
        for s1 = 1:n
            for s2 = 1:n
                Kcd(t1,t2,s1,s2) = lam^(t1+t2+s1+s2)*rho^(abs(t1-s1)+abs(t2-s2));
            end
        end
    end
end

Kcd = c*Kcd;

K1 = Kcd(1,2,:,:);

K2 = Kcd(2,1,:,:);