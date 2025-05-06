function [P,Q] = PoweredGenerators(N, r, m, CM, Cf, indcum, Ua, Va)
%This function returns the generators of (U*V).^m

% C = MultiPermuations(m,r);
C = CM(indcum(m)+1:indcum(m+1),:);
cf = Cf(indcum(m)+1:indcum(m+1));
[num_gen,~] = size(C);

P = zeros(N,num_gen);
Q = zeros(N,num_gen);
tN = 1:N;

for i = 1:num_gen
    %
    p = ones(N,1);
    q = ones(N,1);
    %p = zeros(N,r);
    %q = zeros(N,r);

    for j = 1:r
        row = tN +(j-1)*N;%(j-1)*N+1:j*N;
        col = C(i,j)+1;
        p = p.*Ua(row,col);
        q = q.*Va(row,col);
        %p(:,j) = Ua(row,col);
        %q(:,j) = Va(row,col);
    end
    %p = prod(p,2);
    %q = prod(q,2);

    % cf = multcoef(m,C(i,:));
    % cf = multinomial(m,C(i,:));

    P(:,i) = cf(i)*p;
    Q(:,i) = q;
end

end