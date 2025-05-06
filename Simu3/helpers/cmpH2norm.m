function H2 = cmpH2norm(A,B,C,Q,R)
% Compute the H2 norm of a state-space model with original input (u) and
% innovation seen as inputs (MISO).

nu = size(B,2);

try
    [P,~,K] = dare(A',C',Q,R);
catch me
    P = zeros(size(A));
    K = zeros(size(C,1), size(A,1));
end
s = sqrt(C*P*C' + R); % Compute the standard deviation of the innovations
warning off
H2 = norm(ss(A,[B s*K'],C,[zeros(1,nu) s],1),2);

%H2 = norm(ss(A,[B K'],C,[zeros(1,nu) 1],1),2); % Do not include the innovation variance in the norm
