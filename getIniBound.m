function [X, L, U] = getIniBound(params)

l = params.l;
% Bounds on x and u for each node
Lallnode = [-2*l;-3*l;-1000;-1000;-1000;-1000]; %x y vx vy u l
Uallnode = [ 2*l; 2*l;1000;1000;1000;1000];

Lpernode = [Lallnode;-1000;-5]; %lb gb
Upernode = [Uallnode;1000;5];

L = [];
U = [];
for i = 1:params.N-1
    L = [L;Lpernode];
    U = [U;Upernode];
end

L = [L;Lallnode];
U = [U;Uallnode];

L(5) = 0;
U(5) = 0;
%Mid initial guess
X = 1/2*(L+U);
    