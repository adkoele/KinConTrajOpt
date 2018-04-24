function g = objgrad(X,params)

N = params.N;
nstates = params.nstates;
ncontrols = params.ncontrols;
nvarpernode = params.nvarpernode;

iu = nstates+(1:ncontrols);
g = zeros(size(X));
for i = 1:N
    u = X(iu);
    g(iu) = 2*u/N;
    iu = iu+nvarpernode;
end