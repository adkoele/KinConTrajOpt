function f = objfun(X, params)

N = params.N;
nstates = params.nstates;
ncontrols = params.ncontrols;
nvarpernode = params.nvarpernode;

iu = nstates+(1:ncontrols);
f = 0;
for i = 1:N
    u = X(iu);
    f = f+u^2;
    iu = iu+nvarpernode;
end

f = f/N;