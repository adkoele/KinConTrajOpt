function params = getparams(params)

params.nstates = params.ndof*2;
params.nvarallnode = params.nstates+params.ncontrols+params.nokinconst;
params.nvarpernode = params.nvarallnode+params.nokinconst*2;
params.nvars = params.nvarpernode*(params.N-1)+params.nvarallnode;
params.nconpernode = params.nstates+params.nokinconst*3;
params.ncon = params.nconpernode*(params.N-1)+params.nokinconst*3+params.nstates*2;