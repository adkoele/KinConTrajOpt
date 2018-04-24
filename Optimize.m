function result = Optimize(X0, L, U, params)

if strcmp(params.solver, 'IPOPT')
    funcs.objective = @(X) objfun(X,params);
    funcs.gradient  = @(X) objgrad(X,params);
    funcs.constraints = @(X) confun(X, params);
    funcs.jacobian    = @(X) conjac(X, params);
    funcs.jacobianstructure = @() conjacstructure(L,U,params);
    options.lb = L;
    options.ub = U;
    options.cl = zeros(params.ncon,1);
    options.cu = zeros(params.ncon,1);	

    options.ipopt.max_iter = 10000;
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy = 'adaptive';		% worked better than 'monotone'
    options.ipopt.tol = 1e-4;
    options.ipopt.linear_solver = 'mumps'; %'ma57';% 
    options.ipopt.constr_viol_tol = 1e-3;
    options.ipopt.compl_inf_tol = 1e-3;
    options.ipopt.print_level = 5;
    [X, info] = ipopt(X0,funcs,options);
    disp(['IPOPT status: ' num2str(info.status)]);
    result.info = info.status;
    result.X = X;
    result.params = params;
    result.obj = objfun(X, params);
elseif strcmp(params.solver,'SNOPT')
    testspec.spc = which('testspec.spc');
    snspec ( testspec.spc );
    % Output informative files
    snprint   ([params.snoptname '.out']);
    snsummary ([params.snoptname '.sum']);
    cl = [-inf;zeros(params.ncon,1)];
    cu = [ inf;zeros(params.ncon,1)];
    FL = cl;
    FU = cu;
    xmul = zeros(size(L));
    Fmul = zeros(size(FL));
    xstate = zeros(size(X0));
    Fstate = zeros(size(FL));
    A = [];
    iAfun = [];
    jAvar = [];
    iGfun = params.DGrow;
    jGvar = params.DGcol;
    if params.warmstart == 1;
        snset ('Warm start')
    end
    [X,F,INFO] = snopt(X0,L,U,xmul, xstate,FL,FU,Fmul, Fstate, @objconfunp, ...
        A, iAfun, jAvar, iGfun, jGvar );
    snprint   off;
    snsummary off;
    result.info = INFO;
    result.obj = F(1);
    result.X = X;
    result.params = params;
end
