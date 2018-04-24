 function [f,dfdx,dfdv,dfddx,dfddv,dfdu,dfdlb,dfdgb] = dynconstraint(x,v,dx,dv,u,lb,gb,params)

ndof = params.ndof;
nstates = params.nstates;

[~, jac,~,derivs] = getKinConstraints(x,v,params);
[fdyn,dfdyndx,dfdyndu,dfdyndlb] = dynamics(x,u,lb,params);

f(1:ndof) = v+jac'*gb-dx;
f(ndof+1:nstates) = fdyn-dv;

dfdx = [derivs.djacdx'*gb;dfdyndx];
dfdv = [eye(ndof);zeros(ndof)];
dfddx = [-eye(ndof);zeros(ndof)];
dfddv = [zeros(ndof);-eye(ndof)];
dfdu = [zeros(ndof,1);dfdyndu];
dfdlb = [zeros(ndof,1);dfdyndlb];
dfdgb = [jac';zeros(ndof,1)];