function J = conjac(X, params)

nkinconst = params.nokinconst; %number of kinematic constraints (phi(x)=0)

nstates = params.nstates;
ndof = params.ndof;
ncontrols = params.ncontrols;
N = params.N;
ncon = params.ncon;
nvars = params.nvars;
nvarpernode = params.nvarpernode;
nconpernode = params.nconpernode;
h = params.T/(N-1);

% initialize indices
ix1 = 1:ndof;
iv1 = ndof+1:nstates;
iu1 = nstates+(1:ncontrols);
il1 = nstates+ncontrols+(1:nkinconst); %lambda
ilb = nstates+ncontrols+nkinconst+(1:nkinconst); %lambda bar
igb = nstates+ncontrols+nkinconst*2+(1:nkinconst); %gamma bar

ick = 1:nkinconst; %index for kinematic constraints
ic = 1:nstates; %index for dynamic constraints

% Go through dynamics constraints
J = spalloc(ncon,nvars,params.Jnnz);
for i = 1:N-1
    ix2 = ix1+nvarpernode;
    iv2 = iv1+nvarpernode;
    iu2 = iu1+nvarpernode;
    il2 = il1+nvarpernode;
    
    x1 = X(ix1); % x and y
    x2 = X(ix2);
    v1 = X(iv1); %vx and vy
    v2 = X(iv2);
    u1 = X(iu1);
    u2 = X(iu2);
    l1 = X(il1);
    l2 = X(il2);
    lb = X(ilb);
    gb = X(igb);
    
    % Kinematic constraints
    [~, jac, djac,derivs] = getKinConstraints(x1,v1,params);
    [fbar,dfbardx1,dfbardu1,dfbardl1] = dynamics(x1,u1,l1,params);
    %c((i-1)*nconpernode+ick) = kin;
    J((i-1)*nconpernode+ick,ix1) = derivs.dkindx;
    
    % c((i-1)*nconpernode+nkinconst+ick) = jac*v1;
    J((i-1)*nconpernode+nkinconst+ick,ix1) = derivs.djacdx*v1;
    J((i-1)*nconpernode+nkinconst+ick,iv1) = jac;
    
    %c((i-1)*nconpernode+2*nkinconst+ick) = djac*v1+jac*fbar;
    J((i-1)*nconpernode+2*nkinconst+ick,ix1) = jac*dfbardx1+(derivs.djacdx*fbar)';
    J((i-1)*nconpernode+2*nkinconst+ick,iv1) = djac*eye(1)+(derivs.ddjacdv*v1)';
    J((i-1)*nconpernode+2*nkinconst+ick,iu1) = jac*dfbardu1;
    J((i-1)*nconpernode+2*nkinconst+ick,il1) = jac*dfbardl1;
    
    %USE BE dyns
    %c((i-1)*nconpernode+3*nkinconst+ic) = dynconstraint(x2,v2,(x1+x2)/h,(v1+v2)/h,u2,lb,gb,params);
    [~,dfdx,dfdv,dfddx,dfddv,dfdu,dfdlb,dfdgb] = dynconstraint(x2,v2,(x2-x1)/h,(v2-v1)/h,u2,lb,gb,params);
    J((i-1)*nconpernode+3*nkinconst+ic,ix1) = -dfddx/h;
    J((i-1)*nconpernode+3*nkinconst+ic,iv1) = -dfddv/h;
    J((i-1)*nconpernode+3*nkinconst+ic,ilb) = dfdlb;
    J((i-1)*nconpernode+3*nkinconst+ic,igb) = dfdgb;
    J((i-1)*nconpernode+3*nkinconst+ic,ix2) = dfdx+dfddx/h;
    J((i-1)*nconpernode+3*nkinconst+ic,iv2) = dfdv+dfddv/h;
    J((i-1)*nconpernode+3*nkinconst+ic,iu2) = dfdu;
    
    ix1 = ix2;
    iv1 = iv2;
    iu1 = iu2;
    il1 = il2;
    ilb = ilb+nvarpernode;
    igb = igb+nvarpernode;
end
%Kinematic constraints for time point N
[~, jac, djac,derivs] = getKinConstraints(x2,v2,params);
[fbar,dfbardx2,dfbardu2,dfbardl2] = dynamics(x2,u2,l2,params);
%c(i*nconpernode+ick) = kin;
J(i*nconpernode+ick,ix2) = derivs.dkindx;

% c((i-1)*nconpernode+nkinconst+ick) = jac*v2;
J(i*nconpernode+nkinconst+ick,ix2) = derivs.djacdx*v2;
J(i*nconpernode+nkinconst+ick,iv1) = jac;

%c((i-1)*nconpernode+2*nkinconst+ick) = djac*v1+jac*fbar;
J(i*nconpernode+2*nkinconst+ick,ix2) = jac*dfbardx2+(derivs.djacdx*fbar)';
J(i*nconpernode+2*nkinconst+ick,iv2) = djac*eye(1)+(derivs.ddjacdv*v2)';
J(i*nconpernode+2*nkinconst+ick,iu2) = jac*dfbardu2;
J(i*nconpernode+2*nkinconst+ick,il2) = jac*dfbardl2;

%Task constraints
%c(i*nconpernode+(1:nstates)) = X(1:nstates)-[0; -L;0;0];
J(i*nconpernode+3*nkinconst+(1:nstates),1:nstates) = eye(nstates);

%c(i*nconpernode+nstates+(1:nstates)) = [x2;v2]-[0; L;0;0]; %because the current ix2 is at point N
% J(i*nconpernode+3*nkinconst+nstates+(1:nstates),[ix2 iv2]) = eye(nstates);

%c(i*nconpernode+3*nkinconst+nstates+(1:nstates-1)) = [atan2(x2(2),x2(1));v2]-[atan2(X(2),X(1))+2*pi;0;0];
J(i*nconpernode+3*nkinconst+nstates+(1:nstates-1),ix2) = [1/(1+(x2(2)/x2(1))^2)*(-x2(2)/(x2(1)^2))  1/(1+(x2(2)/x2(1))^2)*(1/x2(1));zeros(ndof)];
J(i*nconpernode+3*nkinconst+nstates+(1:nstates-1),1:2) = -[1/(1+(X(2)/X(1))^2)*(-X(2)/(X(1)^2))  1/(1+(X(2)/X(1))^2)*(1/X(1));zeros(ndof)];
J(i*nconpernode+3*nkinconst+nstates+(1:nstates-1),iv2) = [zeros(1,ndof);eye(ndof)];