function c = confun(X, params)

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
l = params.l;

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
c = zeros(ncon,1);
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
    theta(i) = atan2(x1(2),x1(1));
    
    % Kinematic constraints
    [kin, jac, djac] = getKinConstraints(x1,v1,params);
    fbar = dynamics(x1,u1,l1,params);
    c((i-1)*nconpernode+ick) = kin;
    c((i-1)*nconpernode+nkinconst+ick) = jac*v1;
    c((i-1)*nconpernode+2*nkinconst+ick) = djac*v1+jac*fbar;
    %USE BE dyns
    c((i-1)*nconpernode+3*nkinconst+ic) = dynconstraint(x2,v2,(x2-x1)/h,(v2-v1)/h,u2,lb,gb,params);
    
    ix1 = ix2;
    iv1 = iv2;
    iu1 = iu2;
    il1 = il2;
    ilb = ilb+nvarpernode;
    igb = igb+nvarpernode;
end
%Kinematic constraints for time point N second index corresponds to time
%point N
theta(N) = atan2(x2(2),x2(1));
theta = unwrap(theta);

[kin, jac, djac] = getKinConstraints(x2,v2,params);
fbar = dynamics(x2,u2,l2,params);
c(i*nconpernode+ick) = kin;
c(i*nconpernode+nkinconst+ick) = jac*v2;
c(i*nconpernode+2*nkinconst+ick) = djac*v2+jac*fbar;

%Task constraints
c(i*nconpernode+3*nkinconst+(1:nstates)) = X(1:nstates)-[0; -l;0;0];
% c(i*nconpernode+3*nkinconst+nstates+(1:nstates)) = [x2;v2]-[0; l;0;0]; %because the current ix2 is at point N
c(i*nconpernode+3*nkinconst+nstates+(1:nstates-1)) = [theta(N);v2]-[theta(1)+params.targetangle;0;0];
