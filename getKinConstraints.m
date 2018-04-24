function [kin, jac, djac,derivs] = getKinConstraints(x,v,params)

% This function is problem specific

kin = x(1)^2+x(2)^2-params.l^2;
jac = [2*x(1) 2*x(2)];
djac = [2*v(1) 2*v(2)];

derivs.dkindx = jac';
derivs.djacdx = 2*eye(2);
derivs.ddjacdv = 2*eye(2);