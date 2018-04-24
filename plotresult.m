function [meanu, uall,meanx] = plotresult(name)

%plot kinematic result
filename = [name 'kin.mat'];
load(filename);
plotkinresult(result);

filename = [name 'thet.mat'];
load(filename);

params = result.params;
if ~isfield(params, 'Ks')
    params.Ks = 2;
end
X = result.X;
x1 = reshape(X(1:end-params.Ks),params.nvarpernode, params.N);
x = x1(1:params.nstates/2,:);
u = x1(params.nstates+(1:params.ncontrols),:);

T = result.params.T;
h = T/(result.params.NperSU-1);

subplot(2,1,1)
hold on
plot([0:h:T],x(1,:)/pi*180, 'r')
subplot(2,1,2)
hold on
plot([0:h:T],u(1,:),'r')
legend('Kinematic','Original')
