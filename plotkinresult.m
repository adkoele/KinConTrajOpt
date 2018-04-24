function plotkinresult(result)

params = result.params;
x1 = reshape(result.X(1:params.nvarpernode*(params.N-1)),params.nvarpernode,params.N-1);
x = [x1(1:params.ndof,:) result.X(end-params.nvarallnode+(1:params.ndof))];
v = [x1(params.ndof+1:params.nstates,:) result.X(end-params.nvarallnode+(params.ndof+1:params.nstates))];
u = [x1(params.nstates+(1:params.ncontrols),:) result.X(end-params.nvarallnode+params.nstates+(1:params.ncontrols))];
l = [x1(params.nstates+params.ncontrols+(1:params.nokinconst),:) result.X(end-params.nvarallnode+params.nstates+params.ncontrols+(1:params.nokinconst))];
lb= x1(params.nstates+params.ncontrols+params.nokinconst+(1:params.nokinconst),:);
gb= x1(params.nstates+params.ncontrols+params.nokinconst*2+(1:params.nokinconst),:);

t = linspace(0,params.T,params.N);
% figure
% subplot(2,1,1)
% plot(t,x)
% ylabel('x')
% subplot(2,1,2)
% plot(t,v)
% ylabel('v')
% xlabel('Time')

figure;   
theta = unwrap(atan2(x(2,:),x(1,:)));
subplot(2,1,1)
plot(t,theta/pi*180)
ylabel('\theta')
subplot(2,1,2)
plot(t,u)
xlabel('Time')
ylabel('Input')

% figure
% subplot(2,1,1)
% plot(t,l)
% hold on
% plot(t(1:end-1),lb,'r')
% subplot(2,1,2)
% plot(t(1:end-1),gb)