clear all
close all
clc

Ts = [0.5 1 2 4];
angles = [45 90 180 360];

for k = 1:6
    for i = 1:length(Ts)
        for j = 1:length(angles)
            % global params 

            tic
            params.m = 1;
            params.l = 1;
            params.g = 9.81;
            params.T = Ts(i);
            params.targetangle = angles(j)/180*pi;

            params.N = params.T*100; %collocation nodes
            params.ndof = 2;
            params.ncontrols = 1;
            params.nokinconst = 1;
            params.solver = 'IPOPT';
            params.snoptname = 'test1';
            params.warmstart = 0;

            params = getparams(params);

            [X0, L, U] = getIniBound(params);

            [~, params] = conjacstructure(L, U, params);
            
            %uncomment to check derivatives            
            % [grad,grad_num,cjac, cjac_num] = checkDerivatives(1/2*(L+U)+randn(size(L)), params);
            % keyboard
            
            result = Optimize(L+(U-L).*rand(size(L)), L, U, params);%X0
            result.solvetime = toc;

            % if result.info == 0 %good solution -> plot
                x1 = reshape(result.X(1:params.nvarpernode*(params.N-1)),params.nvarpernode,params.N-1);
                x = [x1(1:params.ndof,:) result.X(end-params.nvarallnode+(1:params.ndof))];
                v = [x1(params.ndof+1:params.nstates,:) result.X(end-params.nvarallnode+(params.ndof+1:params.nstates))];
                u = [x1(params.nstates+(1:params.ncontrols),:) result.X(end-params.nvarallnode+params.nstates+(1:params.ncontrols))];
                l = [x1(params.nstates+params.ncontrols+(1:params.nokinconst),:) result.X(end-params.nvarallnode+params.nstates+params.ncontrols+(1:params.nokinconst))];
                lb= x1(params.nstates+params.ncontrols+params.nokinconst+(1:params.nokinconst),:);
                gb= x1(params.nstates+params.ncontrols+params.nokinconst*2+(1:params.nokinconst),:);

%                 t = linspace(0,params.T,params.N);
%                 figure
%                 subplot(2,1,1)
%                 plot(t,x)
%                 ylabel('x')
%                 subplot(2,1,2)
%                 plot(t,v)
%                 ylabel('v')
%                 xlabel('Time')
% 
%                 figure;   
%                 theta = unwrap(atan2(x(2,:),x(1,:)));
%                 subplot(2,1,1)
%                 plot(t,theta/pi*180)
%                 ylabel('\theta')
%                 subplot(2,1,2)
%                 plot(t,u)
%                 xlabel('Time')
%                 ylabel('Input')
% 
%                 figure
%                 subplot(2,1,1)
%                 plot(t,l)
%                 hold on
%                 plot(t(1:end-1),lb,'r')
%                 subplot(2,1,2)
%                 plot(t(1:end-1),gb)

            % end

            filename = ['Result_' num2str(params.targetangle/pi*180) 'deg_' num2str(params.T) 'sec_' num2str(k) 'kin1.mat'];
            save(filename,'result')
        end
    end
end