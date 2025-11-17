%% Aiman Mushtaq Purra and Danish Rafiq
% Department of Electrical Engineering, University of Kashmir, India.
% Created on 30.10.2025 v1.0
% This file simulates the IEEE-118 Bus System using EN model
% and uses L-SINDy method to obtain reduced sparse model
%% Initialize the IEEE Power System (Load data)
clear;clc;close all
mpc=case2869pegase ;        % (see data file)
mpc.ref_freq=60;  % reference frequency
global data 
[data,details]=EN_model(mpc); %EN, SM (H,D,A,K,gamma,omega_R)
n_oc=length(data.H); %No of oscillators (FOM Size =2*n_oc)
dt=0.001;tf=5;
tspan=dt:dt:tf;
f= @(x) power_func(x);
xdotNL= @(t,x) f(x);
x0= [zeros(n_oc,1); zeros(n_oc,1)];
tic
[tFOM,xFOM]=ode45(xdotNL,tspan,x0);
%[tFOM,xFOM]=implicitEuler(xdotNL,tspan,x0);
FOM_time=toc
delta = xFOM(:,1:n_oc);
omega = xFOM(:,n_oc+1:end);
D=data.D;
H=data.H;
A=data.A;
gamma=data.gamma;
K=data.K;
omega_R=data.omega_R;
Mi=-D./(2*H);
Di=omega_R./(2*H).*A;
disp('Snapshots recorded')
 %% POD (for large data-sets)
 X=xFOM'; Opts.NoSV = 28; %3, 7, 8(15 norm)
 [U,S,V]=svd(X,'econ');
 [xdotr,x0r,Vpod] = POD(X, xdotNL,x0,Opts);
 tic
 [tFOM,xROM]=ode45(xdotr,tspan,x0r);
  ROM_time=toc
 x=xROM;
 dx=zeros(length(tFOM),Opts.NoSV);
 for i=1:length(xROM)
    dx(i,:)=xdotr(0,xROM(i,:)');
 end
disp('POD recorded')
sigma = diag(S);  % Extract singular values
figure(6);
semilogy(sigma, 'ko', 'LineWidth', 2);
xlabel('Index');
ylabel('Singular Values (log scale)');
title('Singular Value Decay Plot');
grid on;
sigma = diag(S);               % Extract singular values
energy = sigma.^2;             % Energy associated with each singular value
cumulative_energy = cumsum(energy) / sum(energy);  % Normalized cumulative energy
figure(7);
plot(cumulative_energy, 'ko', 'LineWidth', 2);
xlabel('Index');
ylabel('Cumulative Energy');
title('Singular Value Energy Plot');
grid on;
%% L-SINDy 
N=Opts.NoSV;
lambda=0.001; % 0.1=4gs
polyorder=1;
usesine=1;
Theta=poolData(x,N,polyorder,usesine);
Xi = sparsifyDynamics(Theta,dx,lambda,N); 
disp('Sindy solution found')
%% Reconstruction
xall2 = [];
dxall2 = [];
    [tB,xB]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,x0r);% approximate
xs=xB;
figure(1);
plot(x,'k','linewidth',2); hold on
plot(xs,'r--','LineWidth',2)
xlabel('$t [s]$','Interpreter','latex')
ylabel('Latent States','Interpreter','latex')
grid on
title(['IEEE', num2str(size(mpc.bus,1)), '-Bus System'],'Interpreter','latex')
set(gca,'ticklabelinterpreter','latex','Fontsize',15)
figure(2);
xFOM_pred=(Vpod*xB')';
plot(tspan,xFOM,'k')
hold on
plot(tspan,xFOM_pred,'r--','linewidth',1.5)
xlabel('$t [s]$','Interpreter','latex')
ylabel(' $\omega$  ','Interpreter','latex')
grid on
set(gca,'ticklabelinterpreter','latex','Fontsize',15)
%% predicted model
figure(3);
subplot(2,1,1)
delta_avg=sum(xFOM(:,1:n_oc),2)./n_oc;
delta_avg_pred=sum(xFOM_pred(:,1:n_oc),2)./n_oc;
plot(tspan,delta_avg,'k','linewidth',2)
hold on
plot(tspan,delta_avg_pred,'r--','linewidth',2)
xlabel('$t$ [s]','Interpreter','latex')
ylabel('Avg. $\delta$ [rad]','Interpreter','latex')
grid on
set(gca,'ticklabelinterpreter','latex','Fontsize',15)
subplot(2,1,2)
omg_avg=sum(xFOM(:,n_oc+1:end),2)./n_oc;
omg_avg_pred=sum(xFOM_pred(:,n_oc+1:end),2)./n_oc;
plot(tspan,omg_avg,'k','linewidth',2)
hold on
plot(tspan,omg_avg_pred,'r--','linewidth',2)
xlabel('$t[s]$','Interpreter','latex')
ylabel('Avg. $\Delta \omega [rad/s]$','Interpreter','latex')
grid on
set(gca,'ticklabelinterpreter','latex','Fontsize',15)
%% Plot Latent States
figure(5)
for i=1:28
    subplot(4,7,i)
    plot(tspan,x(:,i),'k','linewidth',2); hold on
    plot(tspan,xs(:,i),'r--','LineWidth',2)
    xlabel('$t$','Interpreter','latex')
    ylabel(['$ z $ ', num2str(i)],'Interpreter','latex')
    set(gca,'ticklabelinterpreter','latex','Fontsize',15)
    grid on
    xlim([0,5])
end
norm_delta=norm(delta_avg-delta_avg_pred)/norm(delta_avg)
norm_omega=norm(omg_avg-omg_avg_pred)/norm(omg_avg)
norm_Red = norm(x - xs)
 %% Error plot
figure(4)
semilogy(tspan, abs(delta_avg-delta_avg_pred),'r')
hold on
semilogy(tspan, abs(omg_avg-omg_avg_pred),'b')
legend('Delta','Omega')
grid on
grid on
set(gca,'ticklabelinterpreter','latex','Fontsize',15)
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$e(t)$','Interpreter','latex')