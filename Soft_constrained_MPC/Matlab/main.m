%%
clear, close all, clc
Ts = 0.05;
NSim = 200;
TSim = NSim * Ts;
Np = 5; % prediction horizon

%% Define the continuous-time model
Ac = [-.0151 -60.5651 0 -32.174;
      -.0001 -1.3411 .9929 0;
      .00018 43.2541 -.86939 0;
      0      0       1      0];
Bc = [-2.516 -13.136;
      -.1689 -.2514;
      -17.251 -1.5766;
      0        0];
Cc = [0 1 0 0;
      0 0 0 1];
Dc = [0 0;
      0 0];
sys = ss(Ac,Bc,Cc,Dc);
model = c2d(sys,Ts);
[nx,nu] = size(model.B);

A_aug = [model.A, model.B; zeros(nu, nx), eye(nu)];
B_aug = [model.B; eye(nu)];
C_aug = [model.C, zeros(nu,nu)];

%% MPC-to-QP condense construction
Wy = 10*eye(2);
Wdu = 0.1*eye(2);

AiB = B_aug; BB = kron(eye(Np),AiB); 
for i=1:Np-1
    AiB = A_aug*AiB;
    BB = BB+kron(diag(ones(Np-i,1),-i),AiB); 
end

QQ = blkdiag(kron(eye(Np-1),C_aug'*Wy*C_aug),C_aug'*Wy*C_aug); 
RR = kron(eye(Np),Wdu); 
H =(BB'*QQ*BB + RR);

E_x = [0, 1, 0, 0, 0, 0;
       0, 0, 0, 1, 0, 0;
       0, 0, 0, 0, 1, 0;
       0, 0, 0, 0, 0, 1];

G = [kron(eye(Np),E_x)*BB; -kron(eye(Np),E_x)*BB];

%% Closed-loop simulation
ref = [0; 10];
% x = zeros(nx,1);
x = [0; 5; 0; 0]; % No feasible solution
u = zeros(nu,1);

Ref_Records = ref;
U_Records = [];
Y_Records = model.C*x;
Runtime_Records = [];
for k=1:NSim
    if(mod(k,100)==0)
        ref = [0; 0];
    end
    Ref_Records = [Ref_Records, ref];
    %% finish uncompleted MPC-to-QP construction
    x_aug = [x;u];
    ei = A_aug*x_aug;
    ee = ei;
    gg = E_x * ei;
    for i=2:Np
        ei = A_aug*ei;
        ee = [ee; ei];
        gg = [gg; E_x*ei];
    end
    h = BB'*(QQ*ee-repmat(C_aug'*Wy*ref,Np,1));
    g = [repmat([0.5;100;25;25],Np,1)-gg; gg - repmat([-0.5;-100;-25;-25],Np,1)];

    rho = 1000*[repmat([1;1;10;10],Np,1); repmat([1;1;10;10],Np,1)];
    M = G*inv(H)*G'; 
    d = G*inv(H)*h + g;
    H_dual = diag(rho)*M*diag(rho);
    h_dual = diag(rho)*(M*rho+2*d);
    tic
    dual_result = N3_BoxQP(H_dual,h_dual,1e-6,0.3,0.15);
    runtime = toc;
    Runtime_Records = [Runtime_Records; runtime];
    dU = -inv(H)*(h+0.5*G'*(rho.*dual_result+rho));
    
    %%
    u = u+dU(1:nu);
    x = model.A*x + model.B*u;
    U_Records = [U_Records, u];
    Y_Records = [Y_Records, model.C*x];
end
mean(Runtime_Records)
%% plotting
figure(1)
hold on
box on 
set(gca,'linewidth',2,'fontsize',12,'fontweight','bold')
plot(0:Ts:TSim,Y_Records(1,:),'r','linewidth',2),xlabel('time [s]'),ylabel('y_1')
plot(0:Ts:TSim,0.5*ones(NSim+1),'b--','linewidth',2)
plot(0:Ts:TSim,-0.5*ones(NSim+1),'b--','linewidth',2)
plot(0:Ts:TSim,zeros(NSim+1),'b','linewidth',2)


figure(2)
hold on 
box on 
set(gca,'linewidth',2,'fontsize',12,'fontweight','bold')
plot(0:Ts:TSim,Y_Records(2,:),'r','linewidth',2),xlabel('time [s]'),ylabel('y_2'),
plot(0:Ts:TSim,Ref_Records(2,:),'b','linewidth',2)

figure(3)
hold on
box on 
set(gca,'linewidth',2,'fontsize',12,'fontweight','bold')
plot(0:Ts:TSim-Ts,U_Records(1,:),'r','linewidth',2),ylim([-30,30]),xlabel('time [s]'),ylabel('u_1')

plot(0:Ts:TSim,25*ones(NSim+1),'b--','linewidth',2)
plot(0:Ts:TSim,-25*ones(NSim+1),'b--','linewidth',2)

figure(4)
hold on
box on 
set(gca,'linewidth',2,'fontsize',12,'fontweight','bold')
plot(0:Ts:TSim-Ts,U_Records(2,:),'r','linewidth',2),ylim([-30,30]), xlabel('time [s]'),ylabel('u_2')
plot(0:Ts:TSim,25*ones(NSim+1),'b--','linewidth',2)
plot(0:Ts:TSim,-25*ones(NSim+1),'b--','linewidth',2)
