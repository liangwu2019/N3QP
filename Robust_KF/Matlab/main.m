clear, clc, close;
rng(20230914)
N = 1000;
lambda = 2;
x0 = 5*ones(3,1);
P0 = eye(3);
Q = diag([1/40, 1/10, 1/5]);
R = 0.5*eye(2);
%%
v = (randn(N,3)*chol(Q))'; % process noise

% Vector for full state history:
x = zeros(3,N+1);
x(:,1) = x0 + chol(P0)*randn(3,1);

u = zeros(1,N);
u(1,1:500) = 5*ones(1,500);
u(1,501:end) = 2*ones(1,500);

alpha1 = 0.1;
alpha2 = 0.5;
alpha3 = 0.2;
A = [1-alpha1, 0, 0; 
     0, 1-alpha2, 0; 
     alpha1, alpha2, 1-alpha3];
B = [0.5;0.5;0];

% Vector for full measurement history:
z = zeros(2,N+1);
H = [1 0 0; 0 0 1];
w = (randn(N+1,2)*chol(R))'; % measurement noise

% Additional noise term: Simulate sensor failures
y = zeros(N+1,2)';
numberOfOnes = N/20;

indexes = randperm(N);
for i = 1:numberOfOnes
   y(1,indexes(i)) = 25*randn;
end

indexes = randperm(N);
for i = 1:numberOfOnes
   y(2,indexes(i)) = 25*randn;
end
% Simulation
for k=2:N+1
    x(:,k) = A*x(:,k-1) + B*u(k-1) + v(:,k-1);
    z(:,k) = H*x(:,k) + w(:,k) + y(:,k);
end

%% KF
x_hat = zeros(3,N+1);
x_hat(:,1) = x0;
P_hat = P0;
R_inv = inv(R);
for k=2:N+1
    x_p = A*x_hat(:,k-1) + B*u(k-1);
    P_p = A*P_hat*A' + Q;
    
    P_hat = inv(inv(P_p)+H'*R_inv*H);
    x_hat(:,k) = x_p + P_hat*H'*R_inv*(z(:,k)-H*x_p);
end
%% Robust KF
x_hat_RKF = zeros(3,N+1);
x_hat_RKF(:,1) = x0;
P_hat_RKF = P0;
for k=2:N+1
    x_p_RKF = A*x_hat_RKF(:,k-1) + B*u(k-1);
    P_p_RKF = A*P_hat_RKF*A' + Q;

    % P_hat_RKF = P_p_RKF;
    P_hat_RKF = inv(inv(P_p_RKF)+H'*R_inv*H);

    L = P_hat_RKF*H'*inv(H*P_hat_RKF*H'+R);
    myQ = (eye(2)-H*L)'*R_inv*(eye(2)-H*L)+L'*inv(P_hat_RKF)*L;
    e_k = z(:,k)-H*x_p_RKF;

    
    H_inv = inv(myQ); H_inv = 0.5*(H_inv+H_inv');
    f = myQ*e_k;
    sol = N3_BoxQP(H_inv,-H_inv*f,1e-6,0.3,0.15);
    sol = myQ\(f-sol);

    x_hat_RKF(:,k) = x_p_RKF + L*(e_k-sol);
end
%%
RMS_KF = zeros(3,N+1);
for i=1:3
    RMS_KF(i,:) = sqrt((x(i,:)'-x_hat(i,:)').^2/N);
end

RMS_RKF = zeros(3,N+1);
for i=1:3
    RMS_RKF(i,:) = sqrt((x(i,:)'-x_hat_RKF(i,:)').^2/N);
end
imp = zeros(3,1);
for i = 1:3
    imp(i) = 100-100/mean(RMS_KF(i,:))*mean(RMS_RKF(i,:));
    fprintf('\nRMS error of state %i estimate is reduced by %.2f percent.\n',i,imp(i))
end
%%
figure(1)
subplot(311)
plot(0:N, RMS_KF(1,:),'r','LineWidth',1.2), title(['RMS Error state x(',int2str(1),')'], 'FontSize', 12)
hold on; plot(0:N, RMS_RKF(1,:),'b','LineWidth',1.2);
legend('KF','RKF');
xlim([250 500]);
xlabel('Time step t'),ylabel('RMS Error Maginute')
subplot(312)
plot(0:N, RMS_KF(2,:),'r','LineWidth',1.2), title(['RMS Error state x(',int2str(2),')'], 'FontSize', 12)
hold on; plot(0:N, RMS_RKF(2,:),'b','LineWidth',1.2);
legend('KF','RKF');
xlim([250 500]);
xlabel('Time step t'),ylabel('RMS Error Maginute')
subplot(313)
plot(0:N, RMS_KF(3,:),'r','LineWidth',1.2), title(['RMS Error state x(',int2str(3),')'], 'FontSize', 12)
hold on; plot(0:N, RMS_RKF(3,:),'b','LineWidth',1.2);
legend('KF','RKF');
xlim([250 500]);
xlabel('Time step t'),ylabel('RMS Error Maginute')


