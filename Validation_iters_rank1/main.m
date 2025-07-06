clear,close,clc
n = 100;
num_cond = 1e3;
v = [num_cond;(num_cond-1)*rand(n-2,1)+1;1];
U = orth(randn(n,n));
H = U*diag(v)*U';
H = (H+H')/2;
h = 1e3 * (1-2*rand(n,1));

%%
alpha = 0.3;
delta = 0.15;
[z_N3,Iters,Ranks,Gap] = N3QP(H,h,1e-8,alpha,delta);

%%
figure(1)
semilogy(Iters,Gap,'b--','LineWidth',2)
xlabel('Number of Iterations')
ylabel('Duality Gap')
epsilon = 1e-6;
sigma = sqrt(2)*delta*(1+delta)^2*alpha*sqrt((1+alpha)/(1-alpha))+0.5*(1+delta)^2*alpha^2/(1-alpha);
beta = (alpha-sigma)/(1+alpha/sqrt(2*n));
decrease_rate = 1.0-beta/sqrt(2*n);
theory_iters = ceil(-log((2*n+alpha*sqrt(2*n))/epsilon)/log(decrease_rate));
xline(theory_iters,'r--','LineWidth', 2)
%%
hold on
epsilon = 1e-8;
sigma = sqrt(2)*delta*(1+delta)^2*alpha*sqrt((1+alpha)/(1-alpha))+0.5*(1+delta)^2*alpha^2/(1-alpha);
beta = (alpha-sigma)/(1+alpha/sqrt(2*n));
decrease_rate = 1.0-beta/sqrt(2*n);
theory_iters = ceil(-log((2*n+alpha*sqrt(2*n))/epsilon)/log(decrease_rate));
xline(theory_iters,'g--','LineWidth', 2)



%%
figure(2)
semilogy(Ranks,Gap,'b--','LineWidth',2)
xlabel('Number of Rank-1 updates')
ylabel('Duality Gap')
epsilon = 1e-6;
sigma = sqrt(2)*delta*(1+delta)^2*alpha*sqrt((1+alpha)/(1-alpha))+0.5*(1+delta)^2*alpha^2/(1-alpha);
beta = (alpha-sigma)/(1+alpha/sqrt(2*n));
decrease_rate = 1.0-beta/sqrt(2*n);
theory_iters = ceil(-log((2*n+alpha*sqrt(2*n))/epsilon)/log(decrease_rate));
eta = (1+delta)^3*alpha/(1-alpha);
theory_ranks = ceil(4*eta*(theory_iters-1)*sqrt(n)/(1-eta)*log(1+delta));
xline(theory_ranks,'r--','LineWidth', 2)
%%
hold on
epsilon = 1e-8;
sigma = sqrt(2)*delta*(1+delta)^2*alpha*sqrt((1+alpha)/(1-alpha))+0.5*(1+delta)^2*alpha^2/(1-alpha);
beta = (alpha-sigma)/(1+alpha/sqrt(2*n));
decrease_rate = 1.0-beta/sqrt(2*n);
theory_iters = ceil(-log((2*n+alpha*sqrt(2*n))/epsilon)/log(decrease_rate));
eta = (1+delta)^3*alpha/(1-alpha);
theory_ranks = 4*eta*(theory_iters-1)*sqrt(n)/(1-eta)*log(1+delta);
xline(theory_ranks,'g--','LineWidth', 2)

%%
n = 200;
num_cond = 1e3;
v = [num_cond;(num_cond-1)*rand(n-2,1)+1;1];
U = orth(randn(n,n));
H = U*diag(v)*U';
H = (H+H')/2;
h = 1e3 * (1-2*rand(n,1));

%%
alpha = 0.3;
delta = 0.15;
[z_N3,Iters,Ranks,Gap] = N3QP(H,h,1e-8,alpha,delta);

%%
figure(1)
hold on
semilogy(Iters,Gap,'b-','LineWidth',2)
epsilon = 1e-6;
sigma = sqrt(2)*delta*(1+delta)^2*alpha*sqrt((1+alpha)/(1-alpha))+0.5*(1+delta)^2*alpha^2/(1-alpha);
beta = (alpha-sigma)/(1+alpha/sqrt(2*n));
decrease_rate = 1.0-beta/sqrt(2*n);
theory_iters = ceil(-log((2*n+alpha*sqrt(2*n))/epsilon)/log(decrease_rate));
xline(theory_iters,'r-','LineWidth', 2)
%%
hold on
epsilon = 1e-8;
sigma = sqrt(2)*delta*(1+delta)^2*alpha*sqrt((1+alpha)/(1-alpha))+0.5*(1+delta)^2*alpha^2/(1-alpha);
beta = (alpha-sigma)/(1+alpha/sqrt(2*n));
decrease_rate = 1.0-beta/sqrt(2*n);
theory_iters = ceil(-log((2*n+alpha*sqrt(2*n))/epsilon)/log(decrease_rate));
xline(theory_iters,'g-','LineWidth', 2)
yline(1e-6,'k--','LineWidth', 2)
yline(1e-8,'k--','LineWidth', 2)
legend('Practical (n=100)','Theory Bound for Iterations (n=100,\epsilon=10^{-6})','Theory Bound for Iterations (n=100,\epsilon=10^{-8})','Practical (n=200)','Theory Bound for Iterations (n=200,\epsilon=10^{-6})','Theory Bound for Iterations (n=200,\epsilon=10^{-8})')
ylim([1e-9,1e3])
%%
figure(2)
hold on
semilogy(Ranks,Gap,'b-','LineWidth',2)
epsilon = 1e-6;
sigma = sqrt(2)*delta*(1+delta)^2*alpha*sqrt((1+alpha)/(1-alpha))+0.5*(1+delta)^2*alpha^2/(1-alpha);
beta = (alpha-sigma)/(1+alpha/sqrt(2*n));
decrease_rate = 1.0-beta/sqrt(2*n);
theory_iters = ceil(-log((2*n+alpha*sqrt(2*n))/epsilon)/log(decrease_rate));
eta = (1+delta)^3*alpha/(1-alpha);
theory_ranks = 4*eta*(theory_iters-1)*sqrt(n)/(1-eta)*log(1+delta);
xline(theory_ranks,'r-','LineWidth', 2)
%%
hold on
epsilon = 1e-8;
sigma = sqrt(2)*delta*(1+delta)^2*alpha*sqrt((1+alpha)/(1-alpha))+0.5*(1+delta)^2*alpha^2/(1-alpha);
beta = (alpha-sigma)/(1+alpha/sqrt(2*n));
decrease_rate = 1.0-beta/sqrt(2*n);
theory_iters = ceil(-log((2*n+alpha*sqrt(2*n))/epsilon)/log(decrease_rate));
eta = (1+delta)^3*alpha/(1-alpha);
theory_ranks = 4*eta*(theory_iters-1)*sqrt(n)/(1-eta)*log(1+delta);
xline(theory_ranks,'g-','LineWidth', 2)
yline(1e-6,'k--','LineWidth', 2)
yline(1e-8,'k--','LineWidth', 2)
legend('Practical (n=100)','Theory Bound for Rank-1 updates (n=100,\epsilon=10^{-6})','Theory Bound for Rank-1 updates (n=100,\epsilon=10^{-8})','Practical (n=200)','Theory Bound for Rank-1 updates (n=200,\epsilon=10^{-6})','Theory Bound for Rank-1 updates (n=200,\epsilon=10^{-8})')
ylim([1e-9,1e3])
