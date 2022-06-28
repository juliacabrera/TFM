k = [0.8,0.5,0.3,0.2,0.5,0.8];
nk = length(k);

dt = 0.1;
a = 0.35;
c0vect = [0.5,0.45,0.40,0.38,0.34,0.30];
ivect = [0,3];
nl = 6;
iter = 800;

c_vector = zeros(iter,nl);
t_vector = zeros(iter,1);

t_vector(1) = 0.0;
c_vector(1,:) = c0vect;

for n = 1:nl
    
    c0 = c0vect(n);
    kk = k(n);
    
    c = c0;
    
    for i=2:iter
        cp = model0d(dt,c,a,kk);
        c = cp;
        c_vector(i,n) = c;
        t_vector(i) = i*dt;
    end
end

figure(1)
hold on
for xx = 1:nl
    plot(t_vector,c_vector(:,xx),'LineWidth',1.5)
end
yline(a,'LineWidth',1.5)
legend({'K=0.8, {c}_{0}=0.50','K=0.5, {c}_{0}=0.45','k=0.3, {c}_{0}=0.40','k=0.2, {c}_{0}=0.38','k=0.5, {c}_{0}=0.34','k=0.8, {c}_{0}=0.30','a=0.35'},'location','northeast');
hold off
xlabel('t')
ylabel('c')
ylim([0.0 1.0])
set(gca,'xtick',[])

function c_post = model0d(dt,c,a,k)
    c_post = c + dt*(k*c*(1.0-c)*(c-a));
end