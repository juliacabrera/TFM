
d = 50 ; % dimensió de la matriu (quadrada)
dx = 0.2 ; % Interespaiat entre elements de matriu (mida alvèol en mm)
dt = 0.0005 ; % Pas de temps (en hores). 
iter = 100000 ; % Iteracions de temps que volem calcular
D = 0.15e-3 ; % Coeficient de difusió (en mm/h) 
a = 0.16 ; % mateixes unitats que la concentració (g/mm3)
k = 0.5; % S'ha de comprovar amb el valor de l'article (en mm2/h g)
initial_concentration = 1.0 ; %(g/mm3)
tau = 4.5 ; % Proportionality factor for which we multiply k and D to adjust dyamics
xi = 10000;
df = xi*D;

D_les = zeros(d,d) ;
D_fib = zeros(d,d) ;
c_post = zeros(d,d) ;
g = zeros(d,d) ;
temps = zeros(iter) ; 
area = zeros(iter) ;
radi = zeros(iter) ;
dif_lesio_vector = zeros(iter);
dif_fib_vector = zeros(iter);
tfib = zeros(iter);
tc = zeros(iter);

final_radius_eps = 1.0;
final_encaps_eps = 0.2;

% LESIÓ

i0 = 7 ;
j0 = 20 ;

% Initial conditions
c = zeros(d,d) ; 
c0 = 7.0*ones(d,d) ;
eps = 2.0;

% Omplim la matriu c amb les distàncies al centre de la lesió
for i = 1:d
    for j = 1:d
        c(i,j) = (((i-i0)^2.0)+((j-j0)^2.0))^(1.0/2.0) ;
    end
end

c = 0.5*(1-tanh((c-c0)/eps)) ;
c(1,:) = 0.0 ;
c(d,:) = 0.0 ;
c(:,1) = 0.0 ;
c(:,d) = 0.0 ;

%FIBROBLAST

fib = zeros(d,d) ;

% Boundary conditions
fib(1,:) = 1.0 ;
fib(d,:) = 1.0 ;
fib(:,1) = 1.0 ;
fib(:,d) = 1.0 ;


% Diffusion coefficient matrix of the lesion
D_les = D;
% Diffusion coefficient matrix of the fibroblast
[gx,gy] = grad_matlab(c);
gradient_c = ((gx.*gx) + (gy.*gy).^0.5);
D_fib = df*gradient_c;

% Occupation matrix
g = 1.0 -c -fib;

fib0 = fib;

% v = VideoWriter('Encapsulation_2dmodel.mp4');
% open(v);
tic
% Temporal evolution 
for t=1:iter
  
    temps(t) = dt*(t-1);
    
    %Evolució de la lesió
    [c_post,dif_lesio] = lesion_evolution(d,D_les,a,k,dt,dx,tau,c,g);
    c = c_post;
    dif_lesio_vector(t) = max(max(dif_lesio(:)));
    
    % Occupation matrix
    g = 1.0 -c -fib ;
    
    % Diffusion coefficient matrix of the fibroblast
    [gx,gy] = grad_matlab(c);
    gradient_c = ((gx.*gx) + (gy.*gy).^0.5);
    D_fib = df*gradient_c;
   
    dtf = dt/max(max(gradient_c));
    
    %Evolució del fibroblast
    [fib_post,dif_fib] = diffusion(d,D_fib,dtf,dx,fib,g);
    fib = fib_post;
    dif_fib_vector(t) = sum(abs(dif_fib(:)));
    
    %Condicions de contorn per al fibroblast
    fib(1,:) = 1.0 ;
    fib(d,:) = 1.0 ;
    fib(:,1) = 1.0 ;
    fib(:,d) = 1.0 ;
   
    % Occupation matrix
    g = 1.0 -c -fib;
   
    area(t) = calcul_area(c) ;
    radi (t) = sqrt(area(t)/4*3.14) ;

    if t == 10000
        c1 = c;
        fib1 = fib;
    end
    
    if t == 40000
        c2 = c;
        fib2 = fib;
    end
    
    if t == 80000
        c3 = c;
        fib3 = fib;
    end
    
tfib(t) = sum(fib(:));
tc(t) = sum(c(:));

end


figure(5)
subplot(3,2,1)
imagesc(c1)
colorbar
axis equal
title('Lesion')
subplot(3,2,2)
imagesc(fib1)
colorbar
axis equal
title('Fibroblasts')
subplot(3,2,3)
imagesc(c2)
colorbar
axis equal
subplot(3,2,4)
imagesc(fib2)
colorbar
axis equal
subplot(3,2,5)
imagesc(c3)
colorbar
axis equal
subplot(3,2,6)
imagesc(fib3)
colorbar
axis equal


% close(v) ;

toc





function [m_post,DIF] = lesion_evolution(d,D,a,k,dt,dx,tau,c,g)
    m_post = zeros(d,d);
    ix = 1:d-1;
    DIF = zeros(d,d);
    D = D/tau;
    
    DIF(ix,:) = DIF(ix,:) + g(ix,:).*c(ix+1,:) - g(ix+1,:).*c(ix,:);
    DIF(ix+1,:) = DIF(ix+1,:) -g(ix,:).*c(ix+1,:) + g(ix+1,:).*c(ix,:);
    DIF(:,ix) = DIF(:,ix) -g(:,ix+1).*c(:,ix) + g(:,ix).*c(:,ix+1);
    DIF(:,ix+1)= DIF(:,ix+1) +g(:,ix+1).*c(:,ix)  -g(:,ix).*c(:,ix+1);
    
    m_post = c+ dt*((k/tau)*c.*(1-c).*(c-a) + (DIF.*D)/(dx^2.0));
end



function [dif_post,DIF] = diffusion(d,D,dt,dx,fib,g)
    dif_post = zeros(d,d) ;
    ix = 1:d-1;
    DIF = zeros(d,d);
    
    DIF(ix,:) = DIF(ix,:) + g(ix,:).*fib(ix+1,:).*D(ix+1,:) - g(ix+1,:).*fib(ix,:).*D(ix,:);
    DIF(ix+1,:) = DIF(ix+1,:) -g(ix,:).*fib(ix+1,:).*D(ix+1,:) + g(ix+1,:).*fib(ix,:).*D(ix,:);
    DIF(:,ix) = DIF(:,ix) -g(:,ix+1).*fib(:,ix).*D(:,ix) + g(:,ix).*fib(:,ix+1).*D(:,ix+1);
    DIF(:,ix+1)= DIF(:,ix+1) +g(:,ix+1).*fib(:,ix).*D(:,ix)  -g(:,ix).*fib(:,ix+1).*D(:,ix+1);
    
    dif_post = fib +dt*(DIF/(dx^2.0));
end

function grad = gradient_lesio(c,d,dx)
    grad = zeros(d,d) ;
    grad(2:d-1,2:d-1) = ((((c(3:d,2:d-1) - c(1:d-2,2:d-1))/(2*dx)).^2.0) + (((c(2:d-1,3:d) - c(2:d-1,1:d-2))/(2*dx)).^2.0)).^(0.5);
    grad(1,:) = (c(2,:)-c(1,:))/dx;
    grad(d,:) = (c(d-1,:)-c(d,:))/dx;
    grad(:,1) = (c(:,2)-c(:,1))/dx;
    grad(:,d) = (c(:,d-1)-c(:,d))/dx;
end

function [gradx,grady] = grad_matlab(c)
    [gradx,grady] = gradient(c);

end

function A = calcul_area(m)
    A = 200*200*sum(m(m>0.95)); %àrea en micròmetres cúbics
   
end

function [rad] = final_radius(r,radi_lim)
    i = 1;
    while (r(i+1)-r(i))>radi_lim
        i = i+1;
    end
    rad = r(i);
    disp(i)
    disp('radi')
    disp(rad)
end

function [t0,ii] = inici_encapsulacio(mfib,ddt)
   i = 1;
   while mfib == 196
       i = i+1;
   end
   ii =  i;
   t0 = (i-1)*ddt;
end

function d = distancia_fib(x0,y0)
    dd = [x0,y0,50-x0,50-y0];
    d = min(dd)*0.2;
end

function tf = end_encapsulation(diffib, eps_dif_fib,ddt)
    i = 1;
    while diffib(i+1)-diffib(i) > eps_dif_fib
        i = i+1;
    end
    tf = ddt*(i-1);
    disp('radi')
    disp(i)
    disp(tf)
    disp(ddt)
end

function tf = end_encapsulation2(diffib, eps_dif_fib,ddt)
    i = 1;
    while diffib(i) > eps_dif_fib
        i = i+1;
    end
    tf = ddt*(i-1);
    disp('radi')
    disp(i)
    disp(tf)
    disp(ddt)
end

function r = finrad(c,thresh)
    c_filt = zeros(size(c));
    c_filt(c>thresh)=1.0;
    area = sum(c_filt(:));
    r = area*200*200;
%     r = sqrt(area/3.14);
%     r = r/1000.0;
end
