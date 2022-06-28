close all
clear all

% Ranges for the parameters
xi = [10000 30000];
D = [0.15e-3 0.5e-3];
dist0 = [7 10]; 

Ns = 100; % NUmber of points of the LAtin hypercube sample

LIM = [xi' D' dist0'];
DX = (LIM(2,:)-LIM(1,:))/Ns; %matriu de punts de cada parametre
Np = length(DX);
LHS = zeros(Ns,Np);
for i = 1:Ns
    for j = 1:Np
        LHS(i,j) = LIM(1,j) + (i-1+rand)*DX(j);
    end
end
for j=1:Np
    LHS(:,j)=LHS(randperm(Ns),j);
end

fradi = zeros(Ns,1);
tencaps = zeros(Ns,1);
errors = [];

for k = 1:Ns
    tic
    xxi = LHS(k,1);
    DD = LHS(k,2);
    ddist0 = LHS(k,3);
    try
        disp(xxi);
        disp(DD);
        disp(ddist0);
        [fradi(k),tencaps(k)] = model2d(xxi,DD,ddist0);
        
    catch
       [fradi(k),tencaps(k)] = model2d(3000,0.9e-3,randi([5 15],1));
       disp('error')
    end
    toc
end

fradi = fradi(fradi~=0.0);
tencaps = tencaps (tencaps~= 0.0);

mod_LHS = LHS;
final_LHS = LHS;

if isempty(errors) == false
    for i = 1,length(errors)
        mod_LHS(errors(i),:) = [];
        final_LHS = intersect(final_LHS,mod_LHS);
        mod_LHS = LHS;
    end
end



PRCC_fradi = zeros(1,Np);
PRCC_tencaps = zeros(1,Np);

for k = 1:Np
    XX = [fradi LHS(:,k)];
    ZZ = LHS; ZZ(:,k)=[];
    [rho,~]=partialcorr(XX,ZZ,'type','Pearson');
    PRCC_fradi(k) = rho(1,2);
    XX = [tencaps LHS(:,k)];
    ZZ = LHS; ZZ(:,k)=[];
    [rho,~]=partialcorr(XX,ZZ,'type','Pearson');
    PRCC_tencaps(k) = rho(1,2);
end

param = {'\chi','D','d'};

f1 = figure(1);
clf
[s,i] = sort(PRCC_fradi);
bar(1:Np,s)
axis([0 Np+1 -1 1])
ax=gca;
ax.XTick = 1:Np;
ax.XTickLabel = param(i);
ylabel('PRCC with {r}_{f}')
exportgraphics(f1,'FD_SA_rf_provaexp.png')

f2 = figure(2);
clf
[s,i] = sort(PRCC_tencaps);
bar(1:Np,s)
axis([0 Np+1 -1 1])
ax=gca;
ax.XTick = 1:Np;
ax.XTickLabel = param(i);
ylabel('PRCC with {t}_{encaps}')
exportgraphics(f2,'FD_SA_tenc_provaexp.png')

function [r,t_enc] = model2d(xi,D,dist)
    d = 50 ; % dimensió de la matriu (quadrada)
    dx = 0.2 ; % Interespaiat entre elements de matriu (mida alvèol en mm)
    dt = 0.0005 ; % Pas de temps (en hores). Normalment a 0.1
    iter = 100000 ; % Iteracions de temps que volem calcular
    a = 0.16 ; % mateixes unitats que la concentració (g/mm3)
    tau = 4.5 ; % Proportionality factor for which we multiply k and D to adjust dyamics
    df = xi*0.7e-3;
    k = 0.5;

    D_les = zeros(d,d) ;
    D_fib = zeros(d,d) ;
    c_post = zeros(d,d) ;
    g = zeros(d,d) ;
    temps = zeros(iter) ; % Guardem els temps en un vector;
    radi = zeros(iter) ;
    dif_lesio_vector = zeros(iter);
    dif_fib_vector = zeros(iter);
    tfib = zeros(iter);
    
    final_radius_eps = 5.0e-4;
    final_encaps_eps = 0.1;
    
    randim = randi([1 2],1);
    if randim == 1
        i0 = round(dist);
        j0 = 10;
    elseif randim == 2
        j0 = round(dist);
        i0 = 10;
    end
    
    % Initial conditions
    c = zeros(d,d) ; % les concentracions estaran en g/mm^3
    c0 = 10.0*ones(d,d) ;
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

    fib = zeros(d,d) ;

    % condicions de contorn 
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

    

    % figure(1)
    % surf(gradient_c)

    fib0 = fib;

    % v = VideoWriter('Ev_6plots_25_25_5temps.mp4');
    % open(v);
    tic
    % Comencem el bucle per calcular l'evolució temporal
    for t=1:iter

        temps(t) = dt*(t-1);

        %Evolució de la lesió
        [c_post,dif_lesio] = lesion_evolution(d,D_les,a,k,dt,dx,tau,c,g);
        c = c_post;
        dif_lesio_vector(t) = sum(dif_lesio(:));

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
%         if max(max(g)) > 1.0
%             disp(k)
%             disp(xi)
%             disp(D)
%             disp(dist)
%         end

        % Comprovem si la infecció ha tocat la frontera
        % Si no la ha tocat calculem l'àrea de la infecció

        area (t)= calcul_area(c) ;
        radi (t) = sqrt(area(t)/4*3.14) ;
        
        tfib(t) = sum(fib(:));
        
    end

%     r = final_radius(radi,final_radius_eps,dt);
    r = finrad(c,0.001);
    [t0_encaps,index] = inici_encapsulacio(fib,dt);
    tf_encaps = end_encapsulation2(dif_fib_vector(index+1:iter),final_encaps_eps,dt);
    
    t_enc = tf_encaps-t0_encaps;  
end


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

function [gradx,grady] = grad_matlab(c)
    [gradx,grady] = gradient(c);

end

function A = calcul_area(m)
    A = 200*200*sum(m(:)); %àrea en micròmetres cúbics
   
end

function rad = final_radius(r,radi_lim,ddt)
    i = 1;
    while (r(i+1)-r(i))>radi_lim
        i = i+1;
    end
    rad = r(i);
end

function [t0,ii] = inici_encapsulacio(mfib,ddt)
   i = 1;
   while mfib == 196
       i = i+1;
   end
   ii = i;
   t0 = (i-1)*ddt;
end

function d = distancia_fib(x0,y0)
    dd = [x0,y0,50-x0,50-y0];
    d = min(dd)*0.2;
end

function tf = end_encapsulation(diffib, eps_dif_fib,ddt)
    i = 1;
    difference = diffib(i+1)-diffib(i);
    while difference > eps_dif_fib
        i = i+1;
        difference = diffib(i+1)-diffib(i);
        if i == 99997
            difference = 0.5e-4;
        end
    end
    tf = ddt*(i-1);
%     disp('radi')
%     disp(i)
%     disp(tf)
%     disp(ddt)
end


function tf = end_encapsulation2(diffib, eps_dif_fib,ddt)
    i = 1;
    while diffib(i) > eps_dif_fib
        i = i+1;
        if i == 99998
            diffib(i) = 0.6e-4;
        end
    end
%     tf = ddt*(i-1);
    tf = i;
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
    
end
