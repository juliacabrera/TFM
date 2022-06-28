d = 50 ; % dimensió de la matriu (quadrada)
dx = 0.2 ; % Interespaiat entre elements de matriu (mida alvèol en mm)
dt = 0.0005 ; % Pas de temps (en hores). (0.0005)
iter = 15000 ; % Iteracions de temps que volem calcular
D = 5.0e-2; % Coeficient de difusió (en mm2/h) (article baseline value = 1.0e-3)(1.5e-2)
a = 0.12; % mateixes unitats que la concentració (g/mm3) (0.13)
k = 0.5 ; % S'ha de comprovar amb el valor de l'art4cle (en mm3/h g)?
tau = 4.5 ; % Proportionality factor for which we multiply k and D to adjust dyamics (4.5)

D_les = zeros(d,d,d) ;
D_fib = zeros(d,d,d) ;
c_post = zeros(d,d,d) ;
g = zeros(d,d,d) ;
temps = zeros(iter) ; % Guardem els temps en un vector
dif_lesio_vector = zeros(iter);
fif_fib_vector = zeros(iter);
%g_vector = zeros(iter);
radi = zeros(iter);


% LESIÓ
i_0 = 20 ;
j_0 = 4 ;
k_0 = 25 ;

% Initial conditions
c = zeros(d,d,d) ; % les concentracions estaran en g/mm^3
c0 = 5.0*ones(d,d,d) ;
eps = 2.0;

% Omplim la matriu c amb les distàncies al centre de la lesió
for k = 1:d
    for i = 1:d
        for j = 1:d
             c(i,j,k) = (((i-i_0)^2.0)+((j-j_0)^2.0) + ((k-k_0)^2.0))^(1.0/2.0) ;           
        end
    end
end



c = 0.5*(1.0-tanh((c-c0)/eps));




[x,y,z] = meshgrid([1:1:d]); %falta algo aixi i amb xx yy zz



%FIBROBLAST
fib = get_irregular_fibroblast(d);
fib_0 = fib;

c(fib_0 == 1.0) = 0.0;

fib2d = zeros(d,d);
fib2d(:) = c(:,:,25);

figure(1)
imagesc(fib2d)

% Diffusion coefficient matrix of the lesion
D_les = D;
% Diffusion coefficient matrix of the fibroblast
[gx,gy,gz] = grad_matlab(c);
gradient_c = ((gx.*gx) + (gy.*gy) + (gz.*gz)).^0.5;
D_fib = 300*D.*gradient_c;

c2d = zeros(d,d) ;
c2d(:) = c(:,:,25) ;

% figure(1)
% imagesc(c2d)
% colorbar

c00 = zeros(d,d); c00(:) = c(:,:,30);
fib00 = zeros(d,d); fib00(:) = fib(:,:,30);



% grad_dif_zero = zeros(50,50);
% gradc2d = zeros(50,50);
% gradc2d(:) = D_fib(30,:,:);
% grad_dif_zero(gradc2d > 0.5) = 1.0;


% Occupation matrix
g = 1.0 -c -fib;

dt = dt/max(max(max(gradient_c)));

tfib = zeros(1,iter);
tc = zeros(1,iter);

 v = VideoWriter('Ev_lesio_3d_irrfib.mp4');
 open(v);

% Comencem el bucle per calcular l'evolució temporal
for t=1:iter
  
    temps(t) = dt*(t-1);
    
    %Evolució de la lesió
    [c_post,dif_lesio] = lesion_evolution3d(d,D_les,a,k,dt,dx,tau,c,fib,g);
    c = c_post;
    dif_lesio_vector(t) = sum(sum(sum(abs(dif_lesio))));
    radi(t) = radi_lesio(c)/1000.0; %passem el radi a mm
    
%     c(fib_0 == 1.0) = 0.0;
    
    
    % Occupation matrix
    g = 1.0 -c -fib ;
    
    % Diffusion coefficient matrix of the fibroblast
    [gx,gy,gz] = grad_matlab(c);
    gradient_c = ((gx.*gx) + (gy.*gy) + (gz.*gz)).^0.5;
    D_fib = 200.0*D.*gradient_c;
   
    
    %Evolució del fibroblast
    [fib_post,dif_fib] = diffusion3d(d,D_fib,dt,dx,fib,g);
    fib = fib_post;
    dif_fib_vector(t) = sum(sum(sum(abs(dif_fib))));
    
    
    %Condicions de contorn per al fibroblast
    fib(fib_0 == 1.0) = 1.0;
   
    % Occupation matrix
    g = 1.0 -c -fib;
    if max(max(max(g))) > 1.0
        disp('aqui!')
    end
    
    if t == 2000
        c01 = zeros(d,d); c01(:) = c(:,:,30);
        fib01 = zeros(d,d); fib01(:) = fib(:,:,30);
    end   
    
    if t == 8000
       c02 = zeros(d,d); c02(:) = c(:,:,30);
       fib02 = zeros(d,d); fib02(:) = fib(:,:,30); 
    end
    
    if mod(t,500) == 0.0
            disp('iteracions')
            disp(t)
        
%           fsurf = isosurface(x,y,z,fib,0.1);
%           csurf = isosurface(x,y,z,c,0.1);
%           figure(1)
%           clf
%           patch(fsurf,'FaceColor','g','FaceAlpha',0.2,'EdgeColor','none')
%           hold on
%           patch(csurf,'FaceColor','r','FaceAlpha',0.2,'EdgeColor','none')
%           hold off
%           
%           
%           title(num2str(t))
%           axis equal
%           axis([0,50,0,50,0,50])
%           view([45 45])
%           drawnow;

%         c1 = zeros(d,d); c1(:) = c(1,:,:);
%         c10 = zeros(d,d); c10(:) = c(10,:,:);
%         c20 = zeros(d,d); c20(:) = c(20,:,:);
        c30 = zeros(d,d); c30(:) = c(:,:,30);
%         c40 = zeros(d,d); c40(:) = c(40,:,:);
%         c50 = zeros(d,d); c50(:) = c(50,:,:);
        
%         f1 = zeros(d,d); f1(:) = fib(1,:,:);
%         f10 = zeros(d,d); f10(:) = fib(10,:,:);
%         f20 = zeros(d,d); f20(:) = fib(20,:,:);
        f30 = zeros(d,d); f30(:) = fib(:,:,30);
%         f40 = zeros(d,d); f40(:) = fib(40,:,:);
%         f50 = zeros(d,d); f50(:) = fib(50,:,:);
% 
        g30 = zeros(50,50); g30(:) = g(:,:,30);
        dif_les30 = zeros(50,50); dif_les30(:) = dif_lesio(:,:,30);
        dif_fib30 = zeros(50,50); dif_fib30(:) = dif_fib(:,:,30);
        D_fib30 = zeros(50,50); D_fib30(:) = D_fib(:,:,30);
        
        
        
%         grad_dif_zero = zeros(50,50);
%         gradc2d = zeros(50,50);
%         gradc2d(:) = gradient_c(2,:,:);
%         grad_dif_zero(gradc2d ~= 0.0) = 1.0;
%     
% 
%         figure(2)
       
%         imagesc(c30)
%         colorbar
%         drawnow
%         suptitle(num2str(t))
         
%         subplot(2,3,1)
%         imagesc(dif_les30)
%         colorbar
%         title('DIfusio lesio')
%          
%         subplot(2,3,2)
%         imagesc(dif_fib30)
%         colorbar
%         title('difusio fibroblast')
%          
%          
%         subplot(2,3,3)
%         imagesc(D_fib30)
%         colorbar
%         title('coeficient fibroblast')
%        
%         subplot(2,3,4)
%         imagesc(g30)
%         colorbar
%         title('ocupacio')
%          
%         subplot(2,3,5)
%         imagesc(c30)
%         colorbar
%         title('lesio')
%          
%         subplot(2,3,6)
%         imagesc(f30)
%         colorbar
%         title('fibroblast')
%         suptitle(num2str(t))
        
%         figure(3)
%         grad2d = zeros(50,50);
%         grad2d(:) = gradient_c(30,:,:);
%         imagesc(grad2d)
%         colorbar
%         title('gradient a 30')
%         
%         
%         subplot(2,3,6)
%         imagesc(f30)
%         colorbar
%         title('lesio 40')
%         suptitle(num2str(t))
        
        drawnow;
        

%         figure(2)
%         isosurface(xx,yy,zz,fib(2:d-1,2:d-1,2:d-1),0.9);
        
%         isosurface(xx,yy,zz,fib(2:d-1,2:d-1,2:d-1));
%         drawnow;
%         writeVideo(v,getframe(gcf));
        
        
       
    end
% if (sum(c(:)<0))>0
%     a = 1.0;
% end
tfib(t) = sum(fib(:));
tc(t) = sum(c(:));
end

figure(2)

subplot(3,2,1)
imagesc(c00)
colorbar
title('Lesions')

subplot(3,2,2)
imagesc(fib00)
colorbar
title('Fibroblasts')

subplot(3,2,3)
imagesc(c01)
colorbar

subplot(3,2,4)
imagesc(fib01)
colorbar

subplot(3,2,5)
imagesc(c02)
colorbar

subplot(3,2,6)
imagesc(fib02)
colorbar

%  close(v) ;
% 
% figure; plot(temps,tfib(1:length(temps)))
% title('Fibroblast')
% 
% figure; plot(temps,tc(1:length(temps)))
% title('Lesio')
% 
% figure; plot(temps/24.0,radi)
% xlabel('temps(dies)')
% ylabel('radi(mm)')
%  
%  
%  
% figure; semilogy(temps,dif_lesio_vector)
% suptitle('Lesió')
%  
% figure; semilogy(temps,dif_fib_vector(1:length(temps)))
% suptitle('Fibroblast')



function [m_post,DIF] = lesion_evolution3d(d,D,a,k,dt,dx,tau,c,m_fib,g)
    m_post = zeros(d,d,d);
    ix = 1:d-1;
    DIF = zeros(d,d,d);
    D = D/tau;
    
    DIF(ix,:,:) = DIF(ix,:,:) + g(ix,:,:).*c(ix+1,:,:) - g(ix+1,:,:).*c(ix,:,:);
    DIF(ix+1,:,:) = DIF(ix+1,:,:) -g(ix,:,:).*c(ix+1,:,:) + g(ix+1,:,:).*c(ix,:,:);
    DIF(:,ix,:) = DIF(:,ix,:) -g(:,ix+1,:).*c(:,ix,:) + g(:,ix,:).*c(:,ix+1,:);
    DIF(:,ix+1,:)= DIF(:,ix+1,:) +g(:,ix+1,:).*c(:,ix,:)  -g(:,ix,:).*c(:,ix+1,:);
    DIF(:,:,ix) = DIF(:,:,ix) -g(:,:,ix+1).*c(:,:,ix) + g(:,:,ix).*c(:,:,ix+1);
    DIF(:,:,ix+1)= DIF(:,:,ix+1) +g(:,:,ix+1).*c(:,:,ix)  -g(:,:,ix).*c(:,:,ix+1);
    
    m_post = c + dt*((k/tau)*c.*(1.0-m_fib-c).*(c-a) + (DIF.*D)/(dx^2.0));
end



function [dif_post,DIF] = diffusion3d(d,D,dt,dx,fib,g)
    dif_post = zeros(d,d,d) ;
    ix = 1:d-1;
    DIF = zeros(d,d,d);
    
    DIF(ix,:,:) = DIF(ix,:,:) + g(ix,:,:).*fib(ix+1,:,:).*D(ix+1,:,:) - g(ix+1,:,:).*fib(ix,:,:).*D(ix,:,:);
    DIF(ix+1,:,:) = DIF(ix+1,:,:) -g(ix,:,:).*fib(ix+1,:,:).*D(ix+1,:,:) + g(ix+1,:,:).*fib(ix,:,:).*D(ix,:,:);
    DIF(:,ix,:) = DIF(:,ix,:) -g(:,ix+1,:).*fib(:,ix,:).*D(:,ix,:) + g(:,ix,:).*fib(:,ix+1,:).*D(:,ix+1,:);
    DIF(:,ix+1,:)= DIF(:,ix+1,:) +g(:,ix+1,:).*fib(:,ix,:).*D(:,ix,:)  -g(:,ix,:).*fib(:,ix+1,:).*D(:,ix+1,:);
    DIF(:,:,ix) = DIF(:,:,ix) -g(:,:,ix+1).*fib(:,:,ix).*D(:,:,ix) + g(:,:,ix).*fib(:,:,ix+1).*D(:,:,ix+1);
    DIF(:,:,ix+1)= DIF(:,:,ix+1) +g(:,:,ix+1).*fib(:,:,ix).*D(:,:,ix)  -g(:,:,ix).*fib(:,:,ix+1).*D(:,:,ix+1);
     
    dif_post = fib + dt*(DIF/(dx^2.0));
end



function [gradx,grady,gradz] = grad_matlab(c)
    [gradx,grady,gradz] = gradient(c);
end

function r = radi_lesio(c)
    vol = 200*200*200*sum(c(:)); % volum en micrometres cubics
    r = (vol*(3.0/4.0)*(1.0/pi))^(1.0/3.0); %radi en micrometres
    
end

    
    
   
function fibr = get_irregular_fibroblast(d)
    fib = zeros(d,d,d) ;
    
    vect1 = [1,0,0] ;
    vect2 = [vect1;[-1,0,0]] ;
    vect3 = [vect2;[0,1,0]] ;
    vect4 = [vect3;[0,-1,0]];
    vect5 = [vect4;[0,0,1]] ;
    vect6 = [vect5;[0,0,-1]] ;
    vect7 = [vect6;[-1,-1,-1]];
    vect8 = [vect7;[1,-1,-1]];
    vect9 = [vect8;[-1,1,-1]];
    vect10 = [vect9;[1,1,-1]];
    vect11 = [vect10;[1,1,1]];
    vect12 = [vect11;[1,-1,1]];
    vect13 = [vect12;[-1,1,1]];
    vect14 = [vect13;[-1,-1,1]];

    vects = vect14;
    
    p1 = [1,25,25] ;
    p2 = [p1;[50,25,25]] ;
    p3 = [p2;[25,1,25]] ;
    p4 = [p3;[25,50,25]] ;
    p5 = [p4;[25,25,1]] ;
    p6 = [p5;[25,25,50]] ;
    p7 = [p6;[45,45,45]] ; %mirant des de l'eix z punta de baix a la dreta
    p8 = [p7;[2,30,45]] ; % mirant des de l'eix z punta de dalt a la dreta
    p9 = [p8;[30,2,45]] ; % mirant des de l'eix z punta de baix a l'esquerra
    p10 = [p9;[30,2,45]] ; % mirant des de l'eix z punts de dalt a l'esquerra
    p11 = [p10;[2,2,5]] ;
    p12 = [p11;[2,45,5]] ;
    p13 = [p12;[45,2,5]] ;
    p14 = [p13;[45,45,5]] ;


    points = p14;
    
    n = size(vects);

    for i = 1:n(1)
        fib = generar_pla(fib,vects(i,:),points(i,:));
    end

    fib_cavitat = generar_matriu_fib(fib,vects,points);
    
    fibr = zeros(size(fib_cavitat));
    fibr(fib_cavitat ~= 1.0) = 1.0; 
    
end


function mat_pla = generar_pla(mat,vd,p)
    [d1,d2,d3] = size(mat) ;
    d = -(vd(1)*p(1) + vd(2)*p(2) + vd(3)*p(3)) ;
    for i = 1:d1
        for j = 1:d2
            for k = 1:d3
                val_pla = vd(1)*i + vd(2)*j + vd(3)*k + d;
                if val_pla == 0.0
                     mat(i,j,k) = 1.0 ;
                end
            end
        end
    end
    mat_pla = mat ;
end


function mat_fib = generar_matriu_fib(mat,vectors_plans,punts_plans)
    n = size(vectors_plans) ;
    m = size(mat) ;
    p = punts_plans ;
    mat_fib = zeros(m(1),m(2),m(3)) ;
    
    for i = 1:m(1)
        for j = 1:m(2)
            for k = 1:m(3)
                in = [] ;
                for l = 1:n(1)
                    vpp = [i- p(l,1),j-p(l,2),k-p(l,3)];
                    vect = vectors_plans(l,:);
                    sign = dot(vpp,vect) ;
                    if sign > 0.0
                        in(end+1) = 1.0;
                    elseif sign == 0.0
                        in(end +1) = 2.0;
                    elseif sign < 0.0
                        in(end+1) = 3.0;
                    end
                end
                if sum(in(:)==1.0) == n(1)
                    mat_fib(i,j,k) = 1.0 ;
                elseif sum(in(:)==2.0) >= 2.0
                    mat_fib(i,j,k) = 2.0 ;
                else
                    mat_fib(i,j,k) = 3.0 ; 
                end
            end
        end
    end
    
    
    

    
    
end 


function ifib = irregular_fibroblast(d)
    irr_fib = zeros(d,d);
    rand_x = randi([d-10,d-1]);

    irr_fib(1.0,rand_x) = 1.0 ;
    x = rand_x ;
    initial_x = rand_x ;
    y = 1.0;

    while x ~= d
        r1 = randi([0,1]) ;
        r2 = randi([0,1]) ;
    
        x = x + r1 ;
        y = y + r2 ;
    
        irr_fib(y,x) = 1.0 ;
    
    end

    rand_y = randi([d-10,d-1]) ;
    irr_fib(y:rand_y,d) = 1.0 ;
    y = rand_y ;

    while y~= d
    
        r1 = randi([0,1]) ;
        r2 = randi([0,1]) ;
    
        x = x - r1 ;
        y = y + r2 ;
    
        irr_fib(y,x) = 1.0 ;
    
    end

    rand_x = randi([2,10]) ;
    irr_fib(d,rand_x:x) = 1.0 ;
    x = rand_x ;

    while x~= 1
    
        r1 = randi([0,1]) ;
        r2 = randi([0,1]) ;
    
        x = x - r1 ;
        y = y - r2 ;
    
        irr_fib(y,x) = 1.0 ;
    end

    rand_y = randi([2,10]);
    irr_fib(rand_y:y,1) = 1.0 ;
    y = rand_y ;

    while y~= 1
    
        r1 = randi([0,1]) ;
        r2 = randi([0,1]) ;
    
        x = x + r1 ;
        y = y - r2 ;
    
        irr_fib(y,x) = 1.0 ;
    end

    final_x = x ;
    irr_fib(1,final_x:initial_x) = 1.0 ;
   
    for i = 1:d
    yy = [];
    for j = 1:d
       if irr_fib(j,i) == 1.0
           yy(end+1) = j;
       end
    end
    
    y1 = min(yy);
    y2 = max(yy);
    
    irr_fib(1:y1,i) = 1.0;
    irr_fib(y2:d,i) = 1.0;
    
    end
    
    ifib = irr_fib;
        
    end

 
