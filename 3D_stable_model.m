 
d = 50 ; % dimensió de la matriu (quadrada)
dx = 0.2 ; % Interespaiat entre elements de matriu (mida alvèol en mm)
dt = 0.0005 ; % Pas de temps (en hores). (0.0005)
iter = 40001 ; % Iteracions de temps que volem calcular
D = 5.0e-2; % Coeficient de difusió (en mm2/h) (article baseline value = 1.0e-3)(1.5e-2)
a = 0.12; % mateixes unitats que la concentració (g/mm3) (0.13)
k = 0.5 ; % S'ha de comprovar amb el valor de l'article (en mm3/h g)?
initial_concentration = 0.7 ; %(g/mm3)
tau = 4.5 ; % Proportionality factor for which we multiply k and D to adjust dyamics (4.5)


% LESIÓ
i_0 = 2 ;
j_0 = 2 ;
k_0 = 25 ;



    
D_les = zeros(d,d,d) ;
D_fib = zeros(d,d,d) ;
c_post = zeros(d,d,d) ;
g = zeros(d,d,d) ;
temps = zeros(iter) ; % Guardem els temps en un vector
dif_lesio_vector = zeros(iter);
fif_fib_vector = zeros(iter);
%g_vector = zeros(iter);
radi = zeros(iter);
% Initial conditions

c = zeros(d,d,d) ; % les concentracions estaran en g/mm^3
c0 = 7.0*ones(d,d,d) ;
eps = 2.40;

% Omplim la matriu c amb les distàncies al centre de la lesió
for k = 1:d
    for i = 1:d
        for j = 1:d
             c(i,j,k) = (((i-i_0)^2.0)+((j-j_0)^2.0) + ((k-k_0)^2.0))^(1.0/2.0) ;           
        end
    end
end



c = 0.5*(1.0-tanh((c-c0)/eps));
c(1,:,:) = 0.0;
c(d,:,:) = 0.0;
c(:,1,:) = 0.0;
c(:,d,:) = 0.0;
c(:,:,1) = 0.0;
c(:,:,d) = 0.0;

sumc0 = sum(c(:));
disp(sumc0);


[x,y,z] = meshgrid([1:1:d]); %falta algo aixi i amb xx yy zz

% figure(1)
% imagesc(c(:,:,10))
% colorbar


%FIBROBLAST
fib = zeros(d,d,d) ;

% condicions de contorn 
fib(1,:,:) = 1.0 ;
fib(d,:,:) = 1.0 ;
fib(:,1,:) = 1.0 ;
fib(:,d,:) = 1.0 ;
fib(:,:,1) = 1.0 ;
fib(:,:,d) = 1.0 ;



% Diffusion coefficient matrix of the lesion
D_les = D;
% Diffusion coefficient matrix of the fibroblast
[gx,gy,gz] = grad_matlab(c);
gradient_c = ((gx.*gx) + (gy.*gy) + (gz.*gz)).^0.5;
D_fib = 300*D.*gradient_c;



% grad_dif_zero = zeros(50,50);
% gradc2d = zeros(50,50);
% gradc2d(:) = D_fib(30,:,:);
% grad_dif_zero(gradc2d > 0.5) = 1.0;

c00 = zeros(d,d); c00(:) = c(:,:,30);
fib00 = zeros(d,d); fib00(:) = fib(:,:,30);
c003d = c;
fib003d = fib;


% Occupation matrix
g = 1.0 -c -fib;



tfib = zeros(1,iter);
tc = zeros(1,iter);

%  v = VideoWriter('Ev_lesio_3d_costat_encapsulacio.mp4');
%  open(v);

% Comencem el bucle per calcular l'evolució temporal
for t=1:iter

    temps(t) = dt*(t-1);

    %Evolució de la lesió
    [c_post,dif_lesio] = lesion_evolution3d(d,D_les,a,k,dt,dx,tau,c,fib,g);
    c = c_post;
    dif_lesio_vector(t) = sum(sum(sum(abs(dif_lesio))));
    radi(t) = radi_lesio(c)/1000.0; %passem el radi a mm


    % Occupation matrix
    g = 1.0 -c -fib ;

    % Diffusion coefficient matrix of the fibroblast
    [gx,gy,gz] = grad_matlab(c);
    gradient_c = ((gx.*gx) + (gy.*gy) + (gz.*gz)).^0.5;
    D_fib = 200.0*D.*gradient_c;
    
    dtf = dt/max(max(max(gradient_c)));

%     disp(sum(c(:)))

    %Evolució del fibroblast
    [fib_post,dif_fib] = diffusion3d(d,D_fib,dtf,dx,fib,g);
    fib = fib_post;
    dif_fib_vector(t) = sum(sum(sum(abs(dif_fib))));


    %Condicions de contorn per al fibroblast
    fib(1,:,:) = 1.0 ;
    fib(d,:,:) = 1.0 ;
    fib(:,1,:) = 1.0 ;
    fib(:,d,:) = 1.0 ;
    fib(:,:,1) = 1.0 ;
    fib(:,:,d) = 1.0 ;

    % Occupation matrix
    g = 1.0 -c -fib;
    if max(max(max(g))) > 1.0
%         disp('aqui!')
    end
    
    if t == 1000
        c25 = zeros(d,d);
        c25(:) = c(:,:,30);
        fib25 = zeros(d,d);
        fib25(:) = fib(:,:,30);
        c253d = c;
        fib253d = fib;
    end

    if t == 10000
        c49 = zeros(d,d);
        c49(:) = c(:,:,30);
        fib49 = zeros(d,d);
        fib49(:) = fib(:,:,30);
        c493d = c;
        fib493d = fib;
    end
    if mod(t,100) == 0.0
%         disp(t);
%             disp('iteracions')
%             disp(t)

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
        c30 = zeros(d,d); c30(:) = c(:,:,10);
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

%         figure(2)
%         imagesc(c30)
%         colorbar
%         drawnow
%         suptitle(num2str(t))
%     
% 
%         figure(2)
%          
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

% figure(2)
% imagesc(c30)
% colorbar
% drawnow
% suptitle(num2str(t))

%  close(v) ;

%  c30 = zeros(d,d); c30(:) = c(:,:,30);
%  figure(2)
%  imagesc(c30)
%  colorbar
% 
figure(3)
subplot(3,2,1)
imagesc(c00)
title('Lesion')
colorbar
axis equal
subplot(3,2,2)
imagesc(fib00)
title('Fibroblasts')
colorbar
axis equal
subplot(3,2,3)
imagesc(c25)
colorbar
axis equal
subplot(3,2,4)
imagesc(fib25)
colorbar
axis equal
subplot(3,2,5)
imagesc(c49)
colorbar
axis equal
subplot(3,2,6)
imagesc(fib49)
colorbar
axis equal


% figure(7)
% fsurf = isosurface(x,y,z,c003d,0.1);
% csurf = isosurface(x,y,z,fib003d,0.1);
% subplot(1,3,1)
% patch(fsurf,'FaceColor','r','FaceAlpha',0.2,'EdgeColor','none')
% hold on
% patch(csurf,'FaceColor','g','FaceAlpha',0.2,'EdgeColor','none')
% hold off
% axis equal
% axis([0,50,0,50,0,50])
% view([45 45])
% fsurf = isosurface(x,y,z,c253d,0.1);
% csurf = isosurface(x,y,z,fib253d,0.1);
% subplot(1,3,2)
% patch(fsurf,'FaceColor','r','FaceAlpha',0.2,'EdgeColor','none')
% hold on
% patch(csurf,'FaceColor','g','FaceAlpha',0.2,'EdgeColor','none')
% hold off
% axis equal
% axis([0,50,0,50,0,50])
% view([45 45])
% fsurf = isosurface(x,y,z,c493d,0.1);
% csurf = isosurface(x,y,z,fib493d,0.1);
% subplot(1,3,3)
% patch(fsurf,'FaceColor','r','FaceAlpha',0.2,'EdgeColor','none')
% hold on
% patch(csurf,'FaceColor','g','FaceAlpha',0.2,'EdgeColor','none')
% hold off
% axis equal
% axis([0,50,0,50,0,50])
% view([45 45])


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

 
