d = 50 ;
fib2 = get_irregular_fibroblast(d);
disp(sum(fib2(:) == 2.0))
% fib = zeros(d,d,d) ;
% fib2 = zeros(d,d,d);
% 
% vect1 = [1,0,0] ;
% vect2 = [vect1;[-1,0,0]] ;
% vect3 = [vect2;[0,1,0]] ;
% vect4 = [vect3;[0,-1,0]];
% vect5 = [vect4;[0,0,1]] ;
% vect6 = [vect5;[0,0,-1]] ;
% vect7 = [vect6;[-1,-1,-1]];
% vect8 = [vect7;[1,-1,-1]];
% vect9 = [vect8;[-1,1,-1]];
% vect10 = [vect9;[1,1,-1]];
% vect11 = [vect10;[1,1,1]];
% vect12 = [vect11;[1,-1,1]];
% vect13 = [vect12;[-1,1,1]];
% vect14 = [vect13;[-1,-1,1]];
% 
% vects = vect14;

% p1 = [5,25,25] ;
% p2 = [p1;[45,25,25]] ;
% p3 = [p2;[25,5,25]] ;
% p4 = [p3;[25,45,25]] ;
% p5 = [p4;[25,25,5]] ;
% p6 = [p5;[25,25,45]] ;
% p7 = [p6;[40,40,40]] ;
% p8 = [p7;[10,10,40]] ;
% p9 = [p8;[40,10,40]] ;
% p10 = [p9;[10,10,40]] ;
% p11 = [p10;[10,10,5]] ;
% p12 = [p11;[10,40,5]] ;
% p13 = [p12;[40,10,5]] ;
% p14 = [p13;[40,40,5]] ;


% points = p14;
% 
% n = size(vects);

% fib = generar_pla(fib,vect1,p1);

% for i = 1:n(1)
%     fib = generar_pla(fib,vects(i,:),points(i,:));
% end
% 
% fib2 = generar_matriu_fib(fib,vects,points);
% disp(size(fib2))


xx = 1:1:d ;
yy = 1:1:d ;
zz = 1:1:d ;

[x,y,z] = meshgrid(xx,yy,zz) ;

% fsurf = isosurface(x,y,z,fib,0.99);
% figure(1)
% patch(fsurf,'FaceColor','g','FaceAlpha',0.2,'EdgeColor','none')
% axis equal
% axis([0,50,0,50,0,50])
% view([45 45])

f2surf = isosurface(x,y,z,fib2,1.0);
figure(2)
patch(f2surf,'FaceColor','g','FaceAlpha',0.2,'EdgeColor','none')
axis equal
axis([0,50,0,50,0,50])
view([45 45])

% pla2d = zeros(d,d);
% pla2d(:) = fib2(:,25,:);
% 
% figure(2)
% imagesc(pla2d)
% colorbar

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
    
%      p1 = [5,25,25] ;
%     p2 = [p1;[45,25,25]] ;
%     p3 = [p2;[25,5,25]] ;
%     p4 = [p3;[25,45,25]] ;
%     p5 = [p4;[25,25,5]] ;
%     p6 = [p5;[25,25,45]] ;
%     p7 = [p6;[40,40,40]] ;
%     p8 = [p7;[10,10,40]] ;
%     p9 = [p8;[40,10,40]] ;
%     p10 = [p9;[10,10,40]] ;
%     p11 = [p10;[10,10,5]] ;
%     p12 = [p11;[10,40,5]] ;
%     p13 = [p12;[40,10,5]] ;
%     p14 = [p13;[40,40,5]] ;


    points = p14;
    
    n = size(vects);

    for i = 1:n(1)
        fib = generar_pla(fib,vects(i,:),points(i,:));
    end

    fib_cavitat = generar_matriu_fib(fib,vects,points);  
    fibr = fib_cavitat;
    
    dists = distance_planes(14,vects,points);
    disp(dists);
    
    
    
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

function ddd = distance_planes(num,vectors_plans,punts_plans)
    n = size(vectors_plans) ;
    x1 = 25;
    x2 = 25;
    x3 = 25;
    ddd = zeros(num,1);

    for l = 1:num
        vect = vectors_plans(l,:);
        p = punts_plans(l,:);
        a = vect(1);
        b = vect(2);
        c = vect(3);
        D = -a*p(1) -b*p(2) -c*p(3);
        lambda = (a*x1 + b*x2 + c*x3 + D)/((a*a)+(b*b)+(c*c));
        x1f = a*lambda + x1;
        x2f = a*lambda + x2;
        x3f = a*lambda + x3;  
        ddd(l) = (((x1f-x1)^2.0)+((x2f-x2)^2.0)+((x3f-x3)^2.0))^(0.5);
    end 
end
