clear all
close all
load carte.dat
load mesure_accelero

% X = carte(1,:), Y = carte(2,:) Z = carte(3,:) numéro d'amer = numéro de
% colonne
stem3(carte(1,:),carte(2,:),carte(3,:))
axis([0 6000 -1000 1000 -100 1000]);
xlabel('X');
ylabel('Y');
zlabel('Z');
%
% figure
temps = mesure_accelero(:,1);
delta=0.01;
aX = mesure_accelero(:,2);
aY = mesure_accelero(:,3);
aZ = mesure_accelero(:,4);

% plot(mesure_accelero(:,1),mesure_accelero(:,2),'r',mesure_accelero(:,1),mesure_accelero(:,3),'b',mesure_accelero(:,1),mesure_accelero(:,4),'g')
% xlabel('temps (s)');
% ylabel('acceleration (m/s²)');
% legend('axe X','axe Y','axe Z');

k_plot_image = 0; % Pour visualiser une image donnée

figure      
hold on;
for k = 0:1:100

    [filename,err] = sprintf('image%3.3d',k);
    % image(1,:) = numéro d'amer U = image(2,:), V = image(3,:)
    image = load(filename);
    if(k == k_plot_image)
%         figure
%         plot(image(2,:),image(3,:),'.');
%         axis([-512.0 512.0 -512.0 512.0]);
%         xlabel('U');
%         ylabel('V');
    end
    n=size(image,2);
    if(k == 0)
        % initialisez la moyenne et la covariance de l'état en utilisant
        % des couples d'amers et une idée sur la vitesse initiale et les
        % biais des accéleromètres
        tot=n*(n-1)/2;
        %moyenne
        mx=0;
        my=0;
        mz=0;
        for i=1:n-1
            ida=image(1,i);
        	ua=image(2,i);
            va=image(3,i);
            xa=carte(1,ida);
            ya=carte(2,ida);
            za=carte(3,ida);
            for j=i+1:n
                idb=image(1,j);
                ub=image(2,j);
                vb=image(3,j);
                xb=carte(1,idb);
                yb=carte(2,idb);
                zb=carte(3,idb);
                
                if (ua~=ub)
                    z=(512*(xa-xb)-ub*zb+za*ua)/(ua-ub);
                else % (va~=vb)
                    z=(512*(ya-yb)-vb*zb+za*va)/(va-vb);
                end
                
                x=xa-ua*(z-za)/512;
                y=ya-va*(z-za)/512;
                mx=mx+x;
                my=my+y;
                mz=mz+z;
            end
        end 
        mx=mx/tot;
        my=my/tot;
        mz=mz/tot;
        m=[mx my mz];
        
        %covariance
        cx=0;
        cy=0;
        cz=0;
        cxy=0;
        cyz=0;
        czx=0;
        so=[];
        for i=1:n-1
            ida=image(1,i);
        	ua=image(2,i);
            va=image(3,i);
            xa=carte(1,ida);
            ya=carte(2,ida);
            za=carte(3,ida);
            for j=i+1:n
                idb=image(1,j);
                ub=image(2,j);
                vb=image(3,j);
                xb=carte(1,idb);
                yb=carte(2,idb);
                zb=carte(3,idb);
                
                if (ua~=ub)
                    z=(512*(xa-xb)-ub*zb+za*ua)/(ua-ub);
                else % (va~=vb)
                    z=(512*(ya-yb)-vb*zb+za*va)/(va-vb);
                end
                
                x=xa-ua*(z-za)/512;
                y=ya-va*(z-za)/512;
                
                cx=cx+(x-mx)^2;
                cy=cy+(y-my)^2;
                cz=cz+(z-mz)^2;
                cxy=cxy+(x-mx)*(y-my);
                cyz=cyz+(z-mz)*(y-my);
                czx=czx+(x-mx)*(z-mz);
            end
            
        end 
        
        %covariance
        cx=cx/tot;
        cy=cy/tot;
        cz=cz/tot;
        cxy=cxy/tot;
        cyz=cyz/tot;
        czx=czx/tot;
        c=[cx cxy czx; cxy cy cyz; czx cyz cz];
        
        %initialisation vitesse
        Vmoy=[100 0 0];
        Cv=[4 0 0;0 4 0;0 0 4];
        
        %initialisation biais
        Bmoy=[0 0 0];
        Cb=[0.04 0 0; 0 0.04 0 ; 0 0 0.04];

        
        X=[mx my mz Vmoy Bmoy]';
        U=[0;0;0;-2e-5;-2e-5;-2e-5;0;0;0]/delta;
        P=[c zeros(3) zeros(3); zeros(3) Cv zeros(3); zeros(3) zeros(3) Cb];
    else %correction
        % utilisez l'image et la carte pour recaller la prédiction de l'état
  
        s=[];
        ms=[];
        C=[];
        sigma=[];
        V=eye(2*n);

        for i=1:n
            Ui=image(2,i);
            Vi=image(3,i);
            Xi=carte(1,image(1,i));
            Yi=carte(2,image(1,i));
            Zi=carte(3,image(1,i));

            s= [s; Ui ; Vi];
            ms=[ms; 512*(Xi)/(X(3)-Zi); 512*(Yi)/(X(3)-Zi)];
            C=[C; -512*Xi/(X(3)-Zi) 0 -512*(Xi-X(1))/((X(3)-Zi)^2) 0 0 0 0 0 0; 0 -512*Yi/(X(3)-Zi) -512*(Yi-X(2))/((X(3)-Zi)^2) 0 0 0 0 0 0];
        end
        sigma=ms+C*X;
        K=P*C'*inv(C*P*C'+V);
        X=X+K*(s-sigma);
        P=(eye(9)-K*C)*P;
    end
    if(k ~= 100)   %prediction
%         figure      
%         hold on;
        for l=0:1:99
            e=mesure_accelero(100*k+l+1,2:4)' - [0; 0; 1.622];
            Ad=delta*[zeros(3) eye(3) zeros(3);zeros(3) zeros(3) -eye(3); zeros(3) zeros(3) zeros(3)]+eye(9);
            Bd=delta*[zeros(3); eye(3) ; zeros(3)];
            
            X=Ad*X+Bd*e;
        
            P=Ad*P*Ad'+U;

            %evolution des moyennes
            plot(temps(100*k+l+1),X(1),'r+',temps(100*k+l+1),X(2),'b+',temps(100*k+l+1),X(3),'g+')

            %evolution de la matrice de covariance position
%             plot(temps(100*k+l+1),P(1,1),'r+',temps(100*k+l+1),P(2,2),'b+',temps(100*k+l+1),P(3,3),'g+')
            
            %evolution de la matrice de covariance vitesse
%             plot(temps(100*k+l+1),P(4,4),'r+',temps(100*k+l+1),P(5,5),'b+',temps(100*k+l+1),P(6,6),'g+')
            
            %evolution de la matrice de covariance biais
%             plot(temps(100*k+l+1),P(7,7),'r+',temps(100*k+l+1),P(8,8),'b+',temps(100*k+l+1),P(9,9),'g+')
            
        end
        
    xlabel("temps en s")
    title("Evolution de la position")
    legend('x','y','z')

%      title("Evolution de la matrice de covariance des positions")
%      legend('covx','covy','covz')

%       title("Evolution de la matrice de covariance de la vitesse")
%       legend('covvx','covvy','covvz')

%     title("Evolution de la matrice de covariance des biais")
%     legend('covbx','covby','covbz')
    end    

end