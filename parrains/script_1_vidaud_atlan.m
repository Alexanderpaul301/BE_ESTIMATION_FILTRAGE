clear all
close all
load carte.dat
load mesure_accelero

% numéro d'amer = numéro de colonne
% X = carte(1,:)
% Y = carte(2,:)
% Z = carte(3,:)

%stem3(carte(1,:),carte(2,:),carte(3,:))
%axis([0 6000 -1000 1000 -100 1000]);
%xlabel('X');
%ylabel('Y');
%zlabel('Z');

%figure(1)
% temps = mesure_accelero(:,1)
% aX = mesure_accelero(:,2)
% aY = mesure_accelero(:,3)
% aZ = mesure_accelero(:,4)
%plot(mesure_accelero(:,1),mesure_accelero(:,2),'r',mesure_accelero(:,1),mesure_accelero(:,3),'b',mesure_accelero(:,1),mesure_accelero(:,4),'g')
%xlabel('temps (s)');
%ylabel('acceleration (m/s²)');
%legend('axe X','axe Y','axe Z');

k_plot_image = 50; % Pour visualiser une image donnée

for k = 0:1:100
    [filename,err] = sprintf('image%3.3d',k);
    % image(1,:) = numéro d'amer 
    % U = image(2,:)
    % V = image(3,:)
    image = load(filename);
    if(k == k_plot_image)
        %figure(2)
        %plot(image(2,:),image(3,:),'.');
        %axis([-512.0 512.0 -512.0 512.0]);
        %xlabel('U');
        %ylabel('V');
    end
    if(k == 0)
        %positions
        c=0;
        nbr_amer=180; %size(image)
        nbr_couple=nbr_amer*(nbr_amer-1)/2;
        coor=zeros(3,nbr_couple);
        for i=1:1:nbr_amer
            A=image(1,i);
            XA=carte(1,A);
            YA=carte(2,A);
            ZA=carte(3,A);           
            UA=image(2,i);
            VA=image(3,i);
            for j=i+1:1:nbr_amer
                B=image(1,j);
                XB=carte(1,B);
                YB=carte(2,B);
                ZB=carte(3,B); 
                UB=image(2,j);                
                VB=image(3,j);
                if VA==VB
                    z=(512*(XA-XB)+UA*ZA-UB*ZB)/(UA-UB);
                elseif UA==UB
                    z=(512*(YA-YB)+UA*ZA-UB*ZB)/(VA-VB);
                else
                    z=(512*(YA-YB)+UA*ZA-UB*ZB)/(VA-VB);
                end
            x=(512*XA +UA*(ZA-z))/512;
            y=(512*YA+VA*(ZA-z))/512;
            c=c+1;
            coor(1,c)=x;
            coor(2,c)=y;
            coor(3,c)=z;
            end            
        end 

        M_p=[mean(coor(1,:)) mean(coor(2,:)) mean(coor(3,:))];
        
        cov11=1/nbr_couple*sum((coor(1,:)-M_p(1)).^2);
        cov22=1/nbr_couple*sum((coor(2,:)-M_p(2)).^2);
        cov33=1/nbr_couple*sum((coor(3,:)-M_p(3)).^2);
        cov12=1/nbr_couple*sum((coor(1,:)-M_p(1)).*(coor(2,:)-M_p(2)));
        cov13=1/nbr_couple*sum((coor(1,:)-M_p(1)).*(coor(3,:)-M_p(3)));
        cov23=1/nbr_couple*sum((coor(2,:)-M_p(2)).*(coor(3,:)-M_p(3)));
        
        covariance_p=[cov11 cov12 cov13; cov12 cov22 cov23; cov13 cov23 cov33];
        
        %vitesse
        sigma_v=2;
        M_v=[100 0 0];        
        covariance_v=[sigma_v^2 0 0 ; 0 sigma_v^2 0 ; 0 0 sigma_v^2];
        
        %biais
        sigma_b=0.2;
        M_b=[0 0 0];        
        covariance_b=[sigma_b^2 0 0 ; 0 sigma_b^2 0 ; 0 0 sigma_b^2]; 
  
        % initialisez la moyenne et la covariance de l'état en utilisant
        % des couples d'amers et une idée sur la vitesse initiale et les
        % biais des accéleromètres
    else
        % utilisez l'image et la carte pour recaller la prédiction de l'état
    end
    if(k ~= 100)  
        A=[zeros(3) eye(3) zeros(3) ; zeros(3) zeros(3) -eye(3) ; zeros(3) zeros(3) zeros(3) ];
        B=[zeros(3) ; eye(3) ; zeros(3) ];
        delta=0.01;
        g=[0 0 1.622]';
        M=[M_p  M_v  M_b]';
        U=2e-5;
        U_delta=U/delta;
        cov=[covariance_p zeros(3) zeros(3) ; zeros(3) covariance_v zeros(3) ; zeros(3) zeros(3) covariance_b];
        M_tot=M;
        cov_tot=zeros(9,9,100);
        cov_tot(:,:,1)=cov;
        N=[0];
        for n=1:1:100
            a_mesuree=[mesure_accelero(100*k+n,2) ; mesure_accelero(100*k+n,3) ; mesure_accelero(100*k+n,4) ];
            M=(eye(9)+A.*delta)*M+B*delta*(a_mesuree+g);
            cov=(eye(9)+A.*delta)*cov*(eye(9)+A.*delta)'+U_delta;
            M_tot=[M_tot M];
            cov_tot(:,:,n+1)=cov;
            N=[N n];
        end      
    end    
end

%         figure(3)
%         plot(N,M_tot(1,:),N,M_tot(2,:),N,M_tot(3,:))
%         xlabel('temps (s/100)')
%         ylabel('moyenne position (m)')
%         legend('x','y','z')
        
%         figure(4)
%         plot(N,M_tot(4,:),N,M_tot(5,:),N,M_tot(6,:))
%         xlabel('temps (s/100)')
%         ylabel('moyenne vitesse (m/s)')
%         legend('Vx','Vy','Vz')
      
%         figure(5)
%         plot(N,M_tot(7,:),N,M_tot(8,:),N,M_tot(9,:))
%         xlabel('temps (s/100)')
%         ylabel('moyenne biais (m/s2)')
%         legend('bx','by','bz')
        
%         figure(6)
%         plot(N,squeeze(cov_tot(1,1,:)),N,squeeze(cov_tot(2,2,:)),N,squeeze(cov_tot(3,3,:)))
%         xlabel('temps (s/100)')
%         ylabel('covariance position (m)')
%         legend('x','y','z')
        
%         figure(7)
%         plot(N,squeeze(cov_tot(4,4,:)),N,squeeze(cov_tot(5,5,:)),N,squeeze(cov_tot(6,6,:)))
%         xlabel('temps (s/100)')
%         ylabel('covariance vitesse (m/s)')
%         legend('Vx','Vy','Vz')
   
%         figure(8)
%         plot(N,squeeze(cov_tot(7,7,:)),N,squeeze(cov_tot(8,8,:)),N,squeeze(cov_tot(8,8,:)))
%         xlabel('temps (s/100)')
%         ylabel('covariance biais (m/s2)')
%         legend('bx','by','bz')



% Recalage statique


