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
figure
% temps = mesure_accelero(:,1), aX = mesure_accelero(:,2), 
% aY = mesure_accelero(:,3), aZ = mesure_accelero(:,4)
plot(mesure_accelero(:,1),mesure_accelero(:,2),'r',mesure_accelero(:,1),mesure_accelero(:,3),'b',mesure_accelero(:,1),mesure_accelero(:,4),'g')
xlabel('temps (s)');
ylabel('acceleration (m/s²)');
legend('axe X','axe Y','axe Z');

k_plot_image = 50; % Pour visualiser une image donnée

for k = 0:5
    [filename,err] = sprintf('image%3.3d',k);
    % image(1,:) = numéro d'amer U = image(2,:), V = image(3,:)
    image = load(filename);
    if(k == k_plot_image)
        figure
        plot(image(2,:),image(3,:),'.');
        axis([-512.0 512.0 -512.0 512.0]);
        xlabel('U');
        ylabel('V');
    end
    if(k == 0)
        % initialisez la moyenne et la covariance de l'état en utilisant
        % des couples d'amers et une idée sur la vitesse initiale et les
        % biais des accéleromètres
        Xpos=[];
        Ypos=[];
        Zpos=[];
        df=512;
        n=1;
        for j = 1:size(image,2)
            for i = 1:size(image,1)
                if i ~= j
                    XA = carte(1, image(1,j));
                    YA = carte(2, image(1,j));
                    ZA = carte(3, image(1,j));
                    UA= image(2,j);
                    VA=image(3,j);
                    
                    XB = carte(1, image(1,i));
                    YB = carte(2, image(1,i));
                    ZB = carte(3, image(1,i));
                    UB= image(2,i);
                    VB=image(3,i);
                    
                    YE = 1/(UA/VA - UB/VB)*( XB - XA + UA/VA*YA - UB/VB*YB);
                    XE = XA + (YE -YA)* UA/VA;
                    ZE = ZA + df/UA*(XA-XE);

                    
                    Xpos(n)=XE;
                    Ypos(n)=YE;
                    Zpos(n)=ZE;
                end
            end
            n=n+1;
        end
    %Calcul des moyennes
    Xavg=mean(Xpos);
    Yavg=mean(Ypos);
    Zavg=mean(Zpos);

    %Calcul des variances
    % Construire une matrice où chaque colonne est un vecteur
    data = [Xpos(:), Ypos(:), Zpos(:)];
    
    % Calculer la matrice de covariance
    covMatrix = cov(data);
    covMatrix=diag(diag(covMatrix));

    %On a un écart type pour la position en X de 50m environ, la précision
    %du capteur est donc de 50m ce qui est raisonnable étant donné que l'on
    %se situe en moyenne à 1000m des points que l'on observe
    %(mean(Z)=999.98m)
   

    %Construction de la moyenne et de la covariance pour la vitesse.
    Vavg=[100,0,0;
          0,0,0;
          0,0,0];
    covV=[4,0 0;
          0,4,0;
          0,0,4];

    biaisavg=zeros(3,3);
    covbiais=[0.04,0,0;
              0,0.04,0;
              0,0,0.04];
    
    else
        % utilisez l'image et la carte pour recaller la prédiction de l'état
    end
        pas=10^(-2);
        g=1.622;
        X=[Xavg;
            Yavg;
            Zavg;
            100;
            0;
            0;
            0;
            0;
            0];
        
        imagebis=image(2:end,:);
        S=reshape(imagebis, [], 1);

        currentLength = length(S);

        % Desired length
        desiredLength = 360;

        % Check if padding is needed
        if currentLength < desiredLength
            % Calculate the number of zeros to add
            numZerosToAdd = desiredLength - currentLength;
            
            % Append zeros to S
            S = [S; zeros(numZerosToAdd, 1)];
        end

        w=size(S);

        Y=X;

        A=[eye(3,3), eye(3)*pas, zeros(3,3);
            zeros(3,3), eye(3), -eye(3)*pas;
            zeros(3), zeros(3),eye(3)];  % La matrice A est celle que l'on utilise en discret
        
        B=[zeros(3,3);eye(3)*pas;zeros(3,3)]; % La matrice B est celle que l'on utilise en discret
            
        Pk=zeros(9,9);


        COVU=[covMatrix zeros(3) zeros(3);zeros(3) covV zeros(3);zeros(3) zeros(3) covbiais];
    
        traj=zeros(3,99);
    if(k ~= 100)        
        for l=1:1:99
            % utilisez mesure_accelero(100*k+l,2), mesure_accelero(100*k+l,3) et
            % mesure_accelero(100*k+l,4) pour prédire moyenne et
            % covariance de l'état au temps k+0.01*(l+1)

% Prediction %

            %Construction de l'estimation à priori
            Z=A*Y+B*[mesure_accelero(100*k+l,2);mesure_accelero(100*k+l,3);mesure_accelero(100*k+l,4)-g];
            COVY=A*Pk*A'+ COVU;
            f=512;
            n = length(Xpos); % Length of the vectors

            Hz= zeros(2*n,1); 
            dHz= zeros(9,2*n); 
            
            % Calcul de h(z)
            for i = 1:n
                Hz(2*i-1,1) = f*(Xpos(i)-X(1,1))/((Z(3,1))-Zpos(i));
                Hz(2*i,1) = f*(Ypos(i)-X(2,1))/((Z(3,1))-Zpos(i));

                dHz(1,2*i-1) =-f/((Z(3,1))-Zpos(i));
                dHz(2,2*i-1) =0;
                dHz(3,2*i-1) =f*(X(1,1)-Xpos(i))/(((Z(3,1))-Zpos(i))^2);
                dHz(4,2*i-1) =0;
                dHz(5,2*i-1) =0;
                dHz(6,2*i-1) =0;
                dHz(7,2*i-1) =0;
                dHz(8,2*i-1) =0;
                dHz(9,2*i-1) =0;

                dHz(1,2*i) =0;
                dHz(2,2*i) =-f/((Z(3,1))-Zpos(i));
                dHz(3,2*i) =f*(X(2,1)-Ypos(i))/(((Z(3,1))-Zpos(i))^2);
                dHz(4,2*i) =0;
                dHz(5,2*i) =0;
                dHz(6,2*i) =0;
                dHz(7,2*i) =0;
                dHz(8,2*i) =0;
                dHz(9,2*i) =0;
            end 
            dHz=dHz';

% Correction %
            %Mise à jour de l'estimation à posteriori
            COVV=eye(w);
            K= COVY*dHz'*(dHz*COVY*dHz'+ COVV)^(-1);
            Pk_1=(eye(9)-K*dHz)*COVY;
            Y=Z+K*(S-Hz);
            

            %L'accéléromètre va 100 a une fréquence 100 fois plus rapide
            %que le système, ce qui explique pourquoi la double boucle, on
            %doit donc garder la dernière valeur pour mettre à jour l'image.
            traj(1,l)=X(1);
            traj(2,l)=X(2);
            traj(3,l)=X(3);
            figure(1)
            plot(mesure_accelero(100*k+1,1),traj(1,:));
            hold on
            plot(mesure_accelero(100*k+1,1),traj(2,:));
            plot(mesure_accelero(100*k+1,1),traj(3,:));



            figure(1)

            % Subplot for traj(1,:) and X(1,1)
            subplot(3,1,1)
            plot(mesure_accelero(100*k+1,1), traj(1,:), 'DisplayName', 'traj(1,:)');
            hold on
            plot(1:length(X(1,1)), X(1,1), 'DisplayName', 'X(1,1)');
            hold off
            legend show
            title('traj(1,:) and X(1,1)')

            % Subplot for traj(2,:) and X(2,1)
            subplot(3,1,2)
            plot(mesure_accelero(100*k+1,1), traj(2,:), 'DisplayName', 'traj(2,:)');
            hold on
            plot(1:length(X(2,1)), X(2,1), 'DisplayName', 'X(2,1)');
            hold off
            legend show
            title('traj(2,:) and X(2,1)')

            % Subplot for traj(3,:) and X(3,1)
            subplot(3,1,3)
            plot(mesure_accelero(100*k+1,1), traj(3,:), 'DisplayName', 'traj(3,:)');
            hold on
            plot(1:length(X(3,1)), X(3,1), 'DisplayName', 'X(3,1)');
            hold off
            legend show
            title('traj(3,:) and X(3,1)')


        end

    end    
end

