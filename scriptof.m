% Nettoyage
clear all;
close all;
clc;

% Chargement des données
load carte.dat
load mesure_accelero

% Coordonnées des amers dans la carte
X = carte(1, :);
Y = carte(2, :);
Z = carte(3, :); 

% Visualisation de la carte
figure;
stem3(X, Y, Z, 'filled');
axis([0 6000 -1000 1000 -100 1000]);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Carte des amers');

% Affichage des accélérations mesurées
figure;
plot(mesure_accelero(:,1), mesure_accelero(:,2), 'r', ...
     mesure_accelero(:,1), mesure_accelero(:,3), 'b', ...
     mesure_accelero(:,1), mesure_accelero(:,4), 'g');
xlabel('Temps (s)');
ylabel('Accélération (m/s²)');
legend('Axe X', 'Axe Y', 'Axe Z');
title('Accélérations mesurées');

% --- Initialisation des paramètres ---
f = 512; % Distance focale en pixels
g_moon = 1.622; % Gravité lunaire (m/s²)
dt = 1e-2; % Pas de temps (s)
sigma_biais = 0.2; % l'écart-type des biais
sigma_vitesse = 2;
sigma_acc = sqrt(2e-5); % Écart-type du bruit des accéléromètres
[A,B] = compute_transition_matrix(dt); % Transition
% Initialisation de l'état
% mu_pos  Position initiale estimée (X, Y, Z)
mu_vel = [100; 0; 0]; % Vitesse initiale estimée (m/s)
mu_biais = [0; 0; 0]; % Biais initial des accéléromètres
% mu = [mu_pos; mu_vel; mu_biais]; % État initial complet

% Matrices de covariance initiales
% Sigma_pos   incertitude sur la position
Sigma_vel = eye(3) * sigma_vitesse^2;   % Incertitude initiale sur la vitesse
Sigma_biais = eye(3) * sigma_biais^2; % Incertitude sur les biais
%Sigma = blkdiag(Sigma_pos, Sigma_vel, Sigma_biais) Matrice de covariance totale

% Matrice de bruit de processus Q
Q_pos = zeros(3);
Q_vel = eye(3) * (sigma_acc*dt)^2;
Q_biais = zeros(3);
Q = blkdiag(Q_pos, Q_vel, Q_biais);

% --- Boucle principale sur les images ---
num_images = 100; % Nombre total d'images

positions = zeros(num_images + 1, 3); % Tableau pour les positions (X, Y, Z)
U_all = []; % Stockage de tous les points U
V_all = []; % Stockage de tous les points V
K_all = []; % Stockage des indices d'image k
U_est = []; % Stockage de tous les points U
V_est = []; % Stockage de tous les points V
K_est = []; % Stockage des indices d'image k
biais = zeros(num_images + 1, 3); % Tableau pour les biais (b_x, b_y, b_z)

k_plot_image=50;

% Boucle principale sur les images
for k = 0:num_images
    % Chargement de l'image
    filename = sprintf('images/image%3.3d', k);
    if isfile(filename)
        image = load(filename);
    else
        error('Le fichier %s est introuvable.', filename);
    end


    if k == k_plot_image
        figure;
        plot(image(2, :), image(3, :), '.');
        axis([-512.0 512.0 -512.0 512.0]);
        xlabel('U');
        ylabel('V');
        title(['Projection de l''image ', num2str(k)]);
    end


% Recalage statique pour mise à jour de la position
    if k == 0

        Xpos = [];
        Ypos = [];
        Zpos = [];
        n = 1;

        % Calcul initial des positions
        for j = 1:size(image, 2) % Extraction des numéros des amers observés
            XA = carte(1, image(1, j));
            YA = carte(2, image(1, j));
            ZA = carte(3, image(1, j));
            UA = image(2, j);
            VA = image(3, j);

            for i = j:size(image, 2) % Extraction des numéros des amers observés

                if i ~= j              
                    XB = carte(1, image(1, i));
                    YB = carte(2, image(1, i));
                    ZB = carte(3, image(1, i));
                    UB = image(2, i);
                    VB = image(3, i);

                    if (UA ~= UB)
                        ZE = (f * (XA - XB) + UA * ZA - UB * ZB) / (UA - UB);
                    else
                        ZE = (f * (YA - YB) + VA * ZA - VB * ZB) / (VA - VB);
                    end

                    XE = XA - (ZE - ZA) * UA / f;
                    YE = YA - (ZE - ZA) * VA / f;

                    Xpos(n) = XE;
                    Ypos(n) = YE;
                    Zpos(n) = ZE;
                    n = n + 1;
                end

            end

        end
        % Initialisation avec la moyenne des coordonnées des amers visibles
        Xavg = mean(Xpos);
        Yavg = mean(Ypos);
        Zavg = mean(Zpos);
        Sigma_pos = cov([Xpos(:), Ypos(:), Zpos(:)]); % incertitude sur la position
        Sigma = blkdiag(Sigma_pos, Sigma_vel , Sigma_biais); % Matrice de covariance totale


        mu_pos=[Xavg; Yavg; Zavg]; %Position initiale estimée (X, Y, Z)
        mu = [mu_pos; mu_vel; mu_biais]; % Etat initial
        % Récupération des points U et V pour cette image
        U_pred =[]; % Coordonnées U de l'image courante
        V_pred =[]; % Coordonnées V de l'image courante

    else
        % Extraction des numéros des amers observés
        amers_obs = image(1, :); % Numéros des amers visibles
        coord_image = [image(2, :); image(3, :)]; % Coordonnées (U, V)
        coord_3D = [X(amers_obs); Y(amers_obs); Z(amers_obs)]; % Coordonnées réelles des amers
        
        
        %%%%
        U_pred = -f * (coord_3D(1, :) - mu(1)) ./ (coord_3D(3, :) - mu(3));
        V_pred = -f * (coord_3D(2, :) - mu(2)) ./ (coord_3D(3, :) - mu(3));
        z_pred = [U_pred; V_pred];
        z_obs = coord_image(:);

        Zest = A * mu + B * a_real;
        Yest = A * Sigma * A'+ Q*dt;

        % Matrice de Jacobienne H
        H = compute_jacobian(Zest, coord_3D, f); % Taille m x n

        % Kalman gain (avec régularisation pour stabilisation)
        R = eye(size(image(2, :),1)); % Bruit de mesure régularisé
        K = Yest * H' * inv(H * Yest * H' + R); % Gain de Kalman

        % Matrice d'observation S
        S = zeros(size(K,2), 1);
        for i = 1 :size(K,2)/2
            S(2 * i - 1)=image(2, i)-U_pred(i);
            S(2 * i)=image(3, i)-V_pred(i);   
        end

        % Mise à jour de l'état et de la covariance
        mu = Zest + K * S; % Mise à jour de l'estimation a posteriori (y=mu)
        Sigma = (eye(size(Sigma)) - K * H) * Yest; % Mise à jour de la covariance (sigma = P)
       

    end

    % Enregistrement des paramètres

    % Récupération des points U et V pour cette image
    U_current = image(2, :); % Coordonnées U de l'image courante
    V_current = image(3, :); % Coordonnées V de l'image courante
    % Association de l'indice de l'image k à tous les points
    K_current = k * ones(1, length(U_current)); % Même indice k pour tous les points de cette image
    % Concaténation dans les tableaux globaux
    U_all = [U_all, U_current]; % Ajouter les U
    V_all = [V_all, V_current]; % Ajouter les V
    K_all = [K_all, K_current]; % Ajouter les indices d'image
    
    % Association de l'indice de l'image k à tous les points
    K_pred = k * ones(1, length(U_pred)); % Même indice k pour tous les points de cette image
    % Concaténation dans les tableaux globaux
    U_est = [U_est, U_pred]; % Ajouter les U estimé
    V_est = [V_est, V_pred]; % Ajouter les V estimé
    K_est = [K_est, K_pred]; % Ajouter les indices d'image


    positions(k + 1, :) = mu(1:3)'; % Stockage de la position (X, Y, Z)
    biais(k + 1, :) = mu(7:9)'; % Stockage des biais (b_x, b_y, b_z)

    % Intégration dynamique entre les images
    if k ~= num_images
        for l = 0:99
            a_mes = mesure_accelero(100 * k + l + 1, 2:4)'; % Accélérations mesurées
            a_real = a_mes - mu(7:9) + [0; 0; -g_moon]-[sigma_acc; sigma_acc; sigma_acc]; % Accélération réelle (correction des biais)
            mu(1:3) = mu(1:3) + mu(4:6) * dt + 0.5 * a_real * dt^2; % mise a jour de la Position
            Zest = mu(4:6) + a_real * dt; % Vitesse
        end
   end
end

% --- Affichage des trajectoires ---
figure;
subplot(3, 1, 1);
plot(1:size(positions(:, 1), 1), positions(:, 1), '-r', 'DisplayName', 'Position X');
xlabel('Itération');
ylabel('X (m)');
title('Trajectoire - Position X');
legend;

subplot(3, 1, 2);
plot(1:size(positions(:, 2), 1), positions(:, 2), '-g', 'DisplayName', 'Position Y');
xlabel('Itération');
ylabel('Y (m)');
title('Trajectoire - Position Y');
legend;

subplot(3, 1, 3);
plot(1:size(positions(:, 3), 1), positions(:, 3), '-b', 'DisplayName', 'Position Z');
xlabel('Itération');
ylabel('Z (m)');
title('Trajectoire - Position Z');
legend;

% Tracé des positions estimées en 3D
figure;
plot3(positions(:, 1), positions(:, 2), positions(:, 3), 'r-o', 'LineWidth', 2);
grid on;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Trajectoire de l''engin');

% Tracé de l'évolution des biais
time = linspace(0, num_images, num_images + 1); % Temps correspondant aux étapes

figure;
plot(time, biais(:, 1), 'r-', 'LineWidth', 1.5);
hold on;
plot(time, biais(:, 2), 'g-', 'LineWidth', 1.5);
plot(time, biais(:, 3), 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Temps (s)');
ylabel('Biais (m/s²)');
title('Évolution des biais des accéléromètres');
legend('Biais X', 'Biais Y', 'Biais Z');

% --- Tracé 3D ---
figure;
scatter3(K_all, U_all, V_all, 'filled'); % Points (k, U, V)
hold on
scatter3(K_est, U_est, V_est, 'filled'); % Points (k, U, V)
hold off
legend('observée','estimée');
xlabel('Indice de l''image (k)');
ylabel('Coordonnées U');
zlabel('Coordonnées V');
title('Évolution des points (U, V) en fonction des images');
grid on;
