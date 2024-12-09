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
plot(mesure_accelero(:, 1), mesure_accelero(:, 2), 'r', ...
    mesure_accelero(:, 1), mesure_accelero(:, 3), 'b', ...
    mesure_accelero(:, 1), mesure_accelero(:, 4), 'g');
xlabel('Temps (s)');
ylabel('Accélération (m/s²)');
legend('Axe X', 'Axe Y', 'Axe Z');
title('Accélérations mesurées');

% --- Initialisation des paramètres ---
%Donnees
f = 512;
g_moon = 1.622;
dt = 1e-2;

% Ecarts types
sigma_biais = 0.2;
sigma_vitesse = 2;
sigma_acc = sqrt(2e-5);

% Matrices de covariance initiales
Sigma_vel = eye(3) * sigma_vitesse ^ 2;
Sigma_biais = eye(3) * sigma_biais ^ 2;

% Position initiale
mu_vel = [100; 0; 0];
mu_biais = [0; 0; 0];

[A, B] = compute_transition_matrix(dt);

% Matrice de bruit
Q_pos = zeros(3);
Q_vel = eye(3) * sigma_acc ^ 2;
Q_biais = zeros(3);
Q = blkdiag(Q_pos, Q_vel, Q_biais);

% --- Boucle principale sur les images ---
num_images = 100; % Nombre total d'images

% Initialisation des variables
positions = zeros(num_images + 1, 3);
positions0 = zeros(100, 3);
U_all = [];
V_all = [];
K_all = [];
U_est = [];
V_est = [];
K_est = [];
biais = zeros(num_images + 1, 3);

k_plot_image = 50;

% Boucle principale
for k = 0:num_images

    % Chargement de l'image
    filename = sprintf('images/image%3.3d', k);

    if isfile(filename)
        image = load(filename);
    else
        error('Le fichier %s est introuvable.', filename);
    end

    % Extraction des numéros des amers observés
    amers_obs = image(1, :);
    coord_image = [image(2, :); image(3, :)]; % Coordonnées (U, V)
    coord_3D = [X(amers_obs); Y(amers_obs); Z(amers_obs)]; % Coordonnées (X, Y, Z)

    if k == k_plot_image
        figure;
        plot(image(2, :), image(3, :), '.');
        axis([-512.0 512.0 -512.0 512.0]);
        xlabel('U');
        ylabel('V');
        title(['Projection de l''image ', num2str(k)]);
    end

    %% Recalage statique pour mise à jour de la position
    if k == 0
        [mu0, Sigma0] = initialize_filter(image, carte, f, Sigma_vel, Sigma_biais, mu_vel, mu_biais);

        % Récupération des points U et V pour cette image
        U_pred = [];
        V_pred = [];
        Sigma = Sigma0;
        mu = mu0;

    else

        % Prédiction
        U_pred = -f * (coord_3D(1, :) - mu(1)) ./ (coord_3D(3, :) - mu(3));
        V_pred = -f * (coord_3D(2, :) - mu(2)) ./ (coord_3D(3, :) - mu(3));
        z_pred = [U_pred; V_pred];
        z_obs = coord_image;

        % Matrice de Jacobienne H
        H = compute_jacobian(mu, coord_3D, f); % Taille m x n

        % Gain de Kalman
        R = eye(size(H, 1));
        K = Sigma * H' / (H * Sigma * H' + R);

        % Matrice d'observation
        S = zeros(size(K, 2), 1);

        for i = 1:size(K, 2) / 2
            S(2 * i - 1) = z_obs(1, i) - z_pred(1, i); % U
            S(2 * i) = z_obs(2, i) - z_pred(2, i); % V
        end

        % Mise à jour de l'état et de la covariance
        mu = mu + K * S;
        Sigma = (eye(size(Sigma)) - K * H) * Sigma;

        % Vérifiez si l'image correspond
        if k == k_plot_image
            % Calcul des projections avant recalage
            position_avant_recalage = mu(1:3)';
            U_pred_avant = -f * (coord_3D(1, :) - position_avant_recalage(1)) ./ (coord_3D(3, :) - position_avant_recalage(3));
            V_pred_avant = -f * (coord_3D(2, :) - position_avant_recalage(2)) ./ (coord_3D(3, :) - position_avant_recalage(3));
            z_pred_avant = [U_pred_avant; V_pred_avant];

            % Recalage
            H = compute_jacobian(mu, coord_3D, f);
            K = Sigma * H' / (H * Sigma * H' + R);
            S = z_obs(:) - z_pred_avant(:);
            mu = mu + K * S;
            Sigma = (eye(size(Sigma)) - K * H) * Sigma;

            % Calcul des projections après recalage
            position_apres_recalage = mu(1:3)';
            U_pred_apres = -f * (coord_3D(1, :) - position_apres_recalage(1)) ./ (coord_3D(3, :) - position_apres_recalage(3));
            V_pred_apres = -f * (coord_3D(2, :) - position_apres_recalage(2)) ./ (coord_3D(3, :) - position_apres_recalage(3));
            z_pred_apres = [U_pred_apres; V_pred_apres];

            % Calcul des erreurs avant et après recalage
            erreur_avant = sqrt(mean(sum((z_obs - z_pred_avant) .^ 2, 1)));
            erreur_apres = sqrt(mean(sum((z_obs - z_pred_apres) .^ 2, 1)));

            % Affichage démo demandée
            disp("Résultats pour l'image de démonstration :");
            disp(['Image : ', num2str(k_plot_image)]);
            disp(['Erreur moyenne avant recalage : ', num2str(erreur_avant, '%.2f'), ' pixels']);
            disp(['Erreur moyenne après recalage : ', num2str(erreur_apres, '%.2f'), ' pixels']);
            disp(['Position avant recalage : [', num2str(position_avant_recalage, '%.2f '), ']']);
            disp(['Position après recalage : [', num2str(position_apres_recalage, '%.2f '), ']']);
            disp('----------------------------------------');
        end

    end

    %% Enregistrement des paramètres

    % Récupération des points U et V
    U_current = image(2, :);
    V_current = image(3, :);

    K_current = k * ones(1, length(U_current));

    U_all = [U_all, U_current];
    V_all = [V_all, V_current];
    K_all = [K_all, K_current];

    K_pred = k * ones(1, length(U_pred));

    U_est = [U_est, U_pred];
    V_est = [V_est, V_pred];
    K_est = [K_est, K_pred];

    positions(k + 1, :) = mu(1:3)';
    biais(k + 1, :) = mu(7:9)';

    %% Intégration dynamique entre les images
    if k ~= num_images

        for l = 0:99
            % Mesure de l'accélération
            a_mes = mesure_accelero(100 * k + l + 1, 2:4)';

            % Calcul de l'accélération corrigée
            e = a_mes - mu(7:9) + [0; 0; -g_moon];

            bruit = [sigma_acc; sigma_acc; sigma_acc];
            e = e + bruit;

            mu = A * mu + B * e;
            Sigma = A * Sigma * A' + Q * dt;

            if k == 0
                positions0(l + 1, :) = mu(1:3)';
            end

        end

    end

end

% Affichage des trajectoires sur 1s
figure;
time1 = linspace(0, 1, size(positions0(:, 2), 1));
subplot(3, 1, 1);
plot(time1, positions0(:, 1), '-r', 'DisplayName', 'Position X');
xlabel('Itération');
ylabel('X (m)');
title('Trajectoire - Position X sur 1s');
grid on;
legend;

subplot(3, 1, 2);
plot(time1, positions0(:, 2), '-g', 'DisplayName', 'Position Y');
xlabel('Itération');
ylabel('Y (m)');
title('Trajectoire - Position Y sur 1s');
grid on;
legend;

subplot(3, 1, 3);
plot(time1, positions0(:, 3), '-b', 'DisplayName', 'Position Z');
xlabel('Itération');
ylabel('Z (m)');
title('Trajectoire - Position Z sur 1s');
grid on;
legend;

% Tracé des positions estimées en 3D à l'instant 1s
figure;
plot3(positions0(:, 1), positions0(:, 2), positions0(:, 3), 'r-o', 'LineWidth', 2);
grid on;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Trajectoire de l''engin sur 1s');

% Affichage des trajectoires
figure;
subplot(3, 1, 1);
plot(1:size(positions(:, 1), 1), positions(:, 1), '-r', 'DisplayName', 'Position X');
xlabel('Itération');
ylabel('X (m)');
title('Trajectoire - Position X');
grid on;
legend;

subplot(3, 1, 2);
plot(1:size(positions(:, 2), 1), positions(:, 2), '-g', 'DisplayName', 'Position Y');
xlabel('Itération');
ylabel('Y (m)');
title('Trajectoire - Position Y');
grid on;
legend;

subplot(3, 1, 3);
plot(1:size(positions(:, 3), 1), positions(:, 3), '-b', 'DisplayName', 'Position Z');
xlabel('Itération');
ylabel('Z (m)');
title('Trajectoire - Position Z');
grid on;
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
time2 = linspace(0, num_images, num_images + 1); % Temps correspondant aux étapes

figure;
plot(time2, biais(:, 1), 'r-', 'LineWidth', 1.5);
hold on;
plot(time2, biais(:, 2), 'g-', 'LineWidth', 1.5);
plot(time2, biais(:, 3), 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Temps (s)');
ylabel('Biais (m/s²)');
title('Évolution des biais des accéléromètres');
legend('Biais X', 'Biais Y', 'Biais Z');

% Tracé 3D
figure;
scatter3(K_all, U_all, V_all, 'filled'); % Points (k, U, V)
hold on
scatter3(K_est, U_est, V_est, 'filled'); % Points (k, U, V)
hold off
legend('observée', 'estimée');
xlabel('Indice de l''image (k)');
ylabel('Coordonnées U');
zlabel('Coordonnées V');
title('Évolution des points (U, V) en fonction des images');
grid on;
