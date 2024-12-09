%% Nettoyage
clear; clc; close all;

%% Chargement des données
load carte.dat % Carte des amers
load mesure_accelero % Mesures des accéléromètres

% Coordonnées des amers
X = carte(1, :);
Y = carte(2, :);
Z = carte(3, :);

%% Paramètres
f = 512; % Distance focale (pixels)
g_moon = 1.622; % Gravité lunaire (m/s²)
dt = 0.01; % Pas de temps (s)
sigma_biais = 0.2; % Écart-type des biais (m/s²)
sigma_vitesse = 2; % Écart-type de la vitesse initiale (m/s)
sigma_acc = sqrt(2e-5); % Écart-type du bruit accéléromètre (m/s²)

% Matrices de covariance initiales
Sigma_vel = eye(3) * sigma_vitesse^2; % Incertitude initiale sur la vitesse
Sigma_biais = eye(3) * sigma_biais^2; % Incertitude sur les biais

% Matrice de bruit de processus Q
Q_pos = zeros(3);
Q_vel = eye(3) * (sigma_acc * dt)^2;
Q_biais = zeros(3);
Q = blkdiag(Q_pos, Q_vel, Q_biais);

% Initialisation des états
mu_pos = [0; 0; 0]; % Position initiale estimée (m)
mu_vel = [100; 0; 0]; % Vitesse initiale estimée (m/s)
mu_biais = [0; 0; 0]; % Biais initial des accéléromètres
mu = [mu_pos; mu_vel; mu_biais];

% Matrice de covariance initiale
Sigma_pos = eye(3) * 100; % Incertitude initiale sur la position (m)
Sigma = blkdiag(Sigma_pos, Sigma_vel, Sigma_biais);

%% Initialisation des trajectoires
num_images = 100; % Nombre total d'images
positions = zeros(num_images + 1, 3); % Trajectoire estimée
positions(1, :) = mu_pos';

%% Boucle principale sur les images
for k = 1:num_images
    % Chargement de l'image
    filename = sprintf('images/image%03d', k);
    if ~isfile(filename)
        error('Fichier %s introuvable.', filename);
    end
    image_data = load(filename);

    % Extraction des données d'amers observés
    amers_obs = image_data(1, :);
    coord_image = image_data(2:3, :);
    coord_3D = [X(amers_obs); Y(amers_obs); Z(amers_obs)];

    % Prédiction des observations
    U_pred = -f * (coord_3D(1, :) - mu(1)) ./ (coord_3D(3, :) - mu(3));
    V_pred = -f * (coord_3D(2, :) - mu(2)) ./ (coord_3D(3, :) - mu(3));
    z_pred = [U_pred; V_pred];
    z_obs = coord_image;

    % Matrice de Jacobienne H
    H = compute_jacobian(mu, coord_3D, f);

    % Kalman Gain
    R = eye(size(H * Sigma * H', 1)) * 3^2; % Bruit de mesure (pixels)
    K = Sigma * H' / (H * Sigma * H' + R);

    % Mise à jour des états
    innovation = z_obs(:) - z_pred(:);
    mu = mu + K * innovation;
    Sigma = (eye(size(Sigma)) - K * H) * Sigma;

    % Sauvegarde de la position
    positions(k + 1, :) = mu(1:3)';

    % Intégration dynamique
    if k < num_images
        for l = 1:100
            a_mes = mesure_accelero(100 * (k - 1) + l, 2:4)';
            a_corr = a_mes - mu(7:9) + [0; 0; -g_moon];
            mu(1:3) = mu(1:3) + mu(4:6) * dt + 0.5 * a_corr * dt^2;
            mu(4:6) = mu(4:6) + a_corr * dt;
            Sigma = Sigma + Q;
        end
    end
end

%% Tracé des résultats
% Trajectoire estimée
figure;
plot3(positions(:, 1), positions(:, 2), positions(:, 3), 'r-', 'LineWidth', 2);
grid on;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Trajectoire estimée');

% Evolution des biais
figure;
time = (0:num_images) * dt;
plot(time, mu(7:end), 'LineWidth', 1.5);
grid on;
xlabel('Temps (s)');
ylabel('Biais (m/s²)');
title('Évolution des biais des accéléromètres');
legend('Biais X', 'Biais Y', 'Biais Z');

