% Nettoyage
clear all;
close all;
clc;

% Chargement des données
load carte.dat; % Carte des amers
load mesure_accelero; % Données d'accélération

% Coordonnées des amers dans la carte
X = carte(1, :);
Y = carte(2, :);
Z = carte(3, :);

% --- Paramètres de l'engin ---
f = 512; % Distance focale (pixels)
g_moon = 1.622; % Gravité lunaire (m/s²)
dt = 0.01; % Pas de temps (s)

% --- Initialisation des incertitudes ---
sigma_pos = 3; % Écart-type de la position (pixels)
sigma_vel = 2; % Écart-type de la vitesse (m/s)
sigma_biais = 0.2; % Écart-type des biais d'accélération (m/s²)
sigma_acc = sqrt(2e-5); % Bruit des accéléromètres (m/s²)

% Matrices de covariance
Q_pos = zeros(3);
Q_vel = eye(3) * (sigma_acc * dt)^2;
Q_biais = zeros(3);
Q = blkdiag(Q_pos, Q_vel, Q_biais); % Bruit du processus

% Initialisation du filtre
num_images = 100; % Nombre d'images à traiter
positions = zeros(num_images + 1, 3); % Tableau pour les positions
biais = zeros(num_images + 1, 3); % Tableau pour les biais
velocities = zeros(num_images + 1, 3); % Tableau pour les vitesses

% --- Initialisation de l'état ---
mu_pos = [1000; 0; 1000]; % Position initiale (moyenne)
Sigma_pos = eye(3) * sigma_pos^2; % Covariance de la position
mu_vel = [100; 0; 0]; % Vitesse initiale
Sigma_vel = eye(3) * sigma_vel^2; % Covariance de la vitesse
mu_biais = [0; 0; 0]; % Biais initiaux
Sigma_biais = eye(3) * sigma_biais^2; % Covariance des biais
mu = [mu_pos; mu_vel; mu_biais]; % État initial
Sigma = blkdiag(Sigma_pos, Sigma_vel, Sigma_biais); % Covariance totale

% --- Boucle principale ---
for k = 0:num_images
    % Chargement de l'image
    filename = sprintf('images/image%3.3d', k);
    if isfile(filename)
        image = load(filename);
    else
        warning('Image %s introuvable, passage à l\image suivante.', filename);
        continue;
    end
    
    % --- Mesures de l'image ---
    amers_obs = image(1, :); % Numéros des amers visibles
    coord_image = [image(2, :); image(3, :)]; % Coordonnées (U, V)
    coord_3D = [X(amers_obs); Y(amers_obs); Z(amers_obs)]; % Coordonnées 3D des amers observés
    
    % Prédiction des coordonnées dans le plan image
    U_pred = -f * (coord_3D(1, :) - mu(1)) ./ (coord_3D(3, :) - mu(3));
    V_pred = -f * (coord_3D(2, :) - mu(2)) ./ (coord_3D(3, :) - mu(3));
    z_pred = [U_pred; V_pred]; % Prédiction des mesures

    % Observation réelle
    z_obs = coord_image; 

    % Calcul de la Jacobienne H
    H = compute_jacobian(mu, coord_3D, f);

    % Gain de Kalman
    R = eye(size(K)) * sigma_pos^2; % Bruit de mesure
    K = Sigma * H' / (H * Sigma * H' + R);

    % Mise à jour de l'état
    innovation = z_obs(:) - z_pred(:); % Différence entre mesures réelles et prédites
    mu = mu + K * innovation; % Mise à jour de l'état
    Sigma = (eye(size(Sigma)) - K * H) * Sigma; % Mise à jour de la covariance

    % --- Enregistrement ---
    positions(k + 1, :) = mu(1:3)';
    velocities(k + 1, :) = mu(4:6)';
    biais(k + 1, :) = mu(7:9)';

    % --- Intégration dynamique ---
    if k < num_images
        for step = 1:100
            a_mes = mesure_accelero(100 * k + step, 2:4)';
            a_real = a_mes - mu(7:9) + [0; 0; -g_moon]; % Correction des biais
            mu(1:3) = mu(1:3) + mu(4:6) * dt + 0.5 * a_real * dt^2;
            mu(4:6) = mu(4:6) + a_real * dt;
        end
    end
end

% --- Affichage des résultats ---
% Trajectoire
figure;
plot3(positions(:, 1), positions(:, 2), positions(:, 3), 'r', 'LineWidth', 1.5);
grid on;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Trajectoire estimée');

% Évolution des biais
figure;
plot(0:num_images, biais);
xlabel('Temps (s)');
ylabel('Biais (m/s²)');
legend('Biais X', 'Biais Y', 'Biais Z');
title('Évolution des biais des accéléromètres');
grid on;

% --- Fonction pour calculer la Jacobienne ---
function H = compute_jacobian(mu, coord_3D, f)
    n = size(coord_3D, 2);
    H = zeros(2 * n, length(mu));
    for i = 1:n
        dX = coord_3D(1, i) - mu(1);
        dY = coord_3D(2, i) - mu(2);
        dZ = coord_3D(3, i) - mu(3);
        H(2 * i - 1:2 * i, :) = [
            -f / dZ, 0, f * dX / dZ^2, 0, 0, 0, 0, 0, 0;
            0, -f / dZ, f * dY / dZ^2, 0, 0, 0, 0, 0, 0
        ];
    end
end
