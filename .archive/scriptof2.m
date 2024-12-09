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
f = 512; % Distance focale en pixels
g_moon = 1.622; % Gravité lunaire (m/s²)
dt = 1e-2; % Pas de temps (s)
sigma_biais = 0.2; % l'écart-type des biais
sigma_vitesse = 2;
sigma_acc = sqrt(2e-5); % Écart-type du bruit des accéléromètres
[A, B] = compute_transition_matrix(dt); % Transition

% Initialisation de l'état
mu_vel = [100; 0; 0]; % Vitesse initiale estimée (m/s)
mu_biais = [0; 0; 0]; % Biais initial des accéléromètres

% Matrices de covariance initiales
Sigma_vel = eye(3) * sigma_vitesse ^ 2; % Incertitude initiale sur la vitesse
Sigma_biais = eye(3) * sigma_biais ^ 2; % Incertitude sur les biais
Q_pos = zeros(3);
Q_vel = eye(3) * (sigma_acc) ^ 2;
Q_biais = zeros(3);
Q = blkdiag(Q_pos, Q_vel, Q_biais);

% --- Boucle principale sur les images ---
num_images = 100; % Nombre total d'images
positions = zeros(num_images + 1, 3); % Tableau pour les positions (X, Y, Z)
U_all = [];
V_all = [];
K_all = [];
U_est = [];
V_est = [];
K_est = [];
biais = zeros(num_images + 1, 3);
k_plot_image = 50;

for k = 0:num_images
    filename = sprintf('images/image%3.3d', k);

    if isfile(filename)
        image = load(filename);
    else
        error('Le fichier %s est introuvable.', filename);
    end

    amers_obs = image(1, :);
    coord_image = [image(2, :); image(3, :)];
    coord_3D = [X(amers_obs); Y(amers_obs); Z(amers_obs)];

    if k == k_plot_image
        figure;
        plot(image(2, :), image(3, :), '.');
        axis([-512.0 512.0 -512.0 512.0]);
        xlabel('U');
        ylabel('V');
        title(['Projection de l''image ', num2str(k)]);
    end

    if k == 0
        [mu, Sigma] = initialize_filter(image, carte, f, Sigma_vel, Sigma_biais, mu_vel, mu_biais);
        U_pred = [];
        V_pred = [];
    else
        U_pred = -f * (coord_3D(1, :) - mu(1)) ./ (coord_3D(3, :) - mu(3));
        V_pred = -f * (coord_3D(2, :) - mu(2)) ./ (coord_3D(3, :) - mu(3));
        z_pred = [U_pred; V_pred];
        z_obs = coord_image;

        H = compute_jacobian(mu, coord_3D, f);
        R = eye(size(H, 1));
        K = Sigma * H' / (H * Sigma * H' + R);
        %S = z_obs(:) - z_pred(:);
        S = zeros(size(K,2), 1);
        for i = 1 :size(K,2)/2
            S(2 * i - 1) = z_obs(1, i) - z_pred(1, i); % U
            S(2 * i) = z_obs(2, i) - z_pred(2, i); % V 
        end
        mu = mu + K * S;
        Sigma = (eye(size(Sigma)) - K * H) * Sigma;
    end

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

    if k ~= num_images

        for l = 0:99
            a_mes = mesure_accelero(100 * k + l + 1, 2:4)';
            e = a_mes - mu(7:9) + [0; 0; -g_moon];
            bruit = [sigma_acc;sigma_acc;sigma_acc];  % Bruit gaussien ajouté à l'estimation
            e = e + bruit; % Accélération corrigée avec bruit
            mu = A * mu + B * e;
            Sigma = A * Sigma * A' + Q*dt;
        end

    end

end

figure;
subplot(3, 1, 1);
plot(1:size(positions(:, 1), 1), positions(:, 1), '-r');
xlabel('Itération');
ylabel('X (m)');
title('Trajectoire - Position X');

subplot(3, 1, 2);
plot(1:size(positions(:, 2), 1), positions(:, 2), '-g');
xlabel('Itération');
ylabel('Y (m)');
title('Trajectoire - Position Y');

subplot(3, 1, 3);
plot(1:size(positions(:, 3), 1), positions(:, 3), '-b');
xlabel('Itération');
ylabel('Z (m)');
title('Trajectoire - Position Z');

figure;
plot3(positions(:, 1), positions(:, 2), positions(:, 3), 'r-o', 'LineWidth', 2);
grid on;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Trajectoire de l''engin');

figure;
plot(linspace(0, num_images, num_images + 1), biais);
grid on;
xlabel('Temps (s)');
ylabel('Biais (m/s²)');
title('Évolution des biais');
legend('Biais X', 'Biais Y', 'Biais Z');

% --- Tracé 3D ---
figure;
scatter3(K_all, U_all, V_all,'filled'); % Points (k, U, V)
hold on
scatter3(K_est, U_est, V_est,'filled'); % Points (k, U, V)
hold off
legend('observée','estimée');
xlabel('Indice de l''image (k)');
ylabel('Coordonnées U');
zlabel('Coordonnées V');
title('Évolution des points (U, V) en fonction des images');
grid on;