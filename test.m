% Nettoyage
clear all;
close all;
clc;

% Chargement des données
load carte.dat
load mesure_accelero

num_images = 100; % Nombre total d'images


% Initialisation des tableaux globaux
U_all = []; % Stockage de tous les points U
V_all = []; % Stockage de tous les points V
K_all = []; % Stockage des indices d'image k

% Boucle sur les images
for k = 0:num_images
    % Chargement de l'image
    filename = sprintf('images/image%3.3d', k);
    if isfile(filename)
        image = load(filename);
    else
        warning('Le fichier %s est introuvable. Ignoré.', filename);
        continue; % Passer à l'image suivante
    end

    % Récupération des points U et V pour cette image
    U_current = image(2, :); % Coordonnées U de l'image courante
    V_current = image(3, :); % Coordonnées V de l'image courante

    % Association de l'indice de l'image k à tous les points
    K_current = k * ones(1, length(U_current)); % Même indice k pour tous les points de cette image

    % Concaténation dans les tableaux globaux
    U_all = [U_all, U_current]; % Ajouter les U
    V_all = [V_all, V_current]; % Ajouter les V
    K_all = [K_all, K_current]; % Ajouter les indices d'image
end

% --- Tracé 3D ---
figure;
scatter3(K_all, U_all, V_all, 'filled'); % Points (k, U, V)
xlabel('Indice de l''image (k)');
ylabel('Coordonnées U');
zlabel('Coordonnées V');
title('Évolution des points (U, V) en fonction des images');
grid on;
