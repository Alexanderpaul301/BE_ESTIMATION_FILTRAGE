clear all;
close all;
clc;

% Chargement des données
load carte.dat
load mesure_accelero

% LÉGENDES
% Pour carte.dat : numéro d'amer = numéro de colonne et X = carte(1,:), Y = carte(2,:) Z = carte(3,:)
% Pour mesure_accelero : temps = mesure_accelero(:,1), aX = mesure_accelero(:,2), aY = mesure_accelero(:,3), aZ = mesure_accelero(:,4)
% Pour image : image(1,:) = numéro d'amer, U = image(2,:), V = image(3,:)

% --- AFFICHAGES INITIAUX ---

% Visualisation des amers
stem3(carte(1, :), carte(2, :), carte(3, :));
axis([0 6000 -1000 1000 -100 1000]);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Amers');

% Visualisation des mesures d'accélération
figure;
plot(mesure_accelero(:, 1), mesure_accelero(:, 2), 'r', mesure_accelero(:, 1), mesure_accelero(:, 3), 'b', mesure_accelero(:, 1), mesure_accelero(:, 4), 'g');
xlabel('temps (s)');
ylabel('acceleration (m/s²)');
title('Mesures d''accélération');
legend('axe X', 'axe Y', 'axe Z');

% --- PARAMÈTRES GÉNÉRAUX ---
pas = 10 ^ (-2);
g = 1.622;
f = 512; % Distance focale
df = f;

% --- INITIALISATION ---
k_plot_image = 50;

for k = 0:99
    % Chargement de l'image
    filename = sprintf('images/image%3.3d', k);
    image = load(filename);

    if k == k_plot_image
        figure;
        plot(image(2, :), image(3, :), '.');
        axis([-512.0 512.0 -512.0 512.0]);
        xlabel('U');
        ylabel('V');
        title(['Projection de l''image ', num2str(k)]);
    end

    if k == 0
        Xpos = [];
        Ypos = [];
        Zpos = [];
        n = 1;

        % Calcul initial des positions
        for j = 1:size(image, 2)
            XA = carte(1, image(1, j));
            YA = carte(2, image(1, j));
            ZA = carte(3, image(1, j));
            UA = image(2, j);
            VA = image(3, j);

            for i = 1:size(image, 2)

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

        % Calcul des moyennes et covariances
        Xavg = mean(Xpos);
        Yavg = mean(Ypos);
        Zavg = mean(Zpos);
        covMatrix = cov([Xpos(:), Ypos(:), Zpos(:)]);
        covMatrix = diag(diag(covMatrix));

        % Initialisation des vitesses et biais
        Vavg = [100; 0; 0];
        covV = diag([2 ^ 2, 2 ^ 2, 2 ^ 2]);
        biaisavg = [0; 0; 0];
        covbiais = diag([0.2 ^ 2, 0.2 ^ 2, 0.2 ^ 2]);
    end

    % Initialisation des matrices du filtre
    X = [[Xavg; Yavg; Zavg]; Vavg; biaisavg];
    Pk = blkdiag(covMatrix, covV, covbiais);

    traj = zeros(9, 99);

    for l = 1:99
        % Matrices A et B
        A = [eye(3), eye(3) * pas, zeros(3); zeros(3), eye(3), -eye(3) * pas; zeros(3), zeros(3), eye(3)];
        B = [zeros(3); eye(3) * pas; zeros(3)];
        COVU = blkdiag(covMatrix, covV, covbiais);

        % Prédiction
        Z = A * X + B * [mesure_accelero(100 * k + l, 2); mesure_accelero(100 * k + l, 3); mesure_accelero(100 * k + l, 4) - g];
        COVY = A * Pk * A' + COVU;

        % Recalage
        S = reshape(image(2:end, :), [], 1);
        n = size(image, 2);
        Hz = zeros(2 * n, 1);
        dHz = zeros(9, 2 * n);

        for i = 1:n
            Hz(2 * i - 1) = f * (Xpos(i) - Z(1)) / (Z(3) - Zpos(i));
            Hz(2 * i) = f * (Ypos(i) - Z(2)) / (Z(3) - Zpos(i));

            dHz(1, 2 * i - 1) = -f / (Z(3) - Zpos(i));
            dHz(2, 2 * i - 1) = 0;
            dHz(3, 2 * i - 1) = f * (Xpos(i) - Z(1)) / (Z(3) - Zpos(i)) ^ 2;

            dHz(1, 2 * i) = 0;
            dHz(2, 2 * i) = -f / (Z(3) - Zpos(i));
            dHz(3, 2 * i) = f * (Ypos(i) - Z(2)) / (Z(3) - Zpos(i)) ^ 2;
        end

        dHz = dHz';

        % Mise à jour
        COVV = eye(size(Hz, 1));
        K = COVY * dHz' * (dHz * COVY * dHz' + COVV) ^ (-1);
        X = Z + K * (S - dHz*Z);
        Pk = (eye(9) - K * dHz) * COVY;

        traj(:, l) = X;

        %if mod(l, 20) == 0
           % figure;
            %plot(S, 'o-', 'DisplayName', 'S (Observations)');
            %hold on;
            %plot(Hz, 'x-', 'DisplayName', 'Hz (Projections estimées)');
            %legend;
            %xlabel('Index des observations');
            %ylabel('Valeur');
            %title(['Comparaison entre S et Hz à l''itération ', num2str(l)]);
            %hold off;
        %end

    end

end

% --- AFFICHAGE DES TRAJECTOIRES ---
figure;
subplot(3, 1, 1);
plot(1:99, traj(1, :), '-r', 'DisplayName', 'Position X');
xlabel('Itération');
ylabel('X (m)');
title('Trajectoire - Position X');
legend;

subplot(3, 1, 2);
plot(1:99, traj(2, :), '-g', 'DisplayName', 'Position Y');
xlabel('Itération');
ylabel('Y (m)');
title('Trajectoire - Position Y');
legend;

subplot(3, 1, 3);
plot(1:99, traj(3, :), '-b', 'DisplayName', 'Position Z');
xlabel('Itération');
ylabel('Z (m)');
title('Trajectoire - Position Z');
legend;

% --- AFFICHAGE DES BIAIS ---
figure;
plot(1:99, traj(7, :), '-r', 'DisplayName', 'Biais X');
hold on;
plot(1:99, traj(8, :), '-g', 'DisplayName', 'Biais Y');
plot(1:99, traj(9, :), '-b', 'DisplayName', 'Biais Z');
xlabel('Itérations');
ylabel('Valeur des biais');
title('Évolution des biais estimés');
legend;
