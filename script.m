clear all
close all

% Chargement des données
load carte.dat
load mesure_accelero

%LÉGENDES
% Pour carte.dat : numéro d'amer = numéro de colonne et X = carte(1,:), Y = carte(2,:) Z = carte(3,:)
% Pour mesure_accelero : temps = mesure_accelero(:,1), aX = mesure_accelero(:,2), aY = mesure_accelero(:,3), aZ = mesure_accelero(:,4)
% Pour image : image(1,:) = numéro d'amer, U = image(2,:), V = image(3,:)

% --- AFFICHAGES INITIAUX ---

% Visualisation des amers
stem3(carte(1, :), carte(2, :), carte(3, :))
axis([0 6000 -1000 1000 -100 1000]);
xlabel('X');
ylabel('Y');
zlabel('Z');

%Visualisation des mesures d'accélération
figure
plot(mesure_accelero(:, 1), mesure_accelero(:, 2), 'r', mesure_accelero(:, 1), mesure_accelero(:, 3), 'b', mesure_accelero(:, 1), mesure_accelero(:, 4), 'g')
xlabel('temps (s)');
ylabel('acceleration (m/s²)');
title('Mesures d''accélération');
legend('axe X', 'axe Y', 'axe Z');

%Choix d'une image spécifique pour le plot
k_plot_image = 50;

for k = 0:5
    % Chargement de l'image
    [filename, err] = sprintf('images/image%3.3d', k);
    image = load(filename);

    % Affichage de l'image
    if k == k_plot_image
        figure;
        plot(image(2, :), image(3, :), '.');
        axis([-512.0 512.0 -512.0 512.0]);
        xlabel('U');
        ylabel('V');
        title(['Projection de l''image ', num2str(k)]);
    end

    % --- INITIALISATION ---

    % initialisez la moyenne et la covariance de l'état en utilisant
    % des couples d'amers et une idée sur la vitesse initiale et les
    % biais des accéléromètres
    if (k == 0)
        Xpos = [];
        Ypos = [];
        Zpos = [];
        df = 512;
        n = 1;

        for j = 1:size(image, 2)

            for i = 1:size(image, 1)

                if i ~= j
                    XA = carte(1, image(1, j));
                    YA = carte(2, image(1, j));
                    ZA = carte(3, image(1, j));
                    UA = image(2, j);
                    VA = image(3, j);

                    XB = carte(1, image(1, i));
                    YB = carte(2, image(1, i));
                    ZB = carte(3, image(1, i));
                    UB = image(2, i);
                    VB = image(3, i);

                    YE = 1 / (UA / VA - UB / VB) * (XB - XA + UA / VA * YA - UB / VB * YB);
                    XE = XA + (YE -YA) * UA / VA;
                    ZE = ZA + df / UA * (XA - XE);

                    Xpos(n) = XE;
                    Ypos(n) = YE;
                    Zpos(n) = ZE;

                    n = n + 1;
                end

            end

        end

        %Calcul des moyennes
        Xavg = mean(Xpos);
        Yavg = mean(Ypos);
        Zavg = mean(Zpos);

        %Calcul des variances
        data = [Xpos(:), Ypos(:), Zpos(:)];

        % Calculer la matrice de covariance
        covMatrix = cov(data);
        covMatrix = diag(diag(covMatrix));

        % Commentaire
        %On a un écart type pour la position en X de 50m environ, la précision
        %du capteur est donc de 50m ce qui est raisonnable étant donné que l'on
        %se situe en moyenne à 1000m des points que l'on observe
        %(mean(Z)=999.98m)

        %Construction de la moyenne et de la covariance pour la vitesse.
        Vavg = [100, 0, 0]; %j'ai modifié la taille de Vavg (Eric)
        covV = diag([2 ^ 2, 2 ^ 2, 2 ^ 2]);

        biaisavg = [0; 0; 0];
        covbiais = diag([0.2 ^ 2, 0.2 ^ 2, 0.2 ^ 2]);
    end

    % --- PARAMÈTRES DU FILTRE DE KALMAN ---

    pas = 10 ^ (-2);
    g = 1.622;
    X = [[Xavg; Yavg; Zavg]; Vavg'; biaisavg]; % État initial
    Pk = blkdiag(covMatrix, covV, covbiais);

    traj = zeros(9, 99); % 9 états sur 99 itérations

    Y = X; % SUR ?

    if (k ~= 100)

        for l = 1:1:99
            %L'accéléromètre va 100 a une fréquence 100 fois plus rapide
            %que le système, ce qui explique pourquoi la double boucle, on
            %doit donc garder la dernière valeur pour mettre à jour l'image.

            % utilisez mesure_accelero(100*k+l,2), mesure_accelero(100*k+l,3) et
            % mesure_accelero(100*k+l,4) pour prédire moyenne et
            % covariance de l'état au temps k+0.01*(l+1)

            % Prédiction
            A = [eye(3), eye(3) * pas, zeros(3);
                 zeros(3), eye(3), -eye(3) * pas;
                 zeros(3), zeros(3), eye(3)];
            B = [zeros(3); eye(3) * pas; zeros(3)];
            COVU = blkdiag(covMatrix, covV, covbiais);

            %Construction de l'estimation à priori
            Z = A * Y + B * [mesure_accelero(100 * k + l, 2);
                             mesure_accelero(100 * k + l, 3);
                             mesure_accelero(100 * k + l, 4) - g];
            COVY = A * Pk * A' + COVU;

            % Recalage
            S = reshape(image(2:end, :), [], 1);
            f = 512; % Distance focale
            n = size(image, 2); % Nombre d'amers

            % Construction de h(z) et de sa jacobienne
            Hz = zeros(2 * n, 1);
            dHz = zeros(9, 2 * n);

            %Calcul de h(z)
            for i = 1:n
                % J'ai corrigé puisque c'était H(x) et non h(z) qu'on avait mis en séance... (Eric)
                Hz(2 * i - 1) = f * (Xpos(i) - Z(1)) / (Z(3) - Zpos(i));
                Hz(2 * i) = f * (Ypos(i) - Z(2)) / (Z(3) - Zpos(i));

                % Dérivées par rapport aux états (position)
                dHz(1, 2 * i - 1) = -f / (Z(3) - Zpos(i));
                dHz(2, 2 * i - 1) = 0;
                dHz(3, 2 * i - 1) = f * (Xpos(i) - Z(1)) / (Z(3) - Zpos(i)) ^ 2;

                dHz(1, 2 * i) = 0;
                dHz(2, 2 * i) = -f / (Z(3) - Zpos(i));
                dHz(3, 2 * i) = f * (Ypos(i) - Z(2)) / (Z(3) - Zpos(i)) ^ 2;

                % Les autres états (vitesses et biais) ne contribuent pas aux projections
                dHz(4:9, 2 * i - 1) = 0;
                dHz(4:9, 2 * i) = 0;
            end

            % On transpose dHz pour avoir la bonne dimension
            dHz = dHz';

            % Correction %
            %Mise à jour de l'estimation à posteriori
            COVV = eye(size(Hz, 1));
            K = COVY * dHz' * (dHz * COVY * dHz' + COVV) ^ (-1);
            Y = Z + K * (S - Hz);
            Pk_1 = (eye(9) - K * dHz) * COVY;

            if mod(l, 20) == 0 % Enregistrer toutes les 50 itérations
                figure(1);
                plot(S, 'o-', 'DisplayName', 'S (Observations)');
                hold on;
                plot(Hz, 'x-', 'DisplayName', 'Hz (Projections estimées)');
                legend;
                xlabel('Index des observations');
                ylabel('Valeur');
                title(['Comparaison entre S et Hz à l''itération ', num2str(l)]);
                hold off;
                saveas(gcf, ['comparaison_S_Hz_iter_' num2str(l) '.png']); % Enregistrer l'image
            end

            % Calcul de l'erreur
            erreur = abs(S - Hz);
            disp(['Erreur max à l''itération ', num2str(l), ': ', num2str(max(erreur))]);

            % Stockage des trajectoires
            traj(:, l) = Y;

        end

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

% --- AFFICHAGE DES BIAIS ESTIMÉS ---

%On remarque que l'on a notamment un problème sur le biais de X, qui ne commence pas à 0.
%Il faudrait donc revoir la manière dont on initialise les biais.

figure;
plot(1:99, traj(4, :), '-r', 'DisplayName', 'Biais X'); % Biais sur X
hold on;
plot(1:99, traj(5, :), '-g', 'DisplayName', 'Biais Y'); % Biais sur Y
plot(1:99, traj(6, :), '-b', 'DisplayName', 'Biais Z'); % Biais sur Z
xlabel('Itérations');
ylabel('Valeur des biais');
title('Évolution des biais estimés');
legend;
