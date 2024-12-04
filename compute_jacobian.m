function H = compute_jacobian(mu, coord_3D, f)
    X_e = mu(1); % Position estimée X de l'engin
    Y_e = mu(2); % Position estimée Y de l'engin
    Z_e = mu(3); % Position estimée Z de l'engin
    n = size(coord_3D, 2); % Nombre d'amers observés

    % Initialisation de la Jacobienne
    H = zeros(2 * n, length(mu)); % Dimensions : 2n x 9 (état complet)

    % Calcul pour chaque amer observé
    for i = 1:n
        % Coordonnées réelles de l'amer
        X_A = coord_3D(1, i);
        Y_A = coord_3D(2, i);
        Z_A = coord_3D(3, i);

        % Différences
        dZ = Z_A - Z_e;
        dX = X_A - X_e;
        dY = Y_A - Y_e;

        % Dérivées pour U
        H(2 * i - 1, 1) = -f / dZ;            % dU/dx
        H(2 * i - 1, 3) = f * dX / dZ^2;      % dU/dz

        % Dérivées pour V
        H(2 * i, 2) = -f / dZ;                % dV/dy
        H(2 * i, 3) = f * dY / dZ^2;          % dV/dz
    end
end