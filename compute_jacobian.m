% --- Fonctions auxiliaires ---
function H = compute_jacobian(mu, coord_3D, f)
    X_e = mu(1);
    Y_e = mu(2);
    Z_e = mu(3);
    n = size(coord_3D, 2); % Nombre d'amers observés
    H = zeros(2 * n, length(mu));
    for i = 1:n
        X_A = coord_3D(1, i);
        Y_A = coord_3D(2, i);
        Z_A = coord_3D(3, i);
        dZ = Z_A - Z_e;
        dX = X_A - X_e;
        dY = Y_A - Y_e;

        % Jacobienne pour U
        H(2 * i - 1, 1) = -f / dZ;
        H(2 * i - 1, 3) = f * dX / dZ^2;

        % Jacobienne pour V
        H(2 * i, 2) = -f / dZ;
        H(2 * i, 3) = f * dY / dZ^2;
    end
end