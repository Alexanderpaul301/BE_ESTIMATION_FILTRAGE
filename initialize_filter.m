function [mu, Sigma] = initialize_filter(image_data, carte, f,Sigma_vel,Sigma_biais,mu_vel,mu_biais)
    % Calcul initial de la position et des matrices de covariance
    Xpos = [];
    Ypos = [];
    Zpos = [];
    n = 1;

    for j = 1:size(image_data, 2)
        XA = carte(1, image_data(1, j));
        YA = carte(2, image_data(1, j));
        ZA = carte(3, image_data(1, j));
        UA = image_data(2, j);
        VA = image_data(3, j);

        for i = j+1:size(image_data, 2)
            XB = carte(1, image_data(1, i));
            YB = carte(2, image_data(1, i));
            ZB = carte(3, image_data(1, i));
            UB = image_data(2, i);
            VB = image_data(3, i);

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

    mu_pos = [mean(Xpos); mean(Ypos); mean(Zpos)];
    Sigma_pos = cov([Xpos(:), Ypos(:), Zpos(:)]);
    Sigma = blkdiag(Sigma_pos, Sigma_vel, Sigma_biais);
    mu = [mu_pos; mu_vel; mu_biais];
end