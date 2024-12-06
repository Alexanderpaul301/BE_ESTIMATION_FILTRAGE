function [A,B]= compute_transition_matrix(dt)
    A = [eye(3), eye(3) * dt, zeros(3); zeros(3), eye(3), -eye(3) * dt; zeros(3), zeros(3), eye(3)];
    B = [zeros(3), zeros(3), zeros(3); zeros(3),  eye(3) * dt, zeros(3);zeros(3),zeros(3), zeros(3)];
end