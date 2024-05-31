module ss_calcs
    using Symbolics, LinearAlgebra, ControlSystems, ControlSystemsBase
    export ss_calcs

    function ctrb_obsv_gains(A, B, C, D, ζ, tₛ, ζₒ, tₛₒ)
        W = ctrb(A, B);
        Wₒ= obsv(A, C);

        s = Symbolics.variable("s");
        Apoly = expand(det(s*I - A));
        as = [Symbolics.coeff(Apoly, s^(size(A)[1]-i)) for i ∈ 1:size(A)[1]];

        Ã = Matrix{Float64}([-transpose(as);
                            1 0 0 0 0;
                            0 1 0 0 0;
                            0 0 1 0 0;
                            0 0 0 1 0]);

        B̃ = [1;0;0;0;0];

        W̃ = ctrb(Ã, B̃);

        T = W̃*inv(W);

        # ζ = 1.5;
        # tₛ = 2;
        ωₙ = log(50)/(tₛ*real(ζ - sqrt(Complex(ζ^2 - 1))));

        poly = expand((s + 3e2*ζ*ωₙ)*(s + 4e2*ζ*ωₙ)*(s + 5e2*ζ*ωₙ)*(s^2 + 2*ζ*ωₙ*s + ωₙ^2));
        # poly = expand((s + 5*ζ*ωₙ)*(s^2 + (10*ζ)*(1.5*ωₙ)*s + (1.5*ωₙ)^2)*(s^2 + 2*ζ*ωₙ*s + ωₙ^2));
        ps = [Symbolics.coeff(poly, s^(size(A)[1]-i)) for i ∈ 1:size(A)[1]];

        ks = ps - as;

        K̃ = transpose(ks);

        K = K̃*T;

        # ref = [r₀; 0;0;0];
        # Kᵣ = [1 1 0 0 0];

        Wₒᵗ = transpose(obsv(A, C));
        Tₒᵀ = W̃*inv(Wₒᵗ);

        # ζₒ = 2;
        # tₛₒ = 1;
        ωₙₒ = log(50)/(tₛₒ*real(ζₒ - sqrt(Complex(ζₒ^2 - 1))));
        # Opoly = expand((s + 200*ζₒ*ωₙₒ)*(s + 300*ζₒ*ωₙₒ)*(s + 400*ζₒ*ωₙₒ)*(s^2 + 2*ζₒ*ωₙₒ*s + ωₙₒ^2));
        Opoly = expand((s + 10*ζₒ*ωₙₒ)*(s^2 + 2*(ζ)*(4*ωₙ)*s + (4*ωₙ)^2)*(s^2 + 2*ζₒ*ωₙₒ*s + ωₙₒ^2));
        Ops = [Symbolics.coeff(Opoly, s^(size(A)[1]-i)) for i ∈ 1:size(A)[1]];
        ls = Ops - as;
        L̃ᵀ = transpose(ls);
        L = transpose(L̃ᵀ*Tₒᵀ);

        return K, L, T
    end
end