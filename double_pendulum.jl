using Symbolics, LinearAlgebra, Latexify, ControlSystems, ControlSystemsBase, Plots, MAT, LaTeXStrings, .Simulate, .ss_calcs
pgfplotsx()

symbolic = false;
font_size = 18

@variables t;
@variables x(t);
@variables θ(t) ϕ(t) θ̇(t) ϕ̇(t) θ̈(t) ϕ̈(t) θ³(t) ϕ³(t);
# θ̇ = Symbolics.Differential(t)(θ);
# ϕ̇ = Symbolics.Differential(t)(ϕ);
# θ̈ = (Symbolics.Differential(t)^2)(θ);
# ϕ̈ = (Symbolics.Differential(t)^2)(ϕ);
Dt = Symbolics.Differential(t);
Dt2 = Dt^2;
Dt3 = Dt^3;
ψ = θ + ϕ;

v_x = Dt(x);
a_x = Dt2(x);
# ψ̇ = θ̇ + ϕ̇ ;
# ψ̈ = θ̈ + ϕ̈ ;


if symbolic
    @variables g m_cart;
    Ls = Symbolics.variables(:L, 1:2);
    ms = Symbolics.variables(:m, 1:2);
else
    g = 9.8;
    m_cart = 0.5;
    Ls = [0.245; 0.121];
    ms = [0.18; 0.105];
end

Is = (1//12).*ms.*Ls.^2;

D = Symbolics.derivative;

r = [x;0];
v = [v_x; 0];
a = [a_x; 0]

r_O1 = Ls[1]*[cos(θ);sin(θ)];
r_12 = Ls[2]*[cos(ψ);sin(ψ)];

v_c1 = v + (1//2)*D(r_O1, θ)*Dt(θ);
v_c2 = v + D(r_O1, θ)*Dt(θ) + (1//2)*D(r_12, ψ)*Dt(ψ);

a_c1 = a -(1//2)*r_O1*Dt(θ)^2 + (1//2)*D(r_O1, θ)*Dt2(θ)
a_c2 = a -r_O1*(Dt(θ))^2 - (1//2)*r_12*(Dt(ψ))^2 + D(r_O1, θ)*Dt2(θ) + D(r_12, ψ)*Dt2(ψ);

a_c1 = Symbolics.expand(a_c1);
a_c2 = Symbolics.expand(a_c2);

K = (1//2)*ms[1]*dot(v_c1, v_c1) + (1//2)*ms[2]*dot(v_c2, v_c2) + (1//2)*Is[1]*Dt(θ)^2 + (1//2)*Is[2]*Dt(ψ)^2 + (1//2)*m_cart*dot(v, v);
U = ms[1]*g*r_O1[2] + ms[2]*g*r_12[2];

L = K - U;

q = [θ; ϕ];
q̇ = [Dt(θ); Dt(ϕ)];
q̈ = [Dt2(θ); Dt2(ϕ)];

dLdq = Symbolics.gradient(L, q);
dLdq̇ = Symbolics.gradient(L, q̇);
dLdq̇dt = D(dLdq̇, t);

subs = Dict(Dt(θ) => θ̇ , Dt(ϕ) => ϕ̇ , Dt(ψ) => (θ̇ + ϕ̇ ), Dt2(θ) => θ̈ , Dt2(ϕ) => ϕ̈ , Dt2(ψ) => (θ̈ + ϕ̈ ));


dLdq = Vector{Num}(substitute(dLdq, subs));
dLdq̇ = Vector{Num}(substitute(dLdq̇, subs));
dLdq̇dt = Vector{Num}(substitute(dLdq̇dt, subs));

q̇ = Vector{Num}(substitute(q̇, subs));
q̈ = Vector{Num}(substitute(q̈, subs));

# Sub into velocities
v_c1 = Vector{Num}(substitute(v_c1, subs));
v_c2 = Vector{Num}(substitute(v_c2, subs));

# Sub into accelerations
a_c1 = Vector{Num}(substitute(a_c1, subs));
a_c2 = Vector{Num}(substitute(a_c2, subs));


@variables vₓ aₓ;
vel_accel_sub = Dict(Dt(x) => vₓ, Dt2(x) => aₓ);

dLdq = Vector{Num}(substitute(dLdq, vel_accel_sub));
dLdq̇ = Vector{Num}(substitute(dLdq̇, vel_accel_sub));
dLdq̇dt = Vector{Num}(substitute(dLdq̇dt, vel_accel_sub));

# Sub into velocities
v_c1 = Vector{Num}(substitute(v_c1, vel_accel_sub));
v_c2 = Vector{Num}(substitute(v_c2, vel_accel_sub));

# Sub into accelerations
a_c1 = Vector{Num}(substitute(a_c1, vel_accel_sub));
a_c2 = Vector{Num}(substitute(a_c2, vel_accel_sub));

v = Vector{Num}(substitute(v, vel_accel_sub));
a = Vector{Num}(substitute(a, vel_accel_sub));


@variables cₜ cₚ sₜ sₚ ċₜ ċₚ ṡₜ ṡₚ c̈ₜ c̈ₚ s̈ₜ s̈ₚ;
trig_funcs = [cos(θ); cos(θ + ϕ); sin(θ); sin(θ + ϕ); cos(θ̇ ); cos(θ̇ + ϕ̇ ); sin(θ̇ ); sin(θ̇ + ϕ̇ ); cos(θ̈ ); cos(θ̈ + ϕ̈ ); sin(θ̈ ); sin(θ̈ + ϕ̈ )];

trig_stubs = [cₜ; cₚ; sₜ; sₚ; ċₜ; ċₚ; ṡₜ; ṡₚ; c̈ₜ; c̈ₚ; s̈ₜ; s̈ₚ];

trig_sub = Dict([trig_funcs[i] => trig_stubs[i] for i ∈ 1:size(trig_funcs)[1]]);

dLdq = Vector{Num}(substitute(dLdq, trig_sub));
dLdq̇ = Vector{Num}(substitute(dLdq̇, trig_sub));
dLdq̇dt = Vector{Num}(substitute(dLdq̇dt, trig_sub));

# Solve for aₓ in terms of the input force F
@variables F N;
M = m_cart + sum(ms);
F_net = [F; N - M*g];
F_com = m_cart*a + ms[1]*a_c1 + ms[2]*a_c2;

a_cart = Num(Symbolics.solve_for(F_net[1] ~ F_com[1], aₓ));
a_cart = Num(substitute(a_cart, trig_sub));

Asub = Dict(aₓ => a_cart);

dLdq = Vector{Num}(substitute(dLdq, Asub));
dLdq̇ = Vector{Num}(substitute(dLdq̇, Asub));
dLdq̇dt = Vector{Num}(substitute(dLdq̇dt, Asub));


dLdq = expand(dLdq);
dLdq̇ = expand(dLdq̇);
dLdq̇dt = expand(dLdq̇dt);

sol = Vector{Num}(Symbolics.solve_for([dLdq[i] ~ dLdq̇dt[i] for i ∈ 1:size(q)[1]], q̈))

trig_inv_sub = Dict([trig_stubs[i] => trig_funcs[i] for i ∈ 1:size(trig_funcs)[1]]);

q̈ = Vector{Num}(Symbolics.substitute(sol, trig_inv_sub))

dLdq = substitute(dLdq, trig_inv_sub);
dLdq̇ = substitute(dLdq̇, trig_inv_sub);
dLdq̇dt = substitute(dLdq̇dt, trig_inv_sub);
a_cart = expand(Num(substitute(a_cart, trig_inv_sub)));


open("out.txt", "w") do io
    redirect_stdout(io) do 
        println(latexify(dLdq))
        println(latexify(dLdq̇dt))
        print("q̈ = ")
        println(latexify(q̈))
    end
end

a_cart = expand(Num(substitute(a_cart, Dict(θ̈ => q̈[1], ϕ̈ => q̈[2]))));


@variables θᵢ ϕᵢ θᵣ ϕᵣ;
r = [θᵣ; ϕᵣ];
qᵢ = [θᵢ; ϕᵢ];

X = [q; q̇; vₓ];
Ẋ = [q̇; q̈; a_cart];
u = [F];

Aj = Symbolics.jacobian(Ẋ, X);
Bj = Symbolics.jacobian(Ẋ, u);

open("out.txt", "w") do io
    redirect_stdout(io) do 
        println("a_cart")
        println(latexify(a_cart))
        println("dL/dq")
        println(latexify(dLdq))
        println("d/dt(dL/dq̇)")
        println(latexify(dLdq̇dt))
        print("q̈")
        println(latexify(q̈))
        println("Aⱼ")
        println(latexify(Aj))
        println("Bⱼ")
        println(latexify(Bj))
    end
end

if !symbolic
    rs = [[-pi/2; 0], [pi/2, -pi], [pi/2; 0], [-pi/2, -pi], [pi/2; 0], [-pi/2; 0]]
    t_arr = [0, 5, 15, 25, 40, 55, 70]*3;

    # rs = [[-pi/2; 0], [pi/2; -pi]]
    # t_arr = [0, 10, 60];

    X₀ₛ = [[rs[i]; 0;0;0] for i ∈ 1:size(rs)[1]];
    u₀ₛ = [[0] for i ∈ 1:size(rs)[1]];
    
    
    Δt = 1e-3;
    Xₜ₀ = [-pi/2; 0; 0; 0; 0];
    # X̂ₜ₀ = [-pi/2; 0; 0; 0; 0];
    X̂ₜ₀ = [0; 0; 0; 0; 0]
    C = [1 1 0 0 1];
    D = [0];
    ζ = 2;
    tₛ = 30;
    ζₒ = 0.8;
    tₛₒ = 0.5;

    Cₙ = [deg2rad(0.009) deg2rad(0.009) 0 0 1e-9]*10;
    d = [0; 0; 0; 0; 5e-2];

    step = 100;

    plot_obsv = true
    plot_noise = true
    plot_disturbance = false
    plot_speeds = true
    save_plots = true

    if !plot_noise
        Cₙ = Cₙ*0;
    end
    if !plot_disturbance
        d = d*0;
    end

    ts, Xs, X̂s, systems, CLsystems = Simulate.simulate_multiple_references(t_arr, Δt, X₀ₛ, u₀ₛ, Xₜ₀, X̂ₜ₀, X, u, Aj, Bj, C, D, ζ, tₛ, ζₒ=ζₒ, tₛₒ=tₛₒ, Cₙ=Cₙ, d=d)

    
    θᵣₛ = Vector{Float64}()
    ϕᵣₛ = Vector{Float64}()
    tᵣₛ = Vector{Float64}()
    for i ∈ 1:size(X₀ₛ)[1]
        append!(θᵣₛ, X₀ₛ[i][1]);
        append!(ϕᵣₛ, X₀ₛ[i][2]);
        append!(tᵣₛ, t_arr[i]);

        append!(θᵣₛ, X₀ₛ[i][1]);
        append!(ϕᵣₛ, X₀ₛ[i][2]);
        append!(tᵣₛ, t_arr[i+1]);
    end

    if plot_obsv

        p = plot(ts[1:step:end], [X̂s[2,1:step:end] Xs[2,1:step:end]]*180/pi, label=[L"\hat{\phi}(t)" L"\phi(t)"], color=[:peru :dodgerblue3], legend=:topright, tickfontsize=font_size, legendfontsize=font_size, guidefontsize=font_size)
        xlabel!(L"t\,({\rm s})")
        ylabel!(L"\phi\, \left(\circ\right)")
        plot!(tᵣₛ, ϕᵣₛ*180/pi, line=:dash, color="black", label=L"\phi_r")
        if save_plots
            savefig(p, "plots/phi_observer.pdf")
        else
            display(p)
        end

        p = plot(ts[1:step:end], [X̂s[1,1:step:end] Xs[1,1:step:end]]*180/pi, label=[L"\hat{\theta}(t)" L"\theta(t)"], color=[:peru :dodgerblue3], legend=:bottomright, tickfontsize=font_size, legendfontsize=font_size, guidefontsize=font_size)
        xlabel!(L"t\,({\rm s})")
        ylabel!(L"\theta\, (\circ)")
        plot!(tᵣₛ, θᵣₛ*180/pi, line=:dash, color="black", label=L"\theta_r")
        if save_plots
            savefig(p, "plots/theta_observer.pdf")
        else
            display(p)
        end

        if plot_speeds
            p = plot(ts[1:step:end], [X̂s[4,1:step:end] Xs[4,1:step:end]]*180/pi, label=[L"\dot{\hat{\phi}}(t)" L"\dot{\phi}(t)"], color=[:peru :dodgerblue3], legend=:topright, tickfontsize=font_size, legendfontsize=font_size, guidefontsize=font_size)
            xlabel!(L"t\,({\rm s})")
            ylabel!(L"\dot{\phi} (\circ /{\rm s})")
            plot!(tᵣₛ, 0*ϕᵣₛ, line=:dash, color="black", label=L"\dot{\phi}_r")
            if save_plots
                savefig(p, "plots/phi_dot_observer.pdf")
            else
                display(p)
            end

            p = plot(ts[1:step:end], [X̂s[3,1:step:end] Xs[3,1:step:end]]*180/pi, label=[L"\dot{\hat{\theta}}(t)" L"\dot{\theta}(t)"], color=[:peru :dodgerblue3], legend=:bottomright, tickfontsize=font_size, legendfontsize=font_size, guidefontsize=font_size)
            xlabel!(L"t\,({\rm s})")
            ylabel!(L"\dot{\theta} \left(\circ /{\rm s}\right)")
            plot!(tᵣₛ, 0*θᵣₛ, line=:dash, color="black", label=L"\dot{\theta}_r")
            if save_plots
                savefig(p, "plots/theta_dot_observer.pdf")
            else
                display(p)
            end

            p = plot(ts[1:step:end], [X̂s[5,1:step:end] Xs[5,1:step:end]], label=[L"\hat{v}(t)" L"v(t)"], color=[:peru :dodgerblue3], legend=:topright, tickfontsize=font_size, legendfontsize=font_size, guidefontsize=font_size)
            xlabel!(L"t\,({\rm s})")
            ylabel!(L"v ({\rm m/s})")
            plot!(tᵣₛ, 0*θᵣₛ*180/pi, line=:dash, color="black", label=L"v_r")
            if save_plots
                savefig(p, "plots/v_observer.pdf")
            else
                display(p)
            end
        end
    else
        p = plot(ts[1:step:end], Xs[2,1:step:end]*180/pi, label=L"\phi(t)", color=:dodgerblue3, legend=:bottomright, tickfontsize=font_size, legendfontsize=font_size, guidefontsize=font_size)
        xlabel!(L"t\,({\rm s})")
        ylabel!(L"\phi\, \left(\circ\right)")
        plot!(tᵣₛ, ϕᵣₛ*180/pi, line=:dash, color="black", label=L"\phi_r")
        if save_plots
            savefig(p, "plots/phi_no_observer.pdf")
        else
            display(p)
        end

        p = plot(ts[1:step:end], Xs[1,1:step:end]*180/pi, label=L"\theta(t)", color=:dodgerblue3, legend=:topright, tickfontsize=font_size, legendfontsize=font_size, guidefontsize=font_size)
        xlabel!(L"t\,({\rm s})")
        ylabel!(L"\theta\, (\circ)")
        plot!(tᵣₛ, θᵣₛ*180/pi, line=:dash, color="black", label=L"\theta_r")
        if save_plots
            savefig(p, "plots/theta_no_observer.pdf")
        else
            display(p)
        end

        if plot_speeds
            p = plot(ts[1:step:end], Xs[4,1:step:end]*180/pi, label=L"\dot{\phi}(t)", color=:dodgerblue3, legend=:topright, tickfontsize=font_size, legendfontsize=font_size, guidefontsize=font_size)
            xlabel!(L"t\,({\rm s})")
            ylabel!(L"\dot{\phi} (\circ /{\rm s})")
            plot!(tᵣₛ, 0*ϕᵣₛ, line=:dash, color="black", label=L"\dot{\phi}_r")
            if save_plots
                savefig(p, "plots/phi_dot_no_observer.pdf")
            else
                display(p)
            end

            p = plot(ts[1:step:end], Xs[3,1:step:end]*180/pi, label=L"\dot{\theta}(t)", color=:dodgerblue3, legend=:bottomright, tickfontsize=font_size, legendfontsize=font_size, guidefontsize=font_size)
            xlabel!(L"t\,({\rm s})")
            ylabel!(L"\dot{\theta} \left(\circ /{\rm s}\right)")
            plot!(tᵣₛ, 0*θᵣₛ, line=:dash, color="black", label=L"\dot{\theta}_r")
            if save_plots
                savefig(p, "plots/theta_dot_no_observer.pdf")
            else
                display(p)
            end

            p = plot(ts[1:step:end], Xs[5,1:step:end], label=L"v(t)", color=:dodgerblue3, legend=:topright, tickfontsize=font_size, legendfontsize=font_size, guidefontsize=font_size)
            xlabel!(L"t\,({\rm s})")
            ylabel!(L"v ({\rm m/s})")
            plot!(tᵣₛ, 0*θᵣₛ*180/pi, line=:dash, color="black", label=L"v_r")
            if save_plots
                savefig(p, "plots/v_no_observer.pdf")
            else
                display(p)
            end
        end
    end
    
end