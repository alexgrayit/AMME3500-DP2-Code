module Simulate
    
    using Symbolics, LinearAlgebra, Random, Distributions, ControlSystems, Main.ss_calcs
    export simulate



    function simulate(t₀, tₑ, Δt, A, B, C, D, K, X₀, Xᵢ; L=zeros(size(X₀)[1], size(C)[1]),  Cₙ = zero(C), X̂ᵢ=Xᵢ, d = zero(X₀))

        ΔXᵢ = Xᵢ - X₀;
        ΔX̂ᵢ = X̂ᵢ - X₀;
        ts = Vector{Float64}(range(start=t₀, step=Δt, stop=tₑ));

        ΔX = ΔXᵢ
        ΔX̂ = ΔX̂ᵢ
        ΔXₓ = ΔX

        Normal = Distributions.Normal();
        n = zero(X₀);
        # if size(n)[1] == 1
        #     n = Float64(n[1]);
        # end 

        Xs = zeros(size(X₀)[1], size(ts)[1]);
        X̂s = zeros(size(X₀)[1], size(ts)[1]);

        function u(Xᵤ)
            U = -K*Xᵤ;
            # if size(uₒ)[1] == 1
            #     uₒ = Float64(uₒ[1]);
            # end
            return U
        end

        function diffX̂(t, X, X̂)
            return A*X̂ + B*u(X̂) + L*C*X + L*Cₙ*n - L*C*X̂;
            # return A*X̂ + B*u(X̂) + L*C*X + L*Cₙ*n - L*C*X̂;
        end

        function diffX(t, X, X̂)
            return A*X + B*u(X̂) + d
            # return A*X + B*u(X̂) + d
        end

        function rk4step(t, X, X̂, dt, dXfunc, dX̂func)

            hdt = dt/2;
            f1 = dXfunc(t, X, X̂);
            f̂1 = dX̂func(t, X, X̂);
            
            th = t + dt;
            Xtemp = X + hdt*f1;
            X̂temp = X̂ + hdt*f̂1;
            f2 = dXfunc(th, Xtemp, X̂temp);
            f̂2 = dX̂func(t, Xtemp, X̂temp);

            Xtemp = X + hdt*f2;
            X̂temp = X̂ + hdt*f̂2;
            f3 = dXfunc(th, Xtemp, X̂temp);
            f̂3 = dX̂func(th, Xtemp, X̂temp);

            tf =  t + dt;
            Xtemp = X + dt*f3;
            X̂temp = X + dt*f̂3;
            f4 = dXfunc(tf, Xtemp, X̂temp);
            f̂4 = dX̂func(tf, Xtemp, X̂temp);
        
            Xnext = X + dt*(f1 + 2*f2 + 2*f3 + f4)/6;
            X̂next = X̂ + dt*(f̂1 + 2*f̂2 + 2*f̂3 + f̂4)/6;
            return [Xnext, X̂next]
        end

        
        for i ∈ 1:length(ts)
            t = ts[i];

            # Step Function
            n = rand(Normal, size(n)[1])
            ΔX, ΔX̂ = rk4step(t, ΔX, ΔX̂, Δt, diffX, diffX̂);

            Xs[:, i] = ΔX + X₀;
            X̂s[:, i] = ΔX̂ + X₀;

        end


        return [ts, Xs, X̂s]
    end

    function simulate_multiple_references(t_arr, Δt, X₀ₛ, u₀ₛ, Xₜ₀, X̂ₜ₀, X, u, Aⱼ, Bⱼ, C, D, ζ, tₛ; ζₒ=0, tₛₒ=0, d = zero(Xₜ₀), Cₙ = zero(C))
        
        ts = Vector{Float64}();
        for i ∈ 1:(size(t_arr)[1] - 1)
            # println("test")
            tₛᵢ = Vector{Float64}(range(start=t_arr[i], step=Δt, stop=t_arr[i+1]-Δt/2));
            # println(tₛᵢ)
            append!(ts, tₛᵢ);
        end
        # append!(ts, t_arr[size(t_arr)[1]])
        
        Xs = zeros(size(Xₜ₀)[1], size(ts)[1]);
        X̂s = zeros(size(Xₜ₀)[1], size(ts)[1]);
        # print(ts)

        systems = [];
        CLsystems = [];
        
        # Counter for steps
        N = 1
        for i ∈ 1:size(X₀ₛ)[1]
            t₀ = t_arr[i];
            tₑ = t_arr[i+1] - Δt/2;
            X₀ = X₀ₛ[i];
            u₀ = u₀ₛ[i];
            X0ᵢsub = Dict([X[j] => X₀[j] for j ∈ 1:size(X)[1]]);
            u0ᵢsub = Dict([u[j] => u₀[j] for j ∈ 1:size(u)[1]]);
            linear_subᵢ = merge(X0ᵢsub, u0ᵢsub);

            A = Symbolics.substitute(Aⱼ, linear_subᵢ);
            B = Symbolics.substitute(Bⱼ, linear_subᵢ);
            A = Matrix{Float64}(A);
            B = Matrix{Float64}(B);
            if size(B)[2] == 1
                B = B[:,1];
            end
            K, L, T = ss_calcs.ctrb_obsv_gains(A, B, C, D, ζ, tₛ, ζₒ, tₛₒ);

            sys = ss(A, B, C, D);
            CLsys = ss((A - B*K), B, C, D);

            push!(systems, sys);
            push!(CLsystems, CLsys);

            # Simulate
            tsᵢ, Xsᵢ, X̂sᵢ = simulate(t₀, tₑ, Δt, A, B, C, D, K, X₀, Xₜ₀, L=L, X̂ᵢ=X̂ₜ₀, d=d, Cₙ=Cₙ);
            Nᵢ = size(tsᵢ)[1];
            N⁺ = N + Nᵢ;
            # println(size(Xsᵢ))
            # println(size(Xs))
            Xs[:, N:N⁺-1] = Xsᵢ;
            X̂s[:, N:N⁺-1] = X̂sᵢ;
            
            # println((A - B*K)*(Xₜ₀ - X₀))
            # Update initial condition
            Xₜ₀ = Xsᵢ[:, end];
            X̂ₜ₀ = X̂sᵢ[:, end];
            

            # Update N Counter
            N = N⁺;
        end
        return ts, Xs, X̂s, systems, CLsystems
    end

end

