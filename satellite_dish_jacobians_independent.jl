using Symbolics, LinearAlgebra, Latexify, ControlSystems, ControlSystemsBase, MAT, Plots

# using Controlz

cpmatrix(v) = [0    -v[3]  v[2]; 
               v[3]   0   -v[1]; 
              -v[2]  v[1]   0  ];

symbolic = false;
if symbolic
    @variables h r m L g m_arm r_arm;
else
    r = 32;
    m = 300000;
    L = 12;
    g = 9.8;

    m_arm = 700000;
    r_arm = 8;

end

@variables θ ϕ;

I_C_dish = [(5//12)*m*(r^2) 0 0;
            0 (5//12)*m*(r^2) 0;
            0 0 (2//3)*m*(r^2)];

I_C_arm = [(1//12)*m_arm*(3*r_arm^2 + L^2) 0 0;
            0 (1//12)*m_arm*(3*r_arm^2 + L^2) 0;
            0 0 (1//2)*m_arm*r_arm^2];

D = [0;0;-L];
R_D = [0;0;-(r//2)]

I_0_dish = I_C_dish - m*cpmatrix(D+R_D)^2;
I_0_arm = I_C_arm - m_arm*cpmatrix(D/2)^2;
I_0 = I_0_dish + I_0_arm;

Ry = [cos(θ) 0 sin(θ);
      0 1 0;
      -sin(θ) 0 cos(θ)];

Rz = [cos(ϕ) -sin(ϕ) 0;
      sin(ϕ)  cos(ϕ) 0;
      0 0 1];

R = Ry;

I_R = R*I_0*transpose(R);

@variables θ̇  ϕ̇ ;
Γ = Symbolics.variables(:Γ, 1:2); # input torques i.e. t1,t2 (about each axis that you control)

# Angular velocity and torques in 3D (xyz), written in terms of the angular velocities and torques about each of your 2 axes you control
# Basically just decomposing the two axes into the amount they are in each x y z axis.
ω = Symbolics.variables(:ω, 1:3);
τi = Symbolics.variables(:τi, 1:3);
τg = cross(-R*D, [0;0;-m*g]) + cross(-R*D/2, [0;0;-m_arm*g]);
τ = τi + τg;

# Substitute the equations for angular velocity 
Wsub = Dict(ω[1] => 0, ω[2] => θ̇ , ω[3] => ϕ̇ );
ω = Symbolics.substitute(ω, Wsub);

Rsub = Dict(τi[1] => 0, τi[2] => Γ[2], τi[3] => Γ[1]);
τ = Symbolics.substitute(τ, Rsub);

α = inv(I_R)*(τ - cross(ω, I_R*ω));

θ̈ = Num((τg[2] + Γ[2])/I_R[2,2]);
# println(latexify(θ̈ ))
ϕ̈ = Num(Γ[1]/I_R[3,3]);
# println(latexify(ϕ̈ ))

# To choose which angle to control
θ_or_ϕ = θ

if θ_or_ϕ === ϕ
    θrs = [0;0;0;0;0;0;0];
    ϕrs = [-290;-110;-45;30;80;130;360]*pi/180;
    θ_or_ϕ = "p";
elseif θ_or_ϕ === θ
    θrs = [-90;-60;-30;30;60;90;110]*pi/180;
    ϕrs = [0;0;0;0;0;0;0]*pi/180;
    θ_or_ϕ = "t";
else exit()
end

# Iterate over references to linearise the system
for i ∈ 1:size(θrs)[1]

    println(i)

    # MIMO System
    global θr = θrs[i];
    global ϕr = ϕrs[i];

    # global linear_sub = merge(Xsub, usub, rsub);
    global X = [θ; ϕ; θ̇ ; ϕ̇ ];
    global Ẋ = [θ̇ ; ϕ̇ ; θ̈ ; ϕ̈ ];
    global u = Γ

    global r₀ = [θr;ϕr;0;0];
    global X₀ = r₀;
    
    global X0sub = Dict([X[i] => X₀[i] for i ∈ 1:size(X)[1]])

    global Ẋ_intermediate = Vector{Num}(Symbolics.substitute(Ẋ, X0sub));
    global u0 = Vector{Num}(Symbolics.solve_for([Ẋ[i] ~ 0 for i ∈ 3:4], u));
    global u₀ = Vector{Num}(Symbolics.substitute(u0, X0sub));
    global u0sub = Dict([u[i] => u₀[i] for i ∈ 1:size(u)[1]]);

    global linear_sub = merge(X0sub, u0sub);

    # τ₁ (ϕ) - Rotation Control
    global X_t01 = [0;0];
    
    global X₁ = [ϕ; ϕ̇ ];
    global Ẋ₁ = [ϕ̇ ; ϕ̈ ];
    global u₁ = [u[1]]

    global r₁ = [ϕr; 0];
    global X₀₁ = r₁;
    global u₀₁ = u₀[1]

    # Jacobians for linearisation
    global Aj₁ = Symbolics.jacobian(Ẋ₁, X₁);
    global Bj₁ = Symbolics.jacobian(Ẋ₁, u₁);


    # Substitutions are taken from the full system
    global A₁ = Symbolics.substitute(Aj₁, linear_sub);
    global B₁ = Symbolics.substitute(Bj₁, linear_sub);

    global A₁ = Matrix{Float64}(A₁);
    global B₁ = Matrix{Float64}(B₁);
    if size(B₁)[2] == 1
        B₁ = B₁[:, 1]
    end

    global r0₁ = Symbolics.substitute(r₁, linear_sub);
    global X0₁ = Symbolics.substitute(X₁, linear_sub);
    global Ẋ0₁ = Symbolics.substitute(Ẋ₁, linear_sub);

    # Controller parameters
    global ζ₁ = 1.1
    # global tₛ₁ = 10
    global tₛ₁ = abs((X_t01[1]-r0₁[1])*900/(2*pi))
    if tₛ₁ == 0
        tₛ₁ = 10
    end
    global ωₙ₁ = log(50)/(tₛ₁*real(ζ₁ - sqrt(Complex(ζ₁^2 - 1))));

    # Controllability
    global W₁ = ctrb(A₁, B₁);

    global s = Symbolics.variable("s");
    global a_poly₁ = expand(det(s*I-A₁))
    global as₁ = [Symbolics.coeff(a_poly₁, s^(size(A₁)[1]-i)) for i ∈ 1:size(A₁)[1]];

    global poly₁ = Symbolics.expand(s^2 + 2*ζ₁*ωₙ₁*s + ωₙ₁^2);
    global ps₁ = [Symbolics.coeff(poly₁, s^(size(A₁)[1]-i)) for i ∈ 1:size(A₁)[1]];

    global Ã₁ = [-transpose(as₁);
                1  0 ];

    global B̃₁ = [1; 0];

    global W̃₁ = ControlSystems.ctrb(Ã₁, B̃₁);

    global T₁ = W̃₁ * inv(W₁);

    global ks₁ = ps₁ - as₁;
    global K̃₁ = transpose(ks₁);

    # Gain matrix
    global K₁ = K̃₁*T₁;

    print("Nullity of W is ")
    println(size(W₁)[1] - rank(W₁))
    if (size(W₁)[1] - rank(W₁)) == 0;
        println("System is controllable")
    else;
        println("System is not controllable")
    end
    println("f(X0, u0) = ")
    println(Ẋ0₁)

    # Oberservability
    global C₁ = [1 0];

    global Wₒ₁=obsv(A₁, C₁)

    if (size(Wₒ₁)[2] - rank(Wₒ₁)) ==0;
        println("System is observable")
    else;
        println("System is not observable")
    end

    # τ₂ (θ) - Elevation Control
    global X_t02 = [0; 0];
    global r₂ = [θr; 0];
    global X₀₂ = r₂;
    global X₂ = [θ; θ̇ ];
    global Ẋ₂ = [θ̇ ; θ̈ ];
    global u₂ = [u[2]];

    # Jacobians for linearisation
    global Aj₂ = Symbolics.jacobian(Ẋ₂, X₂);
    global Bj₂ = Symbolics.jacobian(Ẋ₂, u₂);

    # Substitutions are taken from the full system
    global A₂ = Symbolics.substitute(Aj₂, linear_sub);
    global B₂ = Symbolics.substitute(Bj₂, linear_sub);

    global A₂ = Matrix{Float64}(A₂);
    global B₂ = Matrix{Float64}(B₂);
    if size(B₂)[2] == 1
        B₂ = B₂[:, 1]
    end

    global X0₂ = Symbolics.substitute(X₂, linear_sub);
    global Ẋ0₂ = Symbolics.substitute(Ẋ₂, linear_sub);
    global r0₂ = Symbolics.substitute(r₂, linear_sub);

    global ζ₂ = 1.1
    global tₛ₂ = abs((X_t02[1] - r0₂[1])*300/(2*pi))
    if tₛ₂ == 0
        tₛ₂ = 10
    end
    global ωₙ₂ = log(50)/(tₛ₂*real(ζ₂ - sqrt(Complex(ζ₂^2 - 1))));

    # Controllability
    global W₂ = ctrb(A₂, B₂);

    global a_poly₂ = expand(det(s*I-A₂))
    global as₂ = [Symbolics.coeff(a_poly₂, s^(size(A₂)[1]-i)) for i ∈ 1:size(A₂)[1]];

    global poly₂ = Symbolics.expand(s^2 + 2*ζ₂*ωₙ₂*s + ωₙ₂^2);
    global ps₂ = [Symbolics.coeff(poly₂, s^(size(A₂)[1]-i)) for i ∈ 1:size(A₂)[1]];

    global Ã₂ = [-transpose(as₂);
                1  0 ];
    global B̃₂ = [1; 0];

    global W̃₂ = ControlSystems.ctrb(Ã₂, B̃₂);

    global T₂ = W̃₂ * inv(W₂);

    global ks₂ = ps₂ - as₂;
    global K̃₂ = transpose(ks₂);

    # Gain matrix
    global K₂ = K̃₂*T₂;
    

    print("Nullity of W is ")
    println(size(W₂)[1] - rank(W₂))
    if (size(W₂)[1] - rank(W₂)) == 0;
        println("System is controllable")
    else;
        println("System is not controllable")
    end
    println("f(X0, u0) = ")
    println(Ẋ0₂)

    # Oberservability
    global C₂ = [1 0];

    global Wₒ₂=obsv(A₂, C₂)

    if (size(Wₒ₂)[2] - rank(Wₒ₂)) ==0;
        println("System is observable")
    else;
        println("System is not observable")
    end
    # Print gain matrix for relevant angle
    if θrs[1,1] == 0 && θrs[2,1] == 0
        println("The gain matrix for ϕ is")
        println(K₁)
    elseif ϕrs[1,1] == 0 && ϕrs[2,1] == 0
        println("The gain matrix for θ is")
        println(K₂) 
    end

    if θ_or_ϕ == "p"
        println("Settle time for ϕ")
        println(tₛ₁)
        
    elseif θ_or_ϕ == "t"
        println("Settle time for θ")
        println(tₛ₂)
    end

    # Write variables to a file that can be accessed in matlab
    matwrite("Code DP2/satellite_dish" * string(i) * ".mat", Dict(
        "X01"=>X0₁,
        "A1"=>A₁,
        "B1" => B₁,
        "C1" => C₁,
        "K1" => K₁,
        "X02" => X0₂,
        "A2" => A₂,
        "B2" => B₂,
        "C2" => C₂,
        "K2" => K₂, 
        "t_or_p" => θ_or_ϕ); compress = false)

end