using LinearAlgebra

const mₑ = 9.109E-31
const c = 2.997E8
const ħ = 1.05457E-34
const e = 1.6022E-19
const a₀ = 5.29177E-11
const ϵ₀ = 8.8542E-12
const au_pol = 1.65E-41
const Eh = ħ^2 / (mₑ * a₀^2)    # Hartree energy
const Ia = 0;                  # Nuclear spin
const mdy=161.926805*1.6605*1e-27;

ω = 2π*404.59745*10^12;  #Hz
λ = c*2π/ω;
k = ω/c;

Γ = 11.2*10^3; #Hz 
Δ_tw = -2π*1.25*10^9;  #Hz

p = 150*10^-3;  # tweezer power in watt
p1 = p; p2 = p;

Er = ħ^2*k^2/(2*mdy);

const Jgs=8;
Jup=zeros(Float64,2Jgs+1,2Jgs+1);
[Jup[i,i+1]=sqrt( Jgs*(Jgs+1) - mJ*(mJ+1) ) for (i,mJ) ∈ enumerate([-Jgs:(Jgs-1);]) ] ; 

Jdown=transpose(Jup);

Jx = (Jup + Jdown)/2;
Jy = (Jup - Jdown)/2/1im;
Jz = diagm([-8:8;]);
Jopr = [Jx, Jy ,Jz];


ϵ1(γ) = [cos(γ); 1im*sin(γ); 0]
ϵ2(γ) = [-sin(γ); 1im*cos(γ); 0]
ϵz = [0, 0, 1];

α0 = 19/51;
α1 = 152/153;
α2 = -40/153;       # Ideal 8->9 transition
Prefac = -1;

Vs(x)=Prefac*diagm(ones(2Jgs+1)) * α0 * (Efield(x)⋅Efield(x))
        
# dot compute the dot product between two vectors. For complex vectors, the first vector is conjugated.

Vv(x)=-Prefac * α1 * 1im / 2 /Jgs * ((conj(cross(conj(Efield(x)),Efield(x)))'Jopr ))

#" ' " operator turns the vector (matrix) before it into its adjoint part. "Adjoint" is a quite different class

# and can contract with vector (matrix) behind it automatically. 

Vt(x)=begin Prefac * α2 /(2 * Jgs * (2*Jgs-1)) * 
         (-2* Jgs* (Jgs + 1 )*diagm(ones(2Jgs+1))* (Efield(x)⋅Efield(x)) +
         3 * ( (Efield(x)'Jopr)*(conj(Efield(x))'Jopr)+(conj(Efield(x))'Jopr)*(Efield(x)'Jopr)) ) end

V(x)= Vs(x) + Vv(x) + Vt(x);

function get_levels(xdum;return_vecs = false)
    levels = zeros(length(xdum), 2Jgs+1)
    expsz = zeros(length(xdum), 2Jgs+1)
    for i ∈ 1:length(xdum)
        levels[i, :] = eigvals( V(xdum[i]) )
    end
    
    if return_vecs
        for i ∈ 1:length(xdum)
            expsz[i,:]=real.(eigvecs(V(xdum[i])).^2)'*[-Jgs:1.:Jgs;]
        end
        return levels,expsz
    end
    return levels
end

w = 1;
δtw = 0.05;  #bilayer distance in the unit of lattice wavelength
Efield(x) = ϵ2(pi/4)*exp(-(x-δtw/2)^2/w^2) + ϵ1(pi/4)*exp(-(x+δtw/2)^2/w^2)
