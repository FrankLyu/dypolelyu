# Dy polarizability

Create: October 3, 2022 2:55 PM
Date: 10/03/2022
Last modify: November 22, 2022 4:11 PM
Month: Oct
Year: 2022
week: 12

The note is according to 
1. Julius’ code (of the bilayer cartoon)
2. Some reference papers
3. PhD thesis 

# Note and revise to the code

### Atomic unit (a.u.)

The database and the code is written in the special unit, but it is widely used in the polarizability measurement.

```julia
linesdat=readdlm("5I8_Spectre_Eexp.txt")
λs = 2*π * 137 ./ linesdat[:, 3]
redMatEls = sqrt.( linesdat[:, 4] )
Γs = linesdat[:, 5]
...
```

The linesdat is a huge csv, recordind transitions between the ground state and 608 excited states. If we define the quantities in SI units with primes, the conversion is defined as:

 

$$
\begin{aligned}
\lambda&=\lambda^{\prime}/a_0\\
\omega&=\omega^{\prime}a_0/\alpha c\\
\lambda&=2\pi/\alpha\omega
\end{aligned}
$$

where $1/\alpha \approx 137, \,a_0 \approx 0.0529\,  \mathrm{nm}.$

One can first verify that the units looks clearer after translating back to SI.

### Linewidth

When dealing with strong D1, D2 lines of alkali, we may assume that we only need to know the linewidth, because we can get the transition strength via 

$$
\frac{1}{\tau}=\frac{\omega_0^3}{3 \pi \epsilon_0 \hbar c^3} \frac{2 J+1}{2 J^{\prime}+1}\left|\left\langle J\|e \mathbf{r}\| J^{\prime}\right\rangle\right|^2 . \tag{\dagger}
$$

Another common convention is to absorb the $(2J+1)$  into the reduced density matrix element square, which is also adapted in the code. For the 421 nm transition, our redefined reduced matrix elements are $\sqrt{2J+1}\left|\left\langle J\|e \mathbf{r}\| J^{\prime}\right\rangle\right|\approx 12 A.U.$ The atomic unit here is $e^2a_0^2/E_h$. Take these into account, we can verify that for the parameters, the 421 nm transition is about 30 MHz, exactly as we know and given by the last column of the csv.

However, some data in the chart does not follow the equation ($\dagger$). This is because the equation only stands if the transition is closed. For an arbitary transition, the excited state population may also leak into other dark state. In general, $1/\tau$ is smaller than the right hand side of the equation, which only considers the transition between the excited state with ground state.

<aside>
⚠️ Some linewidth and reduced matrix elements given by the csv is wrong. The linewidth of the 741 nm transition is 10 times smaller than 1.7 kHz. The 626 nm transition is 2 times larger. I verified that some of the previous paper of Maxence Lepers’ group gave a much smaller prediction for 741 nm transition. So I revised the most important 626 nm and 741 nm transition in the chart.

</aside>

### Verification of the revised data

We can verify the validity of the data by comparing the code result with some special features.

|  | 741 nm tuneout $\lambda_1$ | 741 nm tuneout $\lambda_2$ |  1064 nm scalar | 1064 nm tensor |
| --- | --- | --- | --- | --- |
| code | +6.7 GHz | + 30.5 GHz | 193 a.u. | 1.5 a.u. |
| measure | +7.7 GHz [3] | + 29.9 GHz [3] | 184 a.u. [1] | 1.7 a.u. [1] |

These tuneout wavelengths exist because of a non zero background polarizability around 741 nm transition.

The prediction is not so good as alkali prediciton. But it is already much better than some physical models made around 2010 [6], partly thanks to the collaboration between theory and experiment groups.

### The structure of the main part

The main part of the code is to calculate the potential structure of the energy shift of the ground state $|J=8,m_J\rangle$ in the light field.  However in a future code, we also need to consider magnetic field. Julius used three functions 

```julia
function uu(K, q, pol)
function α_tot(ω, J, mJ, Jp, mJp, pol; Kmax = 2)
function α_tot_precalc(αs, J, mJ, Jp, mJp, pol; Kmax = 2)
```

to rewrite equation (10)(11)(12) in Le Kien [5] (withouting multiplying the amplitude of the electric field). A probable mistake is that he took the real part of the elements of the polarizability matrix. However in general, they could be complex number. So I now turn to a complex number matrix **zeros(ComplexF64,17, 17)** and Hermitian matrix.

$$
\begin{aligned}V_{F M F^{\prime} M^{\prime}}^{E E}=& \frac{1}{4}|\mathcal{E}|^2 \sum_{\substack{K=0,1,2 \\q=-K, \ldots, K}} \alpha_{n J}^{(K)}\left\{\mathbf{u}^* \otimes \mathbf{u}\right\}_{K q} \\& \times(-1)^{J+I+K+q-M} \sqrt{(2 F+1)\left(2 F^{\prime}+1\right)} \\& \times\left(\begin{array}{ccc}F & K & F^{\prime} \\M & q & -M^{\prime}\end{array}\right)\left\{\begin{array}{lll}F & K & F^{\prime} \\J & I & J\end{array}\right\}\end{aligned} \tag{10}
$$

However, I would prefer to use equation (11)(15)(16) in a new code because they seem less mathematical.

$$
\begin{aligned}&V^{E E}=-\frac{1}{4}|\mathcal{E}|^2\left\{\alpha_{n J}^s-i \alpha_{n J}^v \frac{\left[\mathbf{u}^* \times \mathbf{u}\right] \cdot \mathbf{J}}{2 J}\right. \\&\left.\quad+\alpha_{n J}^T \frac{3\left[\left(\mathbf{u}^* \cdot \mathbf{J}\right)(\mathbf{u} \cdot \mathbf{J})+(\mathbf{u} \cdot \mathbf{J})\left(\mathbf{u}^* \cdot \mathbf{J}\right)\right]-2 \mathbf{J}^2}{2 J(2 J-1)}\right\} .\end{aligned} \tag{15}
$$

The following functions are then constructed to get the eigenvalues of matrix $[\alpha]_{ij}$:

```julia
function normalized_polarization_eigenstates(ω, pol; Kmax = 2, offdiag = true, return_vecs = false)
```

The rest of the code—get_levels(…) and so on is specified around the bilayer setup and cartoon generation.

<aside>
⚠️ A probable error in the function get_levels() is using **σp_beam[i] + σm_beam[i]** to calculate the total beam intensity. I changed this line to **σp_beam[i]^2 + σm_beam[i]^2**.
</aside>

# Useful Equation in [1][2][7]

In a simple two-level system, there is only scalar shift.

$$
V^{EE}=\frac{3}{16\pi^2c}\frac{\lambda^3\Gamma}{\Delta\omega}\times I \tag{a}
$$

One can check that equation (a) is the same as equation (1) in [2].

$$
V_i(\mathbf{r}, \omega)=-\frac{2 \pi a_0^3}{c} \tilde{\alpha}_i(\omega) I(\mathbf{r}) \tag{1}
$$

where $\alpha_i$ is in dimensionless polarizability in atomic unit and $a_0$ is the Bohr radius. This equation is easier to use if everything is calculated in atomic units.

Two simple examples in [7]:

1. $\pi$ beam

$$
V\left(m_F\right)=V_0\left(\alpha_F^s+\alpha_F^t \frac{3 m_F^2-F(F+1)}{F(2 F-1)}\right) .
$$

1. $\sigma \pm$ beam

$$
V\left(m_F\right)=V^{EE}\left(\alpha_F^s \pm \alpha_F^v \frac{m_F}{2 F}-\alpha_F^t \frac{3 m_F^2-F(F+1)}{2 F(2 F-1)}\right)
$$

For ideal transition close to 626 (741) nm, $\alpha^s=19/51$, $\alpha^v=152/153$, $\alpha^t=-40/153$. For a σ+ beam , +8 atoms feel exactly the same shift as in (a), but -8 atoms feel 1/153 of it. 
One can also check this in the code.

# Unknown part

What is the AC stark shift for excited states? People don’t know the exact line structure of rare-earth atoms yet. There is some improvement on the theoretical side over the past ten years from [6] to [8][9]. People argued that with the better model, we can calculate the stark shift more accurate and had realized a magic wavelength+polarization transition between GS and 626 intercombination line in 2018[4]. Most theoretical work is done in the group of VV Flambaum, Maxence Lepers and Oliver Dulieu.

# The relevent papers of the last 12 years

### Experiment Side

[1] Rudy Grimm, Accurate Determination of the Dynamical Polarizability of Dysprosium, PRL, 2018

[2] Rudy Grimm, Measurement of the dynamic polarizability of Dy atoms near the 626-nm intercombination line, PRA, 2021
[3] Ben Lev, Anisotropic dependence of tune-out wavelength near Dy 741-nm transition, OE, 2017
[4] Maxence Lepers & Sylvain, Anisotropic light-shift and magic-polarization of the intercombination line of Dysprosium atoms in a far-detuned dipole trap, PRA, 2018

### Theory Side

[5] Le Kien, Dynamical polarizability of atoms in arbitrary light fields: general theory and application to cesium, 2013
[6] V V Flambaum, Theoretical study of some experimentally relevant states of dysprosium, 2010
[7] Davide Dreon @Sylvain, Designing and building an ultracold Dysprosium experiment : a new framework for light-spin interaction, PhD thesis, 2017
[8] Li Hui, Maxence Lepers, Optical trapping of ultracold dysprosium atoms:transition probabilities, dynamic dipole polarizabilities and van der Waals C6 coefficients, 2017
[9] Li Hui, Maxence Lepers, Anisotropic optical trapping as a manifestation of the complex electronic structure of ultracold lanthanide atoms: The example of holmium, 2017
