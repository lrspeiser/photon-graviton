# Radiation-polarized gravity (RPG-1): a nonlinear field response solved on the frozen baryons

Declared before execution, 13 September 2026. Baseline: `main` at 1e9b9d6. This is the second of three branches opened by the project owner. It replaces the companion reservoir by a nonlinear response of the gravitational field. The equation is the AQUAL form of Bekenstein and Milgrom (1984) with the "simple" interpolating function. That equation is established mathematics and an existing modified-gravity proposal, adopted here as a hypothesis for the response of a radiation-polarized medium. Nothing in it is derived here.

No science output (SPARC, Milky Way or lensing) has been computed before this declaration. Only development runs of the solver on analytic test masses were made.

## Postulates

- **R1. Field equation.** ∇·[μ(|∇Φ|/a*)∇Φ] = 4πGρ_b with μ(x) = x/(1+x).
  - ρ_b contains the frozen baryons only; there is no companion reservoir and no dark matter.
  - The owner's radiation-state argument s, in μ(|∇Φ|/a*, s), is not included in this first version.
- **R2. Scale.** a* = 8.563335e-11 m/s². This is the archived simple-μ value, fitted in the algebraic limit to the 89 SPARC training galaxies (`isotropic-galaxy-transfer/model-comparison.py`). It is frozen, not refitted.
- **R3. Link to PF-1.** a* = ξ c|d ln n/dt|.
  - With the archived α, c|d ln n/dt| = c²α = 7.2496e-10 m/s² today, so ξ = 0.1181. With the supernova-preferred cα = 69.76 km/s/Mpc, ξ = 0.1263.
  - A linear index gives a*(z) = (1+z) a*(0) at the epoch of emission.
- **R4. Declared lensing response.** Light is deflected by the same Φ with the relativistic factor two: α̂ = (2/c²)∫∇⊥Φ dl, i.e. Φ = Ψ. This is declared, not derived; a relativistic completion would have to supply it. Lensing by the Newtonian baryonic potential alone is reported as a bracket.
- **R5. Lensing geometry.** PF-1's static Euclidean geometry: D = ln(1+z)/α, D_A = D, and D_ls = D_s − D_l.

## Inputs

These are frozen. Each is labeled as measured or model-inferred.

- **SPARC, 149 galaxies.** 3,150 radii with the original 89/29/31 split.
  - Measured: rotation speeds, 3.6 μm surface brightness, H I mass.
  - Model-inferred, as in `capture-to-orbit/inputs.py`:
    - stellar surface density with Υ_disk = 0.5 and Υ_bulge = 0.7;
    - a spherical bulge built from 0.7 V_bul²;
    - an exponential gas disk of mass 1.33 M_HI whose scale is fitted to V_gas.
  - Declared vertical structure: exponential in |z| with h_z = 0.1 R_d for both stars and gas. A sensitivity run uses h_z = 0.2 R_d.
- **Milky Way.** The archived variants I and II: Miyamoto–Nagai stellar disks, a Plummer bulge, and holed exponential gas disks with a declared h_z = 0.1 kpc.
  - The 38 Eilers circular speeds are Jeans-inferred.
  - The 43 Bovy–Rix values of K_z/2πG at |z| = 1.1 kpc are inferred from stellar kinematics.
- **Six SLACS lenses.**
  - Measured: Einstein radii and redshifts.
  - Model-inferred: population-synthesis stellar masses (Chabrier and Salpeter) and Sérsic light profiles, deprojected as spheres.

## Tests and declared tolerances

- **V1. Solver validation.** It must pass before any science output is interpreted.
  - (a) Plummer sphere under AQUAL against the exact spherical solution: radial field within 1e-6 relative at 0.01–1000 scale radii.
  - (b) Razor-thin Kuzmin disk under AQUAL against the exact solution g = ν(|g_N|/a*) g_N in each half-space (Brada and Milgrom 1995):
    - midplane circular speed within 0.5% at 0.3–60 kpc;
    - K_z at |z| = 1.1 kpc within 0.5% at R = 4–9 kpc.
  - (c) Miyamoto–Nagai disk, Newtonian: midplane speed and K_z(1.1 kpc) within 0.5%.
  - (d) Convergence on three SPARC galaxies, chosen from catalog properties before solving: the gas-richest, the median-mass and the most bulge-dominated.
    - Doubling both grid resolutions changes the predicted speed by less than 0.5 km/s at every observed radius.
    - Multiplying r_out by 10 changes it by less than 0.1 km/s.
    - The nonlinear residual is below 1e-8.
- **T1. SPARC rotation at frozen a*.** Compute AQUAL midplane speeds at all 3,150 radii, and the equal-galaxy RMSE and log RMS for each split. Compare against:
  - the Newtonian solution for the same baryons;
  - the algebraic simple-MOND relation applied to the same baryons' Newtonian midplane field;
  - the archived algebraic rows, built from SPARC's own rotation contributions.

  Report two distributions: AQUAL minus algebraic (the effect of the field equation) and the baryon-model difference. The sensitivity run uses h_z = 0.2 R_d.
- **T2. Milky Way.** For both variants, the rotation RMSE on the 38 Eilers bins and the K_z RMS on the 43 Bovy–Rix values. Compare with the Newtonian baryons and the algebraic relation.
- **T3. Lensing.** The predicted Einstein radius of each SLACS lens under R4–R5 for:
  - Chabrier and Salpeter stellar masses;
  - a* held constant, and a*(z_l) = (1+z_l) a*(0) from the link;
  - the Newtonian-baryon bracket.

  **Decision rule.** The declared lensing response is consistent only if, for at least one of the two IMFs, all six predicted-to-observed ratios lie within 0.9–1.1.
- **T4. Field energy and outer boundary.**
  - E_F(<r) = (a*²/8πG)∫F(|∇Φ|²/a*²)dV, with F(y) = y − 2√y + 2 ln(1+√y).
  - Report dE_F/d ln r at r = 10–100 r_M, where r_M = √(GM/a*). Compare it with the deep-field value M v_f²/3 (v_f⁴ = GMa*) and with the Newtonian gradient energy of the same baryons.
  - State the ln r growth of the potential and its consequence for the total energy.
  - For a future radiation-state dependence, report the radius at which each galaxy's own starlight energy density falls to the microwave-background value of 4.17e-14 J/m³. Use the point-source approximation with L_bol = L[3.6]; this is a declared upper bound, and the radius scales as √(L_bol/L[3.6]). Compare that radius with the outermost measured rotation radius.
- **T5. The link.** Report the values of ξ, the predicted evolution, the data that could test it, and the change in the SLACS predictions between constant and evolving a*.

## Assessment rule

- **Validation first.** RPG-1 counts as a usable gravity response only if V1 passes.
- **SPARC comparison.** It uses the same μ and a* as the archived algebraic simple MOND, so it is judged as the field-equation form of that model. It is set beside the reference and MOND-guided branches, and is not promoted over algebraic MOND unless it improves both validation and test.
- **Lensing.** If the T3 rule fails, the declared lensing response is rejected for RPG-1. Other lensing responses are not thereby excluded.
- **The link.** It counts as supported only by a measured redshift dependence. The value of ξ at z ≈ 0 is a coincidence-level statement, like the familiar a0 ~ cH0/2π, not a test.
- **No fitting.** Nothing in this branch is fitted to the data: not a*, not h_z, not the lensing response. The one-third exponent, spherical geometry and effective-mass mapping remain reference hypotheses of the other branches and are not used here.
