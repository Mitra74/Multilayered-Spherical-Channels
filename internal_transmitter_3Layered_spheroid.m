%% =======================================================================
%  Internal Point-Source Diffusion in a Three-Layer Tumor Spheroid
% ==========================================================================
%
%  Computes the time-domain molecular concentration profile at a chosen
%  observation point, for a point-source transmitter positioned INSIDE
%  the spheroid (in the necrotic core, Layer 1).
%
%  This implements the analytical (Green's function) solution derived in
%  Section IV-A2, "Internal Point Source Delivery," of:
%
%     M. Rezaei, M. Chappell, and A. Noel, "Molecular Communication
%     Model for Drug Delivery in Multi-Layered Spherical Channels."
%
%  Equation numbers referenced in the comments below correspond to that
%  paper. A companion script implements the EXTERNAL transmitter case
%  (Section IV-A1).
%
%  METHOD OVERVIEW
%  ----------------
%   1. The spheroid is modeled as three concentric homogeneous layers
%      (necrotic core / hypoxic layer / proliferating layer), embedded
%      in an infinite outer medium, each with its own effective
%      diffusion coefficient D_i = eps_i^1.5 * D (Eq. 2).
%   2. For every angular frequency omega and every spherical-harmonic
%      degree l = 0,...,N, an 8x8 linear system is assembled from:
%        - concentration/flux continuity at the two inter-layer
%          boundaries R1, R2 and at the spheroid surface R3 (Eqs. 38-43)
%        - continuity/jump conditions at the source location r_tx,
%          which sits inside Layer 1 (Eqs. 44-45)
%      and solved for the unknown modal coefficients.
%   3. The frequency-domain concentration at the observation point
%      (r, theta, phi) is reconstructed by summing the radial solution
%      over spherical harmonics (Eqs. 11-17).
%   4. Sweeping omega and applying an inverse FFT converts the
%      frequency-domain response into a time-domain concentration
%      curve, C(t), at the observation point.
%
%  This reproduces (for r_tx = 45.83 micron) the analytical curves shown
%  in Fig. 3 of the paper, validated there against particle-based
%  simulation.
%
%  OUTPUT
%  -------
%   Conc_dc_Lx : time-domain concentration profile at the observation
%                point, sampled at the time values in vector t.
%
% ==========================================================================

clc
clear

%% 1. Series truncation order ----------------------------------------------
% Number of spherical-harmonic degrees l = 0,...,N kept in the truncated
% modal expansion. Larger N gives higher accuracy at the cost of runtime.
N = 10;

%% 2. Spheroid geometry -------------------------------------------------------
R_SPH = 275e-6;     % spheroid (outer) radius [m]
R_rx  = R_SPH;      % alias used below when defining layer boundaries

%% 3. Layer-3 porosity values used for the paper's sensitivity study --------
% Only the first value (eps_L3 = 0.1697) is used in this script. The other
% two values are used elsewhere to sweep outer-layer porosity (Figs. 4-8).
different_eps_R3 = {0.1697, 0.1, 0.08};

%% 4. Transmitter position: INTERNAL, inside Layer 1 (necrotic core) --------
% r_tx = R_SPH/6 = 45.83 micron -> matches the internal-release case in Fig. 3
r_tx     = (1/3)*R_SPH - (1/6)*R_SPH;
phi_tx   = 3*pi/2;
theta_tx = pi/2;

%% 5. Drug / drug-carrier release parameters (contextual only) ---------------
% Not used in the concentration calculation below; kept to document the
% physical scenario (number of carriers released, drug payload per carrier).
N_dc    = 10e6;   % number of drug carriers released
N_drugs = 100;     % number of drug molecules carried per carrier

%% 6. Spheroid cellular parameters (contextual only) --------------------------
Nc_SPH = 24000;                 % total number of cells in the spheroid
V_SPH  = (4*pi/3)*R_SPH^3;      % spheroid volume [m^3]
vc     = 3.14e-15;              % single-cell volume [m^3] (from [30])

%% 7. Observation points -------------------------------------------------------
% Radius at the mid-point of each layer, plus one point outside the spheroid.
r_L3  = (2/3)*R_SPH + (1/6)*R_SPH;   % mid Layer 3 (outer, proliferating)
r_L2  = (1/3)*R_SPH + (1/6)*R_SPH;   % mid Layer 2 (hypoxic)
r_L1  = (1/3)*R_SPH - (1/6)*R_SPH;   % mid Layer 1 (necrotic core)
r_out = 300e-6;                       % point outside the spheroid

Observation_r = {r_L1, r_L2, r_L3, r_out};
obs_names = {'Layer 1 (necrotic core)', 'Layer 2 (hypoxic)', ...
             'Layer 3 (proliferating)', 'Outside spheroid'};

% ---- Select which observation point to evaluate (1-4, see obs_names) ----
obs_idx = 2;                          % Here: Layer 2
r     = Observation_r{obs_idx};
phi   = 7*pi/4;
theta = pi/2;

%% 8. Layer-specific degradation rates -----------------------------------------
% Direct drug-release scenario: no chemical degradation assumed in any layer
% (all rates are zero here; nonzero values are used only for the pH-triggered
% carrier-release scenario in Section IV-B).
kd_d_L3 = 0;  kd_d_L2 = 0;  kd_d_L1 = 0;                    % drug decay


%% 9. Bulk (free-medium) diffusion coefficients ---------------------------------
D_dc    = 1e-9;    % drug-carrier diffusion coefficient in free medium [m^2/s]
D_drugs = 10e-9;   % drug diffusion coefficient in free medium [m^2/s]

%% 10. Number of concentric spheroid layers --------------------------------------
NR = 3;

%% 11. Layer 3 -- outer, proliferating layer ---------------------------------------
vc_L3  = 3.14e-15;
V_L3   = (4*pi/3)*(((NR-1)*R_SPH/NR)^3 - ((NR-2)*R_SPH/NR)^3);
eps_L3 = different_eps_R3{1};                    % porosity, Eq. (1)
%Nc_L3  = ceil(V_L3*(1-eps_L3)/vc_L3);             % estimated cell count

%% 12. Layer 2 -- hypoxic layer -----------------------------------------------------
vc_L2  = 3.14e-15;
V_L2   = (4*pi/3)*(((NR-2)*R_SPH/NR)^3 - ((NR-3)*R_SPH/NR)^3);
eps_L2 = 0.1196;
%Nc_L2  = ceil(V_L2*(1-eps_L2)/vc_L2);

%% 13. Layer 1 -- necrotic core -----------------------------------------------------
vc_L1  = 3.14e-15;
V_L1   = (4*pi/3)*(((NR-3)*R_SPH/NR)^3 - ((NR-4)*R_SPH/NR)^3);
eps_L1 = 0.2964;
%Nc_L1  = ceil(V_L1*(1-eps_L1)/vc_L1);

%% 14. Effective diffusivities and inter-layer coupling coefficients ---------------
% Effective diffusivity per layer, Eq. (2): D_i = eps_i^1.5 * D
Deff_L1 = eps_L1^1.5 * D_drugs;
Deff_L2 = eps_L2^1.5 * D_drugs;
Deff_L3 = eps_L3^1.5 * D_drugs;

% alpha_Ri = 1/kappa_i, where kappa_i = sqrt(D_{i+1}/D_i) is the interface
% partition coefficient defined below Eq. (4).
alpha_R1 = sqrt(Deff_L1/Deff_L2);
alpha_R2 = sqrt(Deff_L2/Deff_L3);
alpha_R3 = sqrt(Deff_L3/D_drugs);

% Layer boundary radii: 0 < R1 < R2 < R3 = R_SPH
R1 = (1/3)*R_rx;
R2 = (2/3)*R_rx;
R3 = R_rx;

%% 15. Frequency-sweep / inverse-FFT time-axis setup -------------------------------
% The frequency-domain solution is evaluated on a uniform frequency grid and
% later converted to the time domain via an inverse FFT.
ww = 1e-6;      % frequency step [Hz] % decrease this to have smoother results
fw = 0.5;       % maximum frequency [Hz]
omega0 = ww:ww:fw;
domega = omega0(2) - omega0(1);
t = 0:1/max(omega0):1/domega*2 + domega;    % time axis of the reconstructed signal

%% 16. Spherical Bessel / Hankel function handles ------------------------------------
% jn, yn : spherical Bessel functions of the first / second kind
% hn     : spherical Hankel function of the second kind (outgoing-wave solution)
% *d/*dn : their radial derivatives, via the standard recurrence relation
SphBess   = @(x,n) (pi/2./x).^0.5 .* besselj(n+0.5, x);
SphHank   = @(x,n) (pi/2./x).^0.5 .* besselh(n+0.5, 2, x);
SphHankd  = @(x,n) (2*n+1).^-1 .* (n.*SphHank(x,n-1)  - (n+1).*SphHank(x,n+1));
SphBessd  = @(x,n) (2*n+1).^-1 .* (n.*SphBess(x,n-1)  - (n+1).*SphBess(x,n+1));
SphBessn  = @(x,n) (pi/2./x).^0.5 .* bessely(n+0.5, x);
SphBessdn = @(x,n) (2*n+1).^-1 .* (n.*SphBessn(x,n-1) - (n+1).*SphBessn(x,n+1));

%% 17. Main computation: solve the boundary-value problem at every frequency ---------
jj = 1;
tic
for omega = omega0*pi

    % Complex, degradation-adjusted diffusion parameters, Eqs. (6),(9)
    kdp_R1 = -(kd_d_L1 + 1i*omega) / Deff_L1;
    kdp_R2 = -(kd_d_L2 + 1i*omega) / Deff_L2;
    kdp_R3 = -(kd_d_L3 + 1i*omega) / Deff_L3;
    kdpo   = -(0        + 1i*omega) / D_drugs;

    for l = 0:N

        % ----------------------------------------------------------------
        % Coefficient matrix for spherical-harmonic degree l.
        % Unknown vector: [Cn1 An1 Bn1 An2 Bn2 An3 Bn3 Dn]'
        %   Cn1     : Layer-1 solution, region r        < r_tx   (Eq. 36)
        %   An1,Bn1 : Layer-1 solution, region r_tx < r < R1     (Eq. 36)
        %   An2,Bn2 : Layer-2 solution, region R1  < r < R2      (Eq. 37)
        %   An3,Bn3 : Layer-3 solution, region R2  < r < R3      (Eq. 37)
        %   Dn      : outer (infinite) medium, region R3 < r     (Eq. 37)
        %
        %   Rows 1-2 : concentration + flux continuity at R1   (Eqs. 38-39)
        %   Rows 3-4 : concentration + flux continuity at R2   (Eqs. 40-41)
        %   Rows 5-6 : concentration + flux continuity at R3   (Eqs. 42-43)
        %   Rows 7-8 : continuity + jump condition at r_tx     (Eqs. 44-45)
        % ----------------------------------------------------------------
  % Cn1                                                 An1                                               Bn1                                                   An2                                                 Bn2                                                   An3                                                 Bn3                                                 Dn
CoeffMat = [ ...
    0                                                , alpha_R1*SphBess(sqrt(kdp_R1)*R1,l)                , alpha_R1*SphBessn(sqrt(kdp_R1)*R1,l)                , -SphBess(sqrt(kdp_R2)*R1,l)                       , -SphBessn(sqrt(kdp_R2)*R1,l)                      , 0                                                 , 0                                                 , 0 ; ...
    0                                                , Deff_L1*sqrt(kdp_R1)*SphBessd(sqrt(kdp_R1)*R1,l)   , Deff_L1*sqrt(kdp_R1)*SphBessdn(sqrt(kdp_R1)*R1,l)   , -Deff_L2*sqrt(kdp_R2)*SphBessd(sqrt(kdp_R2)*R1,l) , -Deff_L2*sqrt(kdp_R2)*SphBessdn(sqrt(kdp_R2)*R1,l), 0                                                 , 0                                                 , 0 ; ...
    0                                                , 0                                                  , 0                                                   , alpha_R2*SphBess(sqrt(kdp_R2)*R2,l)               , alpha_R2*SphBessn(sqrt(kdp_R2)*R2,l)              , -SphBess(sqrt(kdp_R3)*R2,l)                       , -SphBessn(sqrt(kdp_R3)*R2,l)                      , 0 ; ...
    0                                                , 0                                                  , 0                                                   , Deff_L2*sqrt(kdp_R2)*SphBessd(sqrt(kdp_R2)*R2,l)  , Deff_L2*sqrt(kdp_R2)*SphBessdn(sqrt(kdp_R2)*R2,l) , -Deff_L3*sqrt(kdp_R3)*SphBessd(sqrt(kdp_R3)*R2,l) , -Deff_L3*sqrt(kdp_R3)*SphBessdn(sqrt(kdp_R3)*R2,l), 0 ; ...
    0                                                , 0                                                  , 0                                                   , 0                                                 , 0                                                 , alpha_R3*SphBess(sqrt(kdp_R3)*R3,l)               , alpha_R3*SphBessn(sqrt(kdp_R3)*R3,l)              , -SphHank(sqrt(kdpo)*R3,l) ; ...
    0                                                , 0                                                  , 0                                                   , 0                                                 , 0                                                 , Deff_L3*sqrt(kdp_R3)*SphBessd(sqrt(kdp_R3)*R3,l)  , Deff_L3*sqrt(kdp_R3)*SphBessdn(sqrt(kdp_R3)*R3,l) , -D_drugs*sqrt(kdpo)*SphHankd(sqrt(kdpo)*R3,l) ; ...
    SphBess(sqrt(kdp_R1)*r_tx,l)                     , -SphBess(sqrt(kdp_R1)*r_tx,l)                      , -SphBessn(sqrt(kdp_R1)*r_tx,l)                      , 0                                                 , 0                                                 , 0                                                 , 0                                                 , 0 ; ...
    r_tx^2*sqrt(kdp_R1)*SphBessd(sqrt(kdp_R1)*r_tx,l), -r_tx^2*sqrt(kdp_R1)*SphBessd(sqrt(kdp_R1)*r_tx,l) , -r_tx^2*sqrt(kdp_R1)*SphBessdn(sqrt(kdp_R1)*r_tx,l) , 0                                                 , 0                                                 , 0                                                 , 0                                                 , 0 ];
       
   % Right-hand side: unit-strength source jump condition at r_tx, Eq. (44)
        RHV = [0; 0; 0; 0; 0; 0; 0; 1/Deff_L1];

        % Solve for the modal coefficients. pinv() is used instead of the
        % direct solve (\) for robustness against near-singular matrices
        % that can occur at some (l, omega) combinations.
        AS = pinv(CoeffMat) * RHV;

        Cn1(l+1) = AS(1);
        An1(l+1) = AS(2);
        Bn1(l+1) = AS(3);
        An2(l+1) = AS(4);
        Bn2(l+1) = AS(5);
        An3(l+1) = AS(6);
        Bn3(l+1) = AS(7);
        Dn(l+1)  = AS(8);
    end

    ll = 0:N;

    % --------------------------------------------------------------------
    % Select the radial solution branch matching the observation radius r
    % (Eqs. 36-37).
    % --------------------------------------------------------------------
    if r <= R3
        if r > R2
            CR = An3.*SphBess(sqrt(kdp_R3)*r,ll) + Bn3.*SphBessn(sqrt(kdp_R3)*r,ll);
        elseif (r > R1) && (r <= R2)
            CR = An2.*SphBess(sqrt(kdp_R2)*r,ll) + Bn2.*SphBessn(sqrt(kdp_R2)*r,ll);
        elseif (r > r_tx) && (r <= R1)
            CR = An1.*SphBess(sqrt(kdp_R1)*r,ll) + Bn1.*SphBessn(sqrt(kdp_R1)*r,ll);
        elseif r <= r_tx
            CR = Cn1.*SphBess(sqrt(kdp_R1)*r,ll);
        end
    else
        CR = Dn.*SphHank(sqrt(kdpo)*r,ll);
    end

    % Guard against occasional NaNs (e.g. at the very first, near-zero frequency)
    CR(isnan(CR)) = 0;

    % --------------------------------------------------------------------
    % Angular summation over spherical harmonics (Eqs. 11-17) to obtain the
    % frequency-domain concentration at the observation point (r,theta,phi).
    % --------------------------------------------------------------------
    C(jj) = 0;
    for n = 0:N
        LegP  = legendre(n, cos(theta));       % associated Legendre P_n^m(cos theta)
        LegP0 = legendre(n, cos(theta_tx));    % associated Legendre P_n^m(cos theta_tx)

        for m = 0:n
            Leg  = LegP(m+1);
            Leg0 = LegP0(m+1);

            norm_nm = sqrt((2*n+1)/(4*pi) * factorial(n-m)/factorial(n+m));
            SphHar  = norm_nm * Leg  * cos(m*(phi - phi_tx));
            SphHar0 = norm_nm * Leg0;

            if m == 0
                kk = 1;     % lambda_0 = 1/(2*pi), Eq. (14)
            else
                kk = 2;     % lambda_m = 1/pi for m >= 1
            end

            C(jj) = C(jj) + kk * CR(n+1) * SphHar * SphHar0;
        end
    end

    jj = jj + 1;
end
toc

%% 18. Inverse FFT: frequency-domain response -> time-domain concentration -----
% The one-sided spectrum C(omega) is mirrored into a conjugate-symmetric
% spectrum before the ifft, which guarantees a real-valued time-domain signal.
Conc_dc_Lx = real(ifft([0, C*max(omega0), conj(fliplr(C*max(omega0)))]));

[peak_value, peak_index] = max(Conc_dc_Lx);

%% 19. Plot: concentration vs. time at the chosen observation point -------------
n_samples = 70;                    % number of time samples to display
t_plot = t(1:n_samples);
C_plot = Conc_dc_Lx(1:n_samples);
 
figure('Color', 'w', 'Position', [100 100 720 480]);
plot(t_plot, C_plot, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
grid on
box on
set(gca, 'FontSize', 12, 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 14);
ylabel('Concentration (molecules/m^{3})', 'FontSize', 14);
 
% NOTE on the title string below: MATLAB's default 'tex' text interpreter
% reads backslash-commands greedily, so "\mum" is parsed as one (invalid)
% command "\mum" rather than "\mu" + "m", which throws
% "String ... must have valid interpreter syntax". We avoid the Greek mu
% symbol and the em-dash character entirely and stick to plain ASCII, which
% is also friendlier for terminals, diffs, and non-UTF8 environments.
title_str = sprintf('Internal Point-Source Delivery: r_{tx} = %.2f um, Observation: %s', ...
      r_tx*1e6, obs_names{obs_idx});
title(title_str, 'FontSize', 13, 'Interpreter', 'tex');
 
legend({'Analytical model'}, 'Location', 'best', 'FontSize', 11);
 
% Uncomment to export a publication-ready figure:
% exportgraphics(gcf, 'internal_transmitter_concentration_profile.png', 'Resolution', 300);