%% =======================================================================
%  External Point-Source Diffusion in a Three-Layer Tumor Spheroid
% ==========================================================================
%
%  Computes the time-domain molecular concentration profile at a chosen
%  observation point, for a point-source transmitter positioned OUTSIDE
%  the spheroid, in the unbounded surrounding medium.
%
%  This implements the analytical (Green's function) solution derived in
%  Section IV-A1, "External Point Source Delivery," of:
%
%     M. Rezaei, M. Chappell, and A. Noel, "Molecular Communication
%     Model for Drug Delivery in Multi-Layered Spherical Channels."
%
%  Equation numbers referenced in the comments below correspond to that
%  paper. A companion script implements the INTERNAL transmitter case
%  (Section IV-A2), where the source sits inside the necrotic core.
%
%  METHOD OVERVIEW
%  ----------------
%   1. The spheroid is modeled as three concentric homogeneous layers
%      (necrotic core / hypoxic layer / proliferating layer), embedded
%      in an infinite outer medium, each with its own effective
%      diffusion coefficient D_i = eps_i^1.5 * D (Eq. 2).
%   2. The transmitter sits in the outer medium (region 4) at r_tx > R3.
%      This splits region 4 into two homogeneous sub-regions:
%        - R3 < r < r_tx  (between the spheroid surface and the source)
%        - r_tx < r        (beyond the source, radiating outward)
%   3. For every angular frequency omega and every spherical-harmonic
%      degree l = 0,...,N, an 8x8 linear system is assembled from:
%        - concentration/flux continuity at the two inter-layer
%          boundaries R1, R2 and at the spheroid surface R3 (Eqs. 28-33)
%        - continuity/jump conditions at the source location r_tx,
%          in the outer medium                              (Eqs. 34-35)
%      and solved for the unknown modal coefficients.
%   4. The frequency-domain concentration at the observation point
%      (r, theta, phi) is reconstructed by summing the radial solution
%      over spherical harmonics (Eqs. 11-17).
%   5. Sweeping omega and applying an inverse FFT converts the
%      frequency-domain response into a time-domain concentration
%      curve, C(t), at the observation point.
%
%  This reproduces (for r_tx = 600 micron) the analytical curve shown in
%  Fig. 2 of the paper, validated there against particle-based simulation.
%
%  OUTPUT
%  -------
%   Conc_rx : time-domain concentration profile at the observation
%             point, sampled at the time values in vector t.
%
% ==========================================================================

clc
clear

%% 1. Series truncation order ----------------------------------------------
% Number of spherical-harmonic degrees l = 0,...,N kept in the truncated
% modal expansion. Larger N gives higher accuracy at the cost of runtime.
N = 10;

%% 2. Layer-3 porosity values used for the paper's sensitivity study --------
% Only the first value (eps_R3 = 0.1697) is used in this script. The other
% two values are used elsewhere to sweep outer-layer porosity (Figs. 4-8).
different_eps_R3 = {0.1697, 0.1, 0.08};

%% 3. Transmitter position: EXTERNAL, outside the spheroid ------------------
% r_tx = 600 micron, well beyond the spheroid surface (R_SPH = 275 micron)
% -> matches the external-release case in Fig. 2 of the paper.
r_tx     = 600e-6;
phi_tx   = 3*pi/2;
theta_tx = pi/2;

%% 4. Spheroid geometry and cellular parameters (contextual only) -----------
Nc_SPH = 24000;                 % total number of cells in the spheroid
R_SPH  = 275e-6;                % spheroid (outer) radius [m]
V_SPH  = (4*pi/3)*R_SPH^3;      % spheroid volume [m^3]
vc     = 3.14e-15;              % single-cell volume [m^3] (from [30])

%% 5. Observation points -------------------------------------------------------
% Radius at the mid-point of each layer, one point outside the spheroid,
% and one point essentially at the spheroid center.
r_center = 0.0001e-6;
r_out    = 400e-6;
r_L3     = (2/3)*R_SPH + (1/6)*R_SPH;   % mid Layer 3 (outer, proliferating)
r_L2     = (1/3)*R_SPH + (1/6)*R_SPH;   % mid Layer 2 (hypoxic)
r_L1     = (1/3)*R_SPH - (1/6)*R_SPH;   % mid Layer 1 (necrotic core)

Observation_r = {r_L1, r_L2, r_L3, r_out, r_center};
obs_names = {'Layer 1 (necrotic core)', 'Layer 2 (hypoxic)', ...
             'Layer 3 (proliferating)', 'Outside spheroid', 'Spheroid center'};

% ---- Select which observation point to evaluate (1-5, see obs_names) ----
obs_idx = 2;                          % default: Layer 2, matches Fig. 2
r     = Observation_r{obs_idx};
phi   = 7*pi/4;
theta = pi/2;

%% 6. Bulk (free-medium) diffusion coefficient -------------------------------
D_drugs = 10e-9;   % drug diffusion coefficient in free medium [m^2/s]

%% 7. Number of concentric spheroid layers -----------------------------------
NR = 3;

%% 8. Layer 3 -- outer, proliferating layer ------------------------------------
vc_R3   = 3.14e-15;
V_R3    = (4*pi/3)*(((NR-1)*R_SPH/NR)^3 - ((NR-2)*R_SPH/NR)^3);
eps_R3  = different_eps_R3{1};                     % porosity, Eq. (1)
Nc_R3   = ceil(V_R3*(1-eps_R3)/vc_R3);              % estimated cell count
Deff_R3 = eps_R3^1.5 * D_drugs;                     % effective D, Eq. (2)
alpha_R3 = sqrt(D_drugs/Deff_R3);                   % kappa_3 = sqrt(D4/D3)
R3 = R_SPH;

%% 9. Layer 2 -- hypoxic layer -------------------------------------------------
vc_R2   = 3.14e-15;
V_R2    = (4*pi/3)*(((NR-2)*R_SPH/NR)^3 - ((NR-3)*R_SPH/NR)^3);
eps_R2  = 0.1196;
Nc_R2   = ceil(V_R2*(1-eps_R2)/vc_R2);
Deff_R2 = eps_R2^1.5 * D_drugs;
alpha_R2 = sqrt(Deff_R3/Deff_R2);                   % kappa_2 = sqrt(D3/D2)
R2 = (2/3)*R_SPH;

%% 10. Layer 1 -- necrotic core -------------------------------------------------
vc_R1   = 3.14e-15;
V_R1    = (4*pi/3)*(((NR-3)*R_SPH/NR)^3 - ((NR-4)*R_SPH/NR)^3);
eps_R1  = 0.2964;
Nc_R1   = ceil(V_R1*(1-eps_R1)/vc_R1);
Deff_R1 = eps_R1^1.5 * D_drugs;
alpha_R1 = sqrt(Deff_R2/Deff_R1);                   % kappa_1 = sqrt(D2/D1)
R1 = (1/3)*R_SPH;

%% 11. Layer-specific degradation rates ---------------------------------------
% Direct drug-release scenario: no chemical degradation assumed in any layer
% (nonzero values are used only for the pH-triggered carrier-release
% scenario in Section IV-B).
kd_R3 = 0;
kd_R2 = 0;
kd_R1 = 0;

%% 12. Frequency-sweep / inverse-FFT time-axis setup --------------------------
% The frequency-domain solution is evaluated on a uniform frequency grid and
% later converted to the time domain via an inverse FFT.
ww = 1e-5;      % frequency step [Hz]
fw = 0.5;       % maximum frequency [Hz]
omega0 = ww:ww:fw;
domega = omega0(2) - omega0(1);
t = 0:1/max(omega0):1/domega*2 + domega;    % time axis of the reconstructed signal

%% 13. Spherical Bessel / Hankel function handles ------------------------------
% jn, yn : spherical Bessel functions of the first / second kind
% hn     : spherical Hankel function of the second kind (outgoing-wave solution)
% *d/*dn : their radial derivatives, via the standard recurrence relation
SphBess   = @(x,n) (pi/2./x).^0.5 .* besselj(n+0.5, x);
SphHank   = @(x,n) (pi/2./x).^0.5 .* besselh(n+0.5, 2, x);
SphHankd  = @(x,n) (2*n+1).^-1 .* (n.*SphHank(x,n-1)  - (n+1).*SphHank(x,n+1));
SphBessd  = @(x,n) (2*n+1).^-1 .* (n.*SphBess(x,n-1)  - (n+1).*SphBess(x,n+1));
SphBessn  = @(x,n) (pi/2./x).^0.5 .* bessely(n+0.5, x);
SphBessdn = @(x,n) (2*n+1).^-1 .* (n.*SphBessn(x,n-1) - (n+1).*SphBessn(x,n+1));

%% 14. Main computation: solve the boundary-value problem at every frequency ---
jj = 1;
tic
for omega = omega0*pi

    % Complex, degradation-adjusted diffusion parameters, Eqs. (6),(9)
    kdp_R1 = -(kd_R1 + 1i*omega) / Deff_R1;
    kdp_R2 = -(kd_R2 + 1i*omega) / Deff_R2;
    kdp_R3 = -(kd_R3 + 1i*omega) / Deff_R3;
    kdpo   = -(0     + 1i*omega) / D_drugs;

    for l = 0:N

        % ------------------------------------------------------------------
        % Coefficient matrix for spherical-harmonic degree l.
        % Unknown vector: [ Cn1  An2  Bn2  An3  Bn3  An  Bn  Dn ]'
        %   Cn1     : Layer-1 solution, 0  < r < R1               (Eq. 27)
        %   An2,Bn2 : Layer-2 solution, R1 < r < R2                (Eq. 27)
        %   An3,Bn3 : Layer-3 solution, R2 < r < R3                (Eq. 27)
        %   An,Bn   : outer-medium solution, R3 < r < r_tx         (Eq. 26)
        %   Dn      : outer-medium solution (outgoing wave), r_tx < r (Eq. 26)
        %
        %   Row 1,2 = concentration & flux continuity at R1    (Eqs. 29, 28)
        %   Row 3,4 = concentration & flux continuity at R2    (Eqs. 31, 30)
        %   Row 5,6 = concentration & flux continuity at R3    (Eqs. 33, 32)
        %   Row 7,8 = concentration continuity & flux jump at r_tx (Eqs. 35, 34)
        % ------------------------------------------------------------------
        CoeffMat = [ ...
            SphBess(sqrt(kdp_R1)*R1,l)                      , -alpha_R1*SphBess(sqrt(kdp_R2)*R1,l)             , -alpha_R1*SphBessn(sqrt(kdp_R2)*R1,l)             , 0                                                , 0                                                 , 0                                            , 0                                             , 0 ; ...
            Deff_R1*sqrt(kdp_R1)*SphBessd(sqrt(kdp_R1)*R1,l), -Deff_R2*sqrt(kdp_R2)*SphBessd(sqrt(kdp_R2)*R1,l), -Deff_R2*sqrt(kdp_R2)*SphBessdn(sqrt(kdp_R2)*R1,l), 0                                                , 0                                                 , 0                                            , 0                                             , 0 ; ...
            0                                               , SphBess(sqrt(kdp_R2)*R2,l)                       , SphBessn(sqrt(kdp_R2)*R2,l)                       , -alpha_R2*SphBess(sqrt(kdp_R3)*R2,l)             , -alpha_R2*SphBessn(sqrt(kdp_R3)*R2,l)             , 0                                            , 0                                             , 0 ; ...
            0                                               , Deff_R2*sqrt(kdp_R2)*SphBessd(sqrt(kdp_R2)*R2,l) , Deff_R2*sqrt(kdp_R2)*SphBessdn(sqrt(kdp_R2)*R2,l) , -Deff_R3*sqrt(kdp_R3)*SphBessd(sqrt(kdp_R3)*R2,l), -Deff_R3*sqrt(kdp_R3)*SphBessdn(sqrt(kdp_R3)*R2,l), 0                                            , 0                                             , 0 ; ...
            0                                               , 0                                                , 0                                                 , SphBess(sqrt(kdp_R3)*R3,l)                       , SphBessn(sqrt(kdp_R3)*R3,l)                       , -alpha_R3*SphBess(sqrt(kdpo)*R3,l)           , -alpha_R3*SphBessn(sqrt(kdpo)*R3,l)           , 0 ; ...
            0                                               , 0                                                , 0                                                 , Deff_R3*sqrt(kdp_R3)*SphBessd(sqrt(kdp_R3)*R3,l) , Deff_R3*sqrt(kdp_R3)*SphBessdn(sqrt(kdp_R3)*R3,l) , -D_drugs*sqrt(kdpo)*SphBessd(sqrt(kdpo)*R3,l), -D_drugs*sqrt(kdpo)*SphBessdn(sqrt(kdpo)*R3,l), 0 ; ...
            0                                               , 0                                                , 0                                                 , 0                                                , 0                                                 , SphBess(sqrt(kdpo)*r_tx,l)                   , SphBessn(sqrt(kdpo)*r_tx,l)                   , -SphHank(sqrt(kdpo)*r_tx,l) ; ...
            0                                               , 0                                                , 0                                                 , 0                                                , 0                                                 , r_tx^2*sqrt(kdpo)*SphBessd(sqrt(kdpo)*r_tx,l), r_tx^2*sqrt(kdpo)*SphBessdn(sqrt(kdpo)*r_tx,l), -r_tx^2*sqrt(kdpo)*SphHankd(sqrt(kdpo)*r_tx,l) ];

        % Right-hand side: unit-strength source jump condition at r_tx, Eq. (34)
        RHV = [0; 0; 0; 0; 0; 0; 0; 1/D_drugs];

        % Solve for the modal coefficients. pinv() is used instead of the
        % direct solve (\) for robustness against near-singular matrices
        % that can occur at some (l, omega) combinations.
        AS = pinv(CoeffMat) * RHV;

        Cn1(l+1) = AS(1);
        An2(l+1) = AS(2);
        Bn2(l+1) = AS(3);
        An3(l+1) = AS(4);
        Bn3(l+1) = AS(5);
        An(l+1)  = AS(6);
        Bn(l+1)  = AS(7);
        Dn(l+1)  = AS(8);
    end

    ll = 0:N;

    % --------------------------------------------------------------------
    % Select the radial solution branch matching the observation radius r
    % (Eqs. 26-27).
    % --------------------------------------------------------------------
    if r <= R3
        if r > R2
            CR = An3.*SphBess(sqrt(kdp_R3)*r,ll) + Bn3.*SphBessn(sqrt(kdp_R3)*r,ll);
        elseif (r > R1) && (r <= R2)
            CR = An2.*SphBess(sqrt(kdp_R2)*r,ll) + Bn2.*SphBessn(sqrt(kdp_R2)*r,ll);
        elseif r <= R1
            CR = Cn1.*SphBess(sqrt(kdp_R1)*r,ll);
        end
    elseif r > r_tx
        CR = Dn.*SphHank(sqrt(kdpo)*r,ll);
    else
        CR = An.*SphBess(sqrt(kdpo)*r,ll) + Bn.*SphBessn(sqrt(kdpo)*r,ll);
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

%% 15. Inverse FFT: frequency-domain response -> time-domain concentration ----
% The one-sided spectrum C(omega) is mirrored into a conjugate-symmetric
% spectrum before the ifft, which guarantees a real-valued time-domain signal.
Conc_rx = real(ifft([0, C*max(omega0), conj(fliplr(C*max(omega0)))]));

[peak_value, peak_index] = max(Conc_rx);

% Clip small negative/near-zero numerical noise that can appear before the
% signal has physically arrived (i.e. before the first peak).
Conc_rx(1:peak_index) = max(Conc_rx(1:peak_index), 1e-10);

%% 16. Plot: concentration vs. time at the chosen observation point -----------
n_samples = 50;                    % number of time samples to display
t_plot = t(1:n_samples);
C_plot = Conc_rx(1:n_samples);

figure('Color', 'w', 'Position', [100 100 720 480]);
plot(t_plot, C_plot, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
grid on
box on
set(gca, 'FontSize', 12, 'LineWidth', 1);

xlabel('Time (s)', 'FontSize', 14, 'Interpreter', 'tex');
ylabel('Concentration (molecules/m^{3})', 'FontSize', 14, 'Interpreter', 'tex');

title_str = sprintf('External Point-Source Delivery: r_{tx} = %.2f um, Observation: %s', ...
      r_tx*1e6, obs_names{obs_idx});
title(title_str, 'FontSize', 13, 'Interpreter', 'tex');

legend({'Analytical model'}, 'Location', 'best', 'FontSize', 11, 'Interpreter', 'tex');

% Uncomment to export a publication-ready figure:
% exportgraphics(gcf, 'external_transmitter_concentration_profile.png', 'Resolution', 300);