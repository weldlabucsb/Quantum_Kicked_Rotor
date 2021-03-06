close all; clear;
%This code simulates the atom kicked rotor for a square lattice pulse. The
%final output is a plot of the kinetic energy term in the Hamiltonian
%stroboscopically measured after each kick. The evolution employs diagonal
%propagators over the pulse time and free evolution time, and switches
%bases every iteration from plane waves to Bloch waves and vice versa.
%Here, you must set the actual experimental parameters of V, T and tau.
%% Define physical constants
amu=1.66E-27;       % 1 AMU
m=23*amu;            % Lithium mass
lambda=589E-9;     % Wavelength of light
h=6.626E-34;        % Planck's Const.
kL=2*pi/lambda;     % Wave Vector
hbar=h/(2*pi);        % Reduced planck's constant
Er=hbar^2*kL^2/(2*m); % Recoil Energy
wr = Er/hbar; 
%colororder(lbmap(5,'RedBlue'));
%% Initialize parameters
kicks = 200;  %number of kicks
nk = 201;     %number of plane wave states/bands considered, need more if you expect more momentum
sk = .1;      %momentum width in units of 2kL
k0 = -3*sk:sk/30:3*sk; %k0 =0;  %array of initial momenta centered on 0

% Experimental parameters
tau = .1E-6;   %time of pulse (width for square pulse/approximately characteristic width for other shapes)
T =1.6E-6;     %periodicity of system, time between kicks is T-tau
V = 734;       %lattice depth in recoils

%compute qkr parameters
kbar = 8*wr*T;   %effective Planck constant
talbot = pi/(2*wr);  %Talbot time
K = 2*kL^2*T*V*Er*tau/m;  %Stochastic parameter (kick strength)
D = K^2/(4) * (1-2*besselj(2,K)+2*besselj(2,K)^2);  %Diffusion constant p^2 = 2Dt
tL = D/kbar^2;  %localization time

%% Compute momentum squared
psquare = zeros(length(k0),kicks+1);  %momenta squared matrix over time
weight = exp(- k0.^2/(2*sk^2)); weight = weight/sum(weight);  %gaussian weight
tic
for jj = 1:length(k0)  %can change to parfor if nk big, otherwise for will likely be faster
    k = (-(nk-1)/2:1:(nk-1)/2)+k0(jj);   %plane wave states
    p = 2*hbar*kL*k';                    %unitful momenta
    cvec = zeros(nk,kicks+1);            %coefficient vector
    cvec((nk+1)/2,1) = 1;                %initial momentum state
    [E,btop] = computeBands(V,k,k0(jj)); ptob = btop'; %compute band energies, similarity transformation
    Klat = diag(exp(-1i*E*Er*tau/hbar));  %diagonal propagation in bloch wave basis
    Kfree = diag(exp(-1i*p.^2*(T-tau) / (2*m*hbar)));  %diagonal free propagation in plane wave basis
    propagate = Kfree*btop*Klat*ptob;  %full period propagation matrix
    for ii = 1:kicks  %time evolve
        cii = propagate*cvec(:,ii);  %propagate
        cvec(:,ii+1) = cii/norm(cii);  %normalize
    end
    probs = conj(cvec).*cvec;  %compute probability vector over all kicks
    psquare(jj,:) = (p.^2)'*probs/(2*hbar*kL)^2;  %compute momentum squared over all kicks
end
toc
Etot = weight*psquare; %weight momenta squared
plot(0:kicks,Etot/2,'.','MarkerSize',15); hold on;  %plot kinetic energy term
% plot([0 40],D/kbar^2*[0 40],'handlevisibility','off');  %plot classical diffusion line
% plot([0 25],D/kbar^2*[tL tL],'--k','handlevisibility','off'); %plot localization line
ylim([0 max(Etot/2)*1.2]);
ylabel('$\frac{\langle (p/2\hbar k_L)^2 \rangle}{2}$','interpreter','latex'); xlabel('Kicks');