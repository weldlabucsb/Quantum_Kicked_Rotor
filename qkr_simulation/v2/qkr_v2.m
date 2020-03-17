close all; clear;
%This code simulates the delta kicked rotor for a Gaussian ensemble of
%initial momentum states. The evolution is diagonal over the free
%propagation period and non-diagonal for the instantaneous delta kick. Here
%you can just set K and kbar directly, or instead use
%the experimental parameter section to find K and kbar for the V,T and tau of
%interest.
%% Define physical constants
amu=1.66E-27;       % 1 AMU
m=7*amu;            % Lithium mass
lambda=1064E-9;     % Wavelength of light
h=6.626E-34;        % Planck's Const.
kL=2*pi/lambda;     % Wave Vector
hbar=h/(2*pi);        % Reduced planck's constant
Er=hbar^2*kL^2/(2*m); % Recoil Energy
wr = Er/hbar;       %Recoil Frequency
%colororder(lbmap(4,'RedBlue'));

%% Initialize parameters
kicks = 200;  %number of kicks
nk = 101;     %number of plane wave states/bands considered, need more if you expect more momentum
sk = .1;      %momentum width in units of 2kL
k0 = -3*sk:sk/30:3*sk; %k0 =0;  %array of initial momenta centered on 0

% Set QKR parameters directly
K = 11.6023; kbar = 1;

% Set experimental parameters and compute qkr parameters
% tau = .1E-6;  %time of pulse
% T = .8E-6;     %period on pulse generator, time between kicks is T-tau
% V = 1000;       %lattice depth in recoils
% kbar = 8*wr*T;    %qkr parameters
% talbot = pi/(2*wr);
% K = 2*kL^2*T*V*Er*tau/m;

D = K^2/4 * (1-2*besselj(2,K)+2*besselj(2,K)^2);
tL = D/kbar^2;

%% Compute energy
tic
psquare = zeros(length(k0),kicks+1); %momenta squared matrix
weight = exp(- k0.^2/(2*sk^2)); weight = weight/sum(weight);  %gaussian weights
for jj = 1:length(k0)  %can change to parfor if nk big, otherwise for will likely be faster
    n=(-(nk-1)/2:1:(nk-1)/2)+k0(jj);  %plane wave states
    c = zeros(nk,kicks+1); c((nk+1)/2,1)=1; %initialize coefficient vector
    Tmat = kbar^2*n.^2/2; %kinetic energy matrix for free evolution
    Kickmat = sparse(1:nk-1,2:nk,K/2,nk,nk)+sparse(2:nk,1:nk-1,K/2,nk,nk); %potential energy matrix for kick
    Ukick = expm(-1i*full(Kickmat)/kbar); %kick propagator
    Ufree = diag(exp(-1i*Tmat/kbar));  %free evolution propagator
    propagate = Ufree*Ukick;  %full period propagtor
    for ii = 1:kicks  %time evolve
        cii = propagate*c(:,ii);  %propagate
        c(:,ii+1) = cii/norm(cii);  %normalize
    end
    probs = conj(c).*c;   %compute probability vector at all times
    psquare(jj,:) = (n.^2)*probs;  %compute momentum squared at all times
end
toc
etot = weight*psquare;  %compute weighted kinetic energy
plot(0:1:kicks,etot/2,'.g','MarkerSize',15); hold on;  %plot kinetic energy
% plot([0 50],[0 50]*2*D/kbar^2,'-b','handlevisibility','off');
ylim([0 1.2*max(etot)]);
ylabel('$\langle (p/2\hbar k_L)^2 \rangle$','interpreter','latex'); xlabel('Kicks');

