close all; clear all;
%% Define physical constants in SI
amu=1.66E-27;       % 1 AMU
m=7*amu;            % Lithium mass
lambda=1064E-9;     % Wavelength of light
h=6.626E-34;        % Planck's Const.
r0 = 5.291E-11;     % Bohr radius
kL=2*pi/lambda;     % Wave Vector
d = lambda/2;       % Lattice Spacing
d2 = d/(2*pi);      % Distance units
hbar=h/(2*pi);        % Reduced planck's constant
Er=hbar^2*kL^2/(2*m); % Recoil Energy
Wr = Er/hbar;       % Recoil Frequency

%% Initialize parameters
Nkicks = 100;  %number of kicks

%qkr parameters
kbar = 12.27;  %effective Planck constant
K = 15.4;  %Stochastic parameter (kick strength)
g = 0;

sk = .5; %momentum width in units of 2kL
%%
V = @(x) K*cos(x)/kbar;
%%
Nx = 2^17;
xmin = -5000; xmax = 5000; D = xmax-xmin; dx = D/Nx; dk = 2*pi/D;
tfinal = Nkicks; iter =1; dt = 1/iter; %final time and time spacing

%%
n = -(Nx-1)/2:Nx/2;                              %number of points
x = dx*n';                                       %position mesh
k = dk*n';                                       %momentum mesh
c = zeros(Nx,Nkicks+1);                          %state vector
V = K*cos(x)/kbar; propV = exp(-1i*V);           %delta kick propagator
propK = exp(-1i*kbar*k.^2*dt/2);                 %momentum propagator
c(:,1) = exp(- sk^2 * x.^2/2);                   %initial state
c(:,1) = c(:,1)/sqrt(sum(abs(c(:,1)).^2*dx));    %normalize
tic
for ii = 1:Nkicks                                           %iterate over # of kicks
    cii = propV.*c(:,ii);                                   %kick
    for jj=1:iter                                           %iterate between kick
        cii = exp(-1i*g*abs(cii).^2*dt/2).*cii;             %mean field propagation
        cii = ifft(ifftshift(propK.*fftshift(fft(cii))));   %fft, diffusion, ifft
        cii = exp(-1i*g*abs(cii).^2*dt/2).*cii;             %mean field propagation
    end
    c(:,ii+1) = cii/sqrt(sum(cii'*cii*dx));                 %normalize
end
toc
chat = fftshift(fft(c),1)/sqrt(Nx*dk/dx);         %compute momentum state vector
psquare = k'.^2*abs(chat).^2*dk;                  %compute (p/2hbarkL)^2          

%%
figure(37); hold on;
plot(0:1:Nkicks,psquare,'r.','MarkerSize',15);
xlabel('Kicks'); ylabel('$\langle (p/2 \hbar k_L)^2 \rangle$','interpreter','latex');
ylim([0 max(psquare)*1.2]); xlim([0 Nkicks]);
%%
for ii = 1:2:Nkicks
    f = figure(29); clf; f.Position = [100 200 600 300];
    subplot(121)
    semilogy(x,abs(c(:,ii)).^2); xlim([xmin xmax]); ylim([1E-4 1]);
    xlabel('$x \, (d/2 \pi)$','interpreter','latex'); ylabel('$|\psi(x)|^2 \,(2 k_L)$','interpreter','latex');
    subplot(122)
    semilogy(k,abs(chat(:,ii)).^2); xlim([min(k) max(k)]); ylim([1E-4 1]);
    xlabel('$k \,(2 k_L)$','interpreter','latex'); ylabel('$|\phi(k)|^2 \,(d/2 \pi)$','interpreter','latex')
    
    drawnow;
    del = .01;
    filename = 'dum.gif';
    im = frame2im(getframe(gcf));
    [imind,cm] = rgb2ind(im,256);
    if ii == 1
        imwrite(imind,cm,filename,'gif','Loopcount',Inf,'DelayTime',del);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',del);
    end
end