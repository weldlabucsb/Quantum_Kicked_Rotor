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
Nkicks = 500;  %number of kicks
V0 = 10;
tau = .3E-6;
T=7E-6;

kbar = 4*hbar*kL^2*T/m;
K = 2*V0*Er*tau*kL^2*T/m;

%qkr parameters
% kbar = 1.2;  %effective Planck constant
% K = 10;  %Stochastic parameter (kick strength)
g = 10;

sk = .1; %momentum width in units of 2kL
%%
V = @(x) K*cos(x)/kbar;
Vt = @(t,x) V(x).*(mod(t,1)==0);
%%
Nx = 2^16;
xmin = -15000; xmax = 15000; D = xmax-xmin; dx = D/Nx; dk = 2*pi/D;
tfinal = Nkicks; iter=250; dt = 1/iter; %final time and time spacing

%%
n = -(Nx-1)/2:Nx/2;
x = dx*n';
k = dk*n';
c = zeros(Nx,Nkicks+1);
V = K*cos(x)/(kbar*tau/T);
propK = exp(-1i*kbar*k.^2*dt/2);
c(:,1) = exp(- sk^2 * x.^2/2); c(:,1) = c(:,1)/sqrt(sum(abs(c(:,1)).^2*dx));
tic
for ii = 1:Nkicks
    cii = c(:,ii);
    for jj=1:iter
        if jj<=tau*iter/T
            cii = exp(-1i*(V+g*abs(cii).^2)*dt/2).*cii;
            cii = ifft(ifftshift(propK.*fftshift(fft(cii))));
            cii = exp(-1i*(V+g*abs(cii).^2)*dt/2).*cii;
        else
            cii = exp(-1i*g*abs(cii).^2*dt/2).*cii;
            cii = ifft(ifftshift(propK.*fftshift(fft(cii))));
            cii = exp(-1i*g*abs(cii).^2*dt/2).*cii;
        end
    end
    c(:,ii+1) = cii/sqrt(sum(cii'*cii*dx));
end
toc
chat = fftshift(fft(c),1)/sqrt(Nx*dk/dx); 
psquare = k'.^2*abs(chat).^2*dk;

%%
figure(37); clf;
plot(0:1:Nkicks,psquare,'r.','MarkerSize',15);
xlabel('Kicks'); ylabel('$\langle (p/2 \hbar k_L)^2 \rangle$','interpreter','latex');
ylim([0 max(psquare)*1.2]); xlim([0 Nkicks]);
%%
for ii = 1:5:Nkicks
    f = figure(29); clf; f.Position = [100 200 600 300];
    subplot(121)
    semilogy(x,abs(c(:,ii)).^2*dx); xlim([xmin xmax]); ylim([1E-6 1]);
    xlabel('$x \, (d/2 \pi)$','interpreter','latex'); ylabel('$|\psi(x)|^2 dx \,(2 k_L)$','interpreter','latex');
    subplot(122)
    semilogy(k,abs(chat(:,ii)).^2*dk); xlim([min(k) max(k)]); ylim([1E-6 1]);
    xlabel('$k \,(2 k_L)$','interpreter','latex'); ylabel('$|\phi(k)|^2 dk \,(d/2 \pi)$','interpreter','latex')
    
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