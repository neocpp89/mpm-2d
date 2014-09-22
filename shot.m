#!/usr/bin/env octave

p.rho = 1200;
p.E = 1e6;
p.nu = 0.3;
p.g = -0.8;
p.Np = 60;
p.N = 61;
p.rhop = 3;
p.tol = 1e-12;
p.T = 2;

% shot angle
angle = 0 * pi / 180.0;
mag = 0;

% make an NxN computational grid
N = p.N;
Np = p.Np;
x = linspace(0,1,N);
y = linspace(0,1,N);
h = x(2)-x(1);
dt = 0.4*h*sqrt(p.rho/p.E);
T = p.T;
E = p.E;
nu = p.nu;
tol = p.tol;

Nn = N*N;

% set up material points
%bar_w = 1.0; Nbx = Np; %
bar_w = 0.4; Nbx = Np; %
% Nbx = ceil(p.rhop*bar_w/h);
bar_h = 1.0; Nby = (bar_h/bar_w)*Nbx;
bar_h = 0.4; Nby = (bar_h/bar_w)*Nbx;

bxtmp = linspace(0, bar_w, Nbx+1);
% bar_x = (0.5 - bar_w/2) + bxtmp(1:end-1) + bxtmp(2)/2;
bar_x = 0.5 + bxtmp(1:end-1) + bxtmp(2)/2;
bytmp = linspace(0, bar_h, Nby+1);
bar_y = 0.2 + bytmp(1:end-1) + bytmp(2)/2;

% initialize position
Np = Nbx*Nby;
[xp, yp] = meshgrid(bar_x, bar_y);

cmx = mean(xp(:));
cmy = mean(yp(:));

th = atan2(yp(:) - cmy, xp(:) - cmx);
mm = hypot(yp(:) - cmy, xp(:) - cmx);

xp(:) = mm .* cos(th(:) + angle) + cmx;
yp(:) = mm .* sin(th(:) + angle) + cmy;

% initialize volume
vp = (bar_w/Nbx)*(bar_h/Nby)*ones(size(xp));
% vp = (pi*r^2/Np)*ones(size(xp));

% initialize mass
mp = p.rho*vp;

% initialize velocity
xp_t = mag*cos(angle)*ones(size(xp));
yp_t = mag*sin(angle)*ones(size(yp));

% initialize stress
sxxp = 0*ones(size(xp));
syyp = 0*ones(size(yp));
sxyp = 0*ones(size(yp));

% initialize strain
exxp = zeros(size(xp));
exyp = zeros(size(yp));
eyyp = zeros(size(yp));

% initialize body forces
bxp = zeros(size(xp));
byp = p.g*ones(size(yp));

% plot material points
plot(xp, yp, '.k'); axis equal;
xlim([min(x), max(x)]);
ylim([min(y), max(y)]);

% dump starting config to text file
fp = fopen('generated_particles.txt', 'w');
fprintf(fp, '%d\n', numel(mp));
for i=1:numel(mp)
    fprintf(fp, '%g %g %g %g %g %g %g %g %g\n', mp(i), vp(i), ...
        xp(i), yp(i), xp_t(i), yp_t(i), ...
        sxxp(i), sxyp(i), syyp(i));
end
fclose(fp);

% grid file
fp = fopen('generated_grid.txt', 'w');
fprintf(fp, '%d\n', N);
fprintf(fp, '1\n');
fclose(fp);

% bc file
fp = fopen('bc.cfg', 'w');
fprintf(fp, 'dirichlet-bcs = {\n');
fprintf(fp, '};\n');
fprintf(fp, 'periodic-bcs = {\n');
fprintf(fp, '};\n');
