clc
clear all
%% Loading Data

list = dir('interp/data*.dat');
directory_name='interp/';
M=880; N=400; %Grid points in x and y directions

% Grid data
data = importdata(sprintf('%s%s',directory_name,list(1).name));
xdata = data.data(:,1);
ydata = data.data(:,2);
xgrid = reshape(xdata,M,N); 
ygrid = reshape(ydata,M,N);
xgrid=xgrid';ygrid=ygrid';

X=zeros(2*M*N,length(list));
%% Stacking vectors in time (X) ([u;v])
for k = 1:length(list)
    k
    data = importdata(sprintf('%s%s',directory_name,list(k).name));
    X(1:M*N,k) = data.data(:,3); %uvel
    X(M*N+1:end,k) = data.data(:,4); %uvel
end

Xmean = mean(X,2);
figure;
subplot(1,2,1)
contourf(xgrid,ygrid,reshape(Xmean(1:M*N),M,N)',100,'LineStyle','none')
caxis([-0.1,0.1]);
colormap(bluewhitered)
axis equal;
set(gca,'TickLabelInterpreter','latex','FontSize',16)
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
title('$\bar{u}(x,y)$','Interpreter','latex','FontSize',22)

subplot(1,2,2)
contourf(xgrid,ygrid,reshape(Xmean(M*N+1:end),M,N)',100,'LineStyle','none')
caxis([-0.1,0.1]);
colormap(bluewhitered)
axis equal;
set(gca,'TickLabelInterpreter','latex','FontSize',16)
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
title('$\bar{v}(x,y)$','Interpreter','latex','FontSize',22)

%% Proper Orthogonal Decomposition
Xfluc = X-Xmean; %mean subtraction for POD
[U, S, V] = svd(Xfluc, 'econ'); %computes SVD (same as solving eigenvalue proble for X*X)
Phi_POD = U; %POD Modes
singular_values = diag(S); %singular values
Energy_captured = (singular_values.^2) / sum(singular_values.^2);

% Plot energy captured from each mode
figure;
plot((Energy_captured),'-o','MarkerSize',8,'LineWidth',1);
set(gca,'TickLabelInterpreter','latex','FontSize',16)
xlabel('Number of Modes','Interpreter','latex')
ylabel('Energy','Interpreter','latex')

% Find modes that capture 99.9% energy
Num_modes = find(cumsum(Energy_captured) >= 0.999,1);

% Plotting circle
r = 0.5; % Radius
theta = linspace(0, 2*pi, 100);
xcirc = r*cos(theta);
ycirc = r*sin(theta);
% Plotting 6 POD modes

figure(100);hold on;
for plot_modes = 1:6
    subplot(3,2,plot_modes)
    contourf(xgrid,ygrid,reshape(Phi_POD(1:M*N,plot_modes),M,N)',100,'LineStyle','none')
    hold on
    fill(xcirc,ycirc,[0.5,0.5,0.5])
    caxis([-0.01,0.01]);
    colormap(bluewhitered)
    axis equal;
    set(gca,'TickLabelInterpreter','latex','FontSize',16)
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
end

%% Dynamic Mode Decomposition
X1 = X(:,1:end-1);
X2 = X(:,2:end);

% SVD of X1
[U, S, V] = svd(X1, 'econ');

% Atilde
Atilde = U'*X2*V/S;

%Eigen value decompostion
[W, D] = eig(Atilde); % W: Eigenvectors, D: Eigenvalues
Phi_DMD = X2*V/S*W; % DMD modes

% Frequency and growthrate
lambda = diag(D);
dt = 2; % Time step
sigma = log(lambda) /(2*pi*dt); % Continuous-time frequencies

idx = [5,14,3];
figure;
subplot(1,2,1)
scatter(real(lambda),imag(lambda),15)
hold on
plot(real(lambda(idx)),imag(lambda(idx)),'r.','MarkerSize',20)
set(gca,'TickLabelInterpreter','latex','FontSize',16)
xlabel('$Re(\lambda)$','Interpreter','latex')
ylabel('$Im(\lambda)$','Interpreter','latex')
axis equal
xlim([-1.01,1.01])
ylim([-1.01,1.01])

subplot(1,2,2)
scatter(imag(sigma),real(sigma),15)
hold on
plot(imag(sigma(idx)),real(sigma(idx)),'r.','MarkerSize',20)
ylim()
set(gca,'TickLabelInterpreter','latex','FontSize',16)
xlabel(' Frequency','Interpreter','latex')
ylabel('Growth Rate','Interpreter','latex')
ylim([-0.004,0.001])
%%
%DMD Modes
figure(1000);hold on
for k = 1:length(idx)
    subplot(3,2,2*k-1)
    contourf(xgrid,ygrid,reshape(real(Phi_DMD(1:M*N,idx(k))),M,N)',100,'LineStyle','none')
    hold on
    fill(xcirc,ycirc,[0.5,0.5,0.5])
    caxis([-0.005,0.005]);
    colormap(bluewhitered)
    axis equal;
    set(gca,'TickLabelInterpreter','latex','FontSize',16)
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')

    subplot(3,2,2*k)
    contourf(xgrid,ygrid,reshape(imag(Phi_DMD(1:M*N,idx(k))),M,N)',100,'LineStyle','none')
    hold on
    fill(xcirc,ycirc,[0.5,0.5,0.5])
    caxis([-0.005,0.005]);
    colormap(bluewhitered)
    axis equal;
    set(gca,'TickLabelInterpreter','latex','FontSize',16)
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
end


