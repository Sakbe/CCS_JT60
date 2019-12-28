 
P = Input_struct.p;
T = Input_struct.t;
X = linspace(1, 6, 100);
Y = linspace(-4, 4, 100);
U = x_np(1:length(P));

UXY=tri2grid(P,T,U,X,Y);

[XX, YY] = meshgrid(X, Y);

figure
pdemesh(Input_struct.p,Input_struct.e,[])
hold on
plot3(XX,YY,UXY,'.')


%% Compute field

[psir, psiz] = gradient(UXY, X(2)-X(1), Y(2)-Y(1));
Br = -psiz/2/pi ./ XX; % check signs
Bz = +psir/2/pi ./ XX;

figure
pdemesh(Input_struct.p,Input_struct.e,[])
hold on
quiver(XX,YY,1000*Br,1000*Bz)


%% 


figure(figSens)
iii = find(not(isnan(Input_struct.theta_sens)));
rr = Input_struct.r_sens(iii);
zz = Input_struct.z_sens(iii);
tt = Input_struct.theta_sens(iii);
hold on
quiver(rr,zz,.1*cos(tt),.1*sin(tt))