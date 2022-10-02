 clf;
%% Random Stars
a = -31;
b = 31;
% a = -2;
% b = 2;
% a = -1;
% b = 1;
n = 100;
rnx = a + (b-a).*rand(n,1);
rny = a + (b-a).*rand(n,1);
rnz = a + (b-a).*rand(n,1);

starwall=linspace(b,b,n);
starwall2=linspace(a,a,n);


%% Planet Periods 
Tme=(87.97/365.26); % Mercury period 
Tv=(224.70/365.26); % Venus period 
Te=(365.26/365.26); % Earth period 
Tma=(686.98/365.26); % Mars Period
Tj=(4332.82/365.26); % Jupiter period 
Ts=(10755.70/365.26); % Saturn period 
Tu=(30687.15/365.26); % Uranus period 
Tn=(60190.03/365.26); % Neptune period 

Tc=(75.27); % Comet period 

PH_e=(4/365.26);
AH_e=(184/365.26);

tic
n=1;
%%
npoints=(500*Tn)/n; % Number of points
dt= (0.002*n); % Time step in Years

%% Initialise position of planet in Astronomical Units (AU)

x_me = 0.39; 
y_me = 0;
z_me = 0;


x_v = 0.72; 
y_v = 0;
z_v = 0;


x_e = 1; 
y_e = 0;
z_e = 0;


x_ma = 1.524; 
y_ma = 0;
z_ma = 0;


x_j = 5.20; 
y_j = 0;
z_j = 0;


x_s = 9.54; 
y_s = 0;
z_s = 0;


x_u = 19.19; 
y_u = 0;
z_u = 0;


x_n = 30.06; 
y_n = 0;
z_n = 0;

x_c =  17.80; 
y_c = 0;
z_c = 0;

%% Calculations

% Te=(sqrt(x_e^3)); % Number of orbits/ Period
% Tma=(sqrt(x_ma^3)); % Number of orbits/ Period
% 
% Kc_e=(Te^2)/(x_e^3); % constant across planets, Earth = 0.9997,
% Kc_ma=(Tma^2)/(x_ma^3); % constant across planets, Mars = 1.005

Axme=[x_me];
Ayme=[y_me];

Axv=[x_v];
Ayv=[y_v];

Axe=[x_e];
Aye=[y_e];

Axma=[x_ma];
Ayma=[y_ma];

Axj=[x_j];
Ayj=[y_j];

Axs=[x_s];
Ays=[y_s];

Axu=[x_u];
Ayu=[y_u];

Axn=[x_n];
Ayn=[y_n];

Axc=[x_c];
Ayc=[y_c];





%% Initialise Velocity of planet in AU/year 

v_x_me = 0;
v_y_me = (2*pi*x_me)/Tme;
v_z_me = 0;

v_x_v = 0;
v_y_v = (2*pi*x_v)/Tv;
v_z_v = 0;

v_x_e = 0;
v_y_e = (2*pi*x_e)/Te;
v_z_e = 0;

v_x_ma = 0;
v_y_ma = (2*pi*x_ma)/Tma;
v_z_ma = 0;

v_x_j = 0;
v_y_j = (2*pi*x_j)/Tj;
v_z_j = 0;

v_x_s = 0;
v_y_s = (2*pi*x_s)/Ts;
v_z_s = 0;

v_x_u = 0;
v_y_u = (2*pi*x_u)/Tu;
v_z_u = 0;

v_x_n = 0;
v_y_n = (2*pi*x_n)/Tn;
v_z_n = 0;

v_x_c = 0;
v_y_c = (pi*x_c)/Tc;
v_z_c = 0;

%% Film Figure 
% writerObj=VideoWriter('FinalSolarSystemV012','MPEG-4'); 
% open(writerObj);
fig1=figure(1);
fig1.Position = [10 10 1500 1000];

%% Loop over the timesteps
for step = 1:npoints-1 % for steps 1 to 499


radius_me = sqrt(x_me^2+y_me^2); % Radius: Distance from the sun
radius_v = sqrt(x_v^2+y_v^2); % Radius: Distance from the sun
radius_e = sqrt(x_e^2+y_e^2); % Radius: Distance from the sun
radius_ma = sqrt(x_ma^2+y_ma^2); % Radius: Distance from the sun
radius_j = sqrt(x_j^2+y_j^2); % Radius: Distance from the sun
radius_s = sqrt(x_s^2+y_s^2); % Radius: Distance from the sun
radius_u = sqrt(x_u^2+y_u^2); % Radius: Distance from the sun
radius_n = sqrt(x_n^2+y_n^2); % Radius: Distance from the sun

radius_c = sqrt(x_c^2+y_c^2); % Radius: Distance from the sun

%% Compute the new Velocities in the X and Y directions 
v_x_new_me = v_x_me - (4*pi^2*x_me*dt)/(radius_me^3);   
v_y_new_me = v_y_me - (4*pi^2*y_me*dt)/(radius_me^3);

v_x_new_v = v_x_v - (4*pi^2*x_v*dt)/(radius_v^3);   
v_y_new_v = v_y_v - (4*pi^2*y_v*dt)/(radius_v^3);

v_x_new_e = v_x_e - (4*pi^2*x_e*dt)/(radius_e^3);     
v_y_new_e = v_y_e - (4*pi^2*y_e*dt)/(radius_e^3);

v_x_new_ma = v_x_ma - (4*pi^2*x_ma*dt)/(radius_ma^3);   
v_y_new_ma = v_y_ma - (4*pi^2*y_ma*dt)/(radius_ma^3);

v_x_new_j = v_x_j - (4*pi^2*x_j*dt)/(radius_j^3);   
v_y_new_j = v_y_j - (4*pi^2*y_j*dt)/(radius_j^3);

v_x_new_s = v_x_s - (4*pi^2*x_s*dt)/(radius_s^3);   
v_y_new_s = v_y_s - (4*pi^2*y_s*dt)/(radius_s^3);

v_x_new_u = v_x_u - (4*pi^2*x_u*dt)/(radius_u^3);   
v_y_new_u = v_y_u - (4*pi^2*y_u*dt)/(radius_u^3);

v_x_new_n = v_x_n - (4*pi^2*x_n*dt)/(radius_n^3);   
v_y_new_n = v_y_n - (4*pi^2*y_n*dt)/(radius_n^3);

v_x_new_c = v_x_c - (4*pi^2*x_c*dt)/(radius_c^3);   
v_y_new_c = v_y_c - (4*pi^2*y_c*dt)/(radius_c^3);

%% Euler Cromer Step - Update positions using newly caluclated Velocities
x_new_me = x_me+v_x_new_me*dt;
y_new_me = y_me+v_y_new_me*dt;

x_new_v = x_v+v_x_new_v*dt;
y_new_v = y_v+v_y_new_v*dt;

x_new_e = x_e+v_x_new_e*dt;
y_new_e = y_e+v_y_new_e*dt;

x_new_ma = x_ma+v_x_new_ma*dt;
y_new_ma = y_ma+v_y_new_ma*dt;

x_new_j = x_j+v_x_new_j*dt;
y_new_j = y_j+v_y_new_j*dt;

x_new_s = x_s+v_x_new_s*dt;
y_new_s = y_s+v_y_new_s*dt;

x_new_u = x_u+v_x_new_u*dt;
y_new_u = y_u+v_y_new_u*dt;

x_new_n = x_n+v_x_new_n*dt;
y_new_n = y_n+v_y_new_n*dt;

x_new_c = x_c+v_x_new_c*dt;
y_new_c = y_c+v_y_new_c*dt;



% Create an array for x and y values
Axme(end+1)=x_new_me;
Ayme(end+1)=y_new_me;

Axv(end+1)=x_new_v;
Ayv(end+1)=y_new_v;

Axe(end+1)=x_new_e;
Aye(end+1)=y_new_e;

Axma(end+1)=x_new_ma;
Ayma(end+1)=y_new_ma;

Axj(end+1)=x_new_j;
Ayj(end+1)=y_new_j;

Axs(end+1)=x_new_s;
Ays(end+1)=y_new_s;

Axu(end+1)=x_new_u;
Ayu(end+1)=y_new_u;

Axn(end+1)=x_new_n;
Ayn(end+1)=y_new_n;

Axc(end+1)=x_new_c;
Ayc(end+1)=y_new_c;


%% Update x and y with new positions 

x_me = x_new_me;
y_me = y_new_me;

x_v = x_new_v;
y_v = y_new_v;

x_e = x_new_e;
y_e = y_new_e;

x_ma = x_new_ma;
y_ma = y_new_ma;

x_j = x_new_j;
y_j = y_new_j;

x_s = x_new_s;
y_s = y_new_s;

x_u = x_new_u;
y_u = y_new_u;

x_n = x_new_n;
y_n = y_new_n;

x_c = x_new_c;
y_c = y_new_c;

   	
%% Plot Planet position immediately, rotating around the sun
% set(gcf, 'color', [1 1 1]) % Figure Color
% 
% plot(0,0,"oy","MarkerSize",30,"MarkerFaceColor","yellow");
% Sun=plot(0,0,"oy","MarkerSize",15,"MarkerFaceColor","yellow");
% 
% hold on
% 
% %% Plot Settings 
% 
% 
% plot3(starwall,rny,rnz,".w"); 
% plot3(rnx,starwall,rnz,".w"); 
% plot3(rnx,rny,starwall2,".w");
% 
% ax=gca;
% ax.XColor = '[0.7 0.7 0.7]';
% ax.YColor = '[0.7 0.7 0.7]';
% ax.ZColor = '[0.7 0.7 0.7]';
%  
% set(gca, 'color', [0 0 0]) % Background Color
% 
% xlabel("x(AU)");
% ylabel("y(AU)");
% zlabel("z(AU)");
% axis([a b a b a b]); % sets the axis size for the plot
% 
% 
% %% Inner Planets
% plot(x_new_me,y_new_me,"ok","MarkerSize",2,"MarkerFaceColor","[0.9 0.9 0.9]");
% Mercury=plot(x_new_me,y_new_me,"ok","MarkerSize",2,"MarkerFaceColor","[0.9 0.9 0.9]");
% 
% plot(x_new_v,y_new_v,"ok","MarkerSize",5,"MarkerFaceColor","[0.9290 0.6940 0.1250]");
% Venus=plot(x_new_v,y_new_v,"ok","MarkerSize",5,"MarkerFaceColor","[0.9290 0.6940 0.1250]");
% 
% plot(x_new_e,y_new_e,"ok","MarkerSize",6,"MarkerFaceColor","[0 0 1]");
% Earth=plot(x_new_e,y_new_e,"ok","MarkerSize",6,"MarkerFaceColor","[0 0 1]");
% 
% plot(x_new_ma,y_new_ma,"ok","MarkerSize",3,"MarkerFaceColor","[1 0 0]");
% Mars=plot(x_new_ma,y_new_ma,"ok","MarkerSize",3,"MarkerFaceColor","[1 0 0]");

% %% Inner Traces
% plot(Axme,Ayme,".","MarkerEdgeColor","[0.9 0.9 0.9]");
% plot(Axv,Ayv,".","MarkerEdgeColor","[0.9290 0.6940 0.1250]");
% plot(Axe,Aye,".","MarkerEdgeColor","[0 0 1]");
% plot(Axma,Ayma,".","MarkerEdgeColor","[1 0 0]");
% 
% 
% %% Outer Planets
% plot(x_new_j,y_new_j,"ok","MarkerSize",14,"MarkerFaceColor","[0.9290 0.840 0.1250]");
% Jupiter=plot(x_new_j,y_new_j,"ok","MarkerSize",14,"MarkerFaceColor","[0.9290 0.840 0.1250]");
% 
% plot(x_new_s,y_new_s,"ok","MarkerSize",13,"MarkerFaceColor","[0.8 0.8 0]");
% Saturn=plot(x_new_s,y_new_s,"ok","MarkerSize",13,"MarkerFaceColor","[0.8 0.8 0]");
% 
% plot(x_new_u,y_new_u,"ok","MarkerSize",10,"MarkerFaceColor","[0.6 0.6 1]");
% Uranus=plot(x_new_u,y_new_u,"ok","MarkerSize",10,"MarkerFaceColor","[0.6 0.6 1]");
% 
% plot(x_new_n,y_new_n,"ok","MarkerSize",11,"MarkerFaceColor","[0.2 0.2 1]"); 
% Neptune=plot(x_new_n,y_new_n,"ok","MarkerSize",11,"MarkerFaceColor","[0.2 0.2 1]"); 
% 
% plot(x_new_c,y_new_c,"<w","MarkerSize",8);
% Comet=plot(x_new_c,y_new_c,"<w","MarkerSize",8);
% 
% 
% %% Outer Traces
% plot(Axj,Ayj,".","MarkerEdgeColor","[0.9290 0.840 0.1250]");
% plot(Axs,Ays,".","MarkerEdgeColor","[0.8 0.8 0]");
% plot(Axu,Ayu,".","MarkerEdgeColor","[0.6 0.6 1]");
% plot(Axn,Ayn,".","MarkerEdgeColor","[0.2 0.2 1]");
% 
% plot(Axc,Ayc,".","MarkerEdgeColor","[0.9 0.9 0.9]");
% 
% 
% %% Additional Plot Settings
% lgd=legend([Sun,Mercury,Venus,Earth,Mars,Jupiter,Saturn,Uranus,Neptune,Comet],{'Sun','Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune','Comet'},'TextColor','white','Location','best');
% 
% title(lgd, 'Planets')
% view(-20,40);
% hold off
% 
% drawnow;

%% Update x and y velocities with new velocities 

v_x_me = v_x_new_me;
v_y_me = v_y_new_me;

v_x_v = v_x_new_v;
v_y_v = v_y_new_v;

v_x_e = v_x_new_e;
v_y_e = v_y_new_e;

v_x_ma = v_x_new_ma;
v_y_ma = v_y_new_ma;

v_x_j = v_x_new_j;
v_y_j = v_y_new_j;

v_x_s = v_x_new_s;
v_y_s = v_y_new_s;

v_x_u = v_x_new_u;
v_y_u = v_y_new_u;

v_x_n = v_x_new_n;
v_y_n = v_y_new_n;

v_x_c = v_x_new_c;
v_y_c = v_y_new_c;


% F=getframe(fig1);
% writeVideo(writerObj,F)
















end;
% close(writerObj);
% disp('Video File Written Successfully')



% 
set(gcf, 'color', [1 1 1]) % Figure Color

Sun=plot(0,0,"oy","MarkerSize",15,"MarkerFaceColor","yellow");
hold on

set(gca, 'color', [0 0 0]) % Background Color
xlabel("x(AU)");
ylabel("y(AU)");
zlabel("z(AU)");
axis([a b a b a b]); % sets the axis size for the plot



% plot(Axme,Ayme,"ok","MarkerSize",2,"MarkerEdgeColor","[0.9 0.9 0.9]");
% plot(Axv,Ayv,"ok","MarkerSize",5,"MarkerEdgeColor","[0.9290 0.6940 0.1250]");
% plot(Axe,Aye,"ow","MarkerSize",10,"MarkerFaceColor","[0 0 1]");
% plot(Axma,Ayma,"ok","MarkerSize",10,"MarkerFaceColor","[1 0 0]");
% plot(Axj,Ayj,"ok","MarkerSize",14,"MarkerEdgeColor","[0.9290 0.840 0.1250]");
% plot(Axs,Ays,"ok","MarkerSize",13,"MarkerEdgeColor","[0.8 0.8 0]");
% plot(Axu,Ayu,"ok","MarkerSize",10,"MarkerEdgeColor","[0.6 0.6 1]");
% plot(Axn,Ayn,"ok","MarkerSize",11,"MarkerEdgeColor","[0.2 0.2 1]"); 
% plot(Axc,Ayc,"<w","MarkerSize",8);


% 
% set(gcf, 'color', [1 1 1]) % Figure Color
% 
% plot(0,0,"oy","MarkerSize",15,"MarkerFaceColor","yellow");
% Sun=plot(0,0,"oy","MarkerSize",15,"MarkerFaceColor","yellow");
% 
% hold on

% Plot Settings 


plot3(starwall,rny,rnz,".w"); 
plot3(rnx,starwall,rnz,".w"); 
plot3(rnx,rny,starwall2,".w");

ax=gca;
ax.XColor = '[0.7 0.7 0.7]';
ax.YColor = '[0.7 0.7 0.7]';
ax.ZColor = '[0.7 0.7 0.7]';
 
% set(gca, 'color', [0 0 0]) % Background Color
% 
% xlabel("x(AU)");
% ylabel("y(AU)");
% zlabel("z(AU)");
% axis([a b a b a b]); % sets the axis size for the plot


%Inner Planets
 
% plot(x_new_me,y_new_me,"ok","MarkerSize",2,"MarkerFaceColor","[0.9 0.9 0.9]");
Mercury=plot(x_new_me,y_new_me,"ok","MarkerSize",2,"MarkerFaceColor","[0.9 0.9 0.9]");

% plot(x_new_v,y_new_v,"ok","MarkerSize",5,"MarkerFaceColor","[0.9290 0.6940 0.1250]");
Venus=plot(x_new_v,y_new_v,"ok","MarkerSize",5,"MarkerFaceColor","[0.9290 0.6940 0.1250]");

% plot(Axe,Aye,"ok","MarkerSize",6,"MarkerFaceColor","[0 0 1]");
Earth=plot(x_new_e,y_new_e,"ok","MarkerSize",6,"MarkerFaceColor","[0 0 1]");


% plot(x_new_ma,y_new_ma,"ok","MarkerSize",3,"MarkerFaceColor","[1 0 0]");
Mars=plot(x_new_ma,y_new_ma,"ok","MarkerSize",3,"MarkerFaceColor","[1 0 0]");

plot(Axme,Ayme,".","MarkerEdgeColor","[0.9 0.9 0.9]");
plot(Axv,Ayv,".","MarkerEdgeColor","[0.9290 0.6940 0.1250]");
plot(Axe,Aye,".","MarkerEdgeColor","[0 0 1]");
plot(Axma,Ayma,".","MarkerEdgeColor","[1 0 0]");

% plot(Axme,Ayme,".","MarkerEdgeColor","[1 1 1]");
% plot(Axv,Ayv,".","MarkerEdgeColor","[1 1 1]");
% plot(Axe,Aye,".","MarkerEdgeColor","[1 1 1]");
% plot(Axma,Ayma,".","MarkerEdgeColor","[1 1 1]");


% Outer Planets
% plot(x_new_j,y_new_j,"ok","MarkerSize",14,"MarkerFaceColor","[0.9290 0.840 0.1250]");
Jupiter=plot(x_new_j,y_new_j,"ok","MarkerSize",14,"MarkerFaceColor","[0.9290 0.840 0.1250]");

% plot(x_new_s,y_new_s,"ok","MarkerSize",13,"MarkerFaceColor","[0.8 0.8 0]");
Saturn=plot(x_new_s,y_new_s,"ok","MarkerSize",13,"MarkerFaceColor","[0.8 0.8 0]");

% plot(x_new_u,y_new_u,"ok","MarkerSize",10,"MarkerFaceColor","[0.6 0.6 1]");
Uranus=plot(x_new_u,y_new_u,"ok","MarkerSize",10,"MarkerFaceColor","[0.6 0.6 1]");

% plot(x_new_n,y_new_n,"ok","MarkerSize",11,"MarkerFaceColor","[0.2 0.2 1]"); 
Neptune=plot(x_new_n,y_new_n,"ok","MarkerSize",11,"MarkerFaceColor","[0.2 0.2 1]"); 

% plot(x_new_c,y_new_c,"<w","MarkerSize",8);
Comet=plot(x_new_c,y_new_c,"<w","MarkerSize",8);


plot(Axj,Ayj,".","MarkerEdgeColor","[0.9290 0.840 0.1250]");
plot(Axs,Ays,".","MarkerEdgeColor","[0.8 0.8 0]");
plot(Axu,Ayu,".","MarkerEdgeColor","[0.6 0.6 1]");
plot(Axn,Ayn,".","MarkerEdgeColor","[0.2 0.2 1]");

plot(Axc,Ayc,".","MarkerEdgeColor","[0.9 0.9 0.9]");

lgd=legend([Sun,Mercury,Venus,Earth,Mars,Jupiter,Saturn,Uranus,Neptune,Comet],{'Sun','Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune','Comet'},'TextColor','white','Location','best');

title(lgd, 'Planets')
view(-20,40);
hold off


% toc
% t=toc;

