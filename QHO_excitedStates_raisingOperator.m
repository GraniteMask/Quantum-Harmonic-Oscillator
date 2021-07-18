N=200;  % no of grid points
x=linspace(-5,5,N); % x-coordinate
psi0=exp(-x.*x); % ground state wave function, Griffiths, p. 58, Eq. 2.59
%---forming the raising operator ------------------------------------%
neg_one=-1.*ones(N,1);
pos_one=1.*ones(N,1);
X=spdiags(x,N,N);
DX=0.5.*(spdiags([neg_one pos_one ],[-1 1] ,N,N));
aplus=2.^(-0.5).*(X-DX); % Griffiths, p. 54, Eq. 2.47
%--- Using the raising operator to find excited states------------------%
psi1=(aplus)*psi0'; % first excited state
psi1=psi1-mean(psi1); % normalizing
psi1=(psi1./max(psi1));
psi2=aplus*psi1; % second excited state
psi2=psi2-mean(psi2); 
psi2=psi2./max(psi2);
psi3=aplus*psi2; % third excited state
psi3=psi3-mean(psi3); 
psi3=psi3./max(psi3);
%------rescaling the wavefunction for visualization purpose----------%
psi0=psi0-mean(psi0);  
psi0=psi0./max(psi0);
psi1=psi1+3;
psi2=psi2+7;
psi3=psi3+12;
V=0.5*x.*x; % parabolic oscillator potential
V=V./max(V);
V=100.*V;
%----------------------------------------------------------------------%
figure(1)
plot(x,psi0,'b','linewidth',2)
   text(0, 1.5,'E_0','fontSize',14)
axis([-2 2 0 15])
hold on
plot(x,psi1,'b','linewidth',2)
   line([-2 2], [3 3],'color','k','Linestyle','-.')
   text(-0.5, 3.5,'E_1','fontSize',14)
plot(x,psi2,'b','linewidth',2)
   line([-2 2], [7 7],'color','k','Linestyle','-.')
   text(0.5, 7.5,'E_2','fontSize',14)
plot(x,psi3,'b','linewidth',2)
   line([-2 2], [12 12],'color','k','Linestyle','-.')
   text(0, 12.5,'E_3','fontSize',14)
plot(x,V,'r','linewidth',2)
   h=gca; 
   get(h,'FontSize') 
   set(h,'FontSize',14)
   xlabel('X','fontSize',14);
   ylabel('\psi','fontSize',14);
   title('Stationary States of Harmonic Oscillator','fontsize',14)
   fh = figure(1);
   set(fh, 'color', 'white'); 
   set(gca,'ytick',[]);
hold off