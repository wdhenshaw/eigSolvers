%
%  Solve a PDE by an eigenfunction expansion
%    IN-HOMOGENEOUS BOUNDARY CONDITIONS
%
%    u_tt - c^2 u_xx =0   , 0 < x < 1 , 0 < t <= tFinal
%    u(a,t)=ga(t), u(b,t)=gb(t), 
%    u(x,0) = u_0(x)
%
% USAGE:
%   eve -Nx=200 -Nk=100 -alpha=4.5 -Nt=50
% 
%   eve -alpha=10.5 -Nk=200 -Nt=100
%

function eve(varargin)

  % clear; clf;  
  % Set defaults for plotting 
  fontSize=16; lineWidth=2; markerSize=6; 
  set(0,'DefaultLineMarkerSize',markerSize);
  set(0,'DefaultLineLineWidth',lineWidth);
  set(0,'DefaultAxesFontSize',fontSize);
  set(0,'DefaultLegendFontSize',fontSize);
  % size of matlab figure : default = [560,420] 
  xwidth = 560;
  ywidth = 540; % 570; % 420; % 570;  

  tf=3; 
  Nk= 100;  % number of Fourier modes
  Nx= 500;  % number of space intervals
  Nt= 20;   % number of time-steps 
  a=0; b=pi; 
  c=1; 
  epsu=1e-7; % cutoff for error on the log plot
  savePlots=0; 
  figDir = '../doc/fig/';
  plotName='eveSineBC';

  
  sol = 0; %  g(t)=sin(alpha*t) 

  amp=.5; 
  % sol = 1; amp=.5;  % g(t) = 1 - cos(alpha*t)

  sol = 2;  % Heaviside function, "impulsive start"

  % BC at x=0: g(t)=sin(alpha*t)
  % true u(x,t) = sin(alpha*(t-x/c))
  alpha=2.5; 


  % --- read command line args ---
  for i = 1 : nargin
    line = varargin{i};
    sol       = getInt( line,'-sol',sol );
    Nx        = getInt( line,'-Nx',Nx );
    Nk        = getInt( line,'-Nk',Nk );
    Nt        = getInt( line,'-Nt',Nt );
    savePlots = getInt( line,'-savePlots',savePlots );
    tf        = getReal( line,'-tf',tf );
    alpha     = getReal( line,'-alpha',alpha );
    plotName  = getString( line,'-plotName',plotName );
  end


  if sol==0 
    g = @(t) sin(alpha*t);
    uHat = @(k,t) c^2*(2/pi)*k/(c^2*k^2-alpha^2) * ( sin(alpha*t) - (alpha/k)*sin(k*t));    
  elseif( sol==1 )
    g = @(t) amp*(1 - cos(alpha*t));

    uHat = @(k,t) amp*(2/pi)*c^2*k*(  1./(c^2*k^2) - 1/(c^2*k^2-alpha^2) * cos(alpha*t) ...
                                   - (1./(c^2*k^2) - 1/(c^2*k^2-alpha^2))* cos(k*t) );  
  else
    g = @(t) 1; % implusive start
    % uHat = @(k,t) c^2*(2/pi)*k*( (1/(c*k)^2)*( 1-cos(k*t)) );
    uHat = @(k,t) (2/pi)*(1/k)*( 1-cos(k*t) );
  end

  uTrue = @(x,t) g(t-x/c); 

  dx=(b-a)/Nx;  % grid spacing 
  ia=1; ib=ia+Nx; 
  I=ia:ib;
  x = zeros(Nx+1,1);
  for i=ia:ib
    x(i)=(i-ia)*dx;
  end

  u = zeros(Nx+1,1);
  v = zeros(Nx+1,1);
  ue = zeros(Nx+1,1);
 

  dt = tf/Nt;
  for n=1:Nt
    t=n*dt; 

    for i=ia:ib
      if( x(i)<= c*t )
       % ue(i) = sin(alpha*(t-x(i)/c));
       ue(i) = uTrue(x(i),t);
      else
        ue(i)=0;
      end
    end

    u(:)=0.;
    v(:)=0.; 
    hsign = [0,-1,0,1];
    for k=1:Nk

      % u1 = c^2*(2/pi)*k/(c^2*k^2-alpha^2);
      % uHat = u1*sin(alpha*t) - (alpha/k)*u1*sin(k*t);
      % u(I) = u(I) + uHat*sin(k*x(I));
      u(I) = u(I) + uHat(k,t)*sin(k*x(I)); 

      % v(I) = v(I) + ( uHat - ((2/pi)/k)*g(t) )*sin(k*x(I));
      % k4 = mod(k-1,4)+1;
      % v(I) = v(I) + ( ( (-1)^(k+1) + hsign(k4)  )* ((2/pi)/k)*g(t) )*sin(k*x(I)); % H(x-.5*pi)*g(t)

      % series for v = 1 - x/pi
      vHat =  2/(pi*k);
      v(I) = v(I) + vHat*sin(k*x(I));      
    end

    % subplot(2,1,1);
    figure(1)
    plot(x,u,'r','LineWidth',2); hold on; 
    plot(x,ue,'k-','LineWidth',1.5);
    % plot(x,v,'b-','LineWidth',1.5);
     hold off;
    % legend('u','ue','v');
    % plot(x,ue,'r-', x,u,'b-');
    legend('u','ue');

    title(sprintf('u: t=%6.3f, sol=%d, alpha=%g, Nk=%d',t,sol,alpha,Nk));
    
    grid on; xlabel('x');
    % ylim([-1.25,1.25]);
    if savePlots 
      savePlotFile(sprintf("%s%s",figDir,plotName),'pdf'); 
    end  

    figure(2)
    plot(x, u-ue,'b'); hold on;
    % estimated correction: 
    w = g(t)*( v(I) - (1-x(I)/pi) ); 
    plot(x, w,'k-','Linewidth',1); hold off;
    title(sprintf('u-ue: t=%6.3f, alpha=%g, Nk=%d',t,alpha,Nk));
    legend('u-ue','w'); 
    grid on; xlabel('x');

    if savePlots 
      savePlotFile(sprintf("%s%sErr",figDir,plotName),'pdf'); 
    end  

    figure(3)
    uc = zeros(Nx+1,1); % corrected
    uc(I) = u(I) - w(I); 
    plot(x,uc,'r','LineWidth',2); hold on; 
    plot(x,ue,'k-','LineWidth',1.5);  hold off;
    legend('uc','ue'); grid on; xlabel('x');

    title(sprintf('corrected u: t=%6.3f, sol=%d, alpha=%g, Nk=%d',t,sol,alpha,Nk));

    if savePlots 
      savePlotFile(sprintf("%s%sCorrected",figDir,plotName),'pdf');
    end      


    figure(4)
    plot(x, uc-ue,'b'); hold on;
    title(sprintf('uc-ue: t=%6.3f, alpha=%g, Nk=%d',t,alpha,Nk)); hold off;
    legend('uc-ue'); 
    grid on; xlabel('x');

    if savePlots 
      savePlotFile(sprintf("%s%sCorrectedErr",figDir,plotName),'pdf'); 
    end  


    % subplot(2,1,2);
    if( 1==0 )
      figure(5)
      % plot(x, u-ue,'b');
      semilogy(x, max(abs(u-ue),epsu),'b'); 
      title(sprintf('log error in u: t=%6.3f, alpha=%g, Nk=%d',t,alpha,Nk));
      legend('err');  grid on; xlabel('x');
      if savePlots 
        savePlotFile(sprintf("%s%sLogErr",figDir,plotName),'pdf'); 
      end  
    end

    if n~=Nt 
      pause
    end
  end

  

end



