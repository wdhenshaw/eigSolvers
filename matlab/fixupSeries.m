%
%  Fixup a poorly converging Fourier Series
% USAGE:
% 
%   fixupSeries -Nx=200 -Nk=100 -Nt=50
% 
%   fixupSeries -Nx=1000 -Nk=64 -savePlots=1

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
  Nk= 100;  % number of Fourier modes (EVEN)
  Nx= 5000;  % number of space intervals
  Nt= 20;   % number of time-steps 
  a=-pi; b=pi; 
  c=1; 
  epsu=1e-7; % cutoff for error on the log plot
  savePlots=0; 
  plotName='topHatSeries'; 

  % Hat function interval (scaled by pi below)
  x0 =-1;
  x1 = 0; 

  % BC at x=0: g(t)=sin(alpha*t)
  % true u(x,t) = sin(alpha*(t-x/c))
  % alpha=2.5; 


  % --- read command line args ---
  for i = 1 : nargin
    line = varargin{i};
    Nx        = getInt( line,'-Nx',Nx );
    Nk        = getInt( line,'-Nk',Nk );
    Nt        = getInt( line,'-Nt',Nt );
    savePlots = getInt( line,'-savePlots',savePlots );
    tf        = getReal( line,'-tf',tf );
    x0        = getReal( line,'-x0',x0 );
    x1        = getReal( line,'-x1',x1 );
    plotName  = getString( line,'-plotName',plotName );
  end

  x0 = x0*pi;
  x1 = x1*pi;   

  % Nk = Nk + mod(Nk,2); % make even

  g = @(t) sin(alpha*t);

  dx=(b-a)/Nx;  % grid spacing 
  ia=1; ib=ia+Nx; 
  I=ia:ib;
  x = zeros(Nx+1,1);
  for i=ia:ib
    x(i)= a + (i-ia)*dx;
  end

  u = zeros(Nx+1,1);
  v = zeros(Nx+1,1);
  ue = zeros(Nx+1,1);

  if 1==0 
    for i=ia:ib
      if( x(i)>x0 && x(i)<x1 )
        ue(i) = 1;
      else
        ue(i) = 0;
      end
    end

    % hsign = [0,-1,0,1];
    % There are a total of Nk modes:
    if mod(Nk,2)==0 
      % Nk even:
      kFirst = -Nk/2;
      kLast  =  Nk/2-1;
    else
      % Nk odd:
      kFirst = -(Nk-1)/2;
      kLast  =  (Nk-1)/2;
    end 
    u(:)=0; 
    for k=kFirst:kLast

      if k ~=0 
        uHat = -(1/(2*pi*1i*k))*( exp(-1i*k*x1) - exp(-1i*k*x0) );
      else
        uHat = (x1-x0)/(2*pi); 
      end

      u(I) = u(I) + uHat*exp(1i*k*x(I));

    end

    figure(5)
    plot(x,real(u),'r-',x,imag(u),'b-','LineWidth',2); hold on; 
    plot(x,ue,'k-','LineWidth',1.5);
      % plot(x,v,'b-','LineWidth',1.5);
    hold off;
    legend('Re(u)','Im(u)','ue');
    title(sprintf('u: N_k=%d',Nk));
      
    grid on;
    xlabel('x'); 

    if savePlots 
      savePlotFile(sprintf("%s",plotName),'pdf'); 
    end 

  end

  if( 1==1 )
    % ----- sine series of v= 1 - x/pi,  on [0,pi] ----
    a=0; b=pi;
    dx=(b-a)/Nx;  % grid spacing 
    ia=1; ib=ia+Nx; 
    I=ia:ib;
    x = zeros(Nx+1,1);
    for i=ia:ib
      x(i)= a + (i-ia)*dx;
    end

    v(:)=0; 
    for k=1:Nk
      vHat =  2/(pi*k);
      v(I) = v(I) + vHat*sin(k*x(I));

    end

    figure(5)
    % subplot(2,1,1)
    plot(x,v,'r-'); 
    legend('S_N(v)'); grid on; xlabel('x'); 
    title(sprintf('sine series for v=1 - x/\\pi: N_k=%d',Nk));

    if savePlots 
      savePlotFile('../doc/fig/fourierSeriesOneMinusXoverPi','pdf'); 
    end  

    figure(6)
    % subplot(2,1,2)
    plot(x,v - (1-x/pi),'b-'); 
    legend('S_N - v'); grid on; xlabel('x'); 
    title(sprintf('Difference S_N(1-x/\\pi) - (1-x/\\pi): N_k=%d',Nk));
      
   if savePlots 
      savePlotFile('../doc/fig/errorInFourierSeriesOneMinusXoverPi','pdf'); 
    end      
    

  end

  if( 1==1 )
    % ----- sine series of v=1 on [0,2 pi] ----
    a=0; b=2*pi;
    dx=(b-a)/Nx;  % grid spacing 
    ia=1; ib=ia+Nx; 
    I=ia:ib;
    x = zeros(Nx+1,1);
    for i=ia:ib
      x(i)= a + (i-ia)*dx;
    end

    v(:)=0; 
    for k=1:Nk

      vHat =  1/(2*pi) *(4/k)*(1-cos(k*pi));
      v(I) = v(I) + vHat*sin((k/2)*x(I));

    end
    figure(7)
    plot(x,v,'r-','LineWidth',2); hold on;
    plot( [pi,pi],[0,1],'k-' ); hold off
    % plot(x,ue,'k-','LineWidth',1.5);
      % plot(x,v,'b-','LineWidth',1.5);
    hold off;
    legend('S_N(u)');
    title(sprintf('sine series for u=1: N_k=%d',Nk));
      
    grid on;
    xlabel('x'); 

  end
  % example from "Gibbs phenomenon removal by adding Heaviside functions"

  

  a=-pi; b=pi;
  dx=(b-a)/Nx;  % grid spacing 
  ia=1; ib=ia+Nx; 
  I=ia:ib;
  x = zeros(Nx+1,1);
  for i=ia:ib
    x(i)= a + (i-ia)*dx;
  end


  % -----------------------------------
  function [y] = f1(x)

    n = length(x);
    y = zeros(1,n);
    for i=1:n
      if( x(i) < -pi/3 )
        y(i)= -x(i)*(x(i)/pi^2 + 2/pi);
      elseif( x(i)<2*pi/3 )
        y(i)= x(i)/(3*pi); 
      else
        y(i) = 2 + 2*sin(x(i)-pi/2); 
      end
    end
  end
  % ----------------------------------

  % -----------------------------------
  % Heavi-side function
  function [y] = Heaviside(x,x0)

    n = length(x);
    for i=1:n
      if( x(i) < x0 )
        y(i)= 0;
      else
        y(i) = 1; 
      end
    end
  end
  % ----------------------------------


 relTol=1.e-10; 
 absTol=1.e-10;


  w = f1(x)';
  f1s = zeros(Nx+1,1); % Fourier series approx. to f1 
  h0s = zeros(Nx+1,1); % Fourier series approx to H(x-0)
  h1s = zeros(Nx+1,1); % Fourier series approx to H(x-(-pi/3))
  h2s = zeros(Nx+1,1); % Fourier series approx to H(x-(pi/3))

  

  for k=-Nk:Nk % NOTE range 

    % Compute S_N(f1)
    % integrate : (1/2 pi) int f1(x) exp(- i * k * x ) dx
    f1exp = @(x) f1(x).*exp(-1i*k*x);
    f1Hat =  (1/(2*pi))*integral( f1exp,-pi,pi);

    f1s(I) = f1s(I) + f1Hat*exp(1i*k*x(I));

    % compute S_N(H_0)(x), 
    expk = @(x) exp(-1i*k*x); 

    h0Hat =  (1/(2*pi))*integral( expk,0,pi,'RelTol',relTol,'AbsTol',absTol);
    % h0 = @(x) Heaviside(x,0).*exp(-1i*k*x);
    % h0Hat =  (1/(2*pi))*integral( h0,-pi,pi,'RelTol',relTol,'AbsTol',absTol);
    h0s(I) = h0s(I) + h0Hat*exp(1i*k*x(I));

    h1Hat =  (1/(2*pi))*integral( expk,-pi/3,pi,'RelTol',relTol,'AbsTol',absTol);
    % h1 = @(x) Heaviside(x,-pi/3).*exp(-1i*k*x);
    % h1Hat =  (1/(2*pi))*integral( h1,-pi,pi,'RelTol',relTol,'AbsTol',absTol);
    h1s(I) = h1s(I) + h1Hat*exp(1i*k*x(I));

    h2Hat =  (1/(2*pi))*integral( expk,2*pi/3,pi,'RelTol',relTol,'AbsTol',absTol);
    % h2 = @(x) Heaviside(x,2*pi/3).*exp(-1i*k*x);
    % h2Hat =  (1/(2*pi))*integral( h2,-pi,pi,'RelTol',relTol,'AbsTol',absTol);
    h2s(I) = h2s(I) + h2Hat*exp(1i*k*x(I));

  end

  h02s = zeros(Nx+1,1); % Fourier series approx to S_{2N} H(x/2)
  for k=-2*Nk:2*Nk % NOTE range 

    expk = @(x) exp(-1i*k*x); 
    h0Hat =  (1/(2*pi))*integral( expk,0,pi);
    % h0 = @(x) Heaviside(x,0).*exp(-1i*k*x);
    % h0Hat =  (1/(2*pi))*integral( h0,-pi,pi);

    %  S_{2N} H(x/2)
    h02s(I) = h02s(I) + h0Hat*exp(1i*k*x(I)/2);  % Note evaluate at x/2 

  end  

 

  figure(1)
  plot( x,w,'b-', x,real(f1s),'r-', x,imag(f1s),'g-'); hold on;

  title(sprintf('f1 and series Nk=%d',Nk));
  legend('f1','Re S_N(f)','Im S_N(f)','Location','best'); xlabel('x');
  grid on;

  figure(2)
  plot( x,real(h0s),'b-', x,real(h02s),'r-', x,real(h0s -h02s),'g-' );
  title(sprintf('Example Nk=%d',Nk));
  legend('S_N H(x,0)','S_{2N} H(x/2)',...
         'S_N H(x,0) - S_{2N} H(x/2)',...
         'Location','best'); xlabel('x');
  grid on;

  f1star = zeros(Nx+1,1); % corrected approximation

  heavi1 = zeros(Nx+1,1);
  heavi1(I) = Heaviside(x(I),-pi/3);

  heavi2 = zeros(Nx+1,1);
  heavi2(I) = Heaviside(x(I),2*pi/3);

  f1star(I) = f1s(I) -(8/9)*( h0s(I) - h02s(I) ) ...
            +  (2/3)*( h1s(I) - heavi1(I) ) ...
            - (25/9)*( h2s(I) - heavi2(I) );
  figure(3)
  plot( x,w,'b-', x,real(f1star),'r-' );
  title(sprintf('f1 and corrected FS, Nk=%d',Nk));
  legend('f1','S_N^* (f1)',...
         'Location','best'); xlabel('x');
  grid on;  

  figure(4)
  % size(w)
  % size(f1star)
  % pause
  subplot(2,1,1)
  plot( x, real(w-f1s),'r-');
  title(sprintf('Error in original series, Nk=%d',Nk));
  legend('f1 - S_N (f1)',...
         'Location','best'); xlabel('x');
  grid on;

  subplot(2,1,2);
  plot( x, real(w-f1star),'b-');
  title(sprintf('Error in corrected series, Nk=%d',Nk));
  legend('f1 - S_N^* (f1)',...
         'Location','best'); xlabel('x');
  grid on;




end



