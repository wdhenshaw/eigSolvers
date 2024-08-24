%
%  Compute eigenvalues for PDE BVP's
%
%    u_tt - c^2 u_xx =0   , 0 < x < 1 , 0 < t <= tFinal
%    u(a,t)=ga(t), u(b,t)=gb(t), 
%    u(x,0) = u_0(x)
%
% USAGE:
%     pdeEigs -numRes= -plotOption=
%

function pdeEigs(varargin)

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


numRes=5;
N0=20;     % start with 20 points 
order=2; % order of accuracy 
implicit=1; % 1=implicit scheme
plotOption=0; % 2= movie mode
nplot=10;      % thus many times per time interval

c=1.;            % wave speed
x0=.5; betaPulse=20; % parameters in pulse IC
a=0.; b=1.;      % space interval interval
tFinal=1.;       % final time 
cfl=0.9;         % cfl number 
 
alpha   =.5;     % parameter in implicit scheme
betaBar =.0;     % parameter in implicit scheme

upwind = 0.;      % set to 1 for upwind dissipation

periodic=-1; dirichlet=1;
bc = dirichlet;

  % --- read command line args ---
  for i = 1 : nargin
    line = varargin{i};
    order       =    getInt( line,'-order',order );
    implicit    =    getInt( line,'-implicit',implicit );
    numRes      =    getInt( line,'-numRes',numRes );
    plotOption  =    getInt( line,'-plotOption',plotOption );
    nplot       =    getInt( line,'-nplot',nplot);
    N0          =    getInt( line,'-N0',N0 );

    cfl         =   getReal( line,'-cfl',cfl );
    tFinal      =   getReal( line,'-tf',tFinal );

    upwind      =   getReal( line,'-upwind',upwind );

    alpha       =   getReal( line,'-alpha',alpha );
    betaBar     =   getReal( line,'-betaBar',betaBar );

    betaPulse   =   getReal( line,'-betaPulse',betaPulse );

    bc          = getString( line,'-bc',bc );
    if( strcmp(bc,'d') )
      bc=dirichlet;
    elseif( strcmp(bc,'p'))
      bc=periodic;
    end 
  end

fprintf('--------------------- pdeEigs order=%d ---------------------\n',order);
fprintf(' [a,b]=[%g,%g] c=%g cfl=%g, upwind=%d implicit=%d bc=%s\n',a,b,c,cfl,upwind,implicit,bc);
if( upwind>0 ) upwChar='s'; else upwChar=''; end 
%   fprintf('Wave equation, CFL=%4.2f (%5.3f), tf=%5.3f, order=%d implicit=%d upwind=%g.\n',cfl,tFinal,c*dt/dx,order,implicit,upwind);

%% alpha=.5;  % for implicit scheme 
beta = 1 - 2*alpha; 

% betaBar = alpha - 2*alphaBar - 1/12
% alphaBar = ( alpha - betaBar - 1/12 )/2 
% betaBar=0.; 
alphaBar= ( alpha - betaBar -1./12. )/2. ; % for betaBar=0
if( implicit==1 )
  fprintf('IMPLICIT: alpha=%g, beta=%g, alphaBar=%g, betaBar=%g\n',alpha,beta,alphaBar,betaBar);
  if( betaBar > 2*alphaBar )
    fprintf('WARNING: betaBar > 2*alphaBar -- scheme is not unconditionally stable\n');
  end
end


% parameters for pulse solution:
par.a=a;
par.b=b;
par.c=c;
par.betaPulse=betaPulse;
par.x0=x0; 

% uexact = @(x,t) exp( -( betaPulse*(x-x0-c*t) ).^2 ); % exact solution function
% periodic pulse: 
uexact = @(x,t) pulse(x,t,par);

utexact = @(x,t) (2*betaPulse^2*c)*(x-x0-c*t).*uexact(x,t); 
u0xx = @(x) ( 4.*betaPulse^4*(x-x0).^2 - 2*betaPulse^2).*exp(-betaPulse^2*(x-x0).^2); 
u0 = @(x) uexact(x,0); % initial condition function for u 
u1 = @(x) utexact(x,0); % initial condition function for u_t 
ga = @(t) uexact(a,t); % BC forcing on left
gb = @(t) uexact(b,t); % BC forcing on left



% Loop over grid resolutions: 
for m=1:numRes

  Nx= N0*2^m;        % number of space intervals
  dx=(b-a)/Nx;  % grid spacing 
  
  %       G--G--B--X--X--X--X-- ...    --X--B--G--G
  %             ia                          ia
  ng=1 + order/2;    % number of ghost 
  ia=ng+1;           % index of boundary point at x=a
  ib=ia+Nx;          % index of boundary point at x=b

  i1=ia+1;   % first index to apply equation at.
  i2=ib-1;   % last index to apply equation
  if( bc==periodic ) 
    i1=ia;   % apply EQn at boundary for periodic
  end
  iperiod = ib-ia; % period in index space u(i)=u(i+iperiod)

  I = i1:i2; % apply equation at these points
  Ng=ib+ng;  % number of grid points

  x = linspace(a-ng*dx,b+ng*dx,Ng)';   % spatial grid 

  % allocate space for the solution at two levels
  unm1 = zeros(Ng,1);     % holds U_i^{n-1}
  un   = zeros(Ng,1);     % holds U_i^n
  unp1 = zeros(Ng,1);     % holds U_i^{n+1}
  
  dt = cfl*dx/c;         % time step (adjusted below)
  Nt = round(tFinal/dt); % number of time-steps 
  % fprintf('dt=%9.3e, adjusted=%9.3e, diff=%8.2e\n',dt,tFinal/Nt,abs(dt-tFinal/Nt));
  dt = tFinal/Nt;        % adjust dt to reach tFinal exactly
  
  iplot = round(Nt/nplot); 


  par.Ng        = Ng; 
  par.ng        = ng; 
  par.bc        = bc; 
  par.dirichlet = dirichlet;  
  par.periodic  = periodic;   
  par.ia        = ia;         
  par.ib        = ib;         
  par.i1        = i1;         
  par.i2        = i2;         
  par.iperiod   = iperiod;    
  par.ga        = ga;         
  par.gb        = gb;         
  par.uexact    = uexact;     

  if( order==2 )
     adSosup = c*dt*1./8.;   % CHANGED SIGN FROM cgmx
  elseif( order==4 )
     adSosup = c*dt*5./288.;
  else
     exit(666);
  end 
  uDotFactor=.5; % for D0t 
  adxUpw = upwind*uDotFactor*(adSosup/dx); 

  if( implicit==1 )  
  
    % --- Form the implicit matrix -----
    %    I - alpha*dt^2 *L 
    
    A = sparse(Ng,Ng);
    rhs = zeros(Ng,1);
    
    cdtdx=c*dt/dx;    % CFL parameter
    cdtdxSq=cdtdx^2; 
    cdtdx4=cdtdxSq^2; 
    % Assign interior equations
    if( order==2 )
      for( i=i1:i2 )
        A(i,i-1) = (   cdtdxSq);
        A(i,i  ) = (-2*cdtdxSq); 
        A(i,i+1) = (   cdtdxSq);
        % A(i,i-2) =                           - adxUpw*( -1 );
        % A(i,i-1) =    - alpha*(   cdtdxSq)   - adxUpw*(  4 );
        % A(i,i  ) = 1. - alpha*(-2*cdtdxSq)   - adxUpw*( -6 ); 
        % A(i,i+1) =    - alpha*(   cdtdxSq)   - adxUpw*(  4 );
        % A(i,i+2) =                           - adxUpw*( -1 );        
      end;

    elseif( order==4 )
      for( i=i1:i2 )
        A(i,i-3) =                                                         - adxUpw*(  1 );
        A(i,i-2) =    - alpha*(    -cdtdxSq/12.) + alphaBar*(    cdtdx4 )  - adxUpw*( -6 );
        A(i,i-1) =    - alpha*(  16*cdtdxSq/12.) + alphaBar*( -4*cdtdx4 )  - adxUpw*( 15 );
        A(i,i  ) = 1. - alpha*( -30*cdtdxSq/12.) + alphaBar*(  6*cdtdx4 )  - adxUpw*(-20 ); 
        A(i,i+1) =    - alpha*(  16*cdtdxSq/12.) + alphaBar*( -4*cdtdx4 )  - adxUpw*( 15 );
        A(i,i+2) =    - alpha*(    -cdtdxSq/12.) + alphaBar*(    cdtdx4 )  - adxUpw*( -6 ); 
        A(i,i+3) =                                                         - adxUpw*(  1 );
      end; 
    else
      pause; pause; 
    end
    
    % Boundary conditions - fake Dirichlet for now: 
    if( bc==dirichlet )
      A(ia-2,ia-2)=1;
      A(ia-1,ia-1)=1;
      A(ia  ,ia  )=1;
      A(ib  ,ib  )=1;
      A(ib+1,ib+1)=1; 
      A(ib+2,ib+2)=1; 
    elseif( bc==periodic )

      for( j=1:i1-1 )
        A(j,j)=1; A(j,j+iperiod)=-1;  
      end
      for( j=i2+1:Ng )
        A(j,j)=1; A(j,j-iperiod)=-1;  
      end      

    else
      fprintf('ERROR: bc?\n');
      exit();
    end 


  end ; % end implicit


  numEigs=10;
  [V,D] = eigs(A,numEigs,'smallestabs');
  lam = diag(D);
  m = length(lam);
  for i=1:m
    fprintf('lam(%d) = %12.4e\n',i,lam(i));
  end

  for i=1:m
    plot( x,V(:,i),'-');
    title(sprintf('Eig %d: lam=%12.4e',i,lam(i)));
    xlabel('x'); grid on;
    pause
  end  


  pause; pause; pause;


  t=0.; 
  unm1 = u0(x);    % initial condition at t=0
  t=dt; 
  utt0 = c^2*u0xx(x); 
  %   un = unm1 + dt*u1(x) +.5*dt^2*utt0;  % 3rd order approximation to u at t=dt 
  un = uexact(x,t); % first step

  % un   = uexact(x,t); % u at t=dt  *FIX ME*
  % --- Start time-stepping loop ---
  for( n=2:Nt )
    t = n*dt;        % new time

    % 2nd-order accurate scheme
    if( order==2 )

      if( implicit==0 )

        % --- explicit scheme ---

        unp1(I) = 2.*un(I) - unm1(I) + (c*dt/dx)^2*( un(I+1)-2*un(I)+un(I-1) ); 

        % Testing: corrector -- but only stable to sqrt(2)/2 !!
        numCorrections=0; 
        if( numCorrections>0 )
          % --- corrector ---
          unp1 = applyBoundaryConditions( unp1,t, par );

          unp1(I) = 2.*un(I) - unm1(I) + (.5*(c*dt/dx)^2)*( unp1(I+1)-2*unp1(I)+unp1(I-1) + unm1(I+1)-2*unm1(I)+unm1(I-1) ); 
          % unp1(I) = 2.*un(I) - unm1(I) + ((c*dt/dx)^2)*( .25*(unp1(I+1)-2*unp1(I)+unp1(I-1)) + ...
          %                                                 .5*(un(I+1)  -2*un(I)  +un(I-1)) + ...
          %                                                .25*(unm1(I+1)-2*unm1(I)+unm1(I-1)) ); 
          % unp1(I) = 2.*un(I) - unm1(I) + ( (c*dt/dx)^2 )*( unp1(I+1)-2*unp1(I)+unp1(I-1) ); 

        end

      else

        % --- implicit scheme ----
        %  Dpt Dmt u^n = alpha*L u^{n+1} + alpha*L u^{n-1} 

        rhs(I)=  2.*un(I) - unm1(I) +  beta *(c*dt/dx)^2*( un(I+1)-2*un(I)+un(I-1) ) ...
                                    +  alpha*(c*dt/dx)^2*( unm1(I+1)-2*unm1(I)+unm1(I-1) ) ...
                                    - adxUpw*( -unm1(I+2) +4*unm1(I+1) -6*unm1(I) +4*unm1(I-1) -unm1(I-2) );
        if( bc==dirichlet )
          for j=[ia-2,ia-1,ia,ib,ib+1,ib+2] 
            rhs(j)=uexact(x(j),t);
          end 
        elseif( bc==periodic )
          J = [1:i1-1,i2+1:Ng]; 
          rhs(J)=0; 
        else
          fprintf('ERROR: bc?\n'); exit();
        end 

        unp1 = A\rhs;

      end

    elseif( order==4 )
      % Dpt Dmt u^n = c^2 L_4h u^n + (dt^2/12) L_2h^2 u^n 
      if( implicit==0 )
        unp1(I) = 2.*un(I) - unm1(I) + ((c*dt/dx)^2/12.)*( -un(I+2) +16*un(I+1) -30*un(I)+ 16*un(I-1) -un(I-2) ) ...
                                   + (((c*dt/dx)^4)/12.)*(  un(I+2)  -4*un(I+1)  +6*un(I)  -4*un(I-1) +un(I-2) );
      else
        rhs(I) = 2.*un(I) - unm1(I) + beta *((c*dt/dx)^2/12. )*( -un(I+2) +16*un(I+1) -30*un(I)+ 16*un(I-1) -un(I-2) ) ...
                                 - betaBar *(((c*dt/dx)^4)   )*(  un(I+2)  -4*un(I+1)  +6*un(I)  -4*un(I-1) +un(I-2) ) ...
                                   + alpha *((c*dt/dx)^2/12. )*( -unm1(I+2) +16*unm1(I+1) -30*unm1(I)+ 16*unm1(I-1) -unm1(I-2) ) ...
                                 - alphaBar*(((c*dt/dx)^4)   )*(  unm1(I+2)  -4*unm1(I+1)  +6*unm1(I)  -4*unm1(I-1) +unm1(I-2) ) ...
                                 - adxUpw*( unm1(I+3) -6*unm1(I+2) +15*unm1(I+1) -20*unm1(I) +15*unm1(I-1) -6*unm1(I-2) + unm1(I-3) );

       % rhs(I) = 2.*un(I) - unm1(I) + alpha*((c*dt/dx)^2/12.   )*( -unm1(I+2) +16*unm1(I+1) -30*unm1(I)+ 16*unm1(I-1) -unm1(I-2) ) ...
       %                           - alphaBar*(((c*dt/dx)^4)     )*(  unm1(I+2)  -4*unm1(I+1)  +6*unm1(I)  -4*unm1(I-1) +unm1(I-2) );

        if( bc==dirichlet )
          for j=[ia-2,ia-1,ia,ib,ib+1,ib+2] 
            rhs(j)=uexact(x(j),t);
          end 
        elseif( bc==periodic )
          J = [1:i1-1,i2+1:Ng]; 
          rhs(J)=0; 
        else
          fprintf('ERROR: bc?\n'); exit();
        end         

        unp1 = A\rhs;
        
      end
    else
        pause; pause;
    end;

    if( implicit==0 )
      unp1 = applyBoundaryConditions( unp1,t, par );

      % % --- apply boundary conditions ---
      % if( bc==dirichlet )
      %   unp1(ia)=ga(t); % BC at x=a
      %   unp1(ib)=gb(t); % BC at x=b
      %   if( ng==1 )
      %     % assign ghost points -- do this for now
      %     unp1(ia-ng:ia-1)=uexact(x(ia-ng:ia-1),t);
      %     unp1(ib+1:ib+ng)=uexact(x(ib+1:ib+ng),t);
      %   end;
      % elseif( bc==periodic )
      %   for( j=1:i1-1 )
      %     unp1(j) = unp1(j+iperiod);  
      %   end
      %   for( j=i2+1:Ng )
      %     unp1(j) = unp1(j-iperiod);
      %   end          
      % else
      %   fprintf('ERROR: bc?\n'); exit();
      % end 
    end;

  
    unm1=un;   % set unm1 <- un for next step
    un=unp1;   % Set un <- unp1 for next step

    if( plotOption==2 && mod(n,iplot)==0 )
      ue = uexact(x,t);  % eval exact solution
      subplot(2,1,1);
      plot( x(ia:ib),un(ia:ib),'b-o', x(ia:ib),ue(ia:ib),'k-','Linewidth',2);
      title(sprintf('FD%d%d%s t=%9.2e Nx=%d cfl=%.2f imp=%d',order,order,upwChar,t,Nx,cfl,implicit));
      legend('computed','true');  grid on; 
      ylim([-.1,1.1]); 

      subplot(2,1,2);
      plot( x(ia:ib),un(ia:ib)-ue(ia:ib),'r-','Linewidth',2);
      title(sprintf('FD%d%d%s t=%9.2e Nx=%d cfl=%.2f imp=%d',order,order,upwChar,t,Nx,cfl,implicit));
      legend('error');  grid on; 


      drawnow;     
    end    

  end;
  % --- End time-stepping loop ---
  
  ue = uexact(x,t);  % eval exact solution:
  errMax(m) = max(abs(un(ia:ib)-ue(ia:ib)));  % max-norm error
  fprintf(' t=%10.4e: m=%d Nx=%3d Nt=%4d dt=%9.3e maxErr=%8.2e',t,m,Nx,Nt,dt,errMax(m));
  if( m==1 ) fprintf('\n'); else fprintf(' ratio=%6.2f\n',errMax(m-1)/errMax(m)); end; 

  % plot results 
  if( 1==1 )
    ue = uexact(x,t);  % eval exact solution
    subplot(2,1,1);
    plot( x(ia:ib),un(ia:ib),'b-o', x(ia:ib),ue(ia:ib),'k-','Linewidth',2);
    title(sprintf('FD%d%d%s t=%9.2e Nx=%d cfl=%.2f imp=%d',order,order,upwChar,t,Nx,cfl,implicit));
    legend('computed','true');  grid on; 
    ylim([-.1,1.1]); 
   subplot(2,1,2);
    plot( x(ia:ib),un(ia:ib)-ue(ia:ib),'r-','Linewidth',2);
    title(sprintf('FD%d%d%s t=%9.2e Nx=%d cfl=%.2f imp=%d',order,order,upwChar,t,Nx,cfl,implicit));
    legend('error');  grid on; 
  else

    plot( x(ia:ib),un(ia:ib),'b-o', x(ia:ib),ue(ia:ib),'k-','Linewidth',2);
    legend('computed','true');

    title(sprintf('FD%d%d%s t=%9.3e Nx=%d dt=%8.2e cfl=%5.3f imp=%d',order,order,upwChar,t,Nx,dt,cfl,implicit));
    xlabel('x');  ylabel('u');
    grid on;
  end 

  % set(gca,'FontSize',14);

  % print('-depsc2',sprintf('ps9_WE_Nx_%d.eps',Nx)); % save as an eps file
  pause

end; % end for m 

fprintf('done\n'); 

end

% --- Utility functions ---

% Function getReal: read a command line argument for a real variable
function [ val ] = getReal( line,name,val)
 % fprintf('getReal: val=%g line=[%s] name=[%s]\n',val,line,name);
 if( strncmp(line,strcat(name,'='),length(name)+1) )
   val = sscanf(line,sprintf('%s=%%e',name)); 
   % fprintf('getReal: scan for val=%g\n',val);
 end
end

% Function getInt: read a command line argument for an integer variable
function [ val ] = getInt( line,name,val)
 if( strncmp(line,strcat(name,'='),length(name)+1) )
   val = sscanf(line,sprintf('%s=%%d',name)); 
 end
end

% Function getString: read a command line argument for a string variable
function [ val ] = getString( line,name,val)
 if( strncmp(line,strcat(name,'='),length(name)+1) )
   val = sscanf(line,sprintf('%s=%%s',name)); 
 end
end
