% Structural and Lagrangian growth model
%    of migrating and drifting phytoplankton
%    ** global application **
%
% kai wirtz (HZG) Mar 2021

function gp = glob_prods(offset,season,parsv,pix0)
% the function can be called from multiple processors (e.g., HPC)
%  with arguments
%  offset: defines the first index of the local domain - increments in a parallel simulation
%  season: 1-4 1:'DJF' 2:'MAM' 3:'JJA' 4:'SON'
%  parsv:  string with changed parameter values (e.g., 'bmax=3' for switching off migration)
%  pix0:   starting grid cell in case of continuing a large loop
close all;
%  -------  control parameters
%clear all;
scale=1;
colim0 = 0.5; % degree of co-limitation at NO3(z=0)>0.5
zNel   = 0.0; % relative change in chemocline depth  0: present-day climatology
if exist('offset') & exist('season')
  offset=str2num(offset); season=str2num(season);
else
  offset=0; season=1;
end
pvs='bmax=75';  % default setting including migration
                %  'bmax=2' would discard long-range migration
if exist('parsv')
    pvs=parsv;
end
%  -------  name of I/O file
matfile=['glob_' num2str(scale,'%01d') pvs num2str(offset,'%02d') '_' num2str(season) '.mat']
if exist('pix0') % restart the domain loop from a finite value
  load(matfile);
  pix0=str2num(pix0);
else
  pix0=1;
end

%latsst = -89.5:89.5;lonsst = 0.5:359.5;
tags={'K2';'S1';'HOT';'BATS';'GD_iow';}; %  short names of monitoring sites
ouc=0; ode=0;outc=0;out=-1;tagv=1;di=1;  % initialisation of control parameters
z = [0:0.5:215];zf = z;nzf=length(zf);   % vertical domain and resolution
dz=z(2:end)-z(1:end-1); dz=[dz dz(end)];
%  -------  load coordinates and depths of profile data
load('locs'); load('depthg')
depthg=double(depthg);

%for season=1:4
%‚dx = 1; dy = 1;
fprintf('offset=%d  season=%d\n',offset,season)
% ------- load all boundary setting from file
load(['glob_env_' num2str(season) '.mat'])
dx = 4*scale; dy = 5*scale;

xp_f=double(xp_f);
[nx ny]=size(no3_0_f);
if nx==721, nx=720; end; % use constant grid size
% -------- check grid consistency
maxy=ceil(ny/dy);
if (maxy-1)*dy + 1 + floor(offset/dx)*scale>ny,  maxy=maxy-1; end
maxx=ceil(nx/dx);
if (maxx-1)*dx + 1 + mod(offset,dx)*scale>nx,  maxx=maxx-1; end

% ---------------------------------
%         par settings
% ---------------------------------
param_set; % set all process parameters
omega=omega/1.15; % rescale mortality  %*1.325;%*1.32;

kchl0=kchl; I_al0=I_al; carb0    = 1;
eval([pvs ';']); % reset parameter if given in argument
if bmax<5 % variant without migration
  vsmS0=vs0;
end
vsmax=vsmS0;

% ------- load changed MLD from MPI-ESM, if future scenario is calculated
if zNel >0
  load zmld2
end

fprintf('saving to %s\n',['glob_' num2str(scale,'%01d') pvs num2str(offset,'%02d') '_' num2str(season) '.mat'])
logfile=['log_' pvs num2str(offset,'%02d') '_' num2str(season) '.log'];
fid=fopen(logfile,'w');fclose(fid);

%cvc0=0.1;  % critical cumulative NO3 variation
depthn=depthg;
di=1;li=1;pix=1;piy=1; % initialisation of counters
%size(no3_0_f)
fprintf('nx=%d ny=%d\t px<%d py<%d  x<%d y<%d\tmax %d %d\n',nx,ny,ceil(nx/dx),ceil(ny/dy),(ceil(nx/dx)-1)*dx+1+dx-1,(ceil(ny/dy)-1)*dy+1+dy-1,maxx,maxy);%size(dcm_f)

% ------------------------------------------
% loop over global/recatangular domain
% ------------------------------------------
for pix=pix0:maxx
 for piy=1:maxy
 %for li=1:5%size(loc,1)-2
 %  ix = i_loc(li,2);pix=ix;	   iy = i_loc(li,1) ;piy=iy;

% position on 2D grid, depending on CPU specific counters in parallel mode
  ix = (pix-1)*dx + 1 + mod(offset,dx)*scale;
  iy = (piy-1)*dy + 1 + floor(offset/dx)*scale;

% ------- retrieve boundary condition
% -------   (surface NO3, chemocline depth DCM)
  no3_0=no3_0_f(ix,iy);
  dcm=dcm_f(ix,iy);

  chlnew0= 0;chlnew=[]; vs_v(1)=0;

% ------- is valid marine site?
  if no3_0>-0.5 & dcm>0
   temp=temp_f(ix,iy); % SST
   Tfa = 1;
   Tf   = power(Q10f,(temp-20)/10); % temperature response
   zN = round(max(dcm+20,10)); % correction of DCM

% ------- shift of DCM, e.g. in RCP8.5 scenario
   delzN(pix,piy)=0;
   if zNel>0
     zmld=zmld2(iy,ix,season);
     if ~isnan(zmld)
        delzN(pix,piy)=  sign(zmld) * min(abs(zmld),zNel*zN);
      %  fprintf('zN: %1.1f + %1.1f (%1.1f,%1.1f)\n',zN,delzN(pix,piy),(zmld),zNel*zN)
        zN = zN + delzN(pix,piy);
      end
   end

   Nnew=0;
   par=par_f(ix,iy);    % incident PAR
   chl0=chl0_f(ix,iy);  % surface CHL
   ldif=min(ldif_f(ix,iy),dcm); % mixing depth

   % ---------------------------------
   %         base CHL profile
   % ---------------------------------
   % assume surface contribution of mobiles at low turbulence and co-limitation
   % TODO
   shalfac = 1-1/(1+exp((zN-zNcr)/2));
   if no3_0<2 & ldif<20
      chl0  = chl0*(shalfac*0.8+0.2);
   end

   % ------- uniform surface CHL down to chemcline,
   dcc = max(zN+Tf*23,10); %

   knc = 0.09;
   chlm0  = chl0*(1+exp(-knc*dcc))./(1+exp(knc*(zf-dcc)));
   chl_mod= chlm0;
   % ------- adjust water attenuation due to mixing properties
   dfac = exp(2-120/ldif);

   chl0mi=chlmi(ix,iy);  % annual average surface CHL
   if isnan(chl0mi), chl0mi=0.05; end
   limf = 1.;
   % ------- contribution of CDOM to light attenuation,
   %            estimated based on annual average surface CHL
   fkvi    = min(max(0.84+2.2*chl0mi,0.5),4);

   k0 =k00_0*fkvi*(1+mixatt*dfac); % light attenuation of water
   kchl=kchl0*(1+mixatt*dfac);   %  self-shading: contribution of CHL to light attenuation

   %fprintf('zN=%1.0f chl0=%1.2f\n',zN,chl0);%%1.1f %1.1f\t ,liTf=%1.1f,Tf,watdepth

% ------- exclude sites with very deep chemocline
   if(zN < 225) % heuristic maximum depth
     fprintf('%d %d %d %d\t no3=%1.1f dcm=%1.0f ldif=%1.0f\n',pix,piy,ix,iy,no3_0,dcm_f(ix,iy),ldif);%,liTf=%1.1f,Tf,watdepth

    % cloud=cloud_f(ix,iy);
     Kz  =  0.5*ldif*ldif; % turbulent mixing coefficient
     xp  = xp_f(ix,iy);    % relative seasonal variation of PAR
     zz  = z-zN;           % distance to chemocline center position
     chl00=chl0;

  % ------- temperature dependency
     Tfm = exp(0.12*(temp-20)) *(1+1/(1+exp(-0.5*xp)))*0.5;

  % ------- other growth parameters
     aff_I0 = aff_I00 ;
     aff_I  = aff_I0;
     Pm     = Pm0;
     A_N    = A_N00*Tf*limf;

  % ------- initialisation of auxiliary variables
     subs =0; nppf=0;  zC=310;
     qim=Qpmax; chl2ct0=0.05; chlnew=[];
     b=20; speed=0; carbn=0;
     nppf=zeros(1,nzf); % chl_rect=[];

  % ------- calc light profile
     parz = par*exp(-(k0*zf+kchl*cumsum(chlm0.*dz)));
  % ------- calc light in chemocline
     parN = par*exp(-(k0+kchl*chl00)*max(zN-50,0));

  % ----------------------------------------
  %      photosynthesis & growth model
  % ----------------------------------------
  %   immobile population (drifters)
  % ----------------------------------------
     I_al   = I_al0+3*parN; % photoacclimation
  % ------- N:C of immobile population
     qim    = 0.1*Q0+(Qpmax-0.1*Q0)*z/(zN+10);  % linear increase in quota
  % ------- CHL:C of immobile population
     chl2ct0 = chl2c(qim,parz,I_al,chl2c_min,chl2cf,temp);
   %%fprintf('zN=%1.3f %1.3f k=%1.5f %1.5f %d q=%1.3f parz=%1.1f chl2cf=%1.3f chl2ct0=%1.3f sc=%1.3f\n',zN,shalfac,k0,kchl,length(dz),mean(qim),parz(33),chl2cf,mean(chl2ct0),sum(chlm0.*dz));
  % ------- average C density of immobile population
     carb0 = sum(chlm0.*dz./chl2ct0)/zN;
  % initial guess of RGR (no quota)
     mu  = min(Pm*Tf,aff_I*parN);
    %%fprintf('Pm=%1.3f*%1.3f k=%1.5f %1.5f a*I=%1.3f*%1.3f\n',Pm,Tf,k0,kchl,aff_I,parN);
  %%   mu_pre = mu;

  % ----------------------------------------
  %   mobile population (vertical migrators)
  % ----------------------------------------
  % ------- initial guess of biomass
     mort = omega; % specific mortality rate
    	% TODO check for role/sensitivity of upper bound for numerical safety
     carbn = min(3E3,(mu/(Tfm*mort)-carb0)*2*15);% equilibrium biomass
  % initial guess of CHL:C (no quota)
     c2c0  = chl2c_min+chl2cf*0.4*Qpmax*chl2n(parN,I_al,temp);
  %%  fprintf(' mu=%1.2f\tpar=%1.1f/%1.1f\n',mu,parN,par);
  %%   fprintf('%d I=%1.0f/%1.1f\tC=%1.0f\tzN=%1.0f/%1.0f fZ=%1.2f fkvi=%1.2f NO3=%1.1f\n',di,par,parN,carbn,dcm,zN,shalfac,fkvi,no3_0);
  %%   fprintf('  mu=%1.2f T=%1.2f:%1.2f\tC0=%1.1f %1.1f\tchl0=%1.1f %1.1f\n',mu ,Tf,Tfm,carb0,carbn,chl0,chlm0(1));

  % ----------------------------------------
  %   loop for iterative calculation of
  %        traits and biomass of vertical migrators
  % ----------------------------------------
     iter=1; iv=[];
     % --------------------------------------------------------------------
     % iter<=1: assuming fast convergence, thus low sensitivity of
     %          migration traits on self-shading
     while (carbn > 1E-5 & iter<=1 &  parN >1)
        chln  = carbn*c2c0; % convert carbon to CHL
        % ------- find optimal migration traits
        find_strat

        % ------- viable strategy found?
        if max(vs_v ) <1E-5 | length(iv)==0
            carbn=0;
        else
          % -------  mortality estimate
         mort = omega; % specific mortality rate

   %%     fprintf('  \t   C1=%1.0f\tchln=%1.1f\n',carbn,chln);

         % ------- iterative calculation of mobile biomass
          itj=0;
          chln_0 = chln; chlnv=[];
          while (itj<4)
            % ------- updated self-shading
            %% par_N0/par_N=exp(kcd*min(max(0,zN-(zC-b)),2*b));
            patt  = par_N0/par_N*exp(-0.5*kchl*chln);
            % ------- updated biomass and CHL
            carbn = max((mu_m*patt/(Tfm*mort)-carb0)*2*b,0);
            chlnt = carbn*chl2cmv(ivm);
            chlnv = [chlnv chlnt]; % store in vector
            chln  = mean(chlnv);   % calculate average for better convergence
         %%fprintf(' %d\tmu=%1.3f\t%1.3f new=%1.3f mort=%1.3f C=%1.0f\tchln=%1.1f %1.1f\n',itj,mu_m,patt,10/new,1E3*mort,carbn,chlnt,mean(chlnv));
            itj=itj+1;
          end % while (itj<4)
          % ---------------------------------------------------
          %   resting period during migration
          %      sharp increase in vertical distribution function
          %      and max carbon estimate
          % ---------------------------------------------------
          maxCc  = maxC *exp(4*trest/(2*trest+2*b/speed));
          if max(carbn*xf)>maxCc,
            if out>=1,  text(0.8,190,[num2str(max(carbn*xf),'%1.0f') '->' num2str(maxCc,'%1.0f')],'FontSize',18,'FontWeight','bold'); end
            carbn = carbn*maxCc/max(carbn*xf); % cut unrealistic biomass peaks
          end
          if out>-10,
            fprintf('%d\t%1.1f\tzC=%1.0f b=%1.0f\tT=%1.0f/%1.0f v=%1.0f\t%1.3f %1.3f\t<strong>chl=%1.2f</strong>\t%1.3f\n',di,carbn,zC,b,mean(tr_v(iv)),2*b/mean(vs_v(iv)),mean(vs_v(iv)),sum(cf),max(cf),chln,mu_m);
          end
    %%fprintf('%d\t',di);
          chlnew = carbn*cf; % convert to CHL
          % ------- total chl profile = residual + relative new profile
          chl_mod = chlm0+chlnew;  %
          Nnew    = mean(qm_v(iv))*carbn; % N-conc of migrators
       %   end %if(length(iv)>0)
          subs = 1;
         end %if max(vs_v )if niche (new chl)
        iter=iter+1;
      end % while carbn
    % -------  renormalize base profile of immobiles
    %          depending of surface concentration of mobiles
    if length(chlnew)>0
      chlnew0= chlnew(1);
    end
  end % if zN

  % ------ renormalize if mobiles reach the surface
  fac  = max(chl0-chlnew0,0)/(chl0+0.0001);
  chlm0= chlm0*fac;
  if subs==0 % carbn<0
      chl_mod=chlm0;
    %%     retrieve_sim_chl
  end
  if length(chlnew)>0
     % ------- total chl profile = residual + relative new profile
     chl_mod=chlm0+chlnew;
  end
% ---------------------------------------------------
%  calc NPP of immobile/mixed phytoplankton
  parz = par*exp(-(k0*zf+kchl*cumsum(chl_mod.*dz))); % light profile
  Ilim = aff_Im*parz*Tfa; % light limitation factor

% recalc quota and CHL:C of immobile/mixed phytoplankton
  qim0    = 0.5*Q0+(Qpmax-0.5*Q0)*min(1,zf/zN);  % linear increase in quota
  if ldif <8 % layering under low turbulence
     parc=parz;
     qim1=qim0;
  else
     % -------- average light and nutrient condition in mixed layer
     z1=max(zf-ldif,0);
     z2=min(zf+ldif,200);
     iz1=round(ldif/(zf(2)-zf(1)));
     kk=k0+kchl*chl0;
     parc=par*(exp(-kk*z1)-exp(-kk*z2))./(kk*(z2-z1));

     % -------- step increase in quota by strong up-mixing
     nchemo=(zf+ldif)/zN;
     nchemo=(nchemo>1).*(zf+ldif-zN)./(zf+ldif);
     qim1 =  0.5*(1-nchemo).*(qim0(max(1,(1:nzf)-iz1))+qim0(min(nzf,(1:nzf)+iz1)));
     qim1 = qim1 + nchemo*Qpmax;
     %qim0 = 0.5*(qim0(max(1,(1:nzf)-iz1))+qim0(min(nzf,(1:nzf)+iz2)));  % step increase in quota
  end
  % --------  CHL:C of drifters
  chl2c_Im= chl2c(qim1,parc,I_al,chl2c_min,chl2cf,temp);

  % -------- stronger co-limitation of immobile producers
  %          (no access to deep-water nutrients)
  colim= 1-colim0 + colim0 * exp(-no3_0^2);

  % -------- NPP of drifters in units mol-C/m3.d
  npp_Im = colim * min(Ilim.*qim1./(Q0+qim1),Pm0*Tf).*chlm0./chl2c_Im;%sqrt(limf)*

  %% fprintf(' Il=%1.2f D=%1.3f Pm*Tf=%1.3f*%1.3f\t C_1=%1.2f\n',mean(Ilim),mean(qim1./(Q0+qim1)),Pm0,Tf,chlm0(1)/chl2c_Im(1));
  if(zN < 225) % viable migration strategy found
    % -------- bio-pumping: N-conc of migrators times mortality
    jNv(pix,piy) = Nnew* Tfm*mort*(carbn/(2*b)+mean(chlm0./chl2c_Im));
    % -------- fraction of total NPP by migrators
    npp_frac(pix,piy) = sum(nppf*carbn)/sum(nppf*carbn + npp_Im);
    % --------  total NPP
    nppf   = nppf*carbn + npp_Im;
    %% fprintf('chl0=%1.2f %1.2f %1.2f %1.2f chln=%1.2f c0=%1.2f <strong>npp=%1.1f</strong>\n ',chl00,chlm0(1),fac,shalfac,mean(chlnew),mean(chlm0./chl2c_Im),nppf*dz');
  else % no viable migrators, only drifters at surface
   nppf = npp_Im;
   jNv(pix,piy) = 0;
   npp_frac(pix,piy) = 0;
  end

  % -------- total accumulated water column CHL  (g-CHL/m2)
  chl_modsumv(pix,piy) = chl_mod*dz';
  % -------- total accumulated CHL of migrators (g-CHL/m3)
  chl_modmaxv(pix,piy) = sum(chlnew);
  % -------- depth integrated NPP  (gC/m2.d)
  npptv(pix,piy) = nppf*dz'; % mol-C/gCHL*gC/mol-C
  % -------- chemocline depth
  zNv(pix,piy)   = zN;
  % -------- surface NO3
  no30v(pix,piy) = no3_0;
  % -------- migration amplitude (half distance)
  ampl(pix,piy)  = b;
  % -------- migration center (depth of DCM)
  if max(chl_mod)<(max(chl00,chl_mod(1)))*1.4+0.05
    zCvv(pix,piy) =-1;
  else
    zCvv(pix,piy) = zC;
  end
  % -------- biomass of migrators
  carbn          = max(0,carbn);
  carbnv(pix,piy)= carbn;
  % -------- migration speed
  speedv(pix,piy)= speed;
  % -------- surface CHL of migrators
  if length(chlnew)>0
     chl0mob(pix,piy)= chlnew(1)/(chl_mod(1)+0.0001);
  else
     chl0mob(pix,piy)= 0;
  end
  % -------- incident irradiance
  par0v(pix,piy) = par;
  end % if ~isnan(no3g(ix,iy,1))
 end % for py

 % ------- progress info in log file
 if mod(pix,3)==1
   fid=fopen(logfile,'a');
   fprintf(fid,'%d/%d\n',pix,maxx);
   fclose(fid);
   if pix>maxx*0.7 % save results to file at the end of the loop
     save(matfile,'chl_modsumv','chl_modmaxv','npptv','zCvv','carbnv','npp_frac','zNv','ampl','par0v','no30v','chl0mob','speedv','npp_frac','jNv','dx','dy','delzN');
   end

 end
end %for px
%%

%AreaT=360000000; % total ocean surface area
%fprintf('area=%1.1f %1.1f\t%1.3f Pg-C/yr\n',art*1E-6,AreaT*1E-6,nppglob*1E-12*365);

% ------- output of relevant results to file
save(matfile,'chl_modsumv','chl_modmaxv','npptv','zCvv','carbnv','npp_frac','zNv','ampl','par0v','no30v','chl0mob','speedv','npp_frac','jNv','dx','dy','delzN');
