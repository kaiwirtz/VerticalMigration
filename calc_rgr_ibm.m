% Structural and Lagrangian growth model
%    of migrating and drifting phytoplankton
%
% core model calculating growth and quota changes along a vertical cycle
%
% kai wirtz (HZG) Mar 2021

% ------- reset auxiliary variables
muS = -1;
rgr = -1;
Q0m = Q0;
muS = 0; mur =0; muqS= 0; qm=0; iS = 0;qs=[]; muX=0;

% fast swimming/floating species (dinos, huge diatoms) have lower growth rates
% here to 15% at ve=1 and vsink=60
% ve: strength of pyhsiology-speed trade-off
dv=exp(-0.5*ve);
scalfac = 1./(1+exp((log(vsink)-log(vst)-(1-ve)*2)/dv));

% ------- apply trade-off to growth parameters
aff_I = aff_I0*scalfac;
Pm    = Pm0*scalfac;

% ------- indices defining the vertical migration path
ind = find(z>=zC-b-1E-3 & z<=zC+b+1E-3);
nz=length(ind);

if nz>1 % extending over one vertical layer?
 zzi = zz(ind); % relative layer depth along migration path
                % distance to chemocline center position
 [ff izN] = min(abs(zzi)); % closest layer to center

 if calc==0
  runt = 2*b/vsink;   % migration time
  npp  = zeros(1,nz); % clear storage
  addC = zeros(1,length(z));
  q    = zeros(1,length(zzi));
  clear fxv;

 % ------- light attenuation coefficients
 % ------- self-shading incl resting (mobile) population
  kcd   = kchl*chln/(2*b); %*(1+trest/runt);
% -------  background attenuation + self-shading by immobiles
  k_I   = k0 + kchl*chl00 + kcd; % 1/m

  nn    = kN/k_I; % slope ratio
% ------- light at chemocline
  par_N0 = par*exp(-(k0+kchl*chl00)*zN);

  par_N = par_N0*exp(-kcd*min(max(0,zN-(zC-b)),2*b));
  par_C = par_N;

%    fprintf('k0=%1.4f kc=%1.4f c=%1.4f kcd=%1.4f c=%1.4f\n',k0, kchl,chl0,kcd,chln)
%    fprintf('* par=%1.4f %1.4f k0=%1.6f %1.6f zN=%1.4f kcd=%1.4f\n',par,par_N0 ,k0,kchl*chl0,zN,kcd)

%  aff_IN= aff_I*chl2n(par_N,I_al); % mol-C/mol-N * m2s/µmolE
%  ------- quota-specific light limitation at chemocline
  I_N   = par_N * aff_I * Tfa; % mol-C/mol-N
  rI    = min(par_N/I_al,0.9);
  rCI   = rI;
%  ------- growth and uptake at the chemocline
%  ------- assumes nutrient saturation in central chemocline
  vup   = A_N; % N-uptake rate mol-C/mol-N 1/d
  pmt   = Pm*Tf; % max C-uptake rate
  pN    = I_N; % C-uptake rate at chemocline center
%  ------- critical Q-balance depth
%  zbal=  (log(pN/vup) - rI)/(kN + (1-rI)*k_I);
  zbal  = (log(pN*(1-min(1.5*rI,0.95))*Qpmax/vup))/(4*k_I);
%  fprintf('k_I=%1.4f parN=%1.4f zbal=%1.4f pN=%1.4f %1.4f\tv=%1.4f\n',k_I,par_N,zbal,pN,pmt,vup)

  [fm i3]= min(abs(zbal-zzi));
%  plot(qq,zN+zbal,'g--','LineWidth',2);
%  fprintf('b=%1.1f z:%1.1f-%1.1f\tUp:%1.3f-%1.3f\n',b,zzi(1),zzi(end), vup*exp(kN*zzi(1)), vup*exp(kN*zzi(end)));
%  ------- PAR saturated depth
  xx  = exp(-k_I*zzi); % experienced light level relative to the level at the chemocline center (\zeta in SI)
  parz= par_N*xx;
  nvv = 1/(k_I*vsink); % normalization factor, unit time
  pp2 = pN*nvv*(1-rI); % rescaled C-uptake
  pp2 = pp2*pp2;
  ixi = exp(k_I*zzi);  % inverse of xx
  ixi2= ixi.*ixi;
  %  ------- loop over swimming direction (up, down)
  for ud=[1 -1]
    ui=1+(ud==-1);
    p = -ud*pN*nvv*(1-rCI); % rescaled photosynthesis rate, see Eq(S3) in SI
    %  ------- quota profile along up-/downward migration, see Eq(S4) in SI
    fxv(ui,:)  = -ud*vup*nvv/6*(p*pp2*exp(p*xx).*(-real(expint(p*xx)))+pp2*ixi-p*ixi2+2*ixi2.*ixi);
   end
 end %if calc
 % ------- trajectory starts below chemocline
 if max(isinf(fxv))==0
   % ------- initialize quota
   Qinit = Qpmax*0.15; Qstart=0; %
   Qinit_ode=Qinit;
   t=0;
   tmax=round(16*max(1,2*vsink/b));
   % ------- loop over cycles until convergence of q-profile
   while ( abs(Qstart-Qinit)>9E-5 & t<tmax) %
     Qstart=Qinit;
     trest = 0;
     for ud=[1 -1] % loop over swimming direction (up, down)
       if ud==1
         i0  = nz;  % index of starting point
         i1  = 1;   % index of end point
       else
         i0  = 1;   % index of starting point
         i1  = nz;  % index of end point
       end
       ui=1+(ud==-1);
       % ------- analytical solution of balance equation gives q-profile
       calc_quota;
       qs(ui,:)=q; % store in vector

       % ------- check for resting phase, only at upper turning points
       if ud==1 & q(i1)>Q0m & trm  & zC-b<=16
         rq = pN*xx(i1)*(1-rI*xx(i1)); % rescaled photosynthesis
         if rq>pmt % larger than maximal photosynthesis rate ?
           % ------- remaining nutrient store:
           qc = Q0*(pmt/rq)/(1-pmt/rq);
           % ------- stores allow for unlimited C-uptake?
           if q(i1)>qc  % split into two phases for transition LightSat-Lim (P>Pm)
             p = pmt;   % C-uptake at maximal photosynthesis rate
             trest = log(q(i1)/qc)/p;  % resting time; see Sec.S2 in SI
             p = rq*qc/(Q0+qc); % C-uptake at N-limited photosynthesis rate
  %fprintf('q=%1.3f %1.3f\tp=%1.2f %1.3f\tT:%1.2f %1.2f \n',q(i1),qc,p*exp(-k_I*zzi(i1)), trest, log(qc/Q0m)/p);%*exp(k_I*zzi(i1))
             trest = trest + log(qc/Q0m)/p;  %
           end
         end % if rq>pmt
         if trest<1E-5 % no phase splitting <- N-limitation at upper turning point
            p = rq*q(i1)/(Q0+q(i1));  % photosynthesis rate at resting position
            trest = log(q(i1)/Q0m)/p; % resting time; see Sec.S2 in SI
         end
         qrestf(ui) = Q0m/q(i1); % resting quota factor
       else
         qrestf(ui) = 1; %no rest
       end % if resting

       %  --------------------------------------------------------------------------------
       %  ------- numerical solution of the balance equation for validating the analytical one
       if 0
         dqdt=@(t,x) [-ud*vsink; vup*exp(3*k_I*x(1))*(x(2)<Qpmax)-min(pN*exp(-k_I*x(1))*(1-rI*exp(-k_I*x(1))),pmt*(Q0+x(2))/x(2))*x(2);];
         [tt,x_ode] = ode23s(dqdt,[0 runt],[zzi(i0) Qinit_ode]);
         if(ouc==3)
           plot(x_ode(:,2),zN+x_ode(:,1),'-','Color','r','LineWidth',3+(t>3)+(t>1));
         end
         Qinit_ode = min(x_ode(end,2)*qrestf(ui),Qpmax); %*qrestf(ui);
       end
       %  --------------------------------------------------------------------------------

       Qinit = q(i1)*qrestf(ui); % initialize next round
       % ------- if graphical output: plot quota-hysteresis
       if(ouc==3) plot(q,z(ind),'-','Color','k','LineWidth',1+(t>3)+(t>1)); end%ones(3,1)*max(0.7-t*0.05,0)1+t*0.5
       if (outc|ouc) & 0, fprintf(' %d %d/%d\tQ %1.3f %1.3f .. %1.3f %1.3f\n',ud,t,tmax,Qstart,Qinit,q(i0),q(i1)); end
     end
     t=t+1;
     %%fprintf('+ %d + q T =%1.4f %1.4f %1.4f %1.4f %1.4f \n',t,q(1),q(10),q(20),q(30),q(end))
   end % while

 %% ------- calc and plot RGR
 %% ------- convergence of q-iteration to finite values?
   if Qinit>2E-3 & (t<tmax | Qinit>0.015)
  %fprintf('b,zC = %1.1f %1.1f %d Qinit=%1.4f\n',b,zC,t,Qinit);

    if  qs(1,1)>1.2*Q0m & trm==0 & zC-b<=15 & ldif<20
      trmax=1;
    end
  % ------- store inidices of migration layer for later profile reconstruction
    ii1(j)=ind(1); ii2(j)=ind(end);
    chl2cv   = zeros(1,nz);
  % ------- loop over swimming direction (up, down)
    for ud=[1 -1]
      i0  = (ud==-1) + (ud==1)*nz;
      q=qs(1+(ud==-1),:); % retrieve q-profile
  %    chl2ct = 0.3*q.*chl2n(parz,I_al);
    % ------- calculate CHL:C depending on N:C, PAR (T-dependence is obsolete)
      chl2ct = chl2c(q,parz,I_al,chl2c_min,chl2cf,temp);
      chl2cv = chl2cv+chl2ct;

  % ------- calculate q-dependent gross growth rate
      mu  = min(pN*xx.*q./(Q0+q).*max(1-rI*xx,0),pmt);
      muqS = muqS + mean(mu.*q); % store average of N-assimilation
  %%    fprintf('%1.2f %1.2f\t k=%1.3f %1.3f  pN=%1.2f\t',mu(1),mu(end),k_I ,kchl*chl0,pN);
      muS = muS + mean(mu);      % store average of C-uptake
      muX = max(max(mu),muX);    % maximal realized C-uptake

      npp = npp + mu;
      qm  = qm  + mean(Q0+q);

      % ------- plotting ?
      switch ouc
        case 1 % , 3 plot quota
         %cf=max(min((1-1.5*q/Qmax)*1.3,1),0);
         col = max(ones(3,1)*(1.0-12*mu),0);
         h=colormapline(q,z(ind),[],col');
         set(h,'linewidth',lw,'linestyle','-')
  %       plot(q,z(ind),'-','Color',ones(3,1)*0,'LineWidth',1); end
        case 2
         cf=max(min((1-1.4*q/Qmax)*1.0,1),0);
         col = ones(3,1)*cf;
         %fprintf('%1.2f--%1.2f\n',min(col),max(col))
         h=colormapline(mu,z(ind),[],col');
         set(h,'linewidth',lw,'linestyle','-')
      end
    end % for ud
  %%  fprintf('\n');
  % ------- averages from up/downward sums
    muS = muS*0.5; %
    npp = npp*0.5;
    chl2cv=chl2cv*0.5;

  % ------- adjust average growth in case of resting phase
    if trest>0  %
      dr  = 2*runt/(2*runt+trest);  % ratio of swimming over total (resting+swimming) time
     % TODO: check time integral of resting growth rate
      mur = min(pN*xx(1)*(1-rI*xx(1))*(qs(1,1)/(Q0+qs(1,1))+qs(2,1)/(Q0+qs(2,1))) *0.667,pmt);%Q0m./(Q0+Q0m)+log((qs(1,1)+Q0)/(qs(1,1)+Q0*exp(Pr*trest)))/trest;
  %    fprintf('%1.2f %1.3f\n',Ilim(1),qs(1,1));
    %%  fprintf('T_R=%1.1f dr=%1.2f  mu=%1.3f->%1.3f\tq: %1.3f %1.3f\trI: %1.3f\n',trest,dr,muS,mur,qs(1,1),qs(2,1),rI);
      muS = dr*muS+(1-dr)*mur;  % adjust growth rate: time average
      muqS = dr*muqS+(1-dr)*mur*0.5*(qs(1,1)+qs(2,1));  % approx. quota during resting
  %    addC= kR*exp(-kR*(zzi-zzi(1))) *trest/runt; % relative (C-) contribution of top resting
      bdiff= sqrt(2*Kz); % diffusion length
      addC = exp(-0.5*((zz-zzi(1))/bdiff).^2)/(bdiff*2.5); % relative (C-) contribution of top resting
      npp  = dr*npp+(1-dr)*addC(ind)*mur;    % adjust NPP: vertical profile (time!=space integration)
      muS  = dr*muS+(1-dr)*mur;    % adjust NPP: vertical profile (time!=space integration)
    else
      dr=1;
    end
    % ------- RGR: net growth includes mortality
    rgr = muS - Tfm * omega * carb0; %TODO replace with variable Carbn estimate
   end
  end % if Qinit
end % if ind
