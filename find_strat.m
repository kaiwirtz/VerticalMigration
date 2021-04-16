% Structural and Lagrangian growth model
%    of migrating and drifting phytoplankton
%
% find all near-optimal vertical strategies
%
% kai wirtz (HZG) Mar 2021
% ------- step sizes of strategy screening
db = 5; dzC = 5;

% ------- reference strategy: centered on top of the chemocline
bmaxr=min(bmax,min(ceil(zN)-2*db,bmax+4-ldif));

j=1; % strategy index

% ------- is migration active?
if bmaxr>3
  % ------- estimate RGR for average strategy
  zC = zN*0.85-10;
  %zC = zN*1-0;
  % ------- initial gues of traits
  b=bmaxr/5;  % amplitude
  vsink=2*vs0;% speed
  trest=0;    % resting time at top of cycle

  calc=0;  % switch for full calculation

% -------  fast swimming/floating species (dinos, huge diatoms) have lower growth rates
% -------  strength of pyhsiology-speed trade-off
  ve = shalfac + 2/(1+exp((4-no3_0)));

  % ------- estimate RGR for average strategy
  calc_rgr_ibm;
  calc=0;

  mu_R=0.8*muS; % reference RGR
  rgr_crit = 0.8 * max(rgr,0);
  rgr_m=rgr;rgr_1=rgr;

  % ------- clear storage
  iv=[];
  if exist('mu_v') clear mu_v b_v zC_v tr_v vs_v rgr_v; end
  mu_v(1)=muS; b_v(j)=b; zC_v(j)=zC; tr_v(j)=0;
  vs_v(j)=0;rgr_v(j)=rgr;

  % ---------------------------------------------------------------
  % ------- loop over migration amplitude
  for b = 3:db:bmaxr  %
    mort = omega; % specific mortality rate, without temperature effect
    % ------- biomass of migrants follows from density dependence of mortality and steady state
    carbn = sqrt(1+min(3E3,max(mu_R/(Tfm*mort)-carb0,0)*2*b));%TODO
    chln  = carbn*c2c0; % concomitant CHL of migrants

    % ------- loop over search direction
    for udv=[1 -1]
      zC = max(0.9*zN - 0.9*b - (udv<0)*dzC,b+dzC); % start depth
     %%fprintf('%d %1.1f chln=%1.1f rgr=%1.3f/%1.3f\n',b,zC,chln,rgr_1,rgr_crit);

      % ------- loop over mean migration depth zC (corresponds to SCM center position)
      while ((rgr_1 > max(rgr_crit,0)+0.025 | abs(zC-max(0.9*zN-0.9*b-0,b))<=2*dzC) & zC-b>=0) % viability check
        vsink = vs0;  % initial speed
        vsmaxb= min(4*b,vsmax); % max speed
        muS=0; rgr_1=-1;

      % ------- loop over travel speed
        while vsink<1.001*vsmaxb %
          trmax=0; trm=0;

          % ------- loop over resting times
          while trm<=trmax
           calc_rgr_ibm % estimate RGR for current strategy
           %fprintf('%d zC=%1.0f b=%1.0f v=%1.1f/%1.1f\tmu=%1.4f %1.3f/%1.3f/%1.3f\t%d->%d\t q=%1.3f Par=%1.2f\n ',j,zC,b,vsink,vsmaxb,muS,rgr,rgr_1,rgr_crit,trm,trmax,mean(q),par_N);

           if(rgr>rgr_1) % exceeds maximum RGR so far ?
            % fprintf('%1.1f+-%1.1f/%1.1f \trgr %1.3f %1.3f %1.3f %1.3f\n',zC,b,bmaxr,rgr,rgr_1,rgr_m,rgr_crit )
             rgr_1=rgr;  % store new optimum RGR (net growth=productivity-mortality)
             muM=muS; % productivity
             nppM=npp; vsinkm=vsink; trestm=trest;
             muq=muqS;qmm=qm;drm=dr;qM=q;
           end
           trm=trm+1;
          end
        %  fprintf('%1.1f %1.1f\t%1.1f %d\tmu:%1.3f rgr:%1.3f %1.3f\n',zC,b,vsink,trm,muS,rgr,rgr_1)

          vsink = 2*vsink; % increase in speed
        end %while vsink
    %%   fprintf('%1.1f %1.1f\tchln=%1.1f v=%1.1f %d\tC=%1.3f\t%1.3f %1.3f>%1.3f\n',zC,b,chln,vsinkm,trm,carbn,muM,rgr_1,rgr_crit)
  %%fprintf(' %d b=%d zC=%1.1f\tv=%1.1f\trgr: %1.3f %1.3f\t %1.1f\n',j,b,zC,vsink,rgr_1,rgr_crit,trest);
        % ------- viable strategy found ?
        if(rgr_1 > rgr_crit & ~isnan(mean(q)))
           %  store characteristics
           nppv{j} = nppM;  % NPP profile of migrants
           q_v{j}  = qM;    % N:C profile of migrants
           mu_v(j) = muM;   % mean RGR (of migrants)
           b_v(j)  = b;     % migration amplitude (half distance)
           rgr_v(j)= rgr_1; % RGR (net growth=productivity-mortality)
           muM=muS;         % gross growth rate = productivity
           dr_v(j) = drm;   % ratio of swimming over total (resting+swimming) time
           zC_v(j) = zC;    % mean migration depth = SCM center position
           tr_v(j) = trestm;% resting times
           addCv{j}= addC;  % relative C contribution of top resting
           vs_v(j) = vsinkm;% speed
           muq_v(j)= muq;   % N-assimilation (photosynthesis*quota)
           qm_v(j) = qmm*0.5;% mean N:C quota
           chl2cm(j,ii1(j):ii2(j)) = chl2cv; %CHL:C profile
           chl2cmv(j)=mean(chl2cv); %mean CHL:C
           %  adjust reference RGR if too low
           if(rgr_1>rgr_m) rgr_crit = 0.8 * max(rgr_1,0)-0.; rgr_m=rgr_1; end
           j=j+1;
        end
        zC = zC + udv*dzC; % de/increase center position
    %fprintf('find %d b=%d zC=%1.1f %d\t%1.3f<%1.3f/%1.3f\t%1.1f\n',j,b,zC,udv,rgr_1,rgr_crit,rgr_m,trest.vsinkm);

      end %while zC
    end %for udv
  end % for b
  % ------- calculate statistics for selected strategies
  iv=[];
  if max(vs_v)>0 % ------- viable strategy found ?
    % ------- identify optimal strategy
    [rgr_m ivm] = max(rgr_v);
    mu_m=mu_v(ivm);

    fprintf(' opt %d in %d b=%1.1f zC=%1.1f v=%1.1f ve=%1.2f T=%1.1f/%1.1f mu=%1.2f\tjN=%1.1f %1.1f\n',ivm,length(mu_v),b_v(ivm),zC_v(ivm),vs_v(ivm),ve,tr_v(ivm),2*b_v(ivm)/vs_v(ivm),mu_m,qm_v(ivm),mean(qm_v));

    if rgr_m>-0.11
      % ------- identify near optimal strategy
      if rgr_m>0
        iv=find(rgr_v>rgr_m*crit-0.01); % index of near-optimal strategies
      else
        iv=find(rgr_v>rgr_m/crit);
      end
      e_rgr=exp(rgr_v(iv)*7); %  weighing factor according to RGR
      igf = 1./sum(e_rgr);    %  normalization factor
      e_rgr=e_rgr*igf;        %  normalized weighing factor
    % ------- clean initial profiles ( carbon, chl, npp)
      xf=zeros(1,length(zf)); cf=xf; nppf=xf;qtot=xf;totnew=xf;

    %  ------- construct profiles from best strategies
      b = 0.;zC= 0.; newC=0; speed=0; trest=0;isp=1./sqrt(2*3.1415);bdiff=0;
      for j=1:length(iv) % loop over near-optimal strategies
        ivj=iv(j);
      % ------- calculate mean width for density/new_rat
        b      = b + e_rgr(j)*b_v(ivj);
        zC     = zC + e_rgr(j)*zC_v(ivj);
        newC   = newC + e_rgr(j)*0.5/b_v(ivj);
        speed  = speed+ e_rgr(j)*vs_v(ivj);
      %%fprintf(' %1.0f:%1.0f->%1.3f ',b_v(ivj),e_rgr(j)*1E3,b);
        trest = trest + e_rgr(j)*tr_v(ivj);
        % ------- diffusion length
        bdif = (4.5+ldif)*sqrt(b_v(ivj)/vs_v(ivj))+ 0.5*b_v(ivj); %contributiom from biodiffusion
      %  ------- C density
        new   = (0.25/b_v(ivj)*(erf((zf-zC_v(ivj)+b_v(ivj))/bdif) + erf((-zf+zC_v(ivj)+b_v(ivj))/bdif)))*e_rgr(j);  %*addCv{ivj}
        bdiff = bdiff+e_rgr(j)*bdif;
        totnew= totnew+new;

        %  ------- N profile
        qtotn = new.*[ones(1,ii1(ivj)-1)*q_v{ivj}(1) q_v{ivj} ones(1,length(zf)-ii2(ivj))*q_v{ivj}(end)]; %mol-C -> g-C
        qtot = qtot  + qtotn;

        %Interdiffusion  in two semi-infinite bodies "x=4/sqrt(Dt)"
        newi= ii1(ivj):ii2(ivj);
        %  new(newi)=new(newi)+addCv{ivj}*e_rgr(j);
        %  ------- rescale density with resting population
        new = dr_v(ivj)*new+ (1-dr_v(ivj))*addCv{ivj}*e_rgr(j);
        %  ------- C profile
        xf  = xf + new;
        %  ------- Chl profile
        cf  = cf + new.*[ones(1,ii1(ivj)-1)*chl2cm(ivj,newi(1)) chl2cm(ivj,newi) ones(1,length(zf)-ii2(ivj))*chl2cm(ivj,newi(end))];
        %  ------- NPP profile
        nppfn = new.*[ones(1,ii1(ivj)-1)*nppv{ivj}(1) nppv{ivj} ones(1,length(zf)-ii2(ivj))*nppv{ivj}(end)]; %mol-C -> g-C
        nppf =  nppf + nppfn;
      end
   end
  end
end
