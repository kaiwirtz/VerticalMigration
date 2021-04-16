% Structural and Lagrangian growth model
%    of migrating and drifting phytoplankton
%
%  analytical solution of balance equation during a vertical cycle
%    for the different growth phases (N- vs. PAR limitation)
%
% kai wirtz (HZG) Mar 2021
% ------ initial settings
if0  = i0;      % index of starting point
x0   = xx(i0);  % light level relative to the level at the chemocline center
i2   = 1:nz;    % index until end point
q0   = Qinit;   % initial quota

% ------ distinguishes light-saturation and -limitation
if (pN*x0*Qinit./(Q0+Qinit)*(1-rI*x0)>pmt*0.5 & zzi(i0)>zbal )
% ------------------------------------------------------------
% ------ light-saturation

  % ------ quota profile according to analytical solution Eq.(S6) in SI
   q=calc_lightsat(q,i2,i3,ud,xx,zz,vsink,Qpmax,Q0,pmt,q0,k_I,vup,nz);

  % ------ check for return to light limitation
  % ------     because of reduced N:C and light affinity
   i2 = find(pN*xx.*q./(Q0+q).*(1-rI*xx)<pmt);
   if i2,
     i23  = find(i2(2:end)-i2(1:end-1)>1); % indices of return conditions
     if i23, i2=i2(1):i2(i23); end
     i00 = 1+(ud==1)*(length(i2)-1);
     if0 = i2(i00);
     x0  = xx(if0);
     q0  = q(if0);
     p    = -ud*pN*nvv*(1-rI);
     % ------ light-limited quota profile
     calc_q_lightlim;
   end
else
% ------------------------------------------------------------
% ------ light-limitation
 p    = -ud*pN*nvv*(1-rCI); % rescaled photosynthesis rate
 %calc_q_lightlim
 % ------ light-limited quota profile
 q=calc_lightlim(q,i2,i3,ud,xx,fxv,q0,p,Qpmax,ui,if0,x0);
 %fprintf('LL x0=%1.4f %1.4f %1.4f %1.4f\n%1.4f>%1.4f\tp=%1.4f q=%1.4f\n',x0,pN,Qinit,rI,(pN*x0*Qinit./(Q0+Qinit)*(1-rI*x0)),pmt*0.5,mean(p),mean(q))

  % ------ check for return to light saturation
  % ------     because of increased N:C and light affinity
 i2 = find(pN*xx.*q./(Q0+q).*(1-rI*xx)>pmt);
 if i2,
   i00 = 1+(ud==1)*(length(i2)-1);
   q0 = q(i2(i00)); % N:C at first critical layer
   % ------ light-saturation
   % ------ quota profile according to analytical solution Eq.(S6) in SI
   q=calc_lightsat(q,i2,i3,ud,xx,zz,vsink,Qpmax,Q0,pmt,q0,k_I,vup,nz);
   %calc_q_lightsat;
   % ------ check for relaxation to limitation
   i2 = find(pN*xx.*q./(Q0+q).*(1-rI*xx)<pmt);
   if i2,
     i23  = find(i2(2:end)-i2(1:end-1)>1);
     if i23, i2=i2(1):i2(i23); end
     i00 = 1+(ud==1)*(length(i2)-1);
     if0 = i2(i00);
     x0  = xx(if0);
     q0  = q(if0);
     p    = -ud*pN*nvv*(1-rI);
     calc_q_lightlim;
   end
 end
end
