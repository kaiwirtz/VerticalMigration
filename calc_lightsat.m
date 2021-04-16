% ------ light-saturation
function q=calc_lightsat(q,i2,i3,ud,xx,zz,vsink,Qpmax,Q0,pmt,q0,k_I,vup,nz)

% ------ quota profile according to analytical solution Eq.(S6) in SI
i00  = 1+(ud==1)*(length(i2)-1);
x0   = xx(i2(i00));
tim  = -ud*(zz(i2)-zz(i2(i00)))/vsink;
frat = vup/(pmt -ud*3*k_I*vsink);
arg  = frat*power(xx(i2),-3);

q(i2)= -Q0 +(q0+Q0-arg(i00))*exp(-pmt*tim)+arg;
%fprintf('lightsat %dx0,frat,\tx0=%1.3f frat=%1.3f tim=%1.3f \n',length(i2),x0,frat,tim);q

% ----- re-enforce upper Q-limit (ceasing uptake) by hard clip
i2=find(q>Qpmax+1E-4);
if i2,
  q(i2)=Qpmax;
  if ud==1
    if i3<i2(end)
      i2=1:i3;
      tim  = -ud*(zz(i2)-zz(i3))/vsink;
      arg  = frat*power(xx(i2),-3);
      q(i2) = -Q0 +(Qpmax+Q0-arg(i3))*exp(-pmt*tim)+arg(i2);
    end
  else
    q(i2(end):nz)=Qpmax;
  end
end

end
