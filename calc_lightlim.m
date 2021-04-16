 % ------ light-limited quota profile
function q=calc_lightlim(q,i2,i3,ud,xx,fxv,q0,p,Qpmax,ui,if0,x0)

% ------ quota profile according to analytical solution Eq.(S4) in SI

cc   = (q0 - fxv(ui,if0))*exp(-p*x0);
q(i2)= cc*exp(p*xx(i2)) + fxv(ui,i2);

% ----- re-enforce upper Q-limit (ceasing uptake) by hard clip
i2=find(q>Qpmax+1E-4);
if i2,
 q(i2)=Qpmax;
 if ud==1
   if i3<i2(end)
     x0   = xx(i3);
     cc   = (Qpmax - fxv(ui,i3))*exp(-p*x0);
     q(1:i3) = cc*exp(p*xx(1:i3)) + fxv(ui,1:i3);
   end
 end
end

end
