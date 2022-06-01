function [F,dz,dx]=nearshadowfactor(alpha_f,gamma_v,dframe,tframe,u,dx2,u2,height,width);
%   Calculate shading factor. It is assumed that the distances u and tframe goes all around the window
%   and that the overhang is of infinite width.
%
%   [F,dz,dx]=shadowfactor(alpha_f,gamma_v,dframe,tframe,u,dx2,u2,height,width);
%
%   F: Shading factor
%   dx: Horisontal shading distance on glazing
%       Positive for alpha_f>=0
%       Negative for alpha_f<0
%   dz: Vertical shading distance on glazing
%
%   Input from incidence angle
%   alpha_f: Horisontal shadow angle [deg]
%   gamma_v: Vertical shadow angle [deg]
%
%    
%	dframe: width of frame[m]
%   tframe: frame thickness [m]
%	u: distance from wall surface to frame [m]
%   dx2: distance to overhang from window top [m]
%   u2: distance from wall surface to of overhang [m]
%   height: height of window (vertical length incl. frame) [m]
%   width: width of window (horisontal length incl. frame) [m]

alpha_f=alpha_f*pi/180;
gamma_v=gamma_v*pi/180;

%%%Top of window
if u~=0
    v_k=atan(dframe/u);
else
    v_k=pi/2;
end
logicgamma=(gamma_v>=v_k);
x=logicgamma.*((u+tframe)*tan(gamma_v)-dframe)+(logicgamma==0).*(tframe*tan(gamma_v));
x=(x>0).*x;
logicx=x>(height-2*dframe);
x=(logicx==0).*x+logicx*(height-2*dframe);

%%%Side of window
logicalpha=(abs(alpha_f)>=v_k);
y=logicalpha.*((u+tframe)*tan(abs(alpha_f))-dframe)+(logicalpha==0).*(tframe*tan(abs(alpha_f)));
y=(y>0).*y;
logicy=y>(width-2*dframe);
y=(logicy==0).*y+logicy.*(width-2*dframe);
logicalpha2=alpha_f<0;
dx=-logicalpha2.*y + (logicalpha2==0).*y;

area=x*(width-2*dframe)+y*(height-2*dframe)-x.*y;


%%% SHADE FROM OVERHANG %%%
u=u+tframe+u2;
xover=u*tan(gamma_v);

%Extra area
masx=xover-dx2-dframe-x;
masx=(masx>0).*masx;
masx=min(masx,height-2*dframe-x);
masy=width-2*dframe-y;
masy=max(0,masy);
area0=masx.*masy;
dz=min(height-2*dframe,max(x,xover-dx2-dframe));

%Glass area
Ag=(width-2*dframe)*(height-2*dframe);

%Vector containing shaded areas
F=1-(area+area0)/Ag;
