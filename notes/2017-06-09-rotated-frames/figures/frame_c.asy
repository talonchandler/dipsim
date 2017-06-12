usepackage("amsmath");
settings.outformat = "pdf";
settings.render=16;
settings.embed=true;
settings.prc=true;
import three;
import graph3;
import grid3;
defaultpen(fontsize(20pt));
currentprojection = perspective(5,5,3);

size(12cm,0);

//Draw Axes
pen thinblack = black+1.0;
real ax_len = 1.0;
draw(L=Label("$\mathbf{\hat{x}}$,$\mathbf{\hat{z}'}$", position=Relative(1.1), align=SW), O--ax_len*X,thinblack, Arrow3(emissive(black))); 
draw(L=Label("$\mathbf{\hat{y}}$,$\mathbf{\hat{y}'}$", position=Relative(1.1), align=E), O--ax_len*Y,thinblack, Arrow3(emissive(black))); 
draw(L=Label("$\mathbf{\hat{z}}$", position=Relative(1.1), align=N), O--ax_len*Z,thinblack, Arrow3(emissive(black)));

//Draw Alt. Axes
real psi = 0.5*pi;
triple Xp = cos(psi)*X + -sin(psi)*Z;
triple Yp = Y;
triple Zp = sin(psi)*X + cos(psi)*Z;

draw(L=Label("$\mathbf{\hat{x}'}$", position=Relative(1.1), align=S), O--ax_len*Xp,thinblack, Arrow3(emissive(black))); 
draw(O--ax_len*Yp,thinblack, Arrow3(emissive(black))); 
draw(O--ax_len*Zp,thinblack, Arrow3(emissive(black))); 

// Choose vector
real r = 1;
real q = 0.25*pi; //theta
real f = 0.3*pi; //phi
triple A = r*expi(q,f);

// Find projection of A into Xp, Yp plane
triple P_long = cross(Zp, cross(A, Zp));
triple P = P_long/length(P_long);

// Draw other lines
draw(L=Label("$\mathbf{\hat{r}}$", position=Relative(1.1), align=N), O--A, thinblack, Arrow3);
draw(O--P, thinblack+dashed);
draw(arc(O,A,P), thinblack+dashed);
draw(L=Label("$\psi=\pi/2$", (4,0,4)), arc(O,0.4*Z,0.4*expi(psi,0)), thinblack);
draw(L=Label("$\phi'$", (0,1,-1)), arc(O,0.15*Xp,0.15*P), thinblack);
draw(L=Label("$\theta'$", (1,0,1.5)), arc(O,0.2*length(A)*Zp,0.2*A), thinblack);

shipout(scale(4.0)*currentpicture.fit());
