usepackage("amsmath");
settings.outformat = "pdf";
settings.render=16;
settings.embed=true;
settings.prc=true;
import three;
import graph3;
import grid3;
defaultpen(fontsize(20pt));
currentprojection = perspective(5,3,3);

size(8cm,0);

//Draw Axes
pen thinblack = black+1;
real ax_len = 1.0;
draw(L=Label("$\mathbf{\hat{x}}$", position=Relative(1.1), align=SW), O--ax_len*X,thinblack, Arrow3(emissive(black))); 
draw(L=Label("$\mathbf{\hat{y}}$", position=Relative(1.1), align=E), O--ax_len*Y,thinblack, Arrow3(emissive(black))); 
draw(L=Label("$\mathbf{\hat{z}}$", position=Relative(1.1), align=N), O--ax_len*Z,thinblack, Arrow3(emissive(black)));

//Draw Alt. Axes
real psi = 0.1*pi;
triple Xp = cos(psi)*X + -sin(psi)*Z;
triple Yp = Y;
triple Zp = sin(psi)*X + cos(psi)*Z;

// draw(L=Label("$\mathbf{\hat{x}'}$", position=Relative(1.1), align=SW), O--ax_len*Xp,thinblack, Arrow3(emissive(black))); 
// draw(L=Label("$\mathbf{\hat{y}'}$", position=Relative(1.1), align=E), O--ax_len*Yp,thinblack, Arrow3(emissive(black))); 
// draw(L=Label("$\mathbf{\hat{z}'}$", position=Relative(1.1), align=N), O--ax_len*Zp,thinblack, Arrow3(emissive(black))); 

// Choose vector
real r = 1;
real q = 0.25*pi; //theta
real f = 0.3*pi; //phi
triple A = r*expi(q,f);

// Draw other lines
draw(L=Label("$\mathbf{\hat{p}_{\text{exc}}}$", position=Relative(1.1), align=E), O--expi(pi/2,f), thinblack, Arrow3);
draw(O--expi(pi/2,f), thinblack+dashed);
draw("$\phi_{\text{exc}}$", arc(O,0.15*X,0.15*expi(pi/2,f)), thinblack);
shipout(scale(4.0)*currentpicture.fit());
