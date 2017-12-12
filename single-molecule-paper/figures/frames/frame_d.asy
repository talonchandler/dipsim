usepackage("amsmath");
settings.outformat = "pdf";
settings.render=16;
settings.embed=true;
settings.prc=true;
import three;
import graph3;
import solids;
import grid3;
defaultpen(fontsize(20pt));
currentprojection = orthographic(5,5,1.5);

size(6cm,0);

// Draw the lens.
real r = 1;
real alpha = pi/6;
real f(real x) { return sqrt(r^2 - x^2); }
path s = graph(f, 0, r*sin(alpha), operator..);
path3 b3 = path3(s, plane=YZplane);
surface solidsurface = surface(b3, c=O, axis=Z);
draw(solidsurface, grey+opacity(0.25));
label(L=Label("$\Omega$"), 0.9*r*Z);

// Draw the polarizer
pen thinblack = black+1;
path planeoutline = box((-0.5, -0.5), (0.5, 0.5));
draw(shift(2*r*Z)*surface(planeoutline), surfacepen=grey+opacity(0.25));

//Draw Axes
real ax_len = 0.5;
draw(L=Label("$\mathbf{\hat{x}}$", position=Relative(1.1), align=W), O--ax_len*X,thinblack, Arrow3(emissive(black))); 
draw(L=Label("$\mathbf{\hat{y}}$", position=Relative(1.1), align=E), O--ax_len*Y,thinblack, Arrow3(emissive(black))); 
draw(L=Label("$\mathbf{\hat{z}}$", position=Relative(1.0), align=N), O--ax_len*Z,thinblack, Arrow3(emissive(black)));

// Draw focal length labels with double arrow
triple A = 0.5*X-0.5*Y;
draw(L=Label("$f$", position=Relative(0.5), align=W), A--(A + r*Z), thinblack, Arrow3(emissive(black)));
draw(L=Label("$f$", position=Relative(0.5), align=W), (A+r*Z)--(A + 2*r*Z), thinblack, Arrow3(emissive(black)));
draw(L=Label("$f$", position=Relative(0.5), align=W), (A + r*Z)--A, thinblack, Arrow3(emissive(black)));
draw(L=Label("$f$", position=Relative(0.5), align=W), (A + 2*r*Z)--(A+r*Z), thinblack, Arrow3(emissive(black)));
draw(O--A, dashed);
draw(r*Z--(A+r*Z), dashed);

// Draw alpha
real r = 1;
real q = alpha; //theta
real f = -pi/4; //phi
triple B = r*expi(q,f);
draw(O--B, dashed);
draw(L=Label("$\alpha$", (0.5,0,1.5)), arc(O,0.15*length(B)*Z,0.5*B), thinblack);

// Choose vector
real r = 0.5;
real q = 0.25*pi; //theta
real f = 0.5*pi; //phi
triple A = r*expi(q,f);

// Draw other lines
draw(L=Label("$\mathbf{\hat{\mu}_{\text{abs}}}$", position=Relative(1.1), align=E), O--A, thinblack, Arrow3);
draw(2*Z--(2*Z+0.5*Y), thinblack, Arrow3);
label(L=Label("$\mathbf{\hat{p}_{\text{exc}}}$"), position=1.8*Z+0.5*Y);

shipout(scale(4.0)*currentpicture.fit());
