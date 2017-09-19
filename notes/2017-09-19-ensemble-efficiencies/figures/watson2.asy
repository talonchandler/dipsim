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

size(6cm,6cm);

//Draw Axes
pen thinblack = black+1;
// real ax_len = 1.0;
// draw(L=Label("$\mathbf{\hat{x}}$", position=Relative(1.1), align=SW), O--ax_len*X,thinblack, Arrow3(emissive(black))); 
// draw(L=Label("$\mathbf{\hat{y}}$", position=Relative(1.1), align=E), O--ax_len*Y,thinblack, Arrow3(emissive(black))); 
// draw(L=Label("$\mathbf{\hat{z}}$", position=Relative(1.1), align=N), O--ax_len*Z,thinblack, Arrow3(emissive(black)));

    void circle(real Theta, real Alpha, bool dash, triple color) {
      triple normal = expi(Theta, 0);
      real h = 1 - sqrt(2 - 2*cos(Alpha) - sin(Alpha)^2);
      real radius = sin(Alpha);
      path3 mycircle = circle(c=h*normal, r=radius, normal=normal);
      if (dash) {
        draw(mycircle, p=linetype(new real[] {8,8}, offset=xpart(color))+rgb(xpart(color), ypart(color), zpart(color)));
      } else {
        draw(mycircle, p=rgb(xpart(color), ypart(color), zpart(color)));
      }
    }

void watson(real Theta, real Phi, real kappa, real x, real y, real z) {
     int n_phi = 10;
     int n_theta = 10;

     real max_radius = 0;
     if(kappa > 0){
       max_radius = exp(kappa);
     }
     else{
       max_radius = 1;
     }

     for(int i=0; i <= n_theta; ++i) {
       real Theta_i = pi*i/n_theta;
       real weight = exp(kappa*(cos(Theta_i)**2))/max_radius;     
       path3 mycircle = circle(c=Z*weight*cos(Theta_i), r=weight*sin(Theta_i));
       draw(shift((x, y, z))*rotate(angle=degrees(Phi), u=O, v=Z)*rotate(angle=degrees(Theta), u=O, v=Y)*mycircle);
     }

     triple f(real t) {
       real weight = exp(kappa*(cos(t)**2))/max_radius;
       return (0, weight*sin(t), weight*cos(t));
     }
     path3 phi_path = graph(f, 0, 2pi, operator ..);

     for(int i=0; i <= n_phi; ++i) {
       real Phi_i = 2*pi*i/n_theta;
       draw(shift((x, y, z))*rotate(angle=degrees(Phi), u=O, v=Z)*rotate(angle=degrees(Theta), u=O, v=Y)*rotate(angle=degrees(Phi_i), u=(0,0,0), v=(0,0,1))*phi_path);
     }
}
// watson(0, pi/2, -5, -2, 0, 0);
// watson(pi/8, pi/2, -1, 0, 0, 0);
// watson(pi/4, pi/2, 0, 2, 0, 0);
// watson(3*pi/8, pi/2, 1, 4, 0, 0);
// watson(pi/2, pi/2, 5, 6, 0, 0);
watson(0, 0, 1, 0, 0, 0);

real len = 10;
draw((-len,-len)--(len,-len)--(len,len)--(-len,len)--(-len,-len), white);

shipout(scale(4.0)*currentpicture.fit());
