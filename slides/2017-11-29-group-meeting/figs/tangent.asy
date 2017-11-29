    import three;
    import graph3;
    settings.outformat = "pdf";
    settings.prc = true;
    settings.embed= true;
    settings.render=16;

    size(6cm,6cm);
    currentprojection = orthographic(1, 1, 1);

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

    // void ellipse(real Theta, real Phi, real a, real b, real theta, bool dash, triple color) {
    //   triple normal = expi(Theta, Phi);
    //   real a_scaled = a/max(a, b);
    //   real b_scaled = b/max(a, b);      
    //   path3 mycircle = rotate(degrees(Phi), Z)*rotate(degrees(Theta), Y)*shift(Z)*rotate(degrees(theta), Z)*scale(a_scaled, b_scaled, 1)*circle(c=O, r=0.05, normal=Z);
    //   if (dash) {
    //     draw(mycircle, p=linetype(new real[] {8,8}, offset=xpart(color))+rgb(xpart(color), ypart(color), zpart(color)));
    //   } else {
    //     draw(mycircle, p=rgb(xpart(color), ypart(color), zpart(color)));
    //   }
    // }

    // void mydot(real Theta, triple color) {
    //   triple normal = expi(Theta, 0);
    //   dot(normal, p=rgb(xpart(color), ypart(color), zpart(color)));
    // }

    // void arrow(real Theta, real Phi_Pol, triple color, bool dash) {
    //   if (dash) {
    //     draw(rotate(Theta, Y)*rotate(Phi_Pol, Z)*(Z--(Z+0.2*X)), p=linetype(new real[] {4,4}, offset=xpart(color))+rgb(xpart(color), ypart(color), zpart(color)), arrow=Arrow3(emissive(rgb(xpart(color), ypart(color), zpart(color)))));
    //     draw(rotate(Theta, Y)*rotate(Phi_Pol, Z)*(Z--(Z-0.2*X)), p=linetype(new real[] {4,4}, offset=xpart(color))+rgb(xpart(color), ypart(color), zpart(color)), arrow=Arrow3(emissive(rgb(xpart(color), ypart(color), zpart(color)))));
    //   } else {
    //     draw(rotate(Theta, Y)*rotate(Phi_Pol, Z)*(Z--(Z+0.2*X)), p=rgb(xpart(color), ypart(color), zpart(color)), arrow=Arrow3(emissive(rgb(xpart(color), ypart(color), zpart(color)))));
    //     draw(rotate(Theta, Y)*rotate(Phi_Pol, Z)*(Z--(Z-0.2*X)), p=rgb(xpart(color), ypart(color), zpart(color)), arrow=Arrow3(emissive(rgb(xpart(color), ypart(color), zpart(color)))));
    //   }
    // }

    // void watson(real Theta, real Phi, real kappa, real x, real y, real z) {
    //  int n_phi = 10;
    //  int n_theta = 10;

    //  real max_radius = 0;
    //  if(kappa > 0){
    //    max_radius = exp(kappa);
    //  }
    //  else{
    //    max_radius = 1.0;
    //  }

    //  for(int i=0; i <= n_theta; ++i) {
    //    real Theta_i = pi*i/n_theta;
    //    real weight = exp(kappa*(cos(Theta_i)**2))/max_radius;     
    //    path3 mycircle = circle(c=Z*weight*cos(Theta_i), r=weight*sin(Theta_i));
    //    draw(shift((x, y, z))*rotate(angle=degrees(Phi), u=O, v=Z)*rotate(angle=degrees(Theta), u=O, v=Y)*mycircle);
    //  }

    //  triple f(real t) {
    //    real weight = exp(kappa*(cos(t)**2))/max_radius;
    //    return (0, weight*sin(t), weight*cos(t));
    //  }
    //  path3 phi_path = graph(f, 0, 2pi, operator ..);

    //  for(int i=0; i <= n_phi; ++i) {
    //    real Phi_i = 2*pi*i/n_theta;
    //    draw(shift((x, y, z))*rotate(angle=degrees(Phi), u=O, v=Z)*rotate(angle=degrees(Theta), u=O, v=Y)*rotate(angle=degrees(Phi_i), u=(0,0,0), v=(0,0,1))*phi_path);
    //  }
    // }
    // real len = 10;
    // draw((-len,-len)--(len,-len)--(len,len)--(-len,len)--(-len,-len), white);

draw(unitsphere, surfacepen=material(diffusepen=white+opacity(0.1), emissivepen=grey, specularpen=white));

// Draw points on sphere
dotfactor = 7;
dot(X); 
dot(Y); 
circle(0, pi/2, false, (0, 0, 0));

triple xi = expi(0, 2*pi/3);
triple yi = xi - dot(xi, X)*X;
triple zi = cross(xi, yi);
real d = 0.5;
draw(xi+d*yi+d*zi--xi+d*yi-d*zi--xi-d*yi-d*zi);
dot(xi, green);
dot(xi+0.5*yi, red);
dot(xi+0.5*zi, blue);

dot(Z);shipout(scale(4.0)*currentpicture.fit());
