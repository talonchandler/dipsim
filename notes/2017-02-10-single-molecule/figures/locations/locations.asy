settings.outformat="pdf";
unitsize(1cm);
usepackage("amsmath");


dotfactor=8;
dot((0, 0)); label("o", (-.2, -.2));
draw((0, 0) -- (-2, 5), arrow=Arrow()); label("$\mathbf{r}$", (-1.6, 2.5));
draw((0, 0) -- (2, 0.5), arrow=Arrow()); label("$\mathbf{r'}$", (1, 0));
draw((2, 0.5) -- (-2, 5), arrow=Arrow()); label("$\mathbf{R} \equiv \mathbf{r} - \mathbf{r'}$", (2.1, 2.0));
draw((2, 0.5) -- (2.8, 0.5), arrow=Arrow()); label("$\boldsymbol{\mu}$", (3,0.5));

draw((0,0) -- (0,0.8), arrow=Arrow()); label("$\boldsymbol{\hat{s}}$", (0.3, 0.7));

draw(arc((0,0),6,70,110));
draw(arc((0,11.28),6,70+180,110+180));

draw((0,0) -- (0, 5.2), dashed);


