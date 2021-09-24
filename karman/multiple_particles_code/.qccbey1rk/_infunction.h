
  {
#line 118

    g->x = val_rho(rho,0,0,0)/(val_cm(cm,0,0,0) + 1e-30)*(val_a_x(a.x,0,0,0) + val_a_x(a.x,1,0,0))/2.;
#line 118

    g->y = val_rho(rho,0,0,0)/(val_cm(cm,0,0,0) + 1e-30)*(val_a_y(a.y,0,0,0) + val_a_y(a.y,0,1,0))/2.;}
