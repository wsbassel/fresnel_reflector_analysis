getd ("FRESNEL")

year=2021
    month=3
    day=1
    hour=8
    hr=hour
    mn=0
    sc=0
    
    logitude=-46.6334   
    latit=-23.55
    time_zone=-3
         
   jy_col=-18
   kz_col=16
   ix_col_guess=-3
   j_refl_guess=.4
  [azumith,zenith,i_refl,j_refl,k_refl,ix_col]=fresnel_sun_traking(year,month,day,hour,mn,logitude,latit,time_zone,jy_col,kz_col,ix_col_guess,j_refl_guess)
   disp azumiyh
   disp (azumith)
   disp zenith
 disp (zenith)
 
 disp i_j_k_refl
 disp (i_refl)
 disp (j_refl)
 disp (k_refl)
 disp  ix_col
 disp (ix_col)
 disp  jy_col
 disp (jy_col)
  disp kz_col
  disp (kz_col)
  reflector_angle=asin(j_refl)*180/%pi
  disp 'reflector angle'
  disp (reflector_angle) 
