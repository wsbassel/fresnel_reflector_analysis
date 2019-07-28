function [azumith,zenith,i_refl,j_refl,k_refl,ix_col]=fresnel_sun_traking(year,month,day,hour,mn,logitude,latit,time_zone,jy_col,kz_col,ix_col_guess,j_refl_guess)
   
   //-------------------------------
   ix_col=ix_col_guess
   j_refl= j_refl_guess
    if month==1  then
        N=day
    end
    if month==2 then
        N=30+day
    end
    if  month==3 then

        N=31+28+day
   
    end

    
   
    disp (N)
    if month==4 then
        N=31+28+31+day
    end
    if month==5 then
        N=31+28+31+30+day
    end
    if month==6 then
         N=31+28+31+30+31+day
    end
    if month==7 then
         N=31+28+31+30+31+30+day
    end
    if month==8 then
        N=31+28+31+30+31+30+31+day
    end
     if month==9 then
        N=31+28+31+30+31+30+31+31+day
    end
    
     if month==10 then
        N=31+28+31+30+31+30+31+31+30+day
    end
     if month==11 then
        N=30+28+31+30+31+30+31+31+30+31+day
    end
      if month==12 then
        N=31+28+31+30+31+30+31+31+30+31+30+day
    end
    day_of_year=N
    if  year==2020  | year==2024 then 
        
        day_of_year=N+1
        gama=2*%pi/366*(day_of_year-1+(hour-12)/24)
     else
          day_of_year=N
        gama=2*%pi/365*(day_of_year-1+(hour-12)/24)
    end 
    disp (day_of_year)
    
   //eqtime=229.18*(0.000075 + 0.001868*cos(gama) – 0.032077*sin(gama) – 0.014615*cos(2*gama)-0.040849*sin(2*gamaγ) )
    
    eqtime=229.18*(.000075+.001868*cos(gama)-.032077*sin(gama)-.014615*cos(2*gama)-.040849*sin(2*gama))
    //======================================================================================
    
   // decl = 0.006918 – 0.399912cos(γ) + 0.070257sin(γ) – 0.006758cos(2γ) + 0.000907sin(2γ)
   
//– 0.002697cos(3γ) + 0.00148sin (3γ)

decl=.006918-.399912*cos(gama)+.070257*sin(gama)-.006758*cos(2*gama)-.00907*sin(2*gama)-.002697*cos(3*gama)+.00148*sin(3*gama)
    disp (eqtime)
    disp decl
    disp (decl)
    
    time_offset=eqtime+4*logitude-60*time_zone
    disp (time_offset)
    tst=hr*60+mn+sc/60+time_offset
    ha=tst/4-180
    disp (tst)
    disp (ha)
    ha_rad=ha*%pi/180
    // Fi  zenith
    latit_rad=latit*%pi/180
    //cos() = sin(lat)sin(decl) + cost(lat)cos(decl)cos(ha)
    
    cos_zenith=sin(latit_rad)*sin(decl)+cos(latit_rad)*cos(decl)*cos(ha_rad)  
    zenith=acos(cos_zenith)*180/%pi
   // disp zenith
   // disp(zenith)
   Fi=acos(cos_zenith)
   // eiev=90-zenith
    cos_180_azum=-((sin(latit_rad)*cos(Fi)-sin(decl))/cos(latit_rad)/sin(Fi))
    disp (cos_180_azum)
    azumith=acos(cos_180_azum)*180/%pi
    if hour > 12  then
        azumith=360-azumith
    end
 
   
 
    
    
    
    
    //===========================================================================
    
    
    
      azm_d=azumith
        
              
    elev_d=90-zenith   
                  
    elev_rad=elev_d*%pi/180
      
    delta_azm=0
    azm_c=azm_d
    azm_rad=(90-azm_c)*%pi/180
    zk=sin(elev_rad)          // tan(elev_rad)
    xi=cos(elev_rad)*cos(azm_rad)
    yj=cos (elev_rad)*sin(azm_rad)
    l=sqrt(zk^2+xi^2+yj^2)
    k_sol=zk/l
    i_sol=xi/l
    j_sol=yj/l
    
    i_refl=0
   
  k_refl=sqrt(1-j_refl^2-i_refl^2)
  
   cost=(i_sol*i_refl+j_sol*j_refl+k_sol*k_refl)
  
   
  
   L_col=sqrt(ix_col^2+jy_col^2+kz_col^2)
  
   i_col=ix_col/L_col
   j_col=jy_col/L_col
   k_col=kz_col/L_col
  
   cost=(i_refl*i_col+j_refl*j_col+k_refl*k_col)
   
   a=.5
   b=.5
   f=zeros(7,1)
   f(1)= -cost+(i_sol*i_refl+j_sol*j_refl+k_sol*k_refl)
      f(2)=-cost+(i_refl*ix_col/L_col+j_refl*jy_col/L_col+k_refl*kz_col/L_col)
       f(3)=-1+j_refl^2+k_refl^2+i_refl^2
      f(4)=-L_col+sqrt(ix_col^2+jy_col^2+kz_col^2)
      f(5)=i_sol+a*j_sol+b*k_sol
      f(6)=i_refl+a*j_refl+b*k_refl
      f(7)=ix_col/L_col+a*jy_col/L_col+b*kz_col/L_col
    
      
       
      iter=1
     while max(abs(f))>1e-6
    
      j=zeros(7,7)
       
      j(1,1)=-1
      j(1,2)=j_sol
      j(1,3)=k_sol
     // j(1,8)=i_sol
     //  f(2)=-cost+(i_refl*ix_col/L_col+j_refl*jy_col/L_col+k_refl*kz_col/L_col)
      j(2,1)=-1
      j(2,4)=i_refl/L_col
      j(2,2)=jy_col/L_col
      j(2,5)=-i_refl*ix_col/L_col^2-j_refl*jy_col/L_col^2-k_refl*kz_col/L_col^2
      j(2,3)=kz_col/L_col
      //j(2,8)=ix_col/L_col
    //  f(3)=-1+j_refl^2+k_refl^2+i_refl^2
      j(3,2)=j_refl*2
      j(3,3)=k_refl*2
    // f(4)=-L_col+sqrt(ix_col^2+jy_col^2+kz_col^2)
      j(4,5)=-1
      j(4,4)=1/sqrt(ix_col^2+jy_col^2+kz_col^2)*2*ix_col/2
     
     //    f(5)=i_sol+a*j_sol+b*k_sol
     j(5,6)=j_sol
     j(5,7)=k_sol
   //  f(6)=i_refl+a*j_refl+b*k_refl
     j(6,6)=j_refl
     j(6,2)=a
     j(6,7)=k_refl
     j(6,3)=b
     
     
   j(7,5)=-ix_col/L_col^2-a*jy_col/L_col^2-b*kz_col/L_col^2
   j(7,4)=1/L_col
   j(7,6)=jy_col/L_col
   j(6,2)=a/L_col
   j(7,7)=kz_col/L_col
   
 delta=   j\f                         //ja*f
 
 nf=10
  //     1- cost          2-j_refl     3-k_ref      4  p
      //      5 L_col       6-a      7-b        8  c
      cost=cost-delta(1)/nf
    j_refl=j_refl-delta(2)/nf
    k_refl=k_refl-delta(3)/nf
    ix_col=ix_col-delta(4)/nf
    L_col=L_col-delta(5)/nf
    a=a-delta(6)/nf
    b=b-delta(7)/nf
   // a=a-delta(6)/nf
  
   f=zeros(7,1)
   f(1)= -cost+(i_sol*i_refl+j_sol*j_refl+k_sol*k_refl)
      f(2)=-cost+(i_refl*ix_col/L_col+j_refl*jy_col/L_col+k_refl*kz_col/L_col)
       f(3)=-1+j_refl^2+k_refl^2+i_refl^2
      f(4)=-L_col+sqrt(ix_col^2+jy_col^2+kz_col^2)
      f(5)=i_sol+a*j_sol+b*k_sol
      f(6)=i_refl+a*j_refl+b*k_refl
      f(7)=ix_col/L_col+a*jy_col/L_col+b*kz_col/L_col
    
      iter=iter+1
  end
  
A1=[j_sol k_sol
  j_refl k_refl]
  B1=[-i_sol
  -i_refl]
  R1=A1\B1
  
  A2=[j_refl k_refl
      jy_col/L_col  kz_col/L_col]
      B2=[-i_refl 
           -ix_col/L_col]
           R2=A2\B2
  f_max=max(abs(f))
  disp 'convergence'
  disp (f_max)
 
  
  
  
  acost=acos(cost)*180/%pi
  disp acost
  disp (acost)
  cost2=i_sol*ix_col/L_col+j_sol*jy_col/L_col+k_sol*kz_col/L_col
  acost2=acos(cost2)*180/%pi
  
 
endfunction

