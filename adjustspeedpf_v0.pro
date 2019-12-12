pro adjustspeedpf_v0,jmap,srtotal,stepsize,x0,x1,y0,y1,min=min,arcsec=arcsec,speed=speed,img=img,cdelt1=cdelt1,degree=degree
;pro adjustspeedpf_v0
 ;restore,'C:\sdodata\20110706\171\DeroteData\jmap.sav'
;restore, '/Volumes/Disk1/WorkingOn/Tardelli/speed_automatic/adjustspeed.sav',/verbose
;degree = 2

undefine, xarr
undefine, yarr

  k=[1,1,1]
  jmap2=convol(jmap,k,/norm,/edge_tr)
  for t=0,100 do jmap2=convol(jmap2,k,/norm,/edge_tr)
  jmap_new= smooth(jmap-jmap2,3)                                          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;window,0 & plot_image,hist_equal(jmap_new,per=2),back=255,col=0
  ;cursor,x0,y0 & wait,0.5 & plots,x0,y0,psym=2,color=0
  ;cursor,x1,y1 & plots,x1,y1,psym=2,color=0; & degree=2
  ;x0=6 & x1=20 & y0=3 &  y1=17
  
  ; print,x0,x1,y0,y1
  ; if x0,x1,y0 and y1 are in pixels don't use /min and /arcsec
  
  if keyword_set(arcsec) then begin ; if the y0 and y1 are in arcsec
    y0=y0*(size(jmap,/dim))[1]/srtotal
    y1=y1*(size(jmap,/dim))[1]/srtotal
  endif
  
  if keyword_set(min) then begin ; if the x0 and x1 are in minutes
    x0=x0*60./stepsize
    x1=x1*60./stepsize
  endif

  deslpix=6   ; Desloc. pixels
  
  corrx=x0<deslpix
  corry=x0<deslpix
  ;speedcut=jmap[(x0-deslpix)>0:(x1+deslpix)<(size(jmap,/dim))[0]-1,(y0-deslpix)>0:(y1+deslpix)<(size(jmap,/dim))[1]-1]
  speedcut=jmap_new[(x0-deslpix)>0:(x1+deslpix)<(size(jmap,/dim))[0]-1,(y0-deslpix)>0:(y1+deslpix)<(size(jmap,/dim))[1]-1]
  
  x0in=x0 & x1in=x1
  y0in=y0 & y1in=y1
  
  x1=x1-x0+corrx
  x0=0+corrx
  y1=y1-y0+corry
  y0=0+corry
  
  a=float(y1-y0)/float(x1-x0)
  b=y0-a*x0
  
  x=findgen(50)/49.*(x1-x0)+x0
  y=a*x+b
  
  speedcut1=speedcut
  for j=0,(size(speedcut,/dim))[1]-1 do begin
    for i=0,(size(speedcut,/dim))[0]-1 do begin
    
      ;if  (j gt i*a+b+deslpix) || (j lt i*a+b-deslpix) then begin
      if  (j gt i*a+b+deslpix*1.5) || (j lt i*a+b-deslpix*1.5) then begin   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        speedcut1[i,j]=min(speedcut1)
      endif else begin
        xarr = exist(xarr) ? [xarr,i]:i
        yarr = exist(yarr) ? [yarr,j]:j
      endelse
      
    endfor
  endfor
  
  plot_image,jmap,back=255,col=0
  oplot,xarr+x0in-corrx,yarr+y0in-corry,psym=3
  
  speedcut1=speedcut1>(-0.06)<0.06
;;;  merror=1./(speedcut1(xarr,yarr)+0.0001)
  merror=1./(speedcut1(xarr,yarr)>0.000001<100.)   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  merror = merror-min(merror)+1                                     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  
  if not exist(degree) then degree = 1    
  if degree eq 1 then begin
    y0n=y0 & y1n=y1
    coefs=poly_fit(xarr,yarr,1,MEASURE_ERRORS=sqrt(merror))   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    an=coefs[1]
    bn=coefs[0]
    x0n=float(y0n-bn)/float(an)
    x1n=float(y1n-bn)/float(an)
    xn=findgen(50)/49.*(x1n-x0n)+x0n
    yn=an*xn+bn
    ;window,6
    ;plot,x,y,background=255,color=0
    ;oplot,xn,yn,color=0
  endif
  
  
  if degree eq 2 then begin
    coefs=poly_fit(xarr,yarr,2,MEASURE_ERRORS=merror,yfit=yfit,yband=yband, sigma=sigma)  ;;;;;;;;;;;;;;;;
          ;result = COMFIT(Xarr, Yarr, coefs, /LOGSQUARE, weight=sqrt(merror))
    xx=findgen(1000)*(max(x)-min(x))/999.+min(x)
    y0n=min(xx) & y1n=max(xx)
    x0n=min(x) & x1n=max(x)
    window,1,xs=300,ys=300 & plot_image,speedcut,back=255,color=0
    oplot,xarr,yarr,psym=3
    oplot, xx,coefs[2]*xx^2.+coefs[1]*xx+coefs[0],color=0, thick=3   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;stop
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
    bb = sqrt(n_elements(xarr))
    oplot, xx, (coefs[2]+sigma[2]/bb)*xx^2. + (coefs[1]+sigma[1]/bb)*xx + (coefs[0]+sigma[0]/bb),color=80, thick=2,linestyle=2
    oplot, xx, (coefs[2]-sigma[2]/bb)*xx^2. + (coefs[1]-sigma[1]/bb)*xx + (coefs[0]-sigma[0]/bb),color=80, thick=2,linestyle=2
    
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   img=tvrd()
  
    
  endif
  
  ;srtotal= 35.838534 ;arcsec
  ;stepsize= 12. ;s

  if exist(srtotal) then begin
  if degree eq 1 then begin
    distsun=151816767348.71*0.001
    arcsec=srtotal
    distkm=distsun*tan(arcsec/3600.*!pi/180.)
    
    dist=float(y1n-y0n)*distkm/(size(jmap,/dim))[1]
    tim=float(x1n-x0n)*stepsize
    speed=float(dist)/float(tim)
    speed=[x0n+x0in-corrx,x1n+x0in-corrx,y0n+y0in-corry,y1n+y0in-corry,speed]
    
  endif else begin
    speed=[x0n,x1n,y0n,y1n,x0in-corrx,y0in-corry,coefs[2],coefs[1],coefs[0]]
    distsun=151816767348.71*0.001
    arcsec=cdelt1
    distkmpx=distsun*tan(arcsec/3600.*!pi/180.)
    
    A = 2./(stepsize^2.)*distkmpx*coefs[2]
    Ea = 2./(stepsize^2.)*distkmpx*sigma[2]
    B = 1./stepsize*distkmpx*coefs[1]+2./(stepsize^2.)*distkmpx*coefs[2]*x0n*12.
    Eb = 1./stepsize*distkmpx*sigma[1]+2./(stepsize^2.)*distkmpx*sigma[2]*x0n
    print,'Speed 2nd: V(km/s) = (' + STRTRIM(string(A,FORMAT='(f6.2)'),1) + ' +- ' + STRTRIM(string(Ea,FORMAT='(f6.2)'),1)  + ')*t(s) + (' + STRTRIM(string(B,FORMAT='(f6.2)'),1) + ' +- ' + STRTRIM(string(Eb,FORMAT='(f6.2)'),1) + ')'
    
    ;print,'Speed 2nd: V = (' + STRTRIM(string(2*(coefs[2])),1) + ' +- ' + STRTRIM(string(2*sigma[2]/bb),1)  + ')*x + (' + STRTRIM(string(coefs[1]),1) + ' +- ' + STRTRIM(string(sigma[1]/bb),1) + ')'
  endelse
  endif
  
  ;stop
 

end
