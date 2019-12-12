pro clickjmap,jmap,srtotal,stepsize,npixel,cdelt1,wave=wave,slit_label=slit_label,slit_num=slit_num,spfile=spfile

  a=1
  
  distsun=151816767348.71*0.001
  arcsec=srtotal
  distkm=distsun*tan(arcsec/3600.*!pi/180)
  i=0
  x=fltarr(2)
  y=x
  xp=x & yp=x
  
  wid=!D.Window                                       ; id = "current window"
  WSet, wid
  xsize = float(!D.X_VSize)                           ; X and Y size of the current visible area
  ysize = float(!D.Y_VSize)                           ;
  Window, /Pixmap, /Free, XSize=xsize, YSize=ysize    ; Set a window for displaying.
  pixID = !D.Window
  Device, Copy=[0, 0, xsize, ysize, 0, 0, wid]        ; Copy window contents
  WSet, wid                                           ; Work in the new window, Here we can draw
  wshow,wid                                           ; on the figure without changing the data.
  
  
  
  while a eq 1 do begin
  
  
    cursor,xvar,yvar,/change                        ; Get any cursor movement
    !mouse.button=0                                   ; Clear mouse button flag
    cursor,xvar,yvar,/change                      ; Check any cursor movement
    if !mouse.button EQ 1 THEN BEGIN                  ; If the left mouse button was pressed,
      i++
      xjmap = xvar                                 ; Concatenate X00 array
      yjmap = yvar                                 ; Concatenate Y00 array
      acc_xp = exist(acc_xp) ? [acc_xp,xjmap] : xjmap
      acc_yp = exist(acc_yp) ? [acc_yp,yjmap] : yjmap
      xp[i-1]=xjmap
      yp[i-1]=yjmap
      x[i-1]=ceil(xjmap)*stepsize/60.*1/npixel
      print,'x:',x[i-1]
      sz=size(jmap)
      y[i-1]=ceil(yjmap)*srtotal/sz[2]
      print,'y:',y[i-1]
      wait,0.1
    endif
    
    if exist(acc_xp) && !mouse.button EQ 0 then begin
      acc_xp_tmp = [acc_xp,xvar]
      acc_yp_tmp = [acc_yp,yvar]
    endif else begin
      acc_xp_tmp = xvar
      acc_yp_tmp = yvar
    endelse
    Device, Copy=[ 0, 0, xsize, ysize, 0, 0, pixID]   ; Clear the line drawn.
    plots,acc_xp_tmp,acc_yp_tmp,PSYM=6,SYMSIZE=0.5,color=255
    if exist(ind_xjmap)then begin
    for k=0,n_elements(ind_xjmap)-1 do xyouts,ind_xjmap[k],ind_yjmap[k]+3,STRTRIM(STRING(ind_numb[k]),1),color=255
    endif
    if i eq 2 then begin
      speed=distkm*(y[1]-y[0])/((x[1]-x[0])*60.)/float(srtotal)
      print, 'Speed:',speed
      i=0
      for j=0,n_elements(acc_xp)-2,2 do plots,acc_xp[j:j+1],acc_yp[j:j+1],color=255
      if not exist(speedfile) then begin
        speedfile=float([x,y,speed])       ; Format: X0  X1  Y0  Y1  Speed
      endif else begin
        speedfile=[[speedfile],[x,y,speed]]
      endelse
     
     
      xyouts,xjmap,yjmap+3,STRTRIM(STRING(n_elements(acc_xp)/2),1),color=255
      ind_xjmap= exist(ind_xjmap)?[ind_xjmap,xjmap]:xjmap
      ind_yjmap= exist(ind_yjmap)?[ind_yjmap,yjmap]:yjmap
      ind_numb= exist(ind_numb)?[ind_numb,n_elements(acc_xp)/2]:(n_elements(acc_xp)/2)
      
    endif else begin
      for j=0,n_elements(acc_xp_tmp)-2,2 do plots,acc_xp_tmp[j:j+1],acc_yp_tmp[j:j+1],color=255
    endelse
    
    if !mouse.button eq 4 then a=2
  endwhile
  
  if exist(acc_xp) then begin
    if i ne 0 then begin
      acc_xp=[acc_xp,xvar] & acc_yp=[acc_yp,yvar]
      xj=fltarr(2) & yj=fltarr(2) & x=fltarr(2) & y=fltarr(2)
      xj[1]=acc_xp[n_elements(acc_xp)-1] & xj[0]=acc_xp[n_elements(acc_xp)-2]
      yj[1]=acc_yp[n_elements(acc_yp)-1] & yj[0]=acc_yp[n_elements(acc_yp)-2]
      x[1]=ceil(xj[1])*stepsize/60.*1/npixel
      y[1]=ceil(yj[1])*srtotal/sz[2]
      x[0]=ceil(xj[0])*stepsize/60.*1/npixel
      y[0]=ceil(yj[0])*srtotal/sz[2]
      print,'x:',x[1]
      print,'y:',y[1]
      speed=distkm*(y[1]-y[0])/((x[1]-x[0])*60.)/float(srtotal)
      print, 'Speed:',speed
      speedfile = exist(speedfile) ? [[speedfile],[x,y,speed]] : float([x,y,speed])
    endif
    Device, Copy=[ 0, 0, xsize, ysize, 0, 0, pixID]
    plots,acc_xp,acc_yp,PSYM=6,SYMSIZE=0.5,color=255
    for j=0,n_elements(acc_xp)-2,2 do plots,acc_xp[j:j+1],acc_yp[j:j+1],color=255
    if exist(ind_xjmap)then begin
    for k=0,n_elements(ind_xjmap)-1 do xyouts,ind_xjmap[k],ind_yjmap[k]+3,STRTRIM(STRING(ind_numb[k]),1),color=255
    endif
  endif

if exist(speedfile) then begin   
while 1 do begin 
    delcurve=0
    read,delcurve, prompt='Delete the line number: '

    if delcurve ne 0 then begin
      aux=speedfile
      if (delcurve ne (size(aux))[2]) && (delcurve ne 1) then begin
      speedfile=[[aux[*,0:delcurve-2]],[aux[*,delcurve:*]]]
      aux=acc_xp & acc_xp=[aux[0:delcurve*2-2-1],aux[delcurve*2:*]]
      aux=ind_xjmap & ind_xjmap=[aux[0:delcurve-2],aux[delcurve:*]]
      aux=ind_yjmap & ind_yjmap=[aux[0:delcurve-2],aux[delcurve:*]]
      endif else begin
      if delcurve ne 1 then begin
      speedfile=aux[*,0:delcurve-2]
      aux=acc_xp & acc_xp=aux[0:delcurve*2-2-1]
      aux=ind_xjmap & ind_xjmap=aux[0:delcurve-2]
      aux=ind_yjmap & ind_yjmap=aux[0:delcurve-2]
      endif else begin
      speedfile=aux[*,delcurve:*]
      aux=acc_xp & acc_xp=aux[delcurve*2:*]
      aux=ind_xjmap & ind_xjmap=aux[delcurve:*]
      aux=ind_yjmap & ind_yjmap=aux[delcurve:*]
      endelse
      endelse      
      ind_numb=indgen(n_elements(ind_xjmap))+1
      
      Device, Copy=[ 0, 0, xsize, ysize, 0, 0, pixID]
      plots,acc_xp,acc_yp,PSYM=6,SYMSIZE=0.5,color=255
      for j=0,n_elements(acc_xp)-2,2 do plots,acc_xp[j:j+1],acc_yp[j:j+1],color=255
      if exist(ind_xjmap)then begin
       ; for k=0,n_elements(ind_xjmap)-1 do xyouts,ind_xjmap[k],ind_yjmap[k]+3,STRTRIM(STRING(ind_numb[k]),1),color=255
       for k=0,n_elements(acc_xp)-2,2 do xyouts,max(acc_xp[k:k+1]),max(acc_yp[k:k+1])+3,STRTRIM(STRING(ind_numb[k/2]),1),color=255
      endif
      wait,0.8
    endif else begin
      goto,savefile
    endelse
endwhile
endif
  
  
  
  
  savefile:
  if exist(speedfile) && exist(wave) then begin
  
    save,speedfile,FILENAME = 'speeddata' + wave + slit_label + slit_num + '.sav'
    fname='speeddata' + wave + slit_label + slit_num + '.dat'
    openw,1,fname
    printf,1,speedfile,FORMAT='(F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F8.2)'
    close,1
    spfile=speedfile
  endif
  
end
