;+
; Project     : HT map for many curves with any shape
;
; Name        : HTmap
;
; Purpose     : Height-time map of intensity flutuation for the EUVI/STEREO and AIA/SDO in different wavelenths
;
;
; Syntax      : HTmap
;
; Keyword Parameters     :
;    path - Path of the file list
;    npixel - Width of the curve in pixel
;    set_jmap_type - If '0' plot the npixel HTmap, '1' plot the average of npixel, '2' both
;    slit_label - Label of the series of slits
;    slit_num - Number of each slits in a slit label
;    save_slit - Name of the label of the slit if you want to save it.
;    stepsize - Sample time of analyzed data
;    starttime - Start time
;    endtime - End time
;    running - Set running difference
;    wave - Not necessary any more
;    log - To plot in log scale
;    
; Called routines        :
;    dist_along_curve.pro'
;    timeaxisjmap.pro'
;    clickjmap.pro'
;    adjustspeedpf_v0.pro'
;    equal_dist_ensanchar.pro'
;    ensanchar.pro'
; 
; Calling sequence       :
; 
; HTmap, path=path,npixel=npixel, set_jmap_type=set_jmap_type, $
;    slit_label=slit_label, slit_num=slit_num, save_slit=save_slit, stepsize=stepsize, starttime=starttime,
;    endtime=endtime, running=running, wave=wave, log=log
;    
; Calling Example: 
; HTmap,path='C:\sdodata\20110706\171\DeroteData',npixel=3,save_slit='slittest',stepsize=12
; HTmap,path='C:\sdodata\20110706\171\DeroteData',npixel=3,slit_label='slitb',stepsize=12,slit_num='4',/log
;
;-


;---------------------------------------------------------------
; Create a J-map from a sequence of FITS
;---------------------------------------------------------------


;---------------------------------------------------------------
;       Graphical Interface and definition of slits
;---------------------------------------------------------------
; Open the first image and select region to zoom in.
; - LEFT cursor to move the selected region;
; - MIDDLE cursor to resize region;
; - RIGHT to select.
;
; After zoom in, select the spline curve.
; - LEFT: select points
; - RIGHT: end point
;
; To select more than 1 slit: press spacebar in the command line
; To finish selecting slits: press 's' in the command line
;---------------------------------------------------------------

function graph_intface

  CURSOR, X00, Y00, /DOWN,/DATA                       ; Get the (x,y) coordinates of the first click
  wid=!D.Window                                       ; id = "current window"
  WSet, wid
  xsize = float(!D.X_VSize)                           ; X and Y size of the current visible area
  ysize = float(!D.Y_VSize)                           ;
  Window, /Pixmap, /Free, XSize=xsize, YSize=ysize    ; Set a window for displaying.
  pixID = !D.Window
  Device, Copy=[0, 0, xsize, ysize, 0, 0, wid]        ; Copy window contents
  WSet, wid                                           ; Work in the new window, Here we can draw
  wshow,wid                                           ; on the figure without changing the data.
  plots,X00,Y00,PSYM=6,SYMSIZE=0.5
  
  REPEAT BEGIN                                        ; Repeat until button unpressed
    cursor,px,py,/change,/data                        ; Get any cursor movement
    !mouse.button=0                                   ; Clear mouse button flag
    cursor,pxx,pyy,/change,/data                      ; Check any cursor movement
    X = [X00,pxx]                                     ;
    Y = [Y00,pyy]                                     ;
    if !mouse.button EQ 1 THEN BEGIN                  ; If the left mouse button was pressed,
      X00 = [X00,pxx]                                 ; Concatenate X00 array
      Y00 = [Y00,pyy]                                 ; Concatenate Y00 array
    ENDIF
    if n_elements(X) EQ 2 THEN BEGIN                  ; If n_elements == 2, use linear path
      T0=X                                            ;
      Z0=Y                                            ;
    ENDIF ELSE BEGIN
      SPLINE_P,X,Y,T0,Z0,INTERVAL=0.4                 ; Parametric Cubic Spline Interpolation allow nonmonotonic curve.
    ENDELSE                                           ; The interval was decreased to 0.4 to eliminate the gaps.
    Device, Copy=[ 0, 0, xsize, ysize, 0, 0, pixID]   ; Clear the line drawn.
    plots,T0,Z0,color=255                             ; Draw the new curve.
    for j=0,n_elements(X00)-1 DO BEGIN
      plots,X00[j],Y00[j],PSYM=6,SYMSIZE=0.5,color=255  ; scatter plot for the clicked points
    ENDFOR
  endrep until !mouse.button eq 4                     ; Stop when right mouse button pressed
  
  X00 = [X00,pxx]                                     ; Concatenate X00 array
  Y00 = [Y00,pyy]                                     ; Concatenate Y00 array
  for j=0,n_elements(X00)-1 DO BEGIN
    plots,X00[j],Y00[j],PSYM=6,SYMSIZE=0.2,color=255  ; scatter plot for all the clicked points
  ENDFOR
  slit_pos = {vecX:T0 , vecY:Z0, xclick:x00, yclick:y00}                      ; Structure of the X and Y position vector
  return, slit_pos                                    ; Return Structure of vectors
end


;---------------------------------------------------------------
;       FILTER to eliminate gaps and wrong pixel position
;---------------------------------------------------------------
; - The SPLINE_P interval was decreased to eliminate the gaps.
; - 1st loop to eliminate the overlapping pixels.
; - 2nd loop to satisfy the y(x) or x(y) conditions.
; Result: the correct curve
;---------------------------------------------------------------

function filter_slit,slit_vecX,slit_vecY
  T0=round(slit_vecX)                                    ; Round all the T0
  Z0=round(slit_vecY)                                    ; Round all the Z0
  T00=T0[0]
  Z00=Z0[0]
  i=0
  
  for j=0,n_elements(T0)-3 DO BEGIN                      ;
    dx=T0[j+1]-T0[j]                                     ; Create a T00 and Z00 only with
    dy=Z0[j+1]-Z0[j]                                     ; with non-repeated values
    if (dx || dy) NE 0 THEN BEGIN                        ;
      T00=[T00,T00[i]+dx]                                ;
      Z00=[Z00,Z00[i]+dy]                                ; If T0 or Z0 change the value
      i++                                                ; T00 <- T0 and Z00 <- Z0
    endif                                                ;
  endfor                                                 ;
  
  T0=T00                                                 ;
  Z0=Z00                                                 ;
  T00=T0[0]                                              ;
  Z00=Z0[0]                                              ;
  i=0
  
  for j=0,n_elements(T0)-3 DO BEGIN                      ; Filter to satisfy the y(x) or x(y) conditions.
    dx=T0[j+1]-T0[j]                                     ;
    dy=Z0[j+1]-Z0[j]                                     ;
    dx2=T0[j+2]-T0[j]                                    ; If dx2 and dy2 change both 1 pixel
    dy2=Z0[j+2]-Z0[j]                                    ; there are no horizontal or vertical line,
    if (abs(dx2) EQ 1) && (abs(dy2) EQ 1) THEN BEGIN     ; only diagonal.
      T00=[T00,T00[i]+dx2]                               ;
      Z00=[Z00,Z00[i]+dy2]                               ; Replace this     By this
      j++                                                ;     █               █
    endif else begin                                     ;     █               █
      T00=[T00,T00[i]+dx]                                ;     █               █
      Z00=[Z00,Z00[i]+dy]                                ;    ██              █
    endelse                                              ;    █               █
    i++                                                  ;    █               █
  endfor                                                 ;
  out_filter = {vecX:T00 , vecY:Z00}                     ; Structure of the X and Y position vector
  return,out_filter                                      ; Return Structure of vectors
end


;---------------------------------------------------------------
;                 Enlarge the slit to npixel
;---------------------------------------------------------------
;
; Compare two pixels and determine the new pixel position
; T1 and Z1 Columns:
; 0 - middle curve (guide)
; 1 - 1st external curve
; 2 - 1st internal curve
; 3 - 2nd external curve
; 4 - 2nd internal curve
; ...
;---------------------------------------------------------------

function width_slit,slit_T0,slit_Z0,slit_T1,slit_Z1,npixel
  T0=slit_T0
  Z0=slit_Z0
  T1=slit_T1
  Z1=slit_Z1
  if npixel ne 1 THEN BEGIN                              ; For more than 1 pixel
    FOR j=0,(n_elements(T0)-3) DO BEGIN
      dx=T0[j+1]-T0[j]
      dy=Z0[j+1]-Z0[j]
      if dy EQ 0 THEN BEGIN;
        if abs(dx) EQ 1 THEN BEGIN                       ; If dx=1 and dy=0
          for i=1,npixel-2,2 DO BEGIN                    ; For all the npixels
            T1[j,i]=T0[j+1]
            Z1[j,i]=Z0[j+1]-dx*(i+1)/2
            T1[j,i+1]=T0[j+1]
            Z1[j,i+1]=Z0[j+1]+dx*(i+1)/2
          ENDFOR
        ENDIF
      endif;
      if dx EQ 0 THEN BEGIN
        if abs(dy) EQ 1 THEN BEGIN                        ; If dy=1 and dx=0
          for i=1,npixel-2,2 DO BEGIN                     ; For all the npixels
            T1[j,i]=T0[j+1]+dy*(i+1)/2
            Z1[j,i]=Z0[j+1]
            T1[j,i+1]=T0[j+1]-dy*(i+1)/2
            Z1[j,i+1]=Z0[j+1]
          ENDFOR
        ENDIF
      ENDIF
      if (abs(dx) EQ 1) && (abs(dy) EQ 1) THEN BEGIN      ; If dx=1 and dy=1
        dx2=T0[j+2]-T0[j+1]
        dy2=Z0[j+2]-Z0[j+1]
        if (abs(dy2) EQ 1) && (dx2 EQ 0) THEN BEGIN       ; If dy2=1 and dx2=0
          for i=1,npixel-2,2 DO BEGIN                     ; For all the npixels
            T1[j,i]=T0[j+1]+dy2*(i+1)/2
            Z1[j,i]=Z0[j+1]
            T1[j,i+1]=T0[j+1]-dy2*(i+1)/2
            Z1[j,i+1]=Z0[j+1]
          ENDFOR
        endif
        if (dy2 EQ 0) && (abs(dx2) EQ 1) THEN BEGIN       ; If dy2=0 and dx2=1
          for i=1,npixel-2,2 DO BEGIN                     ; For all the npixels
            T1[j,i]=T0[j+1]
            Z1[j,i]=Z0[j+1]-dx2*(i+1)/2
            T1[j,i+1]=T0[j+1]
            Z1[j,i+1]=Z0[j+1]+dx2*(i+1)/2
          ENDFOR
        endif
        if (abs(dy2) EQ 1) && (abs(dx2) EQ 1) THEN BEGIN  ; If dy2=1 and dx2=1
          for i=1,npixel-2,2 DO BEGIN                     ; For all the npixels
            T1[j,i]=T0[j+1]
            Z1[j,i]=Z0[j+1]-dx2*(i+1)/2
            T1[j,i+1]=T0[j+1]
            Z1[j,i+1]=Z0[j+1]+dx2*(i+1)/2
          ENDFOR
        endif
      ENDIF
    ENDFOR
  ENDIF
  out_width = {arrX:T1 , arrY:Z1}                         ; Structure of the X and Y position arrays
  return,out_width                                        ; Return Structure of arrays
end


;---------------------------------------------------------------
;              Extract the pixels from the Data_cube
;---------------------------------------------------------------
; Make the pixels slit based in the pixels of the discretized curve.
; The number of columns will be: files*npixel
;---------------------------------------------------------------

function extract_slit,slit_vecX,slit_vecY,data_cube,jmap,slit0,n_files,npoint,npixel
  T1=slit_vecX
  Z1=slit_vecY
  FOR j=0,(n_files-1) DO BEGIN                                        ; Make a npixel slit with the pixels of the
    FOR i=0,npoint-1 DO BEGIN                                         ; discretized curve.
      slit0[i,(npixel-1)/2]=data_cube[T1[i,0]:T1[i,0],$               ; Slit0 is the vector of pixels
        Z1[i,0]:Z1[i,0],j]
      FOR k=1,npixel-2,2 DO BEGIN
        slit0[i,(npixel+k)/2]=data_cube[T1[i,k]:T1[i,k],$
          Z1[i,k]:Z1[i,k],j]
        slit0[i,(k-1)/2]=data_cube[T1[i,(npixel-k)]:T1[i,(npixel-k)],$
          Z1[i,(npixel-k)]:Z1[i,(npixel-k)],j]
      ENDFOR
    ENDFOR
    for i=0,npixel-1 DO BEGIN                                          ;
      jmap[npixel*j+i,0:(npoint-1)]=slit0[*,i]                         ; Create the jmap from the slits of all files.
    ENDFOR
  ENDFOR
  ;jmap=jmap[*,0:npoint-3]                                             ; Eliminate the extra pixel of the guide curve
  jmap_avg=make_array(n_files,npoint,/FLOAT)
  
  for m=0,n_files-1 do begin
    for n=0,npoint-1 do begin
      jmap_avg[m,n]=mean(jmap[m*npixel:m*npixel+(npixel-1),n])
    endfor
  endfor
  out_extract_slit = {jmap:jmap , jmap_avg:jmap_avg}
  return,out_extract_slit
end


;---------------------------------------------------------------
;---------------------------------------------------------------
;                        Main Procedure
;---------------------------------------------------------------
;---------------------------------------------------------------


pro HTmap,path=path,npixel=npixel,set_jmap_type=set_jmap_type, $
    slit_label=slit_label,slit_num=slit_num,save_slit=save_slit,stepsize=stepsize,starttime=starttime,endtime=endtime,running=running, $
    wave=wave,log=log
    
  compile_opt strictarr
;    .compile 'C:\codes\dist_along_curve.pro'
;    .compile 'C:\codes\timeaxisjmap.pro'
;    .compile 'C:\codes\clickjmap.pro'
;    .compile 'C:\codes\adjustspeedpf_v0.pro'
;    .compile 'C:\codes\equal_dist_ensanchar.pro'
;    .compile 'C:\codes\ensanchar.pro'
  
  
  
  
  
  ;---------------------------------------------------------------
  ; Configuration Set
  ;---------------------------------------------------------------
  
  ; Define the width of the slit
  if not exist(npixel) then npixel=7                                            ; number of pixels
  
  ; J-map type                                        ; If npixel_type = '0' plot npixel width for image
  if not exist(set_jmap_type) then set_jmap_type=2    ; If npixel_type = '1' plot the mean of npixel width
  ; If npixel_type = '2' plot both
  
  ; Set the path
  if not exist(path) then path='C:\sdodata\20120831\wavelets' ; FITS file directory
  cd, path
  
  filelist=file_list('.','*ts')                   ; Read the file list
  filelist=filelist
  
  
  ;===========================================================
  nfiles=n_elements(filelist)
  header_time=strarr(nfiles)
  time_sec=fltarr(nfiles)
  
  for i=0,nfiles-1 do begin
    header = HEADFITS(filelist[i])
    header_time[i] = SXPAR(header,'DATE-OBS')
    if header_time[i] eq 0 then goto, diffheader
  endfor
  diffheader:
  for i=0,nfiles-1 do begin
    header = HEADFITS(filelist[i])
    header_time[i] = SXPAR(header,'DATE_OBS')
  endfor
  
  hour=long(STRMID(header_time,11,2))
  min=long(STRMID(header_time,14,2))
  sec=long(round(float(STRMID(header_time,17,5))))
  time_sec=long(sec+min*60.+hour*60.*60.)
  
  if not exist(startTime) then begin
    startfile=0 & endfile=nfiles-1
    goto,notimeinput
  endif
  
  ;-----------------------------
  startTimeSec=float(STRMID(startTime,0,2))*60.*60.+float(STRMID(startTime,2,2))*60.+float(STRMID(startTime,4,2))
  
  endTimeSec=float(STRMID(endTime,0,2))*60.*60.+float(STRMID(endTime,2,2))*60.+float(STRMID(endTime,4,2))
  
  startFile=where(time_sec ge startTimeSec)
  startFile=startFile[0]
  endFile=where(time_sec le endTimeSec)
  endFile=endFile[n_elements(endFile)-1]
  
  filelist=filelist[startfile:endfile]
  header_time=header_time[startfile:endfile]
  ;=========================================================================
  
  
  
  
  notimeinput:
  
  
  
  if exist(slit_label) then goto,slitdone
  
  !P.MULTI = [0, 1, 1]
  fits2map,filelist[191],map,header=header              ; Read the first FITS file as a map
  ;fits2map,filelist[422],map,header=header
  
  if keyword_set(log) then map.data = alog(map.data>0.1)
  
  lim_inf = mean(map.data)-2*sigma(map.data)
  lim_sup = mean(map.data)+2*sigma(map.data)
  
  lim_inf = min(map.data)
  lim_sup = max(map.data)/2
  
  size_map = size(map.data,/dim)
  window, 0, xs=size_map[0]<800, ys=size_map[1]<800
  sub_map,map,smap,xrange=xrange, yrange=yrange,irange=irange, $
    color=0, background=255, $
    dmin = lim_inf, dmax=lim_sup                          ; GS:  Select a region to zoom in
    
    
  maxsmap=max(smap.data)
  minsmap=min(smap.data)
  
  
  
  
  
  ;---------------------------------------------------------------
  ;Graphical Interface and definition of slits. Use: graph_intface()
  ;---------------------------------------------------------------
  !P.Background = 255                                 ; White background
  !p.color=0
  
  size_smap = size(smap.data,/dim)
  window, 0, xs=(size_smap[0]*3)<1024, ys=(size_smap[1]*3)<1024
  plot_map, smap, dmin = lim_inf, dmax=lim_sup                   ; plot the new region                        ; plot the new region
  
  kbrd_check=' '                                        ; ' ' (spacebar) is the key to continue making slits.
  T0size=0                                              ; Size of T0
  T0max=0                                               ; Max of T0
  Z0size=0                                              ; Size of Z0
  Z0max=0                                               ; Max of T0
  i=0
  WHILE strcmp(kbrd_check,' ') EQ 1 DO BEGIN            ; While kbrd_check = ' ' continue making slits.
  
    slit_pos=graph_intface()                            ; Make use of the graph_intface function
    
    
    t0=slit_pos.vecx
    z0=slit_pos.vecy
    t0= (t0-xrange[0])/smap.dx       ; GS (generic)
    z0= (z0-yrange[0])/smap.dy
    
    d1 = smap.data
    d1[T0,Z0]=10
    
    
    slitpx=slit_pos
    slitpx.vecx=(slitpx.vecx-xrange[0])/smap.dx
    slitpx.vecy=(slitpx.vecy-yrange[0])/smap.dy
    slitpx.xclick=(slit_pos.xclick-xrange[0])/smap.dx
    slitpx.yclick=(slit_pos.yclick-yrange[0])/smap.dy
    ; ========================================= Sep 5,2012 commented ====================
    ;  sz=size(map.data)
    ;  nx=sz[1]
    ;  ny=sz[2]
    ;
    ;  Xaxismin = min(map.xc-map.dx*(nx-1.)/2.)
    ;  Xaxismax = min(map.xc+map.dx*(nx-1.)/2.)
    ;  Yaxismin = min(map.yc-map.dy*(ny-1.)/2.)
    ;  Yaxismax = max(map.yc+map.dy*(ny-1.)/2.)
    ;
    ;  T0 = round((slit_pos.vecx-Xaxismin)/(Xaxismax-Xaxismin)*(nx-1))
    ;  Z0 = round((slit_pos.vecy-Yaxismin)/(Yaxismax-Yaxismin)*(ny-1))
    ;====================================================================================
    
    ;T0=slit_pos.vecX
    ;Z0=slit_pos.vecY
    
    if i eq 0 then begin                                ; For the 1st slit
      T0size=n_elements(T0)
      Z0size=n_elements(Z0)
    endif else begin                                    ; For the other slits
      T0size=[T0size,n_elements(T0)]                    ; Concatenate all the sizes of the T0s to fix it later
      Z0size=[Z0size,n_elements(Z0)]                    ;
    endelse
    
    T0max = T0max > T0size[i]                           ; Compare with max of T0
    Z0max = Z0max > Z0size[i]
    
    if T0max-T0size[i] gt 0 then T0=[T0,make_array(T0max-T0size[i],1)] ; Complete T0 with 0's if size
    if Z0max-Z0size[i] gt 0 then Z0=[Z0,make_array(Z0max-Z0size[i],1)] ; of T0 is less than the max.
    
    Xout_info = size(Xout)                              ; Size of Xout
    Xout_size = Xout_info[1]                            ; Xout size (number of rows)
    Xout_length = Xout_info[2]                          ; Xout length (number of columns)
    if Xout_info[0] eq 1 then Xout_length = 1
    
    Yout_info = size(Yout)
    Yout_size = Yout_info[1]
    
    if i gt 0 then begin
      if T0max gt Xout_size then Xout=[Xout,make_array(T0max-Xout_size,Xout_length)] ; Complete Xout with 0's if size
      if Z0max gt Yout_size then Yout=[Yout,make_array(Z0max-Yout_size,Xout_length)] ; of T0max is greater than Xout
    endif
    
    if i eq 0 then begin                                ; For the 1st slit
      Xout = T0
      Yout = Z0
    endif else begin                                    ; For the other slits
      Xout = [[Xout],[T0]]
      Yout = [[Yout],[Z0]]
    endelse
    
    kbrd_check = get_kbrd()                             ; Get Keyboard
    
    if strcmp(kbrd_check,'s') THEN BREAK                ; If key = 's' then break
    
    !mouse.button=0
    i++
  ENDWHILE
  
  
  ;Xout=Xout+irange[0]-xrange[0]
  ;Yout=Yout+irange[2]-yrange[0]                         ; from the previously clicked point.
  
  
  ;---------------------------------------------------------------
  ; FILTER to eliminate gaps and wrong pixel position. Use: filter_slit()
  ;---------------------------------------------------------------
  
  Xout_size0=size(Xout)
  Xout_size1=Xout_size0[0]
  Xout_size0=Xout_size0[2]
  if Xout_size1 EQ 1 THEN Xout_size0=1                  ; If Dimension = 1 the third element of SIZE will not be the number of rows
  
  T0max=0
  Z0max=0
  Xout0=Xout
  Yout0=Yout
  T0size0=T0size
  Z0size0=Z0size
  
  FOR k=0,Xout_size0-1 DO BEGIN                         ; Do it for all slits
    IF Xout_size0 EQ 1 THEN BEGIN                       ; If there are only one slit
      T0=Xout0[0:T0size0[k]-1]                          ; One dimension
      Z0=Yout0[0:Z0size0[k]-1]
    ENDIF ELSE BEGIN                                    ; If more than one slit
      T0=Xout0[0:T0size0[k]-1,k]                        ; Two dimensions
      Z0=Yout0[0:Z0size0[k]-1,k]
    ENDELSE
    
    out_filter=filter_slit(T0,Z0)
    
    T0=out_filter.vecX
    Z0=out_filter.vecY
    
    
    if k eq 0 then begin                                ; For the 1st slit
      T0size=n_elements(T0)
      Z0size=n_elements(Z0)
    endif else begin                                    ; For the other slits
      T0size=[T0size,n_elements(T0)]                    ; Concatenate all the sizes of the T0s to fix it later
      Z0size=[Z0size,n_elements(Z0)]
    endelse
    
    T0max = T0max > T0size[k]                           ; Compare with max of T0
    Z0max = Z0max > Z0size[k]                           ; Compare with max of Z0
    
    if T0max-T0size[k] gt 0 then T0=[T0,make_array(T0max-T0size[k],1)] ; Complete T0 with 0's if size
    if Z0max-Z0size[k] gt 0 then Z0=[Z0,make_array(Z0max-Z0size[k],1)] ; of T0 is less than the max.
    
    Xout_info = size(Xout)                              ; Size of Xout
    Xout_size = Xout_info[1]                            ; Xout size (number of rows)
    Xout_length=Xout_info[2]                            ; Xout length (number of columns)
    
    if Xout_info[0] eq 1 then Xout_length = 1           ; If Dimension = 1 the third element of SIZE will not be the number of rows
    
    Yout_info = size(Yout)
    Yout_info = size(Yout)
    Yout_size = Yout_info[1]
    
    if k gt 0 then begin
      if T0max gt Xout_size then Xout=[Xout,make_array(T0max-Xout_size,Xout_length)] ; Complete Xout with 0's if size
      if Z0max gt Yout_size then Yout=[Yout,make_array(Z0max-Yout_size,Xout_length)] ; of T0max is greater than Xout
    endif
    if k eq 0 then begin                                ; For the 1st slit
      Xout = T0
      Yout = Z0
    endif else begin                                    ; For the other slits
      Xout = [[Xout],[T0]]
      Yout = [[Yout],[Z0]]
    endelse
    
    
  ENDFOR
  ;---------------------------------------------------------------
  
  header = HEADFITS(filelist[0])
  cdelt1 = SXPAR(header,'CDELT1')
  
  ;if exist(save_slit) then save, filename=save_slit,xout,yout,irange,cdelt1,t0size
  
  ;======================= 2012-10-03 ==================
  if not exist(save_slit) then save_slit='slit'
  
  
  
  szxout=(size(xout))[2]
  if (size(xout))[0] then szxout=1
  for i=0,szxout-1 do begin
    xout1=xout[*,i] & yout1=yout[*,i]
    t0size1=t0size[i]
    
    save, filename=save_slit+STRTRIM(STRING(i),1)+'.sav',xout1,yout1,irange,cdelt1,t0size1,xrange,yrange
  endfor
  
  if szxout eq 1 then begin
    slit_label=save_slit
    slit_num='0'
    goto,slitdone
  endif else begin
    goto,endjmap
  endelse
  ;stop
  
  slitdone:
  if exist(slit_label) then begin
  
    restore,slit_label+slit_num+'.sav'
    xout=xout1
    yout=yout1
    t0size=t0size1
  endif
  
  ;===========================
  
  T0test=Xout
  Z0test=Yout
  
  
  Xout_info = size(Xout)
  npoint = Xout_info[1]                                  ; Number of points on the linear path
  
  
  
  ;data_cube=sccreadfits(filelist,hdrs)                   ; Read all the de-rotated images into a data cube
  
  n_files=n_elements(filelist)                           ; number of elements in the filelist
  
  
  
  
  
  Xout=float(round(Xout))
  Yout=float(round(Yout))
  
  
  ;---------------------------------------------------------------
  ;         Enlarge the slit to npixel. Use: width_slit()
  ;---------------------------------------------------------------
  
  T1 = MAKE_ARRAY(npoint,npixel)
  Z1 = MAKE_ARRAY(npoint,npixel)
  
  Xout_size0=size(Xout)
  
  Xout_length=Xout_size0[1]
  Xout_size1=Xout_size0[0]
  Xout_size0=Xout_size0[2]
  if Xout_size1 EQ 1 THEN Xout_size0=1
  
  T0max=0
  Z0max=0
  Xout1=Xout
  Yout1=Yout
  ;T0size=Xout_length
  ;Z0size=Xout_length
  T0size1=T0size
  Z0size1=T0size
  
  ;xnew1=make_array(max(t0size),xout_size0)
  ;ynew1=make_array(max(t0size),xout_size0)
  ;xnew2=make_array(max(t0size),xout_size0)
  ;ynew2=make_array(max(t0size),xout_size0)
  
  xnew2=fltarr(max(t0size))
  ynew2=fltarr(max(t0size))
  
  FOR k=0,Xout_size0-1 DO BEGIN                            ; Do it for all slits
    IF Xout_size1 EQ 1 THEN BEGIN                          ; If there are only one slit
      T0=Xout1[0:t0size[k]-1]
      Z0=Yout1[0:t0size[k]-1]
    ;T1[*,0]=Xout1[0:Xout_length-1]                       ; In this case T1 will be 1D
    ;Z1[*,0]=Yout1[0:Xout_length-1]
    ENDIF ELSE BEGIN                                       ; For more than one slit
      T0=Xout1[0:t0size[k]-1,k]
      Z0=Yout1[0:t0size[k]-1,k]
    ;T1[*,0]=Xout1[0:Xout_length-1,k]                     ; In this case T1 will be 2D
    ;Z1[*,0]=Yout1[0:Xout_length-1,k]                     ;
    ENDELSE
    
    
    mask = fltarr( irange[1]-irange[0]+1, irange[3]-irange[2]+1 )
    
    
    ensanchar, T0, Z0, npixel, mask=mask, skeleton=skeleton, cts=cts, xx1=xx1, yy1=yy1, xnew=xnew, ynew=ynew
    xnew2[0:t0size[k]-1]=xnew
    ynew2[0:t0size[k]-1]=ynew
    if k eq 0 then begin
      xnew1=xnew2
      ynew1=ynew2
    endif else begin
      xnew1=[[xnew1],[xnew2]]
      ynew1=[[ynew1],[ynew2]]
    endelse
    
    T1[0:t0size[k]-1,*]=xx1
    Z1[0:t0size[k]-1,*]=yy1
    
    AUXT1=T1
    for i=0,npixel-1 do begin
      T1[*,i]=AUXT1[*,floor(npixel/2.)+ceil(i/2.)*(-1)^(i+1)]
    endfor
    AUXZ1=Z1
    for i=0,npixel-1 do begin
      Z1[*,i]=AUXZ1[*,floor(npixel/2.)+ceil(i/2.)*(-1)^(i+1)]
    endfor
    
    
    ;;=============================================================================================================
    
    
    if k eq 0 then begin                                    ; If only one slit
      Xout2 = T1
      Yout2 = Z1
    endif else begin                                        ; If more than one slit
      if npixel eq 1 then begin                             ; If slits with a 1px width
        Xout2 = [[Xout2],[T1]]                              ; Xout2 will be 2D
        Yout2 = [[Yout2],[Z1]]
      endif else begin                                      ; More than 1px
        Xout2 = [[[Xout2]],[[T1]]]                          ; Xout2 will be 3D
        Yout2 = [[[Yout2]],[[Z1]]]
      endelse
    endelse
  ENDFOR
  
  ;==========================================================================
  jmap=make_array(n_files*(npixel),npoint,/FLOAT)        ; Create array for the J-map
  slit0=make_array(npoint,npixel)                        ; Create Array of the slit
  
  data_cube = fltarr( irange[1]-irange[0]+1, irange[3]-irange[2]+1, n_files)
  FOR j=0,(n_files-1) DO BEGIN
    dum=sccreadfits(filelist[j],hdrs)
    data_cube[*,*,j]=smooth(dum[irange[0]:irange[1], irange[2]:irange[3],*],1,/edge_truncate)
  endfor
  
  ; Filter the wrong images
  ;    value_check1=abs(mean(data_cube)+sigma(data_cube)*3)
  ;    value_check2=abs(mean(data_cube)-sigma(data_cube)*3)
  ;    nerror=0
  ;    meandc=mean(data_cube)
  ;    for i=0,n_files-1 do begin
  ;      if ( abs(mean(data_cube[*,*,i])) gt value_check1 ) || ( abs(mean(data_cube[*,*,i])) lt value_check2 )then begin
  ;        nerror++
  ;        data_cube[*,*,i:n_files-2]=data_cube[*,*,i+1:*]
  ;        filelist[i:n_files-2]=filelist[i+1:*]
  ;        header_time[i:n_files-2]=header_time[i+1:*]
  ;      endif
  ;    endfor
  ;    data_cube=data_cube[*,*,0:n_files-1-nerror]
  ;    filelist=filelist[0:n_files-1-nerror]
  ;    header_time=header_time[0:n_files-1-nerror]
  
  
  if keyword_set(log) and (keyword_set(running) eq 0) then data_cube = alog(data_cube>0.1)
  ;===============================================================================================
  
  
  ;---------------------------------------------------------------
  ;  Running difference
  ;---------------------------------------------------------------
  
  if keyword_set(running) eq 1 then begin
    cube=(data_cube-shift(data_cube,0,0,running))[*,*,running+2:*]
    data_cube=cube
    n_files=(size(data_cube,/dim))[2]
    filelist=filelist[running+2:*]
    header_time=header_time[running+2:*]
  endif
  
  
  
  ;---------------------------------------------------------------
  ; Make the pixels slit based in the pixels of the discretized curve.
  ; The number of columns will be: files*npixel
  ; Use: extract_slit()
  ;---------------------------------------------------------------
  flag_jmaptype=0
  jmapagain:
  jmap=make_array(n_files*(npixel),npoint,/FLOAT)        ; Create array for the J-map
  FOR m=0,Xout_size0-1 DO BEGIN                                         ; Do it for all slits
    IF npixel EQ 1 THEN BEGIN                                           ; If npixel = 1
      T1=Xout2[0:Xout_length-1,m]                                       ; T1 with 2D
      Z1=Yout2[0:Xout_length-1,m]
    ENDIF ELSE BEGIN                                                    ; If npixel > 1
      T1=Xout2[0:Xout_length-1,0:npixel-1,m]                            ; T1 with 3D
      Z1=Yout2[0:Xout_length-1,0:npixel-1,m]
    ENDELSE
    
    
    out_extract_slit=extract_slit(T1,Z1,data_cube,jmap,slit0,n_files,npoint,npixel)
    
    
    ;flag_jmaptype=1
    jmap=out_extract_slit.jmap
    
    
    if set_jmap_type EQ 1 then begin
      jmap=out_extract_slit.jmap_avg
      npixel=1 ; Comment this line you if don't make use of the avg.
    endif
    jmapori=out_extract_slit.jmap_avg
    
    
    if m eq 0 then begin                                                 ; For 1 slit
      jmap_full = jmap
    endif else begin                                                     ; For > 1 slits
      jmap_full = [[[jmap_full]],[[jmap]]]                               ; The final jmap in 3D
    endelse
  ENDFOR
  
  ;  jmapnofilter=jmap
  ;  model = smooth(jmap,5,/edge)
  ;  for iter = 0,100 do model = smooth(model,5,/edge)
  ;  jmap=  (jmap - model) >(-4)<4
  
  
  
  ;---------------------------------------------------------------
  ;                          Plot JMAP
  ;---------------------------------------------------------------
  
  ;totaltime=N_FILES*1.25
  ; 1m15s by file
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  
  
  
  for m=0,Xout_size0-1 DO BEGIN                                          ; For all slits
    fov=1.7                                                              ; Field of View = 1.7 solar radii
    if not exist(stepsize) then stepsize=60
    
    out_timeaxisjmap = timeaxisjmap(jmap_full[*,0:t0size[m]-1,m],set_jmap_type,npixel,filelist,stepsize)
    jmap=out_timeaxisjmap.jmap
    totaltime=out_timeaxisjmap.totaltime/60.     ; in minutes
    xaxis=float([0,totaltime/5.,totaltime/5.*2,totaltime/5.*3,totaltime/5.*4,totaltime])                                  ; X-axis label
    ;npoint=T0size
    srtotal=dist_along_curve(xnew1[0:t0size[m]-1,m],ynew1[0:t0size[m]-1,m])*CDELT1
    yaxis=findgen(7)*(srtotal/6.)                                         ; Y-axis label
    
    ;mapmax=max(map.data)
    ;mapmin=min(map.data)
    
    ;jmap=jmap_full[*,*,m]
    jmapsize=size(jmap)
    
    if jmapsize[1] gt jmapsize[2] then begin
      scalefactor1=100
      scalefactor2=(jmapsize[1]/float(jmapsize[2]))*1/2*100.
    endif else begin
      scalefactor1=(jmapsize[2]/float(jmapsize[1]))*2*100.
      scalefactor2=100
    endelse
    
    str_hour=STRTRIM(string(hour[startfile]),1) & if STRLEN(str_hour) eq 1 then str_hour='0'+str_hour
    str_min=STRTRIM(string(min[startfile]),1) & if STRLEN(str_min) eq 1 then str_min='0'+str_min
    str_sec=STRTRIM(string(sec[startfile]),1) & if STRLEN(str_sec) eq 1 then str_sec='0'+str_sec
    str_starttime=str_hour+str_min+str_sec
    cutcube=0
    if keyword_set(running) eq 1 then cutcube=2 else running=0
    str_hour=STRTRIM(string(hour[endfile-running-cutcube]),1) & if STRLEN(str_hour) eq 1 then str_hour='0'+str_hour
    str_min=STRTRIM(string(min[endfile-running-cutcube]),1) & if STRLEN(str_min) eq 1 then str_min='0'+str_min
    str_sec=STRTRIM(string(sec[endfile-running-cutcube]),1) & if STRLEN(str_sec) eq 1 then str_sec='0'+str_sec
    str_endtime=str_hour+str_min+str_sec
    
    
    !p.color=0
    scalewin=3./npixel

    window,m+flag_jmaptype,XSIZE=scalewin*(size(jmap,/dim))[0], YSIZE=scalewin*(size(jmap,/dim))[1]+70
    ;window,m+flag_jmaptype,XSIZE=1200, YSIZE=500
    ;stop
    
    str_date=header_time[0]
    str_date=STRMID(str_date,0,10)
    
    wave=SXPAR(header,'WAVELNTH')
    wave=STRTRIM(STRING(wave),1)
    
    ;===========2012-10-04==============
    jmap3=jmap
    if flag_jmaptype eq 2 then begin
      k=[1,1,1]
      jmap2=convol(jmap,k,/norm,/edge_tr)
      for t=0,40 do jmap2=convol(jmap2,k,/norm,/edge_tr)
      jmap3=hist_equal(jmap-jmap2,per=5)
    endif
    if keyword_set(running) then jmap3=hist_equal(jmap,per=1)
    
    plot_image,jmap3, $
      charsize=1.2,CHARTHICK=1.4,$;scale=[scalefactor1,scalefactor2]
      xticks=5, xtickname=STRTRIM(string(xaxis,FORMAT='(f6.2)'),1), $
      yticks=6, ytickname=STRTRIM(STRING(yaxis, FORMAT='(f6.2)'), 1), $
      XTITLE = 'Time elapsed start time (min)', YTITLE = 'Distance (arcsec)', /smooth, $
      TITLE = 'wave: '+ wave + '  ' + str_date +'   '+ 'Start: '+ str_starttime + '   ' + 'end: ' + str_endtime +  '   ' + 'slit: '+ slit_num,$
      BACKGROUND = 255, COLOR = 0;,scale=[3,1]
      
      
    ;/smooth, TITLE = 'Slit '+strtrim(m,2),$
    ;max=maxsmap, min=minsmap                                              ; Preserve the same contrast
    fnameout='jmap'+'_'+str_starttime+'_'+str_endtime+'_'+STRTRIM(string(flag_jmaptype),1)+ '_slit_'+slit_num+'.png'
    write_png, fnameout,tvrd()
    clickjmap,jmap,srtotal,stepsize,npixel,cdelt1,wave=wave,slit_label=slit_label,slit_num=slit_num,spfile=spfile
    tvjmap1=tvrd()
  endfor
  if exist(spfile) && ( (size(spfile))[0] gt 1 ) then begin
    window,3 & HISTOGAUSS, spfile[4,*], A,$
      XTITLE= 'Speed (km/s)', YTITLE='Frequency',$
      TITLE = 'wave: '+ wave + '  ' + str_date +'   '+ 'Start: '+ str_starttime + '   ' + 'end: ' + str_endtime +  '   ' + 'slit: '+ slit_num,$
      background=255,color=0,/nofit
    tvhisto1=tvrd()
    tvjmaphisto1=bytarr((size(tvjmap1,/dim))[0],(size(tvjmap1,/dim))[1]+(size(tvhisto1,/dim))[1])+255
    tvjmaphisto1[(size(tvjmap1,/dim))[0]/2-(size(tvhisto1,/dim))[0]/2:(size(tvjmap1,/dim))[0]/2-(size(tvhisto1,/dim))[0]/2+(size(tvhisto1,/dim))[0]-1,0:(size(tvhisto1,/dim))[1]-1]=tvhisto1
    tvjmaphisto1[*,(size(tvhisto1,/dim))[1]:(size(tvjmap1,/dim))[1]+(size(tvhisto1,/dim))[1]-1]=tvjmap1
    
  endif
  
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Speed Determination 1st order
  
  if flag_jmaptype eq 2 && exist(spfile) then begin
    if (size(spfile))[0] ne 1 then newspeed=fltarr((size(spfile,/dim))[0],(size(spfile,/dim))[1]) else $
      newspeed=fltarr((size(spfile,/dim))[0])
    for p=0,n_elements(spfile[0,*])-1 do begin
      x0=spfile[0,p] & x1=spfile[1,p]
      y0=spfile[2,p] & y1=spfile[3,p]
      adjustspeedpf_v0,jmap,srtotal,stepsize,x0,x1,y0,y1,/min,/arcsec,speed=speed,cdelt1=cdelt1,degree=1
      newspeed[*,p]=speed
    endfor
    save,newspeed,FILENAME = 'speeddata_auto' + wave + slit_label + slit_num + '.sav'
    if (size(spfile))[0] ne 1 then begin
    window,4 & HISTOGAUSS, newspeed[4,*], A,$
      XTITLE= 'Speed (km/s)', YTITLE='Frequency',$
      TITLE = 'wave: '+ wave + '  ' + str_date +'   '+ 'Start: '+ str_starttime + '   ' + 'end: ' + str_endtime +  '   ' + 'slit: '+ slit_num,$
      background=255,color=0,/nofit
    tvhisto2=tvrd()
    endif
    
    window,5,XSIZE=scalewin*(size(jmap,/dim))[0], YSIZE=scalewin*(size(jmap,/dim))[1]+70
    plot_image,jmap3, $
      charsize=1.2,CHARTHICK=1.4,$;scale=[scalefactor1,scalefactor2]
      xticks=5, xtickname=STRTRIM(string(xaxis,FORMAT='(f6.2)'),1), $
      yticks=6, ytickname=STRTRIM(STRING(yaxis, FORMAT='(f6.2)'), 1), $
      XTITLE = 'Time elapsed start time (min)', YTITLE = 'Distance (arcsec)', /smooth, $
      TITLE = 'wave: '+ wave + '  ' + str_date +'   '+ 'Start: '+ str_starttime + '   ' + 'end: ' + str_endtime +  '   ' + 'slit: '+ slit_num,$
      BACKGROUND = 255, COLOR = 0
      
      
    if (size(spfile))[0] ne 1 then limsup=(size(newspeed,/dim))[1]-1 else limsup=0
    for r=0,limsup do begin
      plots,newspeed[0,r],newspeed[2,r],PSYM=6,SYMSIZE=0.5,color=255
      plots,newspeed[1,r],newspeed[3,r],PSYM=6,SYMSIZE=0.5,color=255
      oplot,newspeed[0:1,r],newspeed[2:3,r],color=255
      xyouts,newspeed[1,r],newspeed[3,r]+3,STRTRIM(STRING(r+1),1),color=255
    endfor
    tvjmap2=tvrd()
  endif
  if ( exist(newspeed) ) && ( (size(spfile))[0] ne 1 ) then begin
    tvjmaphisto2=bytarr((size(tvjmap2,/dim))[0],(size(tvjmap2,/dim))[1]+(size(tvhisto2,/dim))[1])+255
    tvjmaphisto2[(size(tvjmap2,/dim))[0]/2-(size(tvhisto2,/dim))[0]/2:(size(tvjmap2,/dim))[0]/2-(size(tvhisto2,/dim))[0]/2+(size(tvhisto2,/dim))[0]-1,0:(size(tvhisto2,/dim))[1]-1]=tvhisto2
    tvjmaphisto2[*,(size(tvhisto2,/dim))[1]:(size(tvjmap2,/dim))[1]+(size(tvhisto2,/dim))[1]-1]=tvjmap2
    while !D.WINDOW ne -1 do wdelete
    window,0,xsize=(size(tvjmaphisto1,/dim))[0],ysize=(size(tvjmaphisto1,/dim))[1] & tv,tvjmaphisto1
    fnameout='jmap_histo_manual'+ wave + slit_label + slit_num + '.png'
    write_png, fnameout,tvrd()
    window,1,xsize=(size(tvjmaphisto2,/dim))[0],ysize=(size(tvjmaphisto2,/dim))[1] & tv,tvjmaphisto2
    fnameout='jmap_histo_auto'+ wave + slit_label + slit_num + '.png'
    write_png, fnameout,tvrd()
  endif
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Speed Determination 2nd order
  if exist(spfile) then begin
  if (size(spfile))[0] ne 1 then imgcube=bytarr(300,300,(size(spfile,/dim))[1]) else imgcube=bytarr(300,300,1)
  endif
  
  if flag_jmaptype eq 2 && exist(spfile) then begin
    if (size(spfile))[0] ne 1 then newspeed=fltarr((size(spfile,/dim))[0]+4,(size(spfile,/dim))[1]) else $
    newspeed=fltarr((size(spfile,/dim))[0]+4)
    for p=0,n_elements(spfile[0,*])-1 do begin
      x0=spfile[0,p] & x1=spfile[1,p]
      y0=spfile[2,p] & y1=spfile[3,p]

     
      adjustspeedpf_v0,jmap,srtotal,stepsize,x0,x1,y0,y1,/min,/arcsec,speed=speed,cdelt1=cdelt1,degree=2
      newspeed[*,p]=speed
      if exist(img) then imgcube[*,*,p]=img
    endfor
    save,newspeed,FILENAME = 'speeddata_auto_2nd' + wave + slit_label + slit_num + '.sav'
    
    window,6,XSIZE=scalewin*(size(jmap,/dim))[0], YSIZE=scalewin*(size(jmap,/dim))[1]+70
    plot_image,hist_equal(jmap-jmap2,per=1), $;min=-0.0602702,max=0.06144
      charsize=1.2,CHARTHICK=1.4,$;scale=[scalefactor1,scalefactor2]
      xticks=5, xtickname=STRTRIM(string(xaxis,FORMAT='(f6.2)'),1), $
      yticks=6, ytickname=STRTRIM(STRING(yaxis, FORMAT='(f6.2)'), 1), $
      XTITLE = 'Time elapsed start time (min)', YTITLE = 'Distance (arcsec)', /smooth, $
      TITLE = 'wave: '+ wave + '  ' + str_date +'   '+ 'Start: '+ str_starttime + '   ' + 'end: ' + str_endtime +  '   ' + 'slit: '+ slit_num,$
      BACKGROUND = 255, COLOR = 0
      
      
      if (size(spfile))[0] ne 1 then limsup2=(size(newspeed,/dim))[1]-1 else limsup2=0
    for r=0,limsup2 do begin
      xx = findgen(1000)*(newspeed[1,r]-newspeed[0,r])/999.+newspeed[0,r]
      yy = newspeed[6,r]*xx^2.+newspeed[7,r]*xx+newspeed[8,r]
      
      plots,min(xx+newspeed[4,r]),min(yy+newspeed[5,r])-1.,PSYM=6,SYMSIZE=0.5,color=0
      plots,max(xx+newspeed[4,r]),max(yy+newspeed[5,r])-1.,PSYM=6,SYMSIZE=0.5,color=0
      oplot, xx+newspeed[4,r],yy+newspeed[5,r]-1.,color=0,thick=1
      xyouts,max(xx+newspeed[4,r]),max(yy+newspeed[5,r])+3.,STRTRIM(STRING(r+1),1),color=255
    endfor
  endif
  
  if set_jmap_type eq 2 then begin
    set_jmap_type=1 & m=0 & flag_jmaptype=2 & npixel=1
    goto,jmapagain
  endif
  
  endjmap:
  stop
  return & end
;---------------------------------------------------------------