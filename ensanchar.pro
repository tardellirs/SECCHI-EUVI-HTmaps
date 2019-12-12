 pro ensanchar, T0, Z0, width, mask=mask, skeleton=skeleton, cts=cts, xx1=xx1, yy1=yy1, xnew=xnew, ynew=ynew
   mask1=mask
   ;restore, filename='pruebas.sav'   ;, t0,z0, npixel
   ; width = 31
   npixel=width
   
   x = T0
   y = Z0
 
   npoint = n_elements(x)
   skeleton = lonarr(npoint, width+100)   ; + 100 to be sure
   cts = intarr(npoint)
   
   ns = npoint
   
   ;sg=savgol(19,19,0,2)
   sg=savgol(19,19,0,2)
   if n_elements(x) le 19 then sg=savgol(n_elements(x)/3,n_elements(x)/3,0,2)
   xnew1 =  convol(x,sg,/edge_trun)
   ynew1 =  convol(y,sg,/edge_trun)
   
   out_ede = EQUAL_DIST_ENSANCHAR(xnew1,ynew1)
   xnew=out_ede.xnew
   ynew=out_ede.ynew

   x=xnew
   y=ynew
   x_ = [x,x[ns-1]*2-x[ns-2]]
   y_ = [y,y[ns-1]*2-y[ns-2]]
   dx = x_[1:ns]-x_[0:ns-1]
   dy = y_[1:ns]-y_[0:ns-1]
   der = dy/dx
   
   npoint = n_elements(x)
   ;  idx_der = where(abs(der) le 0.001, ctder)
   ;  if ctder gt 0 then der(idx_der)= 0.001
   ; der = smooth(der,25,/edge,/nan, missing=0)
   ;der = convol(der,sg,/edge_trun)
   
;   window, 3, xs =512, ys=512, $
;     title='Selecting pixels along perpendicular direction to skeleton at each point...'
;;   plot, T0,Z0, psym=3, xrange = [min(T0)-width/2-1, max(T0)+width/2+1], $
;;     yrange = [min(Z0)-width/2-1, max(Z0)+width/2+1]
;   plot, xnew,ynew, psym=3,BACKGROUND = 255, COLOR = 0
     
     
    
;============================================= GS ==============================
     
   for j = 0, npoint-1 do begin
     denuevo:
     if npixel*1000l lt 10 then nn=10000l else nn = 1000l
     if npixel*nn lt 10 then nn=10000000l else nn = 1000l

     xx = xnew[j] + ( (findgen(npixel*nn) - npixel*nn/2) ) /float(nn)
     
     
     ; if (abs(der[j]) le 0.001) then der[j] = 0.001*sgn(der[j])
     if (abs(der[j]) lt 0.05) then begin      ;; CHECK!!! 0.03
       xx = replicate(T0[j], width*10)
       yy = Z0[j] + (findgen(width*10) - width*10/2)/10.
       longitud = max(yy) - min(yy)
       if (longitud gt width+1) then stop
       goto, listo0
     endif
     
     yy = (-1./der[j])* (xx - xnew[j]) + ynew[j]
     longitud = sqrt (  (xx[npixel*nn-1] - xx[0])^2 + (yy[npixel*nn-1] - yy[0])^2   )
     
     if (longitud gt width+1) and (npixel gt 0.1) then begin
       npixel = float(npixel)
       npixel = npixel-0.1
       goto, denuevo
     endif
     if (longitud gt width+1) and (npixel le 0.1) and (npixel gt 0.01) then begin
       npixel = float(npixel)
       npixel = npixel-0.01
       goto, denuevo
     endif
     
     listo0:
     ; help, xx,yy   , der[j]
     ; print, j, min(yy),max(yy), longitud, npixel
     npixel = width
  ;;;;;;;;;   oplot, xx, yy, psym=3,color=50
     if (longitud gt width+1) then stop
     
     
     ;================================ 20120910
     nsample=npixel
     
     xmax=max(xx)
     xmin=min(xx)
     ymax=max(yy)
     ymin=min(yy)
     xmax_pos=where(xx eq xmax)
     xmin_pos=where(xx eq xmin)
    
     x1=xmax & x0=xmin
     y1=yy(xmax_pos(0))
     y0=yy(xmin_pos(0))
        

     m=(y1-y0)/(x1-x0)
     
     if x1 ne x0 then begin
     xx0=findgen(nsample)/(nsample-1)*(x1-x0)+x0
     yy0=m*xx0 + y1-m*x1
     yy0_old=yy0
     
     endif else begin
     xx0=x0+fltarr(npixel)
     yy0=(findgen(nsample)/(nsample-1)*nsample)+y0
     ;m = 9999.
     ;yy0=m*xx0 + y1-m*x1
     ;yy0=yy0_old
    
     endelse
     ;oplot,transpose(xx0), transpose(yy0), psym=-3,color=150
     
     
     if j eq 0 then begin
       yy1=yy0
       xx1=xx0
     endif else begin
       xx1=[[xx1],[xx0]]
       yy1=[[yy1],[yy0]]
     endelse
     ;===============================================

     
;     ; -------------------------- TS ------------------
;     if j eq 0 then begin
;       xx1=xx
;       yy1=yy
;       
;     endif else begin
;       n_row_exc=n_elements(xx1)/j-n_elements(xx) ; number of row that exceed
;       
;       ;lower_xx=lower_xx<n_elements(xx)
;       
;       if n_row_exc gt 0 then begin
;         xx1=[[xx1],[xx,fltarr(n_row_exc)*!VALUES.F_NAN]]
;         yy1=[[yy1],[yy,fltarr(n_row_exc)*!VALUES.F_NAN]]
;       endif else begin
;         if n_row_exc lt 0 then begin
;           xx1=[[xx1,fltarr(-n_row_exc)*!VALUES.F_NAN],[xx]]
;           yy1=[[yy1,fltarr(-n_row_exc)*!VALUES.F_NAN],[yy]]
;         endif else begin
;           xx1=[[xx1],[xx]]
;           yy1=[[yy1],[yy]]
;         endelse
;       endelse
;     endelse
     
     
     
     
   ;;;;;;;;;;;;;;;;;;;;  mask[xx,yy] = 1
   ;;;;;;;;;;;;;;;;;;;;  mask1[xx,yy]=1
   ;;;;;;;;;;;;;;;;;;;;  idx = where(mask eq 1,ct)
     
   ;;;;;;;;;;;;;;;;;;;;  skeleton[j,0:ct<100-1] = idx[0:ct<100-1]
   ;;;;;;;;;;;;;;;;;;;;  cts[j] = ct<100
     ;   mask1[reform(skeleton[j,0:ct-1])] = 1
   ;;;;;;;;;;;;;;;;;;;;  mask[*,*]=0
   endfor
;==========================================================================

;xx1size=size(xx1)
; for i=0,xx1size(2)-2 do begin
;   if yy1[0,i+1]-yy1[0,i] gt 1 then begin
;     yy1_1=yy1[*,0:i]
;     yy1_2=reverse(yy1[*,i+1:xx1size(2)-1],2)
;     xx1_1=xx1[*,0:i]
;     xx1_2=reverse(xx1[*,i+1:xx1size(2)-1],2)
;   endif
; endfor
; ayy1=[[yy1_1],[yy1_2]]
; axx1=[[xx1_1],[xx1_2]]



 
;     xx1size=size(xx1)
;     pxx1=fltarr(npixel,xx1size(2))
;     pyy1=pxx1
;     
;     stepsize=xx1size[1]/npixel
;     for i= 0,xx1size[1]-1,(npixel-1)*stepsize do begin
;     pxx1[i/(stepsize-1),*]=xx1[i,*]
;     pyy1[i/(stepsize-1),*]=yy1[i,*]
;     endfor

     

; for i=0,npoint-2 do begin
;   if abs(yy1[0,i+1]-yy1[0,i]) gt 3 then begin
;     yy1_1=yy1[*,0:i]
;     yy1_2=reverse(yy1[*,i+1:npoint-1],2)
;     xx1_1=xx1[*,0:i]
;     xx1_2=reverse(xx1[*,i+1:npoint-1],2)
;   endif
; endfor
; yy1=[[yy1_1],[yy1_2]]
; xx1=[[xx1_1],[xx1_2]]



xx1=transpose(xx1)
yy1=transpose(yy1)



;   S = REPLICATE(1, 3, 3)
;   mask1 = erode(dilate(mask1, s),s)
;   tvscl, mask1
   
 ; for t=0,npoint-1,2 do begin
 ;   mask1[skeleton[t,*]]=2
 ;   tvscl, mask1
 ;   wait, 1
 ;endfor
   
 ;stop
 END