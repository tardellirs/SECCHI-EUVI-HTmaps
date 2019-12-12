 function equal_dist_ensanchar,xnew,ynew
 ;restore,'C:\sdodata\20110706\171\DeroteData\testensanchar.sav'
 SPLINE_P, xnew, ynew, Xr, Yr,interval=0.001
 distcurve=dist_along_curve(xnew,ynew)
 distpoints=distcurve/(n_elements(xnew)-1)
 xout=xnew[0]
 yout=ynew[0]
 offset=float(0)
 
 for j=0,n_elements(xnew)-1 do begin
   for i=900,1100 do begin
     if i+offset ge n_elements(xr) then goto,end_of_xr
     lastxout=n_elements(xout)-1
     rest=abs( sqrt( (xout[lastxout] - xr[i+offset]) * (xout[lastxout] - xr[i+offset]) + (yout[lastxout] - yr[i+offset]) * (yout[lastxout] - yr[i+offset]) ) - distpoints )
     if exist(minrest) then minrest_old=minrest
     minrest = exist(minrest) ? rest < minrest : rest
     if exist(minrest_old) && (minrest_old ne minrest) then begin
     i_offset=i+offset
     endif
   endfor
   end_of_xr:
   offset=i_offset
   ;if exist(offset_old)then print,minrest,offset,offset-offset_old else print,minrest,offset
   offset_old=offset
   ;stop
   xout=[xout,xr[i_offset]]
   yout=[yout,yr[i_offset]]
   minrest=20
   minrest_old=20
 endfor
 if n_elements(xout) gt n_elements(xnew) then begin
 xout=xout[0:n_elements(xout)-2]
 yout=yout[0:n_elements(yout)-2]
 endif
out_ede={xnew:xout, ynew:yout}
;stop
return,out_ede
 END