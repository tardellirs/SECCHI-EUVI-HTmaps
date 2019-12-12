;if (x[i] eq x[i-1]) || (y[i] eq y[i-1]) then distcurve++ else distcurve= distcurve + 1.41421

function dist_along_curve,x,y

; === 2012/09/14 ===
  x1 = [shift(x, -1), 0.]
  y1 = [shift(y, -1), 0.]
  d = sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1))
  distcurve=total(d[0:n_elements(d) - 2L])
  goto,fin
 
; ==== TARDELLI =====
  curve_size=n_elements(x)
  distcurve=0
  for i=1,curve_size-1 do $
    distcurve = (x[i] eq x[i-1]) || (y[i] eq y[i-1]) ? $
    (distcurve+1) : (distcurve+1.41421)
   
goto, fin
    
 ; ==== GUILLERMO   
    
origin = [x(0),y(0)]
n = n_elements(x)
dist_int = fltarr(n)
dist_lin = fltarr(n)


; ---------------------------------------------------------------
x_ori = x
y_ori = y
sg=savgol(33<(n_elements(x)/2-1), 33<(n_elements(x)/2-1), 0, 2)         ; SOLO PARA MAYOR ACCURACY EN LA DERIVADA              
                  ; 11 para dummy1

xnew = convol(x,sg,/edge_trun)
ynew = convol(y,sg,/edge_trun)
x=xnew
y=ynew

  ns = (size(x,/dim))(0)
 ; derivada = deriv(x, y)   
  
  x_  = [x,x(ns-1)*2-x(ns-2)] ;extrapolate 1 additional value 
  y_  = [y,y(ns-1)*2-y(ns-2)]
 dx = x_(1:ns)-x_(0:ns-1)
 dy = y_(1:ns)-y_(0:ns-1)
 derivada = dy/dx
; -----------------------------------------------------------

; =============================================================
; === this is algorithm to compute the  "distance along a curve" ================
; =============================================================
for j = 5, n-1 do begin    ; start in 5 to avoid border effects
;;;;    der = deriv( x(0:j), y(0:j) )
    der = derivada(0:j)       ; it's more efficient and accurate than using the pre-defined routine "deriv"
    F = (1 +(der^2.) ) ^0.5
    dist_int (j) = int_tabulated(x(0:j),F)
    dist_lin(j) = ( (origin(1) - y(j) )^2. + ( origin(0) - x(j) )^2. )^0.5
endfor
;=============================================================

window, 5, xs = 300, ys = 200
plot, x, dist_lin
oplot,x, dist_int, linestyle = 1
    
    
fin:
    
    
    
  return,distcurve
  stop
end