pro jmap_periods2, jmap, xaxis,yaxis,header_time,row=row,width=width,step=step,full=full,aver=aver, $
    xout=xout, yout=yout,jmapout=jmapout,wave=wave,slit_label=slit_label,slit_num=slit_num
    
    
  hour=float(STRMID(header_time,11,2))
  min=float(STRMID(header_time,14,2))
  sec=float(round(float(STRMID(header_time,17,5))))
  time_min=sec/60.+min+hour*60.
  x = time_min
  
  if keyword_set(step) && keyword_set(row) then begin
    print,'Select only one: STEP or ROW' & STOP
  endif
  
  if keyword_set(step) then begin
    row=findgen( ((size(jmap))[2]-width/2) / step)*step+width/2
    if row[n_elements(row)-1]+width/2 ge (size(jmap))[2] then row=row[0:n_elements(row)-2]
  endif
  
  if keyword_set(full) then begin
    width=(size(jmap))[2]-1
    row=(size(jmap))[2]/2
  endif
  
  for j=0,n_elements(row)-1 do begin
  
    widthname=width
    if keyword_set(aver) then begin
      if width eq 1 then y=jmap[*,row[j]] else $
      y=total(jmap[*,(-width/2)+row[j]:(width/2)+row[j]],2)/width
      width=1
    endif
    
    for i = (-width/2),(width/2) do begin
      if not keyword_set(aver) then y=jmap[*,row[j]+i]
      r=lnp_test(x,y,wk1=wk1,wk2=wk2,ofac=50)
      

        maxposi= where(r[0] eq wk2)
        maxposifreq=wk1[maxposi]
        maxposiSet= exist(maxposiSet) ? [maxposiSet,maxposifreq] : maxposifreq
        sigset= exist(sigset) ? [sigset,r[1]] : r[1]
        freqSet = wk1
        amplSet = exist(amplSet) ? [[amplSet],[wk2]] : wk2

    endfor
    
    if exist(amplset) then begin
    
      if (size(amplset))[0] ne 1 then amplSet=avg(amplSet,1)
      
      window,j & plot,1./freqSet,amplSet, xrange=[0,10], $
        XTITLE='Period (min)', YTITLE='Lomb Normalized Periodogram', $
        TITLE = 'Wave: '+ wave + '  slit: '+ slit_num, $
        ;TITLE='Range: ' + STRTRIM(STRING(row[j]-widthname/2,FORMAT='(i6)'),1) + ' - ' + STRTRIM(STRING(row[j]+widthname/2,FORMAT='(i6)'),1), $  ; Periodicities in minuts
        BACKGROUND = 255, COLOR = 0
       if exist(wave) && exist(slit_label) && exist(slit_num) then begin
       fnameout='Periodogram'+ wave + slit_label + slit_num + '_range_'+STRTRIM(STRING(row[j]-widthname/2,FORMAT='(i6)'),1) +'_'+ STRTRIM(STRING(row[j]+widthname/2,FORMAT='(i6)'),1) + '.png'
    write_png, fnameout,tvrd()
    endif
      if (n_elements(row) eq 1) then begin
        if not keyword_set(aver) && width ne 1 then begin
          window,n_elements(row)
          plot,sigset,psym=1,Title='Significance of the maximum peak',BACKGROUND = 255, COLOR = 0,yrange=[min(sigset)-0.1,max(sigset)+0.1],$
          xtitle='Row of the J-map',ytitle='Significance of the maximum peak'
          window,n_elements(row)+1
          plot,1./maxposiSet,psym=1,Title='Position of the maximum peak',BACKGROUND = 255, COLOR = 0,$
          xtitle='Row of the J-map',ytitle='Period (min) of the maximum peak'
        endif
      endif else begin
        tmp=temporary(amplset)
        tmp=temporary(maxposiSet)
      endelse
    endif
  endfor
  
;  xout = x
;  yout = y
  jmapout=jmap
  
END