function timeaxisjmap,jmap,set_jmap_type,npixel,filelist,stepsize

; stepsize in sec.

  ;cd, 'C:\sdodata\20120831\wavelets'                 ; FITS file directory
  ;filelist=file_list('.','*.fts')
  ;filelist = filelist[0:100]
  
  nfiles=n_elements(filelist)
  header_time=strarr(nfiles)  
 
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
    time_sec=sec+min*60.+hour*60.*60.
  
  
;  small_interv=time_sec[0,1]-time_sec[0,0]
;  for i= 1, (nfiles-2) do $
;    small_interv = small_interv < (time_sec[0,i+1] - time_sec[0,i])
    
  
  
  time_sec_old=time_sec[*,*]-time_sec[0,0]
  time_sec[*,*]=round((time_sec[*,*]-time_sec[0,0])/stepsize)*stepsize
  RMSD=sqrt(total((time_sec_old-time_sec)^2.)/nfiles)
  
  xaxis_time=time_sec[nfiles-1]         ; Time in seconds
  xaxis_size=xaxis_time/stepsize+1          ; Quantity of columns in the new J-map
  xaxis_position=time_sec/stepsize        ; Position of each columm of the old J-map in the new J-map
  
  sizejmap=size(jmap)
  slit_length=sizejmap[2]
  
  if set_jmap_type eq 1 then npixel=1
  
  jmap_new=make_array(xaxis_size*npixel,slit_length)
  ;jmap_new[*,*]=!Values.F_NAN
  
  for i=0,nfiles-1 do begin
    jmap_new[xaxis_position[i]*npixel:xaxis_position[i]*npixel+npixel-1,*]=jmap[i*npixel:i*npixel+npixel-1,*]
  endfor

  out_timeaxisjmap = {jmap:jmap_new , totaltime:time_sec_old[nfiles-1]}
  
 
  return,out_timeaxisjmap
  
end