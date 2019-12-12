;+
; Project     :
;
; Name        : DEROTATE_TEST
;
; Purpose     :
;
; Category    :
;
; Syntax      : DEROTATE_TEST
;
; Inputs      :
;
; Calling     : derotate_test,path='c:/sdodata',timeRef=053010,startTime=025020,endTime=0830010
;-

;---------------------------------------------------------
; De-rotate all the images in a path and save in a new one
;---------------------------------------------------------

pro derotate_test,path=path,timeRef=timeRef,startTime=startTime,endTime=endTime ;time= hhmmss
  compile_opt strictarr
  
  
  cd, path                          ; Define the path of the original files
  
  new_dir='.\DeroteData\'
  file_mkdir, new_dir                          ; The path to save the de-rotated images (in this case C:\StereoDatatest\DeroteData\)
  
  
  filelist=file_list('.','*.fts')                   ; Create a list of all the FITS file in the path '.'

  
  nfiles=n_elements(filelist)
  header_time=strarr(nfiles)
  time_sec=fltarr(nfiles)
  
  for i=0,nfiles-1 do begin
    header = HEADFITS(filelist[i])
    header_time[i] = SXPAR(header,'DATE-OBS')
  endfor
    hour=long(STRMID(header_time,11,2))
    min=long(STRMID(header_time,14,2))
    sec=long(round(float(STRMID(header_time,17,5))))
    time_sec=sec+min*60.+hour*60.*60.

  ;-----------------------------
  if not exist(startTime) then begin
  startfile=0 & endfile=n_elements(filelist)-1
  goto,cutfilelist
  endif
  startTimeSec=float(STRMID(startTime,0,2))*60.*60.+float(STRMID(startTime,2,2))*60.+float(STRMID(startTime,4,2))
  endTimeSec=float(STRMID(endTime,0,2))*60.*60.+float(STRMID(endTime,2,2))*60.+float(STRMID(endTime,4,2))
  
  startFile=where(time_sec ge startTimeSec)
  startFile=startFile[0]
  endFile=where(time_sec le endTimeSec)
  endFile=endFile[n_elements(endFile)-1]
  
  cutfilelist:
  filelist=filelist[startfile:endfile]
  ;---------------------------------
  if not exist(timeRef) then begin
    fileRef='000000'
  endif else begin
    timeRefSec=float(STRMID(timeRef,0,2))*60.*60.+float(STRMID(timeRef,2,2))*60.+float(STRMID(timeRef,4,2))
    nearRef=abs(time_sec-timeRefSec)
    fileRef=where(nearRef eq min(nearRef))
    fileRef=fileRef[0]
  endelse
    fileRef=fileRef-startfile


filenames=STRMID(filelist, 2)                    ; Create a list of only the filenames (without path) '.'



fits2map,filelist[fileRef],map1,header=header           ; Make a map from the 1st file

map1.time=strmid(map1.time,0,12)+strmid(timeref,0,2)+':'+strmid(timeref,2,2)+':'+strmid(timeref,4,2)+'.000'

if ~valid_map(map1) then begin                     ; Test if map is ok
 pr_syntax,'DEROTATE_TEST'
 return
endif


hdr_time_ref=STRTRIM(get_map_time(map1),1)                        ; Get the time of the 1st file that will be used as reference
hdr_time_ref=map1.time


;xc=map1.xc
;yc=map1.yc


FOR i=0,(n_elements(filelist)-1) DO BEGIN         ; For 1 to the n_element of the filelist1

  fits2map,filelist[i],map,header=header          ; Make a map for each file in the list
  
  ;new_time=get_map_time(map)
  
 ;drot_xy,xc,yc,new_time,cur_time,xr,yr
  
  
  ;drot_test=drot_map(map,time=cur_time,CENTER=[xr,yr],/keep_limb)            ; Differentially rotate the image in the time used as reference (1st file)
  drot_test=drot_map(map,ref_map=map1,/keep_limb)
  

  new_filename=new_dir+filenames[i]       ; Filename + path

  ;map2fits,drot_test,new_filename                  ; Save the Differentially rotate image in the 'new_filename'
  
  FXADDPAR, HEADER, 'DROTTIME', hdr_time_ref,AFTER='DATE-OBS'
  writefits,new_filename,drot_test.data,header

ENDFOR



;------------------------------------------------------------------------
; Read all the de-rotated images in a data cube and save the Median model
;------------------------------------------------------------------------

;cd, 'C:\StereoDatatest1\DeroteData'                 ; Path of de-rotated images
;
;filelist2=file_list('.','*.fts')                   ; Read the file list of de-rotated images
;
;data_cube=sccreadfits(filelist2,hdrs)              ; Read all the de-rotated images into a data cube
;
;data_median=median(data_cube, dimension=3)         ; Median model of the data cube
;
;filename_median='file_median.fts'                  ; Name of the FITS file within the Median model
;
;writefits,filename_median,data_median              ; Save Median model


return & end