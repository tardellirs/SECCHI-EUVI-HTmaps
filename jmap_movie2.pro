pro jmap_movie2, jmap, cube, header_time, xrange, yrange, xaxis, yaxis, xout1,yout1,t0size1, format=format

  ;jmap=hist_equal(jmap,per=1)
  nfiles=n_elements(header_time)
  ;nfiles=10
  xs=800 & ys=300
  if not exist(format) then format=1
  if format eq 2 then begin
    xs=ys*1.2
  endif
  
  yscube=ys
  cube_im=bytarr(xs,yscube,nfiles)
  
  xaxiscube=findgen(5)*(xrange[1]-xrange[0])/4+xrange[0]
  yaxiscube=findgen(5)*(yrange[1]-yrange[0])/4+yrange[0]
  
  
  
  minmag=(mean(cube)-5*sigma(cube)) > min(cube)
  maxmag=(mean(cube)+5*sigma(cube)) < max(cube)
  ;minmag=min(cube)
  ;maxmag=max(cube)/2
  
  titlecube=STRMID(header_time,0,10) + ' ' + STRMID(header_time,11,11)
  xtickname=STRTRIM(string(xaxiscube,FORMAT='(i6)'),1)
  ytickname=STRTRIM(string(yaxiscube,FORMAT='(i6)'),1)
  
  ;slitpx.xclick=[slitpx.xclick[0],slitpx.xclick[n_elements(slitpx.xclick)-1]]
  ;slitpx.yclick=[slitpx.yclick[0],slitpx.yclick[n_elements(slitpx.yclick)-1]]
  
  for i=0,nfiles-1 do begin
    window,0,xs=xs,ys=yscube
    PLOT_IMAGE,cube[*,*,i], $
      MIN=minmag, MAX=maxmag, $
      BACKGROUND = 255, COLOR = 0, XTICKS=4, YTICKS=4, $
      XTICKNAME=xtickname, YTICKNAME=ytickname, $
      TITLE=titlecube[i], $
      XTITLE='X (arcsecs)', YTITLE='Y (arcsecs)',/smooth
    plots,xout1[0:t0size1-1],yout1[0:t0size1-1],color=255,thick=1,linestyle=2
    ;plots,slitpx.xclick,slitpx.yclick,color=255,psym=7,symsize=0.7
    cube_im[*,*,i]=tvrd()
  endfor
  
  
  
  
  hour=long(STRMID(header_time,11,2))
  min=long(STRMID(header_time,14,2))
  sec=long(round(float(STRMID(header_time,17,5))))
  time_sec=long(sec+min*60.+hour*60.*60.)
  time_sec=time_sec-min(time_sec)
  szjmap=size(jmap)
  jmapin=time_sec*(szjmap[1]-1)/max(time_sec)
  ;jmapin=[jmapin,jmapin[n_elements(jmapin)-1]]
  
  
  jmapnew=fltarr(szjmap[1],szjmap[2],nfiles)
  for i=0,nfiles-1 do begin
    jmapnew[*,*,i]=jmap
    jmapnew[jmapin[i],*,i]=jmapnew[jmapin[i],*]*1.8+mean(jmap)/5
    ;jmapnew[jmapin[i+1],*,i]=jmapnew[jmapin[i],*]
  endfor
  
  jmapsize=size(jmap)
  proportion=3
  if format eq 2 then proportion=2
  if jmapsize[1] gt jmapsize[2] then begin
    scalefactor1=100
    scalefactor2=(jmapsize[1]/float(jmapsize[2]))*1/proportion*100.
  endif else begin
    scalefactor1=(jmapsize[2]/float(jmapsize[1]))*proportion*100.
    scalefactor2=100
  endelse
  
  
  
  ysjmap=xs/3
  if format eq 2 then begin
    xs_old=xs
    xs=ys*2
    ysjmap=ys
  endif
  
  xtickname=STRTRIM(string(xaxis,FORMAT='(f6.2)'),1)
  ytickname=STRTRIM(STRING(yaxis,FORMAT='(f6.2)'),1)
  
  ;dotline=findgen((size(jmap))[2])*scalefactor2
  jmap_im=bytarr(xs,ysjmap,nfiles)
  for i=0,nfiles-1 do begin
    window,0,xs=xs,ys=ysjmap
    PLOT_IMAGE,jmapnew[*,*,i], $
      MIN=min(jmap), MAX=max(jmap), $
      BACKGROUND = 255, COLOR = 0, $
      xticks=5, xtickname=xtickname, $
      yticks=6, ytickname=ytickname, $
      XTITLE = 'Time (min)', YTITLE = 'Distance (arcsec)', /smooth, $
      scale=[scalefactor1,scalefactor2]
      ;plots,jmapin[i+1]*scalefactor1,dotline,color=255,thick=2,linestyle=1
    jmap_im[*,*,i]=tvrd()
  endfor
  
  
  if format eq 1 then begin
    fullfig=bytarr(xs,yscube+ysjmap,nfiles)
    for i=0,nfiles-1 do begin
      fullfig[*,*,i]=[[jmap_im[*,*,i]],[cube_im[*,*,i]]]
    endfor
  endif else begin
    fullfig=bytarr(xs_old+xs,ys,nfiles)
    for i=0,nfiles-1 do begin
      fullfig[*,*,i]=[cube_im[*,*,i],jmap_im[*,*,i]]
    endfor
  endelse
  
  generic_movie,fullfig
  ;stop
END