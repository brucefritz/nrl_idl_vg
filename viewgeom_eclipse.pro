;***************************************************************************************
;+
; NAME: viewgeom_eclipse.pro
;
; PURPOSE:
;             This program computes some key variables for computing
;             viewing geometry information useful to a variety of
;             programs and analysis techniques.
;
;             This requires inputs for
;             (1) the HRT (high rate telemetry) structure of VVIPRE science data to get scan angle
;             (2) the ISS GSE packet APID 674, format 6 with USGNC info
;             (3) an approximate TLE set that will be fitted to the USGNC data
;
; CALLING SEQUENCE:
;   vg = viewgeom_eclipse(hrt, gse, tle)
;
; INPUT ARGUMENTS: Variable Name		Type		Description
;                 hrt    		structure array  relevant information of VVIPRE science data, esp time and mirror angle
;                 gse       structure array  ISS location and attitude
;                 tle           structure    TLE appropriate for the ISS ephmeris epoch
;
; OUTPUTS: Variable Name         Type                    Description
;		returns structure of s/c and instrument location, look vectors, tangent points,
;       slit orientation, time and solar and all kinds of viewing information.
;	-- The possible values copied into the structure are given by the string values,
;	-- viewing geometry names, as follows:
;		''			NO ACTION TAKEN
;		'met'		MET (mission elapsed time in ticks)
;		'met_hp'	MET high precision double
;		'year'
;		'month'
;	    'day'
;    	'ut_sec'	universal time second of day
;	  	'look_vec'	as array of x,y,z vector
;	 	  'look_vec_x'
;		  'look_vec_y'
;    	'look_vec_z'
;	    'look_az_deg'
;    	'look_za_deg'
;    	'look_ra_deg'
;    	'look_dec_deg'
;	    'lookqual'        ;4-bit quality of look direction Y,P,R: 0=>best
;    	'sc_x'
;    	'sc_y'
;    	'sc_z'
;    	'sc_lat'
;    	'sc_lon'
;    	'sc_alt'
;    	'sc_re'
;	    'sc_loc_vert'	as array of x,y,z vector
;	    'sc_loc_vert_x'
;	    'sc_loc_vert_y'
;	    'sc_loc_vert_z'
;	    'sc_loc_zen_dot_radial'
;	    'tanpt_vec'		as array of x,y,z vector
;	    'tanpt_vec_x'
;	    'tanpt_vec_y'
;	    'tanpt_vec_z'
;	    'tanpt_lat'
;	    'tanpt_lon'
;   	'tanpt_alt'
;	  	'tanpt_re'
;	 	'tanpt_loc_vert'		as array of x,y,z vector
;	    'tanpt_loc_vert_x'
;	    'tanpt_loc_vert_y'
;	    'tanpt_loc_vert_z'		as array of x,y,z vector
;    	'sun_vec_x'
;    	'sun_vec_y'
;    	'sun_vec_z'
;    	'gmst_hrs'
;    	'sc_sza'
;    	'tanpt_sza'
;    	'ss_lat'
;    	'ss_lon'
;    	'slit_az_deg'
;    	'slit_roll_deg'
;
; INPUT KEYWORDS:
;
; OUTPUT KEYWORDS:
;
; COMMON BLOCKS:
;
; CALLED ROUTINES:
;	met2ut
;
; SIDE EFFECTS:
;
; RESTRICTIONS: This IDL routine is to be distributed only by the Naval
;				Research Laboratory, Code 7607.
;
; NOTES:
;
; MODIFICATION HISTORY:
;		Version 0.1: Scott Alan Budzien, NRL Code 7607, 12/04/97
;		Version 0.9: Scott Alan Budzien, NRL Code 7607, 12/13/99
;		Version 1.0: Clyde Fortna, SFA inc. @NRL Code 7607, 2002/06/26
;                        split into generic module with specific wrapper functions
;		Version 1.1: Clyde Fortna, SFA inc. @NRL Code 7607, 07/10/2002
;                        IDL5.5 array/structure update put reform() around structure
;       	Version 1.2: Scott Budzien, NRL Code 7607
;                        replaced argos references, added velocity estimation section,
;                        added ssuli x,y,z to output structure
;               Version 1.3: Scott Budzien, NRL Code 7607, 07/06/2005
;                        Implemented new methods for velocity vector estimates as
;                        helper routines.  Handle case of unsorted input structure.
;-
;***************************************************************************************
;*************************************************************************
;
; This function interpolates the USGNC data from the STP-H9 GSE, which has time gaps,
; onto the time grid of the VVIPRE HRT data.
;
; We could try this a different way to improve the positional errors...
; Possibly in chunks of 0.25 or 0.1 day.  However, we would need to figure 
; out a way to merge the segments smoothly.  For now, we simply use one
; TLE for the whole day, optimized for a best fit.  (This will not handle
; maneuvers.)
;
; Here there are three approaches to generating the ISS position, which is
; for computing the ISS local LVLH axes.  The problem is that USGNC times do
; not align with measurement times and there are time gaps in USGNC.
;
; /use_sgp4
; (1) The simplest is using SGP4 modeled positions only.  Although this
;     is roughly correct and continuous, there are positional differences 
;     with the USGNC reported values.
; /scaled_merge
; (2) The second approach is interpolating USGNC values and patching the 
;     USGNC missing values directly with scaled SGP4 values.  This shows 
;     evidence of errors at the edges of the patches and errors up to 10 km.
; /additive_merge
; (3) The third approach is interpolating the USGNC values, but the patching
;     of USGNC missing values use interpolated additive differences from
;     the SGP4 values, not multiplicatively scaled values.
;
; OUTPUT QUALITY FLAGS:
;     spv_flag= 0 USGNC interpolated values
;               1 fitted SGP4 ephemeris
;               2 additively adjusted SGP4 ephemeris
;               4 multiplicatively adjusted SGP4 ephemeris
;
;*************************************************************************
;
FUNCTION VGE_CHECK_LRT_ISS_TIME_OVERLAP, time_lrt, time_iss
; Mapping where the USGNC and LRT time arrays overlap
  tmin_ = MIN([MIN(time_iss), MIN(time_lrt)])
  tmax_ = MIN([MAX(time_iss), MAX(time_lrt)])  ; ??? MAX ???
  t_ = LINDGEN(tmax_-tmin_+1)+tmin_        ;sets up a grid 1/sec for maximum time range
  
  f_ = REPLICATE(0L,tmax_-tmin_+1)         ;flag array for USGNC and HRT times
  i = wave2bin(t_,time_iss)                ;indices for USGNC times
  f_[i] = 1                                ;mark USGNC times
  i = wave2bin(t_,time_lrt)                ;indices for HRT times
  f_[i] += 2                               ;mark HRT times
  w_ = WHERE(f_ EQ 3,nw_)
  return, nw_ LT .8*N_ELEMENTS(time_lrt)
END

FUNCTION VGE_CONVERT_GPS_TIME_TO_SGP4_EPH, gps_time
  ; T_EPH for SGP4 --> yyyyddd.ffff
  jd_time = gpscumulative2julian(gps_time)
  ut      = julian2ut(jd_time)
  return, ut.year*1d3 + ut.dayofyear + ut.secondofday/86400d0
END

FUNCTION VGE_CHECK_LRT_TLE_OVERLAP, time_lrt, tle
  ; Convert the time in GPS sec to what SGP4 expects: yyyyddd.ffff and check the TLE epoch
  t_eph_out = VGE_CONVERT_GPS_TIME_TO_SGP4_EPH(time_lrt)
  t_eph_median =  ABS(MEDIAN(t_eph_out))
  t_eph_mcentury = DOUBLE(FIX(t_eph_median/1d5))

  return, ((t_eph_median-t_eph_mcentury*1d5)-tle.epoch_time) GT 7
END

FUNCTION VGE_FIND_COARSE_GRID_INDEX, vec, resolution
  ; Create a grid at specified resolution
  vec_min = MIN(vec, MAX=vec_max)
  nt = FIX((vec_max-vec_min)/DOUBLE(resolution)+1, TYPE=3) ; Updated to use LONG integer
  vec_grid_res = vec_min + INDGEN(nt)*DOUBLE(resolution)
  ;
  ; Find the nearest index to each time
  iv_res = wave2bin(vec, vec_grid_res)
  vec_res = vec[iv_res]
  ; Eliminate duplicate indices that may be associate with gaps
  dv_res = vec_res - vec_grid_res
  w = where(abs(dv_res) le .5, nres)
  
  return, iv_res[w]
END

FUNCTION VGE_GENERATE_INITIAL_ORBIT_MODEL, time_iss, iss_posn, tle
  ;----------------------------------------------------------------------------
  ; Function to create an initial model (M) of orbit parameters based on the
  ; ISS time / position. Use SGP4 for use in fitting to fill in gaps where
  ; ISS is bad or missing.
  ;   Inputs: ISS Time and Postion
  ;           TLE
  ;   Output: M [struct] containing orbit information
  ;----------------------------------------------------------------------------
  ;
  ; Pulling out the distinct variables from the structure
  it30a= VGE_FIND_COARSE_GRID_INDEX(time_iss, 30d0)
  n30a = N_ELEMENTS(it30a)
  t30a = REFORM(time_iss[it30a])
  x30a = REFORM(DOUBLE(iss_posn[0, it30a]))
  y30a = REFORM(DOUBLE(iss_posn[1, it30a]))
  z30a = REFORM(DOUBLE(iss_posn[2, it30a]))
  ;
  ; Convert the time in GPS sec to what SGP4 expects: 
  t_eph_30a = VGE_CONVERT_GPS_TIME_TO_SGP4_EPH(t30a)
  ;
  ; Determine orbital phase of the data relative to the TLE orbit.
  ; If the TLE epoch is very close to the data time, the shift will be minimal
  ; Otherwse the shift could be substantial.
  ; Perform a coarse grid check for mean_anomaly for just a couple of orbits.
  ; Does the epoch fall within the dataset?
  ; Identify the epoch index, and orbit befoe the epoch and an orbit after the epoch
  Ptle=2d0*!dpi*60/tle.mean_motion                 ;Period in seconds
  iepoch = wave2bin(t_eph_30a-2d6,tle.epoch_time+[-Ptle,0,Ptle]/86400d0)
  ilowhigh = wave2bin(t30a,[min(t30a)+Ptle*2,max(t30a)-Ptle*2])
  if (iepoch[1] lt ilowhigh[0]) then begin
    imin = 0
    imax = ilowhigh[0]
  endif else if (iepoch[2] ge n30a-1) then begin
    imin = ilowhigh[1]
    imax = n30a-1
  endif else begin
    imin = iepoch[0]
    imax = iepoch[2]
  endelse
  ;
  ; In SGP4 all the times are computed relative to the epoch, which can mean
  ; that a change in the period or mean anomaly at the epoch is easily computed,
  ; but becomes complicated at start and end times that may be far from the
  ; anomaly.  The following fit only adjusts the mean anomaly and period.
  ; The mean anomaly adjusts where along the orbit the spacecraft is at the epoch
  ; and the period adjusts both the altitude and time period.
  ; Note that the ascending node will change with each orbit due to precession,
  ; but this affects x and y, not z.
  ;
  ; Parameter  Scalar              Profile         Variable
  ;-----------------------------------------------------------------------
  M = {model, $
    descriptor: ['Mean Anomaly (rad)','Period(sec)','Inclination (rad)',  'Asc Node (rad)',         'Eccen',       'Arg Perigee'],$
    minit:      [    tle.mean_anomaly,         Ptle,    tle.inclination,tle.ascending_node,tle.eccentricity,tle.argument_perigee],$
    m:          dblarr(6),$
    lowerlimit: dblarr(6),$
    upperlimit: dblarr(6),$
    ineq_flag:  [                 'b',          'b',                'b',               'b',             'b',                 'b'],$
    change:     [                   1,            1,                  1,                 1,               1,                   1],$
    covest:     0d0,$
    covfit:     dblarr(6,6),$
    tlefit:     tle,$
    tleinit:    tle}
  M.m = M.minit
  M.lowerlimit[0]  = 0.
  M.upperlimit[0]  = 2*!dpi
  M.lowerlimit[1]  = Ptle*.98
  M.upperlimit[1]  = Ptle*1.02
  M.lowerlimit[2]  = M.minit[2]*0.97
  M.upperlimit[2]  = M.minit[2]*1.03
  M.lowerlimit[3]  = 0.
  M.upperlimit[3]  = 2*!dpi
  M.lowerlimit[4]  = (M.minit[4]-0.02)>0d0
  M.upperlimit[4]  = (M.minit[4]+0.02)
  M.lowerlimit[5]  = (M.minit[5]-0.1)>0.
  M.upperlimit[5]  = (M.minit[5]+0.1)<(2*!dpi)
  ;
  ; Initial DINV fit ignores X,Y and looks only at Z to tweak inclination and period
  ;
  nw =imax-imin+1
  D0 = {d  :[x30a[imin:imax],y30a[imin:imax],z30a[imin:imax]],$    ;x, y, and z position
    din:[x30a[imin:imax],y30a[imin:imax],z30a[imin:imax]],$
    t  :t_eph_30a[imin:imax],$                                 ;time for each position
    covdata: [replicate(1d10,nw*2),replicate(1d,nw)]}          ;1-km uncertainty in z (1 km2 covariance)
  Mz = M
  ;
  ; Perform a coarse grid check for the mean anomaly
  nshift = 300
  ma = dindgen(nshift)/double(nshift-1d0)*2d0*!dpi
  chi2 = ma*0.
  for i=0,nshift-1 do begin
    Mz.m[0] = ma[i]
    sgp4_forward,Mz,D0,/nostop
    diff = (d0.d-d0.din)
    chi2[i] = total(diff*diff/d0.covdata)/nw
    ;    print,i,ma[i],chi2[i]
  endfor
  w = where(chi2 eq min(chi2))
  Mz.minit[0] = ma[w[0]]
  Mz.lowerlimit[0]  = (Mz.minit[0]-0.04)>0d0
  Mz.upperlimit[0]  = (Mz.minit[0]+0.04)<(2d0*!dpi)
  Mz.m = Mz.minit
  Mz.change[0] = 0          ; Ref to M.m, unsure this does anything since it gets changed shortly after
  sgp4_forward,Mz,D0,/nostop
  ;
  ; Now that an approximate mean anomaly is estimated, now do the full z-fit
  ; for mean anomaly, period, and inclination for the full z-dataset
  D1 = {d  :[x30a,y30a,z30a],$    ;x, y, and z position
    din:[x30a,y30a,z30a],$
    t  :t_eph_30a,$                                         ;time for each position
    covdata: [replicate(1d10,n30a*2),replicate(1d,n30a)]}   ;1-km uncertainty in z (1 km2 covariance)

  Mz.change = [1,1,1,0,0,0]
  dinv,'sgp4_forward',Mz,D1,/no_n,deriv_eps=0.0003d0
  ;
  ; From the z-fit we can assume that the inclination and period are constrained.
  ; The argument of the ascending node w (which affects the x-y position most) is
  ; probably already pretty close, since the inclination and precession are almost
  ; always constant.
  ;
  M.minit = Mz.m                ;initialize M based upon latest fit
  M.change = [1,0,0,1,1,1]      ;holding period and inclination constant
  D2 = D1
  D2.covdata = replicate(1d0,3*n30a)                             ;1-km uncertainty in x,y,z (1 km2 covariance)
  dinv,'sgp4_forward',M,D2,/no_n,deriv_eps=0.0003d0
  ;
  ; Some diagnostics - how good the coarse fit matches the input
  ;ds = sqrt(total(reform((d2.din[0:3*n30a-1]-d2.d[0:3*n30a-1])^2,n30a,3),2,/double))
  ;rin = sqrt(total(reform((d2.din[0:3*n30a-1])^2,n30a,3),2,/double))
  ;r = sqrt(total(reform((d2.d[0:3*n30a-1])^2,n30a,3),2,/double))
  ;dr = r-rin

  RETURN, M
END

FUNCTION VGE_GENERATE_SPV, time_hrt, time_iss, sc_posn, timepos
;-----------------------------------------------------------------------
;
; Spacecraft position vector
;
; Here there are three approaches to generating the SPV, which is essential
; for computing the local LVLH axes.  The problem is that USGNC times do
; not align with measurement times and there are time gaps in USGNC.
;
; (1) The simplest is using SGP4 modeled positions only.  Although this
;     is roughly correct and continuous, there are positional differences
;     with the USGNC reported values.  This is not preferred.
; (2) The second approach is interpolating USGNC values and patching the
;     USGNC missing values directly with scaled SGP4 values.  This shows
;     evidence of errors at the edges of the patches and errors up to 10 km.
;     Found not to have sufficient accuracy.
; (3) The third approach is interpolating the USGNC values, but the patching
;     of USGNC missing values use interpolated additive differences from
;     the SGP4 values, not multiplicatively scaled values.  This seems to
;     be best and is the default.
;
  print, 'Generating SPV ... '
  xref = REFORM(DOUBLE(sc_posn[0,*]))
  yref = REFORM(DOUBLE(sc_posn[1,*]))
  zref = REFORM(DOUBLE(sc_posn[2,*]))
  rref = SQRT(xref*xref + yref*yref + zref*zref)
  
  nt_hrt = N_ELEMENTS(time_hrt)
  nt_iss = N_ELEMENTS(time_iss)
  spv_flag = REPLICATE(0B,nt_hrt)
  ;
  ; Here is where we need to carefully calculate the primary GEI coordinates
  ; from USGNC values, not use what we estimated earlier from the SGP4 model.
  ; These accurate SC position estimates are needed for both the additive and
  ; multiplicative fitted patches from SGP4
  ;
  ; First define where there is good coverage between time_hrt and time_iss (USGNC_sec)
  ; Denote gaps with 1 flags
  ;
;  svm_coarse = time_hrt[VGE_FIND_COARSE_GRID_INDEX(time_hrt, 1)]
;  nt_crs = N_ELEMENTS(svm_coarse)
  mask = REPLICATE(1.,nt_hrt)
  
  w = where(time_iss GE time_hrt[0] AND time_iss LE time_hrt[nt_hrt-1])
;  w = where(time_iss ge (svm_coarse[0]-0.5) and time_iss le (svm_coarse[nt_crs-1]))
;  i_mask = wave2bin(svm_coarse,time_iss[w])
;  mask[i_mask] = 0
  
  gapl = []
  gapr = []
  gap_time_count = 0.
  FOR i=0,nt_iss-2 DO BEGIN
    IF time_iss[i+1]-time_iss[i] GT 4.1 THEN BEGIN 
      gapl = [gapl, LONG(i)]
      gapr = [gapr, LONG(i+1)]
      gap_time_count += time_iss[i+1]-time_iss[i] 
    ENDIF
  ENDFOR
  print, ' ... found total of ' + STRTRIM(gap_time_count, 2) + ' seconds missing'
  ngapl = N_ELEMENTS(gapl)
  ngapr = N_ELEMENTS(gapr)
  
;  IF ngapl EQ 0 THEN gapl = [0]
;  IF ngapr EQ 0 THEN gapr = [-1]
  ; reconcile the number of gaps
  csum = TOTAL(mask,/cumul)
  n = N_ELEMENTS(csum)
  gapl = where(csum[0:n-3] eq csum[1:n-2] and csum[2:n-1] gt csum[1:n-2],ngapl)
  gapr = where(csum[0:n-3] lt csum[1:n-2] and csum[2:n-1] eq csum[1:n-2],ngapr)
  if ((ngapr eq 1) and (ngapl eq 0)) then begin
    ngapl = 1
    gapl = [0]
  endif else if ((ngapr ge 1) and (gapr[0] lt gapl[0])) then begin
    ngapl++
    gapl = [0,gapl]
  endif
  if ((ngapl eq 1) and (ngapr eq 0)) then begin
    ngapr = 1
    gapl = [n-3]
  endif else if ((ngapl ge 1) and (gapr[ngapr-1] lt gapl[ngapl-1])) then begin
    ngapr++
    gapr = [gapr,n-2]
  endif
  gapl += 1     ;these adjust position, since we looked at the derivatives/changes
  gapr += 2     ; this region is made 1 element larger on each side, to accomodate
  ; differences going from time_iss to time_hrt
  gapl = gapl>0
  gapr = gapr<(nt_hrt-1)
  mask[gapl] = 1            ;marking the mask for the expanded size of gaps
  mask[gapr] = 1
  
  ;
  ; Clear out any small gaps--these will just be interpolated over
;  smallgap = where((gapr-gapl) le 4,nsg,complement=biggap)
;  if (nsg gt 0) then begin
;    for i=0,nsg-1 do mask[gapl[smallgap[i]]:gapr[smallgap[i]]] = 0
;    gapl = gapl[biggap]
;    gapr = gapr[biggap]
;  endif
;  
  ; Expand the gaps by one element on each side
  ngaps = N_ELEMENTS(gapl)
  
  ;
  ; Simple linear interpolation. Where data are present, this does a pretty good
  ; job, even where small gaps occur.  Big gaps are the problem.
  ; First interpolate the whole array.  The gaps will be handled differently later.
  ;
  time_iss2 = time_iss-time_iss[0]
  time_hrt2 = time_hrt-time_iss[0]
  
  radial_x_all = intrpl(time_iss2,xref,time_hrt2,/double)
  radial_y_all = intrpl(time_iss2,yref,time_hrt2,/double)
  radial_z_all = intrpl(time_iss2,zref,time_hrt2,/double)
  
  ;
  ; The SGP4 fit is pretty good, on average, but particular sections may be
  ; slightly offset from USGNC.  Look at the average difference between SGP4
  ; and USGNC where both are present, at the edges of gaps, then infer the
  ; predicted difference between SGP4 and USGNC within the gaps.
  ;
  ; If we look at the difference between the interpolated values and the SGP4 model
  ; the difference is fairly smooth (especially on a +/-10 km scale), but with gaps and
  ; spikes.  A median filter can eliminate the small spikes.
  ;
  diff_x = radial_x_all - timepos.x_km      ;difference between true and calculated positions
  diff_y = radial_y_all - timepos.y_km
  diff_z = radial_z_all - timepos.z_km
  diffm_x = median(diff_x,61)               ;get rid of spikes
  diffm_y = median(diff_y,61)
  diffm_z = median(diff_z,61)

  wbad = where(mask eq 1,complement=wgood)  ;only places with bad USGNC values
  diff_x[wbad] = diffm_x[wbad]              ;substitute median-filtered values at bad spots
  diff_y[wbad] = diffm_y[wbad]
  diff_z[wbad] = diffm_z[wbad]
  spv_flag[wbad] = 2                        ;denote substitution of median
  
  unc = 0.1 + mask*1e7                      ;deemphasize the gaps by increasing uncertainty
  ;
  ; Make output arrays based on USGNC interpolations
  ;
  radial2_x_all = radial_x_all
  radial2_y_all = radial_y_all
  radial2_z_all = radial_z_all
  ;
  ; For each gap, fit the background difference between SGP4 and USGNC
  ; linearly on either side of the gap.
  ;
;  for i=0,ngaps-1 do begin
;    i2 = (gapl[i]-1)>0
;    i3 = (gapr[i]+1)<(nt_hrt-1)
;    i0 = (i2-180)>0            ;refers to time_hrt2, where index is
;    i1 = (i3+180)<(nt_hrt-1)
;    w = where(mask[i0:i1] eq 0,nw)
;    if (nw eq 0) then begin
;      i0 = (i2-250)>0            ;refers to time_hrt2, where index is
;      i1 = (i3+250)<(nt_hrt-1)
;      index_ = [lindgen(i2-i0+1)+i0,lindgen(i1-i3+1)+i3] ;indices outside the gap
;    endif else $
;      index_ = w + i0
;    t_ = time_hrt2[index_]
;    dx_ = diff_x[index_]
;    dy_ = diff_y[index_]
;    dz_ = diff_z[index_]
;    u_ = unc[index_]
;    cx = linfit_nrl(t_,dx_,u_,1)
;    cy = linfit_nrl(t_,dy_,u_,1)
;    cz = linfit_nrl(t_,dz_,u_,1)
;    ;
;    ; The baseline USGNC-SGP4 background in the gap
;    ;
;    t_ = time_hrt2[i2:i3]                                  ;times inside the gap
;    gapdx = cx[0] + cx[1]*t_
;    gapdy = cy[0] + cy[1]*t_
;    gapdz = cz[0] + cz[1]*t_
;    ;
;    ; Estimating the baseline USGNC-SGP4 within the gap
;    ;
;    radial2_x_all[i2:i3] = timepos[i2:i3].x_km + gapdx
;    radial2_y_all[i2:i3] = timepos[i2:i3].y_km + gapdy
;    radial2_z_all[i2:i3] = timepos[i2:i3].z_km + gapdz
;    spv_flag[wbad] += 4
;
;    ; temporary output for diagnostics
;            !p.multi=[0,3,2]
;            plot,time_hrt2[i0:i1],radial_x_all[i0:i1],xtit='Second of Day',ytit='X position (km)'
;            oplot,time_hrt2[i0:i1],timepos[i0:i1].x_km,col=150
;            oplot,time_hrt2[i0:i1],radial2_x_all[i0:i1],col=60
;            plot,time_hrt2[i0:i1],radial_y_all[i0:i1],xtit='Second of Day',ytit='Y position (km)'
;            oplot,time_hrt2[i0:i1],timepos[i0:i1].y_km,col=150
;            oplot,time_hrt2[i0:i1],radial2_y_all[i0:i1],col=60
;            plot,time_hrt2[i0:i1],radial_z_all[i0:i1],xtit='Second of Day',ytit='Z position (km)'
;            oplot,time_hrt2[i0:i1],timepos[i0:i1].z_km,col=150
;            oplot,time_hrt2[i0:i1],radial2_z_all[i0:i1],col=60
;            plot,time_hrt2[i0:i1],radial2_x_all[i0:i1]-timepos[i0:i1].x_km,yr=[-4,4],$
;                xtit='Second of Day',ytit='X position difference (km)'
;            plot,time_hrt2[i0:i1],radial2_y_all[i0:i1]-timepos[i0:i1].y_km,yr=[-4,4],$
;                xtit='Second of Day',ytit='Y position difference (km)'
;            plot,time_hrt2[i0:i1],radial2_z_all[i0:i1]-timepos[i0:i1].z_km,yr=[-4,4],$
;                xtit='Second of Day',ytit='Z position difference (km)'
;  endfor
  spv_gei = transpose([[radial2_x_all],[radial2_y_all],[radial2_z_all]])
  RETURN, spv_gei
END

FUNCTION VGE_GENERATE_LVLH_AXIS_GEI, time_iss, time_hrt, spv_gei, qinert, qLVLH, M
  PRINT, 'Generating LVLH Axes ... '
  ;-------------------------------------------------------------------------------------
  ;   The spacecraft position radial direction defines the negative LVLH z-axis
  ;   (in GEI).  The x-axis of LVLH is roughly in the velocity direction, but it
  ;   changes very rapidly due to orbital motion.  Evidence from 2011 suggests
  ;   that the LVLH quaternions are calculated from the inertial quaternions and
  ;   position, even though they correspond to different times.
  ;
  ;   The y-axis in LVLH corresponds to the negative orbit normal, which will
  ;   precess slightly throughout the day.  This vector will precess in x and y
  ;   linearly in azimutal angle versus time, whereas the z vector should be
  ;   roughly constant for a constant inclination.  Since the vector is nearly
  ;   constant second-to-second, any discrepancy in reporting times can be ignored.
  ;
  ;   The y-LVLH axis in GEI is determined by reported LVLH and INERTIAL quaternions.  
  ;   Although the times of quaternion reports does not exactly correspond to 
  ;   measurement times and the positional time, since y-LVLH is perpendicular
  ;   the radial direction and changes slowly, this should introduce little error.
  ;                 
  ;   x-LVLH == y-LVLH (cross) z-LVLH
  ;
  ;--------------------------------------------------------------------------------------
  q_lvlh_to_inert = QMULTIPLY(qinert, QCONJ(qLVLH))
  ylvlh_gei = APPLY_QUATERNION(q_lvlh_to_inert,[0,1d,0])
  PRINT, '    Find Coarse Grid Index ... '
  it30a= VGE_FIND_COARSE_GRID_INDEX(time_iss, 30d0)
  n30a = N_ELEMENTS(it30a)
  t30a = REFORM(time_iss[it30a])
  PRINT, '    Fit Orbit Normal ... '
  FIT_ORBIT_NORMAL,t30a,ylvlh_gei[0,it30a],M.m[1],M1a
  FIT_ORBIT_NORMAL,t30a,ylvlh_gei[1,it30a],M.m[1],M1b
  FIT_ORBIT_NORMAL,t30a,ylvlh_gei[2,it30a],M.m[1],M1c,/zfit
  PRINT, '    Forward model sinudoid model ... '
  D_all = {d  :REPLICATE(0d, N_ELEMENTS(time_hrt)),$             ;x, y, or z position
    t  :time_hrt - t30a[0]     }
  SINUSOID_FORWARD,M1a,D_all         ;compute based on x
  ylvlh_gei_fitx = D_all.d
  SINUSOID_FORWARD,M1b,D_all         ;compute based on y
  ylvlh_gei_fity = D_all.d
  SINUSOID_FORWARD,M1c,D_all         ;compute based on z
  ylvlh_gei_fitz = D_all.d
  ylvlh_gei_fitr = SQRT(ylvlh_gei_fitx*ylvlh_gei_fitx + $
                        ylvlh_gei_fity*ylvlh_gei_fity + $
                        ylvlh_gei_fitz*ylvlh_gei_fitz   )
  
  ;
  ; A potential flaw here:  We have a good expression for the orbit normal.  However,
  ; if we use the SGP4 fit to get the radial positions, we can introduce on the order
  ; of +/-0.10 deg errors or +/-1.7 mrad between the position vector and the orbit
  ; normal.  (Maybe this relates to slight differences in fitted vs actual position,
  ; but this would require 12 km of error, and there is no evidence of this.)  This error
  ; oscillates with orbit, but is not in phase with radius; it is most closely anti-
  ; correlated with Z actuals.  Is there an orbital effect that causes an instantaneous
  ; orbit normal to oscillate about the mean?  Oblateness?
  ;
  ; We will continue to use this reference to interpolate pitch, roll, and yaw, but we
  ; for obtaining the instantaneous LVLH vectors in GEI we will try to use the USGNC
  ; values later.
  ;

  ;
  ;; The LVLH coordinate system is defined with +Xa in roughly the velocity vector direction,
  ;; +Za in the nadir, and +Ya in the negative orbit normal.  Expressing these changing
  ;; unit axis vector (Xa, Ya, Za) in GEI coordinates are
  ;;

  ; Make a unit vector
  radial2_r_all = SQRT(spv_gei[0,*]*spv_gei[0,*]+$
    spv_gei[1,*]*spv_gei[1,*]+$
    spv_gei[2,*]*spv_gei[2,*])
  radial_x_all = spv_gei[0,*]/radial2_r_all
  radial_y_all = spv_gei[1,*]/radial2_r_all
  radial_z_all = spv_gei[2,*]/radial2_r_all

  ya_lvlh_x = ylvlh_gei_fitx/ylvlh_gei_fitr
  ya_lvlh_y = ylvlh_gei_fity/ylvlh_gei_fitr
  ya_lvlh_z = ylvlh_gei_fitz/ylvlh_gei_fitr
  ;This should be a pretty good value that matches USGNC
;  za_lvlh = {x:-radial_x_all, y:-radial_y_all, z:-radial_z_all}
  za_lvlh_x = -radial_x_all            ;This should be a pretty good value that matches USGNC
  za_lvlh_y = -radial_y_all
  za_lvlh_z = -radial_z_all

  ;pslaser,'Fig39 LVLH-y data and fit',/col
  ;!p.multi = [0,1,3]
  ;plot,t30a,ylvlh_gei[0,it30a],xs=1,ps=4,symsize=.4,$
  ;    ytit='LVLH Y-axis GEI-x component',xtit='GPS Time (s)',charsize=1.5
  ;oplot,time_hrt,ya_lvlh_x,col=150,thick=2
  ;plot,t30a,ylvlh_gei[1,it30a],xs=1,ps=4,symsize=.4,$
  ;    ytit='LVLH Y-axis GEI-y component',xtit='GPS Time (s)',charsize=1.5
  ;oplot,time_hrt,ya_lvlh_y,col=150,thick=2
  ;plot,t30a,ylvlh_gei[2,it30a],xs=1,ps=4,symsize=.4,$
  ;    ytit='LVLH Y-axis GEI-z component',xtit='GPS Time (s)',charsize=1.5
  ;oplot,time_hrt,ya_lvlh_z,col=150,thick=2
  ;store
  ;!p.multi=0

  ;;
  ;; The instantaneous Xa axis is best defined by the cross product of Ya x Za
  ;;
  xa_lvlh_x = ya_lvlh_y*za_lvlh_z - za_lvlh_y*ya_lvlh_z
  xa_lvlh_y = ya_lvlh_z*za_lvlh_x - za_lvlh_z*ya_lvlh_x
  xa_lvlh_z = ya_lvlh_x*za_lvlh_y - za_lvlh_x*ya_lvlh_y
  xa_lvlh_r = SQRT(xa_lvlh_x*xa_lvlh_x+$
    xa_lvlh_y*xa_lvlh_y+$
    xa_lvlh_z*xa_lvlh_z)
  xa_lvlh_x /= xa_lvlh_r
  xa_lvlh_y /= xa_lvlh_r
  xa_lvlh_z /= xa_lvlh_r
  ;;
  ;; Ya should already be perpendicular to Xa, but might not be perpendicular
  ;; to Za, depending on the accuracy of the orbit normal fit and spacecraft
  ;; position fit.
  ;;
  err=ASIN(za_lvlh_x*ya_lvlh_x+za_lvlh_y*ya_lvlh_y+za_lvlh_z*ya_lvlh_z)*180d0/!dpi
  maxerr=MAX(err,min=minerr)
  rmserr=STDDEV(err)
  PRINT,"$('Orbit normal-Position error from 90=',f6.3,' to ',f6.3,' deg  rms=',f6.3,'deg')",$
    minerr,maxerr,rmserr

  ;;
  ;; Redefine to a true orthonormal set using Xa and Ya.  This
  ;; assumes that Xa may be slightly in error.  Given how gradually
  ;; Ya moves, this may not be a bad assumption.
  ;;
  za_lvlh_x = xa_lvlh_y*ya_lvlh_z - ya_lvlh_y*xa_lvlh_z
  za_lvlh_y = xa_lvlh_z*ya_lvlh_x - ya_lvlh_z*xa_lvlh_x
  za_lvlh_z = xa_lvlh_x*ya_lvlh_y - ya_lvlh_x*xa_lvlh_y
  za_lvlh_r = SQRT(za_lvlh_x*za_lvlh_x+$
    za_lvlh_y*za_lvlh_y+$
    za_lvlh_z*za_lvlh_z)
  za_lvlh_x /= za_lvlh_r
  za_lvlh_y /= za_lvlh_r
  za_lvlh_z /= za_lvlh_r

  ; Older alternative: Redefine to a true orthonormal set using Za and Xa
  ;
  ;ya_lvlh_x = za_lvlh.y*xa_lvlh_z - xa_lvlh_y*za_lvlh.z
  ;ya_lvlh_y = za_lvlh.z*xa_lvlh_x - xa_lvlh_z*za_lvlh.x
  ;ya_lvlh_z = za_lvlh.x*xa_lvlh_y - xa_lvlh_x*za_lvlh.y
  ;ya_lvlh_r = sqrt(ya_lvlh_x*ya_lvlh_x+$
  ;                 ya_lvlh_y*ya_lvlh_y+$
  ;                 ya_lvlh_z*ya_lvlh_z)
  ;ya_lvlh_x /= ya_lvlh_r
  ;ya_lvlh_y /= ya_lvlh_r
  ;ya_lvlh_z /= ya_lvlh_r


  ; Now put this back in terms of the GEI representation of the instantaneous LVLH axes
  axis_gei = {x:TRANSPOSE([[xa_lvlh_x],[xa_lvlh_y],[xa_lvlh_z]]),$
              y:TRANSPOSE([[ya_lvlh_x],[ya_lvlh_y],[ya_lvlh_z]]),$
              z:TRANSPOSE([[za_lvlh_x],[za_lvlh_y],[za_lvlh_z]])}
  PRINT, '           ... Complete'
  
  RETURN, axis_gei
END

FUNCTION XIP_M1_TIME_INTERP, hrt, xip=xip
  if STRUPCASE(xip) ne 'MIP' and STRUPCASE(xip) ne 'TIP' then begin
    print, 'Invalid instrument type! Must specify MIP or TIP'
    return, -1
  endif
  
  if STRUPCASE(xip) eq 'TIP' then begin
    txip_raw = hrt.tip_m1_runtime
    tgps_raw = hrt.tip_m1_ccsds_gps_time
  endif

  if STRUPCASE(xip) eq 'MIP' then begin
    txip_raw = hrt.mip_m1_runtime
    tgps_raw = hrt.mip_m1_ccsds_gps_time
  endif
  
  txip = DBLARR(N_ELEMENTS(txip_raw))
  tgps = DBLARR(N_ELEMENTS(tgps_raw))
  tgps0 = tgps_raw[0]
  txip0 = txip_raw[0]
  uidx = ULONG64(N_ELEMENTS(txip_raw)/10.)
  
  for i=0L, uidx-1 do begin
    for j=0,9 do begin
      txip[i*10 + j] = txip_raw[i*10 + j]
      tgps[i*10 + j] = tgps0 - txip0 + txip[i*10 + j]
;      if i gt 30000 and i lt 30005 then print, txip_raw[i*10 + j], 0.1*j
;      if i gt 30000 and i lt 30005 then print, txip_raw[i*10 + j] + 0.1*j
    endfor
  endfor
;  STOP
  return, {xip_m1_tgps:tgps, xip_m1_txip:txip}
END

FUNCTION SUVM_TIME_INTERP, hrt, xip=xip
;;; Viewer arg may be unnecessary
  if STRUPCASE(xip) ne 'MIP' and STRUPCASE(xip) ne 'TIP' then begin
    print, 'Invalid instrument type! Must specify MIP or TIP'
    return, -1
  endif
  
  if STRUPCASE(xip) eq 'TIP' then begin
    tsvm_raw = hrt.tip_suvm_system_counter
    tgps_raw = hrt.tip_suvm_ccsds_gps_time
  endif

  if STRUPCASE(xip) eq 'MIP' then begin
    tsvm_raw = hrt.mip_suvm_system_counter
    tgps_raw = hrt.mip_suvm_ccsds_gps_time
  endif
  
  tgps = DBLARR(N_ELEMENTS(tsvm_raw))
  tsvm0 = tsvm_raw[0]
  tgps0 = tgps_raw[0]
  
  for i=0,N_ELEMENTS(tgps)-1 do tgps[i] = tgps0 - tsvm0 + tsvm_raw[i]
  
  return, tgps
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;   MAIN FUNCTION
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
FUNCTION VIEWGEOM_ECLIPSE,hrt,iss,tle,xip,viewer, $
  stopit=stopit,plot=plotfl,deltaya=deltaya,deltaza=deltaza,latencysec=latencysec

;______________________________________________________________________________
; Parse the parameters
if STRUPCASE(xip) NE 'MIP' AND STRUPCASE(xip) NE 'TIP' THEN BEGIN
  print, 'Invalid instrument type! Must specify MIP or TIP'
  return, -1
endif
if STRUPCASE(viewer) NE 'LIMB' AND STRUPCASE(viewer) NE 'DISK' THEN BEGIN
  print, 'Invalid instrument type! Must specify LIMB or DISK'
  return, -1
endif
if KEYWORD_SET(latencysec) THEN latencysec = FLOAT(latencysec) ELSE latencysec = 0.
if KEYWORD_SET(deltaya) THEN deltaya = DOUBLE(deltaya) ELSE deltaya = 0.
if KEYWORD_SET(deltaza) THEN deltaza = DOUBLE(deltaza) ELSE deltaza = 0.
;______________________________________________________________________________
; Define basic time ranges for the input
if ISA(hrt) eq 0 then begin
  message,'No HRT data for '+viewer,/info
  return,-1
endif

if TAG_EXIST(iss,'usgnc_sec') eq 0 then begin
  message,'No valid ISS data',/info
  return,-1
endif
time_iss = iss.usgnc_sec          ; Formerly tref *** Updated for eclipse array
nt_iss = N_ELEMENTS(time_iss)     ; Formerly nref

if STRUPCASE(xip) eq 'TIP' then begin
  if TAG_EXIST(hrt,'TIP_M1_CCSDS_GPS_TIME') eq 0 then begin
    print, TAG_EXIST(hrt,'TIP_M1_CCSDS_GPS_TIME')
    print, ISA(hrt)
    message,'No data for TIP '+viewer,/info
    return,-1
  endif
  time_lrt = hrt.tip_m0_ccsds_gps_time
  hrt_struct = xip_m1_time_interp(hrt, xip=xip)
  time_svm = suvm_time_interp(hrt, xip=xip)
  IF STRUPCASE(viewer) EQ 'DISK' THEN angle_svm = hrt.tip_suvm_encoder_angle>0.1<45.
  IF STRUPCASE(viewer) EQ 'LIMB' THEN angle_svm = hrt.tip_suvm_encoder_angle>1.<25.

ENDIF ELSE BEGIN
  if TAG_EXIST(hrt,'MIP_M1_CCSDS_GPS_TIME') eq 0 then begin
    message,'No data for MIP '+viewer,/info
    return,-1
  endif
  time_lrt = hrt.mip_m0_ccsds_gps_time
  hrt_struct = xip_m1_time_interp(hrt, xip=xip)
  time_svm = suvm_time_interp(hrt, xip=xip)
  
;  time_hrt = hrt.mip_m1_ccsds_gps_time   ; Formerly mip_m1_gps_time   ; Formerly t_all
;  offset_i = UNIQ(hrt.mip_suvm_time)
;  IF offset_i[0] LT 9 THEN BEGIN
;    offset_t = FLOAT(offset_i[0] + 1)/10.
;    time_svm = hrt.mip_suvm_system_counter - hrt.mip_suvm_system_counter[0] + hrt.mip_suvm_time[0] - offset_t
;  ENDIF ELSE BEGIN
;    time_svm = hrt.mip_suvm_system_counter - hrt.mip_suvm_system_counter[0] + hrt.mip_suvm_time[0]
;  ENDELSE
  IF STRUPCASE(viewer) EQ 'DISK' THEN angle_svm = hrt.mip_suvm_encoder_angle>0.1<45.
  IF STRUPCASE(viewer) EQ 'LIMB' THEN angle_svm = hrt.mip_suvm_encoder_angle>1.<25.
ENDELSE

time_hrt = hrt_struct.xip_m1_tgps
nt_hrt = N_ELEMENTS(time_hrt)             ; *** Updated for eclipse array

IF latencysec GT 0 THEN FOR i=0,nt_hrt-1 DO time_hrt[i] -= latencysec

; Interpolate view_angle onto the XIP time basis
angle_hrt = INTRPL(time_svm, angle_svm, time_hrt)

;______________________________________________________________________________
; Placeholader for "Plot time example
; 
; Find a fix bad data points
; ;;; SHOULD BE UNNECESSARY NOW
;FOR i=0,nt_hrt-2 DO IF (ABS(time_hrt[i+1] - time_hrt[i]) GT 10.) THEN time_hrt[i+1] = time_hrt[i] + 0.1

; Check the time ranges  
IF VGE_CHECK_LRT_ISS_TIME_OVERLAP(time_lrt, time_iss) THEN BEGIN
    message,"The HRT and GSE time overlap is less than 80%.",/info
    return,-1
ENDIF
IF VGE_CHECK_LRT_TLE_OVERLAP(time_lrt, tle) THEN BEGIN
    message,"The epoch of the TLE differs by more than one week from the data.  Find a closer TLE.",/info
    return,-1
ENDIF

M = VGE_GENERATE_INITIAL_ORBIT_MODEL(time_iss, iss.usgnc_posn_inert, tle)

; Use forward model to calculate ISS position for each ECLIPSE time
D_data = {d :REPLICATE(0d,nt_hrt*3),$                    ;x, y, or z position
          t :VGE_CONVERT_GPS_TIME_TO_SGP4_EPH(time_hrt)} ;time for each position
sgp4_forward,M,D_data,vout=v_data

; Compute time since ascending node
launchtime = H9_LAUNCH_INFO()
jd_launch = FJULIAN(launchtime.year,launchtime.month,launchtime.day+$
    launchtime.hour_ut/24d0+launchtime.minute/1440d0+launchtime.second/86400d0)
;
timepos_init = {met_stph9:0d,x_km:0d,y_km:0d,z_km:0d}
timepos = REPLICATE(timepos_init, nt_hrt)
;
timepos.met_stph9 = (GPSCUMULATIVE2JULIAN(time_hrt) - jd_launch)*86400d0  ; high-precision MET
timepos.x_km = D_data.d[0*nt_hrt:1*nt_hrt-1]
timepos.y_km = D_data.d[1*nt_hrt:2*nt_hrt-1]
timepos.z_km = D_data.d[2*nt_hrt:3*nt_hrt-1]
;
; Use TSAN to calculate orbits numbers.  Period here is not used
tsan = TIME_SINCE_ASCENDING_NODE(timepos, period=p, asc_node=metasc, orbit=orb)
utinfo = MET2UT(timepos.met_stph9, /stph9)
spv_gei = VGE_GENERATE_SPV(time_hrt, time_iss, iss.usgnc_posn_inert, timepos)

;______________________________________________________________________________
; Placeholder for "plotfl"

; Ensure quaternions are in correct order
qinert = QSXYZ2XYZS(iss.usgnc_quat_inert)   ; Re-order ISS quaternions
qLVLH = QSXYZ2XYZS(iss.usgnc_quat_lvlh)

ypr = QUAT2YAWPITCHROLL(qLVLH)

print,"$('Mean yaw=',f7.2,', pitch=',f7.2,', roll=',f7.2,'.')",MEAN(ypr.yaw),MEAN(ypr.pitch),MEAN(ypr.roll)
print,"If these values don't look correct, verify that quaternions are in (X,Y,Z,S) format."

; Generate LVLH coordinates in GEI on instrument time basis
axis_gei = VGE_GENERATE_LVLH_AXIS_GEI(time_iss, time_hrt, spv_gei, qinert, qLVLH, M)

; Identify gaps in ISS GSE information
i_all = WAVE2BIN(iss.gps_cumul, time_hrt)

dt_all = ABS(time_hrt - iss.gps_cumul[i_all])
weight = REPLICATE(-1.,nt_hrt)
w = WHERE(dt_all LE 5.,nw)            ;all gaps 10 s or less
IF (nw GT 0) THEN weight[w] = 0.
w = WHERE(dt_all GE 15.,nw)
IF (nw GT 0) THEN weight[w] = 1.
w = WHERE(dt_all GT 5. AND dt_all LT 15.,nw)
IF (nw GT 0) THEN $
    weight[w] = STEPFUNCTION(dt_all[w],x0=10.,delta=1.2,/stepup) 
;______________________________________________________________________________
; Fix gaps in ROLL
FIT_ROLLYAW, time_iss, ypr.roll, M.m[1], M2r
Dr_all = {d  :replicate(0d,nt_hrt),$                 ;roll
          t  :time_hrt-time_iss[0]}                  ;time for each position
SINUSOID_FORWARD, M2r, Dr_all
; Combine interpolated and fitted functions
roll_all = (1d0-weight)*INTRPL(time_iss-time_iss[0],ypr.roll,Dr_all.t)+weight*Dr_all.d
roll_all = MEDIAN(roll_all,181) ; Median filter smooths out most spikes
;
; quality is an integer, the lowest 12 bits are
; bits 11-8 yaw quality, bits 7-4 pitch, bits 3-0 roll
qual = ROUND(weight*15)                             ;0=interpolated,15=fitted
PRINT, 'Finished calculating roll ...'
;______________________________________________________________________________
; Fix gaps in YAW
FIT_ROLLYAW,time_iss,ypr.yaw,M.m[1],M2y
Dy_all = {d  :replicate(0d,nt_hrt),$          ;yaw
          t  :time_hrt-time_iss[0]}           ;time for each position
SINUSOID_FORWARD, M2y, Dy_all
; Combine interpolated and fitted functions
yaw_all = (1d0-weight)*INTRPL(time_iss-time_iss[0],ypr.yaw,Dy_all.t)+weight*Dy_all.d
yaw_all = MEDIAN(yaw_all,181) ; Median filter smooths out most spikes
;
; quality is an integer, the lowest 12 bits are
; bits 11-8 yaw quality, bits 7-4 pitch, bits 3-0 roll
qual += ROUND(weight*15)*256                             ;0=interpolated,15=fitted
PRINT, 'Finished calculating yaw ...'
;______________________________________________________________________________
; Fix gaps in PITCH
n = N_ELEMENTS(ypr)
; Identify anomalous pitch values
dpdt = (ypr[1:*].pitch-ypr[0:*].pitch)/(time_iss[1:*]-time_iss[0:*])
dpdt31 = dpdt
dpdt31[15] = (ypr[31:*].pitch-ypr[0:*].pitch)/(time_iss[31:*]-time_iss[0:*])
; Identify 10-sec chunks with anomalous behavior
dphifl = FLOAT(ABS(dpdt gt 0.007))      ;empirically determined threshold on day 2023140
dphith = SMOOTH(dphifl,11)
wosc = WHERE(dphith ge 0.1,nosc)
if (nosc gt 0 ) then begin ; Reconcile # of oscilalatory sections (threshold = 1/11)
  th = 0.0909
  oscl = where(dphith[0:n-3] lt th and dphith[1:n-2] ge th,noscl)
  oscr = where(dphith[0:n-3] ge th and dphith[1:n-2] lt th,noscr)
  if ((noscr eq 1) and (noscl eq 0)) then begin
      noscl = 1
      oscl = [0]
  endif else if ((noscr ge 1) and (oscr[0] lt oscl[0])) then begin
      noscl++
      oscl = [0,oscl]
  endif
  if ((noscl eq 1) and (noscr eq 0)) then begin
      noscr = 1
      goscr = [n-2]
  endif else if ((noscl ge 1) and (oscr[noscr-1] lt oscl[noscl-1])) then begin
      noscr++
      oscr = [oscr,n-2]
  endif
  oscl += 2     ;these adjust position, since we looked at the derivatives/changes
  oscr += 1
  ;
  abflag = bytarr(n)
  for j=0,n_elements(oscl)-1 do begin
      print,"$('Oscillatory segment ',i2,' bins ',i5,'-',i5)",j,oscl[j],oscr[j]
      t0 = time_iss[oscl[j]-2]
      p0 = ypr[oscl[j]-2].pitch
      t1 = time_iss[oscl[j]-1]
      p1 = ypr[oscl[j]-1].pitch
      for i=oscl[j],oscr[j] do begin
          if (abs(dpdt[i-1]) gt 0.007) then $
              abflag[i] = 1-abflag[i-1] $
          else $
              abflag[i] = abflag[i-1]
       ; double-check based on the last couple of normal points as to
       ; whether this is really abnormal
          if (abflag[i] eq 1) then begin
              pest = p1 + (p1-p0)/(t1-t0)*(time_iss[i]-t1)
              if (abs(pest-ypr[i].pitch) lt 0.007*(time_iss[i]-t1)) then $
                  abflag[i] = 0
          endif
          if (abflag[i] eq 0) then begin
              t0 = t1
              p0 = p1
              t1 = time_iss[i]
              p1 = ypr[i].pitch
          endif
      endfor
  endfor
  wpnormal = where(abflag eq 0)
endif else begin
    abflag = intarr(nt_iss)
    wpnormal = lindgen(nt_iss)
endelse
;
; quality is an integer, the lowest 12 bits are
; bits 11-8 yaw quality, bits 7-4 pitch, bits 3-0 roll
qual += round(intrpl(time_iss-time_iss[0],float(abflag),time_hrt-time_iss[0]))*16    ;mark values that look oscillatory
PRINT, 'Finished calculating pitch ...'


;------------------------------------------------------------------
; The look direction vector depends upon yaw and pitch.  Roll affects the slit orientation vector.
;
; The pitch should not change too much, at most a few 10ths of a degree up or down,
; unless there is a maneuver
; 1. Make an interpolation from where the repaired USGNC pitch data (p2) exists
; 3. Finally, perform a 30-sec median to get rid of spikes
;
pitch_all = INTRPL(time_iss[wpnormal]-time_iss[0],ypr[wpnormal].pitch,time_hrt-time_iss[0])
pitch_all = MEDIAN(pitch_all,31)
          
          ;!p.multi=[0,1,2]
          ;plot,time_iss,ypr.pitch,ps=3,xs=1,ys=3,charsize=1.5,tit='Pitch',xr=[min(time_hrt),max(time_hrt)]
          ;oplot,time_hrt,pitch_all,col=60 
          ;p2 = intrpl(time_hrt-time_iss[0],pitch_all,time_iss-time_iss[0])
          ;plot,time_iss,ypr.pitch-p2,xs=1,ys=3,charsize=1.5,tit='Pitch Diff',xr=[min(time_hrt),max(time_hrt)]


;------------------------------------------------------------------
; Now that we have yaw, pitch, and roll we can generate LVLH quaternions for the
; full data set at specified times
;
;______________________________________________________________________________
; Recalculate quaternions based on corrected ROLL/PITCH/YAW
eulerang = TRANSPOSE([[roll_all],[pitch_all],[yaw_all]])
quat_lvlh_all = EULER2QUAT(eulerang*!dpi/180d0)

;______________________________________________________________________________
; Since the rotations of the LVLH quaternions are about the instatntaneous 
; X, Y, and Z axes for LVLH, we need to compute the look direction in
; the LVLH frame.  This is pretty easy, since the body axes are nominally
; aligned with LVLH
; 
; Thus far, we have an angle which is positive toward zenith, and negative away from zenith
; The look direction vector is a combination of the za body (roughly negative radial) and 
; the xa body (roughly velocity) vectors, according to the mirror scan angle
; 
; VVIPRE is angle-wedged down from -xa axis of ISS, 48 deg
; ECLIPSE is mounted flush to the STP-H9 body, aligned with -xa (LIMB) or +za (DISK)
;
;         ISS Body Frame
;
;   -Za (Zenith) ^
;                |  / +Xa (RAM)
;                | /
;                |/
;    -Ya <-------|-------> +Ya
;               /|
;              / |
;             /  |
;       - Xa /   v +Za (Nadir)
;
; ECLIPSE LIMB points WAKE, in the +Za / -Xa quadrant
; ECLIPSE DISK points NADIR, in the +Za / +/- Ya Quadrants
;
; Note that the LVLH Za vector is directed downward toward nadir
;
mounting_angle = 0.0
view_angle = mounting_angle + 45. - 2.0*angle_hrt +deltaza
;______________________________________________________________________________
; Calculate vectors for viewing geometry calculations
IF STRUPCASE(viewer) EQ 'LIMB' THEN BEGIN
  ldv_body_x = -COS(view_angle*!dpi/180)
  ldv_body_y = REPLICATE(0d,nt_hrt)
  ldv_body_z = SIN(view_angle*!dpi/180)
  ;____________________________________________________________________________
  ; Account for any y-axis adjustment -- likely 0
  cy = COS(deltaya*!dpi/180d0)
  sy = SIN(deltaya*!dpi/180d0)
  ldv_body_x2 = cy*ldv_body_x - sy*ldv_body_y
  ldv_body_y2 = sy*ldv_body_x + cy*ldv_body_y
  ldv_body = TRANSPOSE([[ldv_body_x2],[ldv_body_y2],[ldv_body_z]])
ENDIF
IF STRUPCASE(viewer) EQ 'DISK' THEN BEGIN
  ldv_body_x = REPLICATE(0d,nt_hrt)
  ldv_body_y = SIN(view_angle*!dpi/180)
  ldv_body_z = COS(view_angle*!dpi/180)
  ldv_body = TRANSPOSE([[ldv_body_x],[ldv_body_y],[ldv_body_z]])
  zmin=100.
ENDIF
;______________________________________________________________________________
; Rotate this look direction according to LVLH quaternion
ldv_body_lvlh = APPLY_QUATERNION(quat_lvlh_all,ldv_body)
ldv_body_lvlh = ldv_body_lvlh[0:2,*]
lvdbr = SQRT(TOTAL(ldv_body_lvlh*ldv_body_lvlh,1))
;______________________________________________________________________________
; Compute look direction and slit orientation based on GEI coordinates
ldv_gei = ([1d0,1d0,1d0]#(ldv_body_lvlh[0,*]/lvdbr))*axis_gei.x +$
          ([1d0,1d0,1d0]#(ldv_body_lvlh[1,*]/lvdbr))*axis_gei.y +$
          ([1d0,1d0,1d0]#(ldv_body_lvlh[2,*]/lvdbr))*axis_gei.z
;ldv_gei = ([1d0,1d0,1d0]#(ldv_body_lvlh[0,*]/lvdbr))*TRANSPOSE(axis_gei.x) +$
;          ([1d0,1d0,1d0]#(ldv_body_lvlh[1,*]/lvdbr))*TRANSPOSE(axis_gei.y) +$
;          ([1d0,1d0,1d0]#(ldv_body_lvlh[2,*]/lvdbr))*TRANSPOSE(axis_gei.z)
ldvgr = SQRT(TOTAL(ldv_gei*ldv_gei,1))
ldv_gei = ldv_gei/([1d0,1d0,1d0]#ldvgr)
;______________________________________________________________________________
; Now do the slit orientation vector--Not sure if this is +ya or -ya
sov_body2 = APPLY_QUATERNION(quat_lvlh_all,[0,1,0]#REPLICATE(1,nt_hrt))        
sov_gei = ([1d0,1d0,1d0]#sov_body2[0,*])*axis_gei.x +$
          ([1d0,1d0,1d0]#sov_body2[1,*])*axis_gei.y +$
          ([1d0,1d0,1d0]#sov_body2[2,*])*axis_gei.z
sovr = SQRT(TOTAL(sov_gei*sov_gei,1))
sov_gei = sov_gei/([1d0,1d0,1d0]#sovr)
;______________________________________________________________________________
; Now feed to the viewing_geometery_vector_engine
zmin = 0.
viewing_geometry_vector_engine, spv_gei,ldv_gei,sov_gei,utinfo, $ ; INPUTS
	look_azi,lza,look_ra,look_dec,sc_lat,sc_lon,sc_alt, $ ;OUTPUT VARIABLES
	sc_zen,sc_radial,tp,tp_lat,tp_lon,tp_alt, $ ;OUTPUT VARIABLES
    tp_zen,sc_sza,tp_sza,suninfo,slit_azimuth,slit_roll, $ ;OUTPUT VARIABLES
    zmin=zmin ; INPUT KEYWORD

out_str = create_struct(name=strname,'met',0L,'met_hp',0d,'gps_sec',0d,$
                    'year',0,'month',0,'day',0,$
                    'ut_sec',0.,'look_vec',DBLARR(3),$
                    'look_az_deg',0.,'look_za_deg',0.,$
                    'look_ra_deg',0.,'look_dec_deg',0.,'lookqual',0,$
                    'sc_xkm',0d,'sc_ykm',0d,'sc_zkm',0d,$
                    'sc_lat',0.,'sc_lon',0.,$
                    'sc_alt',0.,'sc_re',0.,$
                    'sc_loc_vert',DBLARR(3),$
                    'sc_loc_zen_dot_radial',0D,$
                    'tanpt_vec',DBLARR(3),$
                    'tanpt_lat',0.,'tanpt_lon',0.,$
                    'tanpt_alt',0.,'tanpt_re',0d,$
                    'tanpt_loc_vert',DBLARR(3),$
                    'sun_vec',dblarr(3),'gmst_hrs',0.,$
                    'sc_sza',0.,'tanpt_sza',0.,$
                    'ss_lat',0.,'ss_lon',0.,'orbit',0,$
                    'slit_az_deg',0.,'slit_roll_deg',0.,$
                    'sc_yaw',0.,'sc_pitch',0.,'sc_roll',0.,'view_angle',0.)

; transfer the results to instrument specific structure
vg_names = TAG_NAMES(out_str)
vg = REPLICATE(out_str, nt_hrt)
for i=0,N_ELEMENTS(vg_names)-1 do begin
	CASE vg_names[i] OF
		'': dummy=0 ; SKIP THIS ONE
		'MET': vg.(i)=utinfo.met
		'MET_HP': vg.(i)=timepos.met_stph9
    'GPS_SEC': vg.(i) = time_hrt
    'YEAR': vg.(i) = utinfo.year
    'MONTH': vg.(i) = utinfo.month
	  'DAY': vg.(i)=utinfo.day
    'UT_SEC': vg.(i) = LONG(utinfo.secondofday)
    'LOOK_VEC': vg.(i)=ldv_gei ; array[3]
    'LOOK_VEC_X': vg.(i)=reform(ldv_gei[0,*]) ; x,y,z
    'LOOK_VEC_Y': vg.(i)=reform(ldv_gei[1,*]) ; x,y,z
    'LOOK_VEC_Z': vg.(i)=reform(ldv_gei[2,*]) ; x,y,z
	  'LOOK_AZ_DEG': vg.(i)=look_azi
    'LOOK_ZA_DEG': vg.(i)=lza
    'LOOK_RA_DEG': vg.(i)=look_ra
    'LOOK_DEC_DEG': vg.(i)=look_dec
    'LOOKQUAL': vg.(i)=qual
    'SC_VEC': vg.(i)=spv_gei
    'SC_XKM': vg.(i)=REFORM(spv_gei[0,*])
    'SC_YKM': vg.(i)=REFORM(spv_gei[1,*])
    'SC_ZKM': vg.(i)=REFORM(spv_gei[2,*])
    'SC_LAT': vg.(i)=sc_lat
    'SC_LON': vg.(i)=sc_lon
    'SC_ALT': vg.(i)=sc_alt
    'SC_RE': vg.(i)=REFGEOID(GEODGEOC(sc_lat))
	  'SC_LOC_VERT': vg.(i)=TRANSPOSE(sc_zen) ; array[3]
	  'SC_LOC_VERT_X': vg.(i)=sc_zen[*,0] ; x,y,z
	  'SC_LOC_VERT_Y': vg.(i)=sc_zen[*,1] ; x,y,z
	  'SC_LOC_VERT_Z': vg.(i)=sc_zen[*,2] ; x,y,z
	  'SC_LOC_ZEN_DOT_RADIAL': vg.(i)=TOTAL(sc_radial*sc_zen,2)
	  'TANPT_VEC': vg.(i)=TRANSPOSE(tp) ; array[3]
	  'TANPT_VEC_X': vg.(i)=tp[*,0] ; x,y,z
	  'TANPT_VEC_Y': vg.(i)=tp[*,1] ; x,y,z
	  'TANPT_VEC_Z': vg.(i)=tp[*,2] ; x,y,z
	  'TANPT_LAT': vg.(i)=tp_lat
	  'TANPT_LON': vg.(i)=tp_lon
    'TANPT_ALT': vg.(i)=tp_alt
    'TANPT_RE': vg.(i)=REFGEOID(GEODGEOC(tp_lat))
    'TANPT_LOC_VERT': vg.(i)=TRANSPOSE(tp_zen) ; array[3]
	  'TANPT_LOC_VERT_X': vg.(i)=tp_zen[*,0] ; x,y,z
	  'TANPT_LOC_VERT_Y': vg.(i)=tp_zen[*,1] ; x,y,z
	  'TANPT_LOC_VERT_Z': vg.(i)=tp_zen[*,2] ; x,y,z
    'SUN_VEC': vg.(i)=TRANSPOSE([[suninfo.x_gei2000],[suninfo.y_gei2000],[suninfo.z_gei2000]])  ; array[3]
    'SUN_VEC_X': vg.(i)=suninfo.x_gei2000 ; x,y,z
    'SUN_VEC_Y': vg.(i)=suninfo.y_gei2000 ; x,y,z
    'SUN_VEC_Z': vg.(i)=suninfo.z_gei2000 ; x,y,z
    'GMST_HRS': vg.(i)=suninfo.gmst_hrs
    'SC_SZA': vg.(i)=sc_sza
    'TANPT_SZA': vg.(i)=tp_sza
    'SS_LAT': vg.(i)=suninfo.sslat_geo
    'SS_LON': vg.(i)=suninfo.sslon_geo
    'ORBIT': vg.(i)=orb
    'SLIT_AZ_DEG': vg.(i)=REFORM(slit_azimuth)
    'SLIT_ROLL_DEG': vg.(i)=REFORM(slit_roll)
    'SC_YAW':   vg.(i) = yaw_all
    'SC_PITCH': vg.(i) = pitch_all
    'SC_ROLL':  vg.(i) = roll_all
    'VIEW_ANGLE': vg.(i) = view_angle
    	
	ELSE : begin
		   print,'vg value tag name=',vg_names(i),' does not match options.'
		   STOP
		   return,-1
	       endelse
	ENDCASE ; vg_names(i)
ENDFOR ; i

IF KEYWORD_SET(stopit) then STOP

RETURN, vg
END

;;; Plot time example
;plot, time_iss, time_iss, $
;  PSYM=4, THICK=0, $
;  xrange=[min(time_iss), min(time_iss)+8], yrange=[min(time_iss), min(time_iss)+8]
;oplot, time_hrt, time_hrt, THICK=0, PSYM=1, COLOR=150
;oplot, time_hrt, time_hrt, THICK=0, PSYM=6, COLOR=150

;;; plotfl routine
;if keyword_set(plotfl) then begin
;    if (sgp4fl) then begin
;       message,"Differences for SGP4 model vs SGP4 model not shown.",/info
;    endif else begin
;        if (addmergfl) then $
;            tit = 'Additive merge vs SGP4 model' $
;        else $
;            tit = 'Scaled merge vs SGP4 model'
;        pmulti = !p.multi
;        !p.multi = [0,2,2]
;        plot,time_hrt,spv_gei[0,*]-timepos.x_km,yr=[-abs(plotfl),abs(plotfl)],$
;            ys=1,xs=1,ytit='X Position Difference',tit=tit
;        plot,time_hrt,spv_gei[1,*]-timepos.y_km,yr=[-abs(plotfl),abs(plotfl)],$
;            ys=1,xs=1,ytit='Y Position Difference',tit=tit
;        plot,time_hrt,spv_gei[2,*]-timepos.z_km,yr=[-abs(plotfl),abs(plotfl)],$
;            ys=1,xs=1,ytit='Z Position Difference',tit=tit
;        spv0_gei = transpose([[timepos.x_km],[timepos.y_km],[timepos.z_km]])
;        plot,time_hrt,sqrt(total(spv_gei*spv_gei,1))-sqrt(total(spv0_gei*spv0_gei,1)),$
;            yr=[-abs(plotfl),abs(plotfl)],ys=1,xs=1,ytit='Radial Difference',tit=tit
;        !p.multi = pmulti
;    endelse
;endif

