;;GEN_SMESPECS
;;
;;Wrapper procedure for running Spectroscopy Made Easy at the IDL command
;;line instead of through the GUI. Useful for running large batches of
;;spectra instead of doing them one by one through the GUI. 

;;All inputs are string variables!!!

;;NOTE: Now TEFFS and GRAVS are *paired* so one spectrum is generated for
;;each *pair* of values. If only a single gravity is specified, gravity
;;is assumed to be the same for all TEFFS.

pro gen_SMEspecs,teffs,gravs,wints,zvals_in=zvals_in,sres=sres,lfls_in=lfls_in,$
	sme_dir=sme_dir,line_dir=line_dir

	if not keyword_set(lfls_in) then begin

		print,'No line file(s) specified! Returning...'
		retall

	endif

	;;Setup some directories. If nothing specified then assume working directory
	if not keyword_set(sme_dir) then sme_dir=''
	if not keyword_set(line_dir) then line_dir=''

	;;Define wints
	nints=n_elements(wints[0,*])

	;;Define temperatures to synthesize spectra for
	temps=teffs
	ntemps=n_elements(temps)

        ;;Define gravities to synthesize spectra for
        gvals=gravs
        ngvals=n_elements(gvals)
	if ngvals ne ntemps then begin

		print,'Assuming same gravity value for all Teff...'
		gvals=strarr(ntemps)+gvals[0]

	endif

	;;Define line files
	if n_elements(lfls_in) ne ntemps then begin

		print,'Assuming same VALD line file for all spectra...'
		line_files = line_dir + strarr(ntemps) + lfls_in[0]

	endif else line_files = line_dir + lfls_in

	;;In rare case of metallicity input...
	if keyword_set(zvals_in) then zvals=zvals_in else zvals=strarr(ntemps)+'0.0'
        if n_elements(zvals) ne ntemps then begin

                print,'Assuming same Z value for all Teff...'
                zvals=strarr(ntemps)+zvals[0]

        endif

	;;Make new mu-values
	muvals=[1.,.9,.8,.7,.6,.5,.4,.3,.2,.1,.05,.03,.01,.005]
	nmu=n_elements(muvals)

	;;Restore template SME structure
	restore,sme_dir+'sme_template.sav'
	sme0=sme
	struct_replace_field,sme0,'NMU',nmu
	struct_replace_field,sme0,'MU',muvals
	if keyword_set(sres) then struct_replace_field,sme0,'IPRES',long(sres)	
	sme_orig=sme0

	;;Loop through gravities and temperatures doing a spectrum for each PAIR
	for k=0,ntemps-1 do begin

		gval=gvals[k]
		teff=temps[k]
		zval = zvals[0]
		line_file = line_files[k]

		;;Loop through wavelength intervals and run SME_MAIN	
		for j=0,nints-1 do begin

			sme0=sme_orig

	                struct_replace_field,sme0,'GRAV',double(gval)
                        struct_replace_field,sme0,'TEFF',float(teff)

			;;Define line file
			wint=wints[*,j]
			lfl = line_file

			;;Restore line file data
			readcol,lfl,spec,ion,wl,excit,vmic,loggf,rad,stark,waals,fac,depth,format='A,D,D,D,D,D,D,D,D,D',/silent

			;;Have to sort b/c some lists from VALD are not sorted...
			ww=sort(wl)
			spec=spec[ww] & ion=ion[ww] & wl=wl[ww] & excit=excit[ww] & vmic=vmic[ww] & loggf=loggf[ww]
			rad=rad[ww] & stark=stark[ww] & waals=waals[ww] & fac=fac[ww] & depth=depth[ww]
			ion=strmid(strtrim(ion,2),0,1)
			spec=strmid(spec,1)
			elemparse,spec+' '+ion, anum, ionstage
			atom0=[transpose(anum),transpose(dblarr(n_elements(spec))+1.),transpose(wl),transpose(excit),$
				transpose(loggf),transpose(rad),transpose(stark),transpose(waals)]
			struct_replace_field,sme0,'SPECIES',spec+' '+ion
			struct_replace_field,sme0,'ATOMIC',atom0
			struct_replace_field,sme0,'LANDE',fac
			struct_replace_field,sme0,'DEPTH',depth
			struct_replace_field,sme0,'WRAN',[float(wint[0]),float(wint[1])]

			;;Save file
			if teff ge 1d4 then nts=5 else nts=4
			sfl=sme_dir+'s'+strmid(strtrim(teff,2),0,nts)+'_'+gval+'_'+zval+'_'+wint[0]+'_'+wint[1]+'.out'
		
			;;Run SME_MAIN
			sme_main,sme0
			sme=sme0
			save,sme,filename=sfl

		endfor

	endfor

end
