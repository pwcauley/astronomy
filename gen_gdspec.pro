;;GEN_GDSPEC
;;
;;Procedure to use a gravity darkening solution and previously generated
;;SME (Spectroscopy Made Easy) spectra to produce a spectrum for the
;;gravity darkened stellar surface at some inclination angle to the 
;;viewer.  
;;
;;INPUTS:
;;	gdfile - string specifying the location and filename of the
;;		 Theta-R-T grav dark solution.
;;	sme_files - string *array* with the filenames of the pre-computed
;;		    SME spectra, can contain full path or just the file
;;                  names which will be combined with SPECDIR keyword
;;      tvec - array containing the temperatures at which SME_FILES
;;	       have been calculated. NOTE: TVEC and SME_FILES must be
;;	       the same length!
;;	wvec - input wavelength vector in ANGSTROMS. 
;;	gres - grid resolution in stellar radii (0.01 is a good choice)
;;	phi - inclination angle of the stellar rotation axis to the *viewer*
;;	      in degrees
;;	prot - stellar rotation period in hours
;;	
;;OUTPUTS:
;;	spec_out - the summed spectrum over the gravity darkened stellar
;;		   surface, same number of elements as WVEC
;;		   NOTE: WVEC and SPEC_OUT are also saved as columns in
;;		   an ascii table 
;;	 	
;;OPTIONAL KEYWORDS:
;;	specdir - a string containing the file location of the SME spectra,
;;		  assumed to be in current directory if not supplied or
;;		  full path can be specified using the SME_FILES input.
;;	one_temp - a floating point or double precision value specifying
;;		   a single temperature for the stellar surface if that's
;;		   desired
;;	save_str - a string containing the desired file name and path for
;;		   the output spectrum ascii file, filename is automatically
;;		   generated and saved in the current directory if not
;;		   specified
;;	rpow - resolving power. Note that the spectrum is NOT resampled
;;	       so significantly lowering the resolution from the original SME
;;	       spectrum will leave the spectrum highly oversampled

pro gen_gdspec,gdfile,sme_files,tvec,wvec,gres,phi,prot,spec_out,$
	save_str=save_str,specdir=specdir,one_temp=one_temp,$
	rpow=rpow

	;;Check to make sure TVEC and sme_files are same length
	if n_elements(tvec) ne n_elements(sme_files) then begin

		print,'TVEC and SME_FILES do not have same length! Check '
		print,'the number of elements in each vector and try again...'
		retall

	endif

        ;;Constant (cgs units except clight)
        msun=1.988d33
        rsun=6.96d10
        g=6.67d-8
        clight=299792.458d0

	;;If no SPEC_DIR defined, assume SME spectra are in the
	;;current directory.
	if not keyword_set(specdir) then specdir=''
	sme_files = specdir+sme_files

	;;Derive some important variables from the Theta-R-T file.
	;;Note that logg isn't actually used in any calculations
	;;since there is a one-to-one relationship between T and logg
	;;and the SME spectra at each T,logg pair are generated before
	;;running GEN_GDSPEC. logg is used in some file name stuff
	;;so we're keeping it in here.
	readcol,gdfile,colat,rad,temp,logg,/silent
	tpole=max(temp)
	logg_pole = max(logg) & logg_eq = min(logg)
        rpole=min(rad) & req=max(rad)
        f=(req-rpole)/req ;;oblateness
        phi0=(!Pi*(90.-phi))/180.
        ieq=where(rad eq req)
        colat=colat[0:ieq] & tlat=temp[0:ieq]
	logg=logg[0:ieq]

        ;;Derive vsini from phi and Prot
        vsini=(cos(phi0)*2.*!Pi*req*rsun)/(prot*3600.)/1d5

        ;;Build stellar grid
        gres_str=req*gres
        xstr=findgen(2.*req/gres_str+1)*gres_str-req
        ystr=xstr
        nx=n_elements(xstr) & ny=n_elements(ystr)
        xarr=xstr#(dblarr(ny)+1.)
        yarr=(dblarr(nx)+1.)#ystr
        vstarr=vsini*(xarr/req)
        dstarr=sqrt(xarr^2.+yarr^2.)	

        ;;Now calculate z for every x-y position on disk
        ;;FOUND: Typo in Eq. 14 from Barnes 2009! Third line, should be
        ;;(1-f)^2 not (1-f^2)
        d=4.*yarr^2.*(1.-(1.-f^2.))^2.*sin(phi0)^2.*cos(phi0)^2.$
                -4.*((cos(phi0)^2.*(1.-f)^2.+sin(phi0)^2.)*$
                ((yarr^2.*sin(phi0)^2.-req^2.+xarr^2.)*(1.-f)^2.+yarr^2.*cos(phi0)^2.))
        z=(-2.*yarr*(1.-(1.-f)^2.)*sin(phi0)*cos(phi0)+sqrt(d))/$
                (2.*((1.-f)^2.*cos(phi0)^2.+sin(phi0)^2.))

        ;;Rotated coordinates
        x0=xarr*rsun
        y0=(yarr*cos(phi0)+z*sin(phi0))*rsun
        z0=(-yarr*sin(phi0)+z*cos(phi0))*rsun
        isurf=where(finite(z0) eq 1,nsurf,complement=inosurf)

        ;;Co-latitude from poles (this is the variable we get from John)
        theta=!Pi/2.-abs(asin(y0/sqrt(x0^2.+y0^2.+z0^2.)))

        ;;Mu-grid using R(theta)
        mugrid=sqrt(1.-(dstarr/sqrt((x0/rsun)^2.+(y0/rsun)^2.+(z0/rsun)^2.))^2.)
        mugrid[inosurf]=0.

        ;;Now loop over surface and interpolate onto Ahler's vector
	;;for both TGRID. Note that T and logg have a one-to-one
	;;relationship so we use T as a proxy for logg and interpolate
	;;between pre-synthesized spectra
        tgrid=fltarr(nx,ny)
        for i=0,nsurf-1 do tgrid[isurf[i]]=interpol(tlat,colat,theta[isurf[i]])

	xstr=xstr/req & ystr=ystr/req

	if keyword_set(one_temp) then begin
	
		tgrid[where(tgrid ne 0.)]=one_temp
		tpole=one_temp

	endif

        if tpole le min(tvec) or tpole gt max(tvec) then begin

                print,'Max temperature is outside range of temperature grid!'
                print,'Check TVEC or make sure GDFILE is correct'
                retall

        endif

        ;;Define some other variables
        gstr=strmid(strtrim(gres,2),0,6)
        vstr=strmid(strtrim(vsini,2),0,3)
	nx=n_elements(xstr) & ny=nx
	w0=mean(wvec)
	vvec=(wvec-w0)*clight/w0
	nlam=n_elements(vvec)
	ares=gres^2.  ;;area element in stellar radii

	;;Begin loop over stellar surface and add to SOUT every time
	;;a new pixel is finished.
	igd=isurf & ngd=nsurf
	sout=0.
	for i=0,ngd-1 do begin

		jj=array_indices(mugrid,igd[i])
		temp=tgrid[igd[i]]
		tdiff=tvec-temp
		tpos=where(tdiff gt 0.)
		tneg=reverse(where(tdiff lt 0.))
		muval=mugrid[igd[i]]

		;;If temperature matches one of the TGRID values just use that spectrum
		if where(tdiff eq 0.) ne [-1] then begin

			;;Find the index where it matches 
			itemp=where(tdiff eq 0.)
			specfl = sme_files[itemp] 

			restore,specfl
			vsme=(sme.wint-w0)*clight/w0
			musme=sme.mu
			mudiff=musme-muval

			;;If equal to a mu value in MUSME just use that spectrum
			if where(mudiff eq 0.) ne [-1] then begin

				spec_sme=sme.sint[*,where(mudiff eq 0.)]
				speci=interpol(spec_sme,vsme+vstarr[igd[i]],vvec)
				sout=sout+speci	

			endif else begin

				;;If not, interpolate between spectra
                        	imulow=min(where(mudiff lt 0.))
                        	imuhigh=max(where(mudiff gt 0.))
                        	if muval le min(musme) then begin

                        	        mu1=0.
                        	        spec1=fltarr(n_elements(vsme))

                        	endif else begin

                                	mu1=musme[imulow]
                                	spec1=sme.sint[*,imulow]

                        	endelse
                        	mu2=musme[imuhigh]
                        	spec2=sme.sint[*,imuhigh]
                        	spec1=interpol(spec1,vsme+vstarr[igd[i]],vvec)
                        	spec2=interpol(spec2,vsme+vstarr[igd[i]],vvec)
                        	wt1=1.-(muval-mu1)/(mu2-mu1) & wt2=1.-(mu2-muval)/(mu2-mu1)
			
				speci=spec1*wt1+spec2*wt2
				sout=sout+speci

			endelse

		endif else begin 

			;;If no match, interpolate between TGRID spectra
			;;REad in SME spectra and interpolate onto correct Mu-value
			t1=tvec[tneg[0]] & t2=tvec[tpos[0]]
			twt1=1.-abs(temp-t1)/abs(t2-t1) & twt2=1.-abs(t2-temp)/abs(t2-t1)
			specfl1 = sme_files[tneg[0]] & specfl2 = sme_files[tpos[0]]

			restore,specfl1
			sme1=sme & musme1=sme1.mu
			vsme1=(sme1.wint-w0)*clight/w0
			mudiff=musme1-muval

			if where(mudiff eq 0.) ne [-1] then begin

				izz=where(mudiff eq 0.)
				spec_one=sme1.sint[*,izz]
				spec_one=interpol(spec_one,vsme1+vstarr[igd[i]],vvec)

			endif else begin

				imulow=min(where(mudiff lt 0.))
				imuhigh=max(where(mudiff gt 0.))
				if muval le min(musme1) then begin

					mu1=0.
					spec1=fltarr(n_elements(vsme1))

				endif else begin

					mu1=musme1[imulow]
					spec1=sme1.sint[*,imulow]

				endelse
				mu2=musme1[imuhigh]
				spec2=sme1.sint[*,imuhigh]
				spec1=interpol(spec1,vsme1+vstarr[igd[i]],vvec)
				spec2=interpol(spec2,vsme1+vstarr[igd[i]],vvec)
				wt1=1.-(muval-mu1)/(mu2-mu1) & wt2=1.-(mu2-muval)/(mu2-mu1)
				spec_one=spec1*wt1+spec2*wt2

			endelse

			;;Do second temperature
                        restore,specfl2
                        sme2=sme & musme2=sme2.mu
                        vsme2=(sme2.wint-w0)*clight/w0		
                        mudiff=musme2-muval

                        if where(mudiff eq 0.) ne [-1] then begin

                                izz=where(mudiff eq 0.)
                                spec_two=sme2.sint[*,izz]
                                spec_two=interpol(spec_two,vsme2+vstarr[igd[i]],vvec)

                        endif else begin

	                        imulow=min(where(mudiff lt 0.))
	                        imuhigh=max(where(mudiff gt 0.))
	                        if muval le min(musme1) then begin

	                                mu1=0.
	                                spec1=fltarr(n_elements(vsme2))

        	                endif else begin

        	                        mu1=musme2[imulow]
        	                        spec1=sme2.sint[*,imulow]

                	        endelse
                	        mu2=musme2[imuhigh]
                        	spec2=sme2.sint[*,imuhigh]
                        	spec1=interpol(spec1,vsme2+vstarr[igd[i]],vvec)
                        	spec2=interpol(spec2,vsme2+vstarr[igd[i]],vvec)
                        	wt1=1.-(muval-mu1)/(mu2-mu1) & wt2=1.-(mu2-muval)/(mu2-mu1)
                        	spec_two=spec1*wt1+spec2*wt2

			endelse

			;;Create final weighted spectrum for SCUBE
			speci=spec_one*twt1+spec_two*twt2
			sout=sout+speci

		endelse

	endfor

	;;Finally scale by grid pixel area
	sout=sout*ares

	;;Convolve to lower resolution if desired
	if keyword_set(rpow) then begin

                dv=vvec[1]-vvec[0]
                res0=dv   ;;resolution of synthetic spectrum in km/s
                res_new=clight/rpow       ;;desired resolution
                kwid=sqrt(abs(res_new^2.-res0^2.))  ;;width of Gaussian kernal in km/s
                nwid=(kwid/dv)/(2.*sqrt(2.*alog(2.)))   ;;number of elements for kernel FWHM
                nkern=10.*nwid ;;number of elements for full kernel
                if nkern mod 2 ne 0 then nkern=nkern+1
                xkern=findgen(nkern)-float(nkern)/2.
                gkern=gaussian(xkern,[1.,0.,nwid])

	        sout=convol(sout,gkern,/edge_mirror,/normalize)

	endif

	;;Write the spectrum to a csv file 
	if not keyword_set(save_str) then begin

	        if tpole ge 10000. then nts=5 else nts=4
	        grav = strmid(strtrim(mean(logg),2),0,3)
	        teff_str=strmid(strtrim(tpole,2),0,nts)
	        if min(wvec) lt 10000. then nws1 = 4 else nws1 = 5
	        if max(wvec) lt 10000. then nws2 = 4 else nws2 = 5
	        specstr = strmid(strtrim(min(wvec),2),0,nws1)+'_'+strmid(strtrim(max(wvec),2),0,nws2)
		fnam='specgd_'+specstr+'_'+teff_str+'_'+grav+'_'+gstr+'.dat'

	endif else fnam=save_str
	write_csv,fnam,wvec,sout
	spec_out = sout


end
