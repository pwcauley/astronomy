pro sme_gdspec_examp

	gfl_in='~/Dropbox/atmot_work/data/gdark_files/KELT9_RadiusTempProfile.txt'
	wtest=findgen(10000)*.1+2000.  ;;0.1 angstrom sampling from 2000 - 3000 angstroms
	temps=[9900.,10000.,10100.,10200.,10300.,10400.,10500.]
	sme_fls = ['s9900_4.1_0.0_2000_3000.out','s10000_4.1_0.0_2000_3000.out',$
		's10100_4.1_0.0_2000_3000.out','s10200_4.1_0.0_2000_3000.out',$
		's10300_4.1_0.0_2000_3000.out','s10400_4.1_0.0_2000_3000.out',$
		's10500_4.1_0.0_2000_3000.out']
	sdir='~/Dropbox/sme_tool/'
	save_str0='~/Dropbox/atmot_work/data/gdark_files/kelt9_tests.dat'
	res0=500

	;;Now run gen_gdspec with all of the keywords included, should generate SME
	;;spectra for temperatures in TVEC and the wavelength interval in SSTR
	gen_gdspec,gfl_in,sme_fls,temps,wtest,0.025,30.,16.,spec_test,$
        	save_str=save_str0,specdir=sdir,rpow=res0

end
