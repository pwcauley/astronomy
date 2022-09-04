;;Example of how to use GEN_SMESPECS
;;The biggest changes are:
;;
;; -Same number of temperatures and gravities are now required and the
;;  spectra are computed in pairs, i.e., T1,g1 then T2,g2 then T3,g3 etc.
;;  The same is true for metallicity values but you probably won't use
;;  this option as much. If the number of gravities given to the program
;;  is 1, then it will assume the same gravity for all T
;; 
;; -You can input an array of VALD line file names to be used with each
;;  T,g pair. If only a single VALD line file is given then the program
;;  assumes the same line file for all T,g pairs. The line file is now
;;  a *required* input, although it is programmed as a keyword. If you
;;  don't specify one, the program will exit with a warning.
;;
;; -GEN_SMESPECS is no longer called in GEN_GDSPEC! All SME spectra
;;  must be generated before running GEN_GDSPEC and their file names
;;  are now used as an input to GEN_GDSPEC. This simplifies things
;;  and allows more customization regarding which SME spectra get
;;  used to create the gravity darkened SED.

pro gen_smespec_examp

	;;Specify the line file to be used. Here we use a single file for all
	;;T,g pairs but you can make this an array with a different line file
	;;for each pair
	lfls0='lines_10000_4.1_0.0_1000_20000.lin'

	;;Specify the directory where the SME spectra will get saved
	smedir0='~/Dropbox/sme_tool/'

	;;Specify the directory where the line files live
	ldir0='~/Dropbox/idl_files/sme_522/line_data/'

	;;Create array of temperature values for which to compute
	;;spectra. Can now be string or float, doesn't matter.
	temps_in = ['9900','10000','10100']

        ;;Create a *string* array of temperature values for which to compute
        ;;spectra. 
        gravs_in = ['4.1']

	;;Specify which spectral ranges across which to compute the spectra.
	;;Remember to keep these small, on the order of 2000 angstroms or less.
	;;SME does not like when they get large since the spectra are calculated
	;;at high resolution and the arrays can get unmanageable for huge
	;;wavelength regions.
	wints_in = [['5000','6000'],['8000','9000']]

	;;Run gen_smespecs to synthesize the spectra
	gen_smespecs,temps_in,gravs_in,wints_in,$
		lfls_in=lfls0,line_dir=ldir0,sme_dir=smedir0


end
