# gravdark_specs
IDL procedures for generating gravity-darkened spectral energy distributions. See
individual procedures for details related to running the programs. 
Note that the entire Spectrsocopy Made Easy library for IDL is required to
run gen_smespecs.pro. You can get it here: https://www.stsci.edu/~valenti/sme.html.

There are two additional files that are not IDL procedures: sme_template.sav
and lines_10000_4.1_0.0_1000_20000.lin. The sme_template.sav is a template
SME structure that gets read in by gen_smespecs.pro and modified to create
new spectra with the specified parameters. You don't have to use this one; any 
previously existing SME structure will do. The other file, lines*, is an atomic
line file from the VALD database. Check out the SME website above for information
on the format of these line files. They are necessary to generate spectra with
SME.

Also be aware that I did not put a lot of effort into making these procedures
easily usable by the public so they are not thoroughly documented in some cases.
If you are interested in generating a gravity-darkened SED with this code, feel
free to reach out to me at pwcauley at gmail.com if you can't get it to work.
