

some useful information:

LHAPDF location:
/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/lhapdf/6.1.6-giojec/share/LHAPDF

to use, compile

make

to run,

in bin directory.

for m50:

./local_run.sh z z_m50_nnlo_inclusive_NNPDF30_nlo_as_0118 input_z_m50.txt my_histograms.txt results_z_m50_nnlo_inclusive_NNPDF30_nlo_as_0118 .. 3

where 
"z" means run on Zll, 
z_m50_nnlo_inclusive_NNPDF30_nlo_as_0118: output directory
input_z_m50.txt : configuration
my_histograms.txt : hist format
results_z_m50_nnlo_inclusive_NNPDF30_nlo_as_0118 : the results tag
.. : the PDF directory parent directory
3 : 3 threads, can be any number

to submit to lsf on lxplus:
./sub_lsf.sh z z_m50_nnlo_inclusive_NNPDF30_nlo_as_0118 input_z_m50.txt my_histograms.txt results_z_m50_nnlo_inclusive_NNPDF30_nlo_as_0118 .


collect computing results:

./finish.sh z_m50_nnlo_inclusive_NNPDF30_nlo_as_0118 NNLO.results_z_m50_nnlo_inclusive_NNPDF30_nlo_as_0118

.
