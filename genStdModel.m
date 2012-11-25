function [src_map,bg_model]=genStdModel(num_sigma,exptime)
%GENSTDMODEL Generate standard model for calibrating the resolution of DDM.
% num_sigma: significance of point source in the generated model, in number
% of sigma.
%
% stand source: crab (1000 mcrab)
% effective area of the detector: 300 cm^2 per module, 5 narrow FOV modules
% for each position angle.

global config
NUMPIX=config.SAMPLING_FREQ;
clmt_cm2=300*5;
photon_keV=20;
mCrab2keV=1.43*6.2415e-3;

src_cts=(1000*mCrab2keV*clmt_cm2/photon_keV)*exptime;
bg_std=src_cts/num_sigma;
bg_obs=bg_std^2;

src_map=zeros(NUMPIX);
src_map(round(NUMPIX/2),round(NUMPIX/2))=src_cts/exptime;

PSF=genPSFClmt(0);
bg_model=bg_obs/sum(PSF(:))/exptime;

src_map=src_map+bg_model;
return
