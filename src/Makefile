# Makefile for timbreIDLib

PDDIR = /Applications/Pd-0.48-2.app/Contents/Resources

ldlibs = -lfftw3f
  
lib.name = timbreIDLib

# add this if I decide to put tIDLib.h into a separate include directory
# cflags = -I../include

lib.setup.sources = timbreIDLib.c

class.sources = attackTime.c attackTime~.c bark.c bark2freq.c barkSpec.c barkSpecBrightness.c barkSpecBrightness~.c barkSpecCentroid.c barkSpecCentroid~.c barkSpecFlatness.c barkSpecFlatness~.c barkSpecFlux.c barkSpecFlux~.c barkSpecIrregularity.c barkSpecIrregularity~.c barkSpecKurtosis.c barkSpecKurtosis~.c barkSpecRolloff.c barkSpecRolloff~.c barkSpecSkewness.c barkSpecSkewness~.c barkSpecSlope.c barkSpecSlope~.c barkSpecSpread.c barkSpecSpread~.c barkSpec~.c bark~.c bfcc.c bfcc~.c bin2freq.c binWrangler.c cepstrum.c cepstrum~.c cepstrumPitch.c cepstrumPitch~.c dct.c dct~.c featureAccum.c freq2bark.c freq2bin.c freq2mel.c magSpec.c magSpec~.c maxSample.c maxSample~.c maxSampleDelta.c maxSampleDelta~.c mel2freq.c melSpec.c melSpec~.c mfcc.c mfcc~.c minSample.c minSample~.c minSampleDelta.c minSampleDelta~.c nearestPoint.c peakSample.c peakSample~.c phaseSpec.c phaseSpec~.c sampleBuffer~.c specBrightness.c specBrightness~.c specCentroid.c specCentroid~.c specFlatness.c specFlatness~.c specFlux.c specFlux~.c specHarmonicity.c specHarmonicity~.c specIrregularity.c specIrregularity~.c specKurtosis.c specKurtosis~.c specRolloff.c specRolloff~.c specSkewness.c specSkewness~.c specSlope.c specSlope~.c specSpread.c specSpread~.c tempo~.c tID_fft.c tID_fft~.c tID_mean.c tID_std.c tIDLib.c timbreID.c tabletool.c waveSlope.c waveSlope~.c waveDirChange.c waveDirChange~.c zeroCrossing.c zeroCrossing~.c

datafiles = attackTime-help.pd attackTime~-help.pd bark-help.pd barkSpec-help.pd barkSpecBrightness-help.pd barkSpecBrightness~-help.pd barkSpecCentroid-help.pd barkSpecCentroid~-help.pd barkSpecFlatness-help.pd barkSpecFlatness~-help.pd barkSpecFlux-help.pd barkSpecFlux~-help.pd barkSpecIrregularity-help.pd barkSpecIrregularity~-help.pd barkSpecKurtosis-help.pd barkSpecKurtosis~-help.pd barkSpecRolloff-help.pd barkSpecRolloff~-help.pd barkSpecSkewness-help.pd barkSpecSkewness~-help.pd barkSpecSlope-help.pd barkSpecSlope~-help.pd barkSpecSpread-help.pd barkSpecSpread~-help.pd barkSpec~-help.pd bark~-help.pd bfcc-help.pd bfcc~-help.pd binWrangler-help.pd cepstrum-help.pd cepstrum~-help.pd cepstrumPitch-help.pd cepstrumPitch~-help.pd dct-help.pd dct~-help.pd featureAccum-help.pd  magSpec-help.pd magSpec~-help.pd maxSample-help.pd maxSample~-help.pd maxSampleDelta-help.pd maxSampleDelta~-help.pd melSpec-help.pd melSpec~-help.pd mfcc-help.pd mfcc~-help.pd minSample-help.pd minSample~-help.pd minSampleDelta-help.pd minSampleDelta~-help.pd nearestPoint-help.pd peakSample-help.pd peakSample~-help.pd phaseSpec-help.pd phaseSpec~-help.pd sampleBuffer~-help.pd specBrightness-help.pd specBrightness~-help.pd specCentroid-help.pd specCentroid~-help.pd specFlatness-help.pd specFlatness~-help.pd specFlux-help.pd specFlux~-help.pd specHarmonicity-help.pd specHarmonicity~-help.pd specIrregularity-help.pd specIrregularity~-help.pd specKurtosis-help.pd specKurtosis~-help.pd specRolloff-help.pd specRolloff~-help.pd specSkewness-help.pd specSkewness~-help.pd specSlope-help.pd specSlope~-help.pd specSpread-help.pd specSpread~-help.pd tempo~-help.pd tID_fft-help.pd tID_fft~-help.pd tID_mean-help.pd tID_std-help.pd tID-conversion-help.pd timbreID-help.pd tabletool-help.pd waveSlope-help.pd waveSlope~-help.pd waveDirChange-help.pd waveDirChange~-help.pd zeroCrossing-help.pd zeroCrossing~-help.pd

# build a multi-object library
make-lib-executable=yes

# include Makefile.pdlibbuilder from submodule directory 'pd-lib-builder'
PDLIBBUILDER_DIR=/Users/williambrent/dev/pd-lib-builder/
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder