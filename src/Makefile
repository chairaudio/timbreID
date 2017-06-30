SOURCES = attackTime.c attackTime~.c bark.c bark2freq.c barkSpec.c barkSpecBrightness.c barkSpecBrightness~.c barkSpecCentroid.c barkSpecCentroid~.c barkSpecFlatness.c barkSpecFlatness~.c barkSpecFlux.c barkSpecFlux~.c barkSpecIrregularity.c barkSpecIrregularity~.c barkSpecKurtosis.c barkSpecKurtosis~.c barkSpecRolloff.c barkSpecRolloff~.c barkSpecSkewness.c barkSpecSkewness~.c barkSpecSlope.c barkSpecSlope~.c barkSpecSpread.c barkSpecSpread~.c barkSpec~.c bark~.c bfcc.c bfcc~.c bin2freq.c binWrangler.c cepstrum.c cepstrum~.c cepstrumPitch.c cepstrumPitch~.c dct.c dct~.c featureAccum.c freq2bark.c freq2bin.c freq2mel.c magSpec.c magSpec~.c maxSample.c maxSample~.c maxSampleDelta.c maxSampleDelta~.c mel2freq.c melSpec.c melSpec~.c mfcc.c mfcc~.c minSample.c minSample~.c minSampleDelta.c minSampleDelta~.c nearestPoint.c peakSample.c peakSample~.c sampleBuffer~.c specBrightness.c specBrightness~.c specCentroid.c specCentroid~.c specFlatness.c specFlatness~.c specFlux.c specFlux~.c specHarmonicity.c specHarmonicity~.c specIrregularity.c specIrregularity~.c specKurtosis.c specKurtosis~.c specRolloff.c specRolloff~.c specSkewness.c specSkewness~.c specSlope.c specSlope~.c specSpread.c specSpread~.c tID_mean.c tID_std.c tabletool.c timbreID.c waveSlope.c waveSlope~.c zeroCrossing.c zeroCrossing~.c


# ****SUPPLY THE LOCATION OF PD SOURCE****
pd_src = /Applications/Pd-0.47-1.app/Contents/Resources
# pd_src = /home/yourname/somedirectory/pd-0.47-1


CFLAGS = -DPD -I$(pd_src)/src -I/usr/local/include -Wall -W -g
LDFLAGS =
LIBS = -lfftw3f
LIBDIR = /usr/local/lib

UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
  SOURCES +=
  EXTENSION = pd_darwin
  OS = macosx
  OPT_CFLAGS = -O3 -ftree-vectorize
  FAT_FLAGS = -arch i386 -arch x86_64
  CFLAGS += -fPIC $(FAT_FLAGS)
  LDFLAGS += -bundle -undefined dynamic_lookup $(FAT_FLAGS)
  LIBS += -lc -L$(LIBDIR)
  STRIP = strip -x
 endif
ifeq ($(UNAME),Linux)
  EXTENSION = pd_linux
  OS = linux
  OPT_CFLAGS = -O6 -funroll-loops -fomit-frame-pointer
  CFLAGS += -fPIC
  LDFLAGS += -Wl,--export-dynamic -shared -fPIC
  LIBS += -lc
  STRIP = strip --strip-unneeded -R .note -R .comment
endif
ifeq (MINGW,$(findstring MINGW,$(UNAME)))
  SOURCES +=
  EXTENSION = dll
  OS = windows
  OPT_CFLAGS = -O3 -funroll-loops -fomit-frame-pointer -march=i686 -mtune=pentium4
  WINDOWS_HACKS = -D'O_NONBLOCK=1'
  CFLAGS += -mms-bitfields $(WINDOWS_HACKS)
  LDFLAGS += -s -shared -Wl,--enable-auto-import
  LIBS += -L$(pd_src)/bin -L$(pd_src)/obj -lpd -lwsock32 -lkernel32 -luser32 -lgdi32
  STRIP = strip --strip-unneeded -R .note -R .comment
endif

CFLAGS += $(OPT_CFLAGS)



all: $(SOURCES:.c=.$(EXTENSION))

%.o: %.c tIDLib.c tIDLib.h
	$(CC) $(CFLAGS) -o "$*.o" -c "$*.c"

%.$(EXTENSION): %.o
	$(CC) $(LDFLAGS) -o "$*.$(EXTENSION)" "$*.o" tIDLib.o $(LIBS)
	chmod a-x "$*.$(EXTENSION)"
	$(STRIP) $*.$(EXTENSION)
	rm -f -- $*.o


tIDLib.o: tIDLib.c tIDLib.h
	$(CC) $(CFLAGS) -o tIDLib.o -c tIDLib.c


.PHONY: clean

clean:
	-rm -f -- tIDLib.o
	-rm -f -- $(SOURCES:.c=.o)
	-rm -f -- $(SOURCES:.c=.$(EXTENSION))