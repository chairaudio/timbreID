SOURCES = attackTime.c attackTime~.c bark.c bark2freq.c barkSpec.c barkSpecBrightness.c barkSpecBrightness~.c barkSpecCentroid.c barkSpecCentroid~.c barkSpecFlatness.c barkSpecFlatness~.c barkSpecFlux.c barkSpecFlux~.c barkSpecIrregularity.c barkSpecIrregularity~.c barkSpecKurtosis.c barkSpecKurtosis~.c barkSpecRolloff.c barkSpecRolloff~.c barkSpecSkewness.c barkSpecSkewness~.c barkSpecSlope.c barkSpecSlope~.c barkSpecSpread.c barkSpecSpread~.c barkSpec~.c bark~.c bfcc.c bfcc~.c bin2freq.c binWrangler.c cepstrum.c cepstrum~.c cepstrumPitch.c cepstrumPitch~.c dct.c dct~.c featureAccum.c freq2bark.c freq2bin.c freq2mel.c magSpec.c magSpec~.c maxSample.c maxSample~.c maxSampleDelta.c maxSampleDelta~.c mel2freq.c melSpec.c melSpec~.c mfcc.c mfcc~.c minSample.c minSample~.c minSampleDelta.c minSampleDelta~.c nearestPoint.c peakSample.c peakSample~.c phaseSpec.c phaseSpec~.c sampleBuffer~.c specBrightness.c specBrightness~.c specCentroid.c specCentroid~.c specFlatness.c specFlatness~.c specFlux.c specFlux~.c specHarmonicity.c specHarmonicity~.c specIrregularity.c specIrregularity~.c specKurtosis.c specKurtosis~.c specRolloff.c specRolloff~.c specSkewness.c specSkewness~.c specSlope.c specSlope~.c specSpread.c specSpread~.c tempo~.c tID_fft.c tID_fft~.c tID_mean.c tID_std.c tIDLib.c timbreID.c timbreIDLib.c tabletool.c waveSlope.c waveSlope~.c waveDirChange.c waveDirChange~.c zeroCrossing.c zeroCrossing~.c

# ****SUPPLY THE LOCATION OF PD SOURCE****
# pd_src = /home/yourname/somedirectory/pd-0.48-1
pd_src = /Users/williambrent/Dropbox/PUREDATA/pd-0.48-1
# pd_src = "C:\Program Files (x86)\Pd"

CFLAGS = -DPD -I$(pd_src)/src -I/usr/local/include -Wall -W -g


UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
  EXTENSION = pd_darwin
  OS = macosx
  LIBDIR = /usr/local/lib
  OPT_CFLAGS = -O3 -ftree-vectorize
  FAT_FLAGS = -arch i386 -arch x86_64
  CFLAGS += -fPIC $(FAT_FLAGS)
  LDFLAGS += -bundle -undefined dynamic_lookup $(FAT_FLAGS)
  LIBS += -lc -L$(LIBDIR) -lfftw3f
  STRIP = strip -x
 endif
ifeq ($(UNAME),Linux)
  EXTENSION = pd_linux
  OS = linux
  LIBDIR = /usr/local/lib
  OPT_CFLAGS = -O6 -funroll-loops -fomit-frame-pointer
  CFLAGS += -fPIC
  LDFLAGS += -Wl,--export-dynamic -shared -fPIC
  LIBS += -lc -lfftw3f
  STRIP = strip --strip-unneeded -R .note -R .comment
endif
ifeq (MINGW,$(findstring MINGW,$(UNAME)))
  CC = gcc
  EXTENSION = dll
  OS = windows
  LIBDIR = "C:\MinGW\lib\fftw-3.3.5-dll32"
  OPT_CFLAGS = -O3 -funroll-loops -fomit-frame-pointer \
    -march=pentium4 -mfpmath=sse -msse -msse2
  WINDOWS_HACKS = -D'O_NONBLOCK=1'
  CFLAGS += -mms-bitfields $(WINDOWS_HACKS) \
    -I"C:\MinGW\lib\fftw-3.3.5-dll32"
  LDFLAGS += -static-libgcc -s -shared \
    -Wl,--enable-auto-import $(pd_src)/bin/pd.dll
  LIBS += -L$(LIBDIR) -L$(pd_src)/bin -lpd -lfftw3f-3 \
    -lwsock32 -lkernel32 -luser32 -lgdi32
  STRIP = strip --strip-unneeded -R .note -R .comment
endif

CFLAGS += $(OPT_CFLAGS)


all: $(SOURCES:.c=.o)
	$(CC) $(LDFLAGS) -o timbreIDLib.$(EXTENSION) $(SOURCES:.c=.o) $(LIBS)
	chmod a-x timbreIDLib.$(EXTENSION)
	$(STRIP) timbreIDLib.$(EXTENSION)
	rm -f -- $(SOURCES:.c=.o)

%.o: %.c tIDLib.h
	$(CC) $(CFLAGS) -o "$*.o" -c "$*.c"


.PHONY: clean

clean:
	-rm -f -- $(SOURCES:.c=.o)
	-rm -f -- timbreIDLib.$(EXTENSION)
