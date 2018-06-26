/*

tempo~

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

static t_class *tempo_tilde_class;

typedef struct _tempo_tilde
{
    t_object x_obj;
    t_symbol *x_objSymbol;
    t_float x_sr;
    t_float x_n;
	t_uInt x_hop;
    t_sampIdx x_window;
    t_sampIdx x_windowHalf;
    t_windowFunction x_windowFunction;
    t_uShortInt x_overlap;
    t_bool x_normalize;
    t_bool x_powerSpectrum;
    t_bool x_alignFirstPeak;
    t_float x_growthThresh;
    t_float x_belowThreshDefault;
    t_binIdx x_loBin;
    t_binIdx x_hiBin;
    t_uInt x_onsetsBufSize;
    t_uInt x_numHarm;
    t_float x_loTempo;
    t_float x_hiTempo;
    t_uInt x_loTempoIdx;
    t_uInt x_hiTempoIdx;
    double x_lastDspTime;
    t_uInt x_dspTicks;
    t_sample *x_signalBuffer;
    t_float *x_fftwInForwardWindow;
    t_float *x_fftwInBackWindow;
    t_float *x_onsetsBuffer;
    fftwf_complex *x_fftwOutForwardWindow;
    fftwf_complex *x_fftwOutBackWindow;
	fftwf_plan x_fftwPlanForwardWindow;
	fftwf_plan x_fftwPlanBackWindow;
    t_float *x_blackman;
    t_float *x_cosine;
    t_float *x_hamming;
    t_float *x_hann;
    t_bool x_squaredDiff;
	t_atom *x_listOut;
    t_outlet *x_onsetList;
    t_outlet *x_tempo;
    t_float x_f;

} t_tempo_tilde;



// this needs to be included as part of tIDLib.c, and used in tabletool.c as well as here
static t_float tempo_tilde_hps(t_tempo_tilde *x, t_float loIdx, t_float hiIdx, t_float numHarm)
{
	t_sampIdx i, j, numIndices, maxIdx;
	t_float *yValues, maxVal;

	// if the highest harmonic of the hiIdx is beyond the end of the array data, we'll reduce the requested number of harmonics and post a warning
	while(hiIdx*numHarm > x->x_onsetsBufSize)
	{
		numHarm--;
		post("%s: WARNING: reducing numHarm to %f", x->x_objSymbol->s_name, numHarm);

		if(numHarm<=1)
		{
			pd_error(x, "%s: HPS function: second harmonic of hiIdx is out of table bounds. Aborting.", x->x_objSymbol->s_name);
			return(-1);
		}	
	}

	numIndices = hiIdx - loIdx + 1;
	// need a safety check that numIndices is at least 1
	numIndices = (numIndices<1)?1:numIndices;

	yValues = (t_float *)t_getbytes(numIndices*sizeof(t_float));

	// init yValues arrays to zero
	for(i=0; i<numIndices; i++)
		yValues[i] = 0.0;

	for(i=0; i<numIndices; i++)
	{
		t_float thisProduct;
		
		thisProduct = 1.0;
		
		for(j=0; j<numHarm; j++)
		{
			t_sampIdx thisIdx;
			
			thisIdx = (loIdx+i)*j;
			thisProduct = thisProduct * atom_getfloat(x->x_listOut + thisIdx);
			//post("i: %i, thisProduct: %f",  i, thisProduct);
		}
		
		yValues[i] = thisProduct;
	}

	maxVal = -1.0;
	maxIdx = UINT_MAX;
	
	for(i=0; i<numIndices; i++)
	{
		if(yValues[i]>maxVal)
		{
			maxVal = yValues[i];
			maxIdx = i;
		}
	}

	// free local memory
	t_freebytes(yValues, numIndices * sizeof(t_float));
	
	// if maxIdx is somehow not updated, output -1 to indicate failure.
	// otherwise, add the loIdx offset back in to get the array index of the HPS peak
	if(maxIdx==UINT_MAX)
		return(-1.0);
	else
		return(loIdx + maxIdx);
}
	

/* ------------------------ tempo~ -------------------------------- */

static void tempo_tilde_analyze(t_tempo_tilde *x)
{
    t_sampIdx i, j, startIdx, periodIdx, window, windowHalf;
    t_float growth, tempo, *windowFuncPtr;

    window = x->x_window;
    windowHalf = x->x_windowHalf;

	// construct forward analysis window
	for(i=0, j=0; i<window; i++, j++)
		x->x_fftwInForwardWindow[i] = x->x_signalBuffer[window + j];

	// construct back analysis window x->x_hop samples earlier
	for(i=0, j=0; i<window; i++, j++)
		x->x_fftwInBackWindow[i] = x->x_signalBuffer[window - x->x_hop + j];
	
	switch(x->x_windowFunction)
	{
		case rectangular:
			break;
		case blackman:
			windowFuncPtr = x->x_blackman;
			break;
		case cosine:
			windowFuncPtr = x->x_cosine;
			break;
		case hamming:
			windowFuncPtr = x->x_hamming;
			break;
		case hann:
			windowFuncPtr = x->x_hann;
			break;
		default:
			windowFuncPtr = x->x_blackman;
			break;
	};

	// if x_windowFunction == 0, skip the windowing (rectangular)
	if(x->x_windowFunction!=rectangular)
		for(i=0; i<window; i++, windowFuncPtr++)
			x->x_fftwInForwardWindow[i] *= *windowFuncPtr;

	fftwf_execute(x->x_fftwPlanForwardWindow);

	// put the result of power calc back in x_fftwIn
	tIDLib_power(windowHalf+1, x->x_fftwOutForwardWindow, x->x_fftwInForwardWindow);
	
	if(!x->x_powerSpectrum)
		tIDLib_mag(windowHalf+1, x->x_fftwInForwardWindow);

	if(x->x_normalize)
		tIDLib_normal(windowHalf+1, x->x_fftwInForwardWindow);
			
	switch(x->x_windowFunction)
	{
		case rectangular:
			break;
		case blackman:
			windowFuncPtr = x->x_blackman;
			break;
		case cosine:
			windowFuncPtr = x->x_cosine;
			break;
		case hamming:
			windowFuncPtr = x->x_hamming;
			break;
		case hann:
			windowFuncPtr = x->x_hann;
			break;
		default:
			windowFuncPtr = x->x_blackman;
			break;
	};

	// if x_windowFunction == 0, skip the windowing (rectangular)
	if(x->x_windowFunction!=rectangular)
		for(i=0; i<window; i++, windowFuncPtr++)
			x->x_fftwInBackWindow[i] *= *windowFuncPtr;

	fftwf_execute(x->x_fftwPlanBackWindow);

	// put the result of power calc back in x_fftwIn
	tIDLib_power(windowHalf+1, x->x_fftwOutBackWindow, x->x_fftwInBackWindow);
	
	if(!x->x_powerSpectrum)
		tIDLib_mag(windowHalf+1, x->x_fftwInBackWindow);

	growth=0;

	for(i=x->x_loBin; i<=x->x_hiBin; i++)
	{
		t_float diff;

		diff = x->x_fftwInForwardWindow[i] - x->x_fftwInBackWindow[i];
		
		// only look at growth (positive flux)
		if(diff>0)
		{
			if(x->x_squaredDiff)
				diff = diff*diff;
		}
		else
			diff = 0.0;
			
		growth += diff;
	}

	for(i=0; i<x->x_onsetsBufSize-1; i++)
		x->x_onsetsBuffer[i] = x->x_onsetsBuffer[i+1];
	
	if(growth>x->x_growthThresh)
		x->x_onsetsBuffer[x->x_onsetsBufSize-1] = growth;
	else
		x->x_onsetsBuffer[x->x_onsetsBufSize-1] = x->x_belowThreshDefault; //set below thresh growth to 1.0 and NOT ZERO, which would ruin the running product in the HPS algo


	startIdx=0;
			
	// align first peak to beginning?
	if(x->x_alignFirstPeak)
	{
		while(x->x_onsetsBuffer[startIdx]==x->x_belowThreshDefault)
		{
			startIdx++;

			if(startIdx>=x->x_onsetsBufSize-2)
				break;
		}
	}

	for(i=startIdx, j=0; i<x->x_onsetsBufSize; i++, j++)
		SETFLOAT(x->x_listOut+j, x->x_onsetsBuffer[i]);

	for(; j<x->x_onsetsBufSize; j++)
		SETFLOAT(x->x_listOut+j, x->x_belowThreshDefault);

	periodIdx = tempo_tilde_hps(x, x->x_hiTempoIdx, x->x_loTempoIdx, x->x_numHarm);

	tempo = (t_float)(periodIdx*x->x_hop); // frames to samples
	tempo /= x->x_sr; // samples to seconds
	tempo = 60.0f/tempo; // seconds to BPM

 	outlet_list(x->x_onsetList, 0, x->x_onsetsBufSize, x->x_listOut);
	outlet_float(x->x_tempo, tempo);
}


static void tempo_tilde_print(t_tempo_tilde *x)
{
	post("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr/x->x_overlap));
	post("%s block size: %i", x->x_objSymbol->s_name, (t_uShortInt)x->x_n);
	post("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
	post("%s window: %i", x->x_objSymbol->s_name, x->x_window);
	post("%s normalize: %i", x->x_objSymbol->s_name, x->x_normalize);
	post("%s power spectrum: %i", x->x_objSymbol->s_name, x->x_powerSpectrum);
	post("%s window function: %i", x->x_objSymbol->s_name, x->x_windowFunction);
	post("%s hop: %i", x->x_objSymbol->s_name, x->x_hop);
	post("%s squared difference: %i", x->x_objSymbol->s_name, x->x_squaredDiff);
	post("%s align first peak: %i", x->x_objSymbol->s_name, x->x_alignFirstPeak);
	post("%s onset buffer size: %i", x->x_objSymbol->s_name, x->x_onsetsBufSize);
	post("%s growth thresh: %f.", x->x_objSymbol->s_name, x->x_growthThresh);
	post("%s below threshold default value: %f.", x->x_objSymbol->s_name, x->x_belowThreshDefault);
	post("%s frequency range: bin %i through %i.", x->x_objSymbol->s_name, x->x_loBin, x->x_hiBin);
	post("%s tempo range: frame %i through %i.", x->x_objSymbol->s_name, x->x_loTempoIdx, x->x_hiTempoIdx);
	post("%s harmonics: %i.", x->x_objSymbol->s_name, x->x_numHarm);
}


static void tempo_tilde_window(t_tempo_tilde *x, t_floatarg w)
{
	t_sampIdx i, window, windowHalf;
	
	window = w;
	
	if(window<MINWINDOWSIZE)
	{
		window = WINDOWSIZEDEFAULT;
		post("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, MINWINDOWSIZE, WINDOWSIZEDEFAULT);
	}

	if(x->x_hop > window)
	{
		post("%s change in window size caused frame hop to be less than %i samples apart. Setting frame hop to half of current window size instead.", x->x_objSymbol->s_name, window);
		x->x_hop = window*0.5;
	}
	
	windowHalf = window*0.5;
	
	x->x_signalBuffer = (t_sample *)t_resizebytes(x->x_signalBuffer, (x->x_window*2) * sizeof(t_sample), (window*2) * sizeof(t_sample));
	x->x_fftwInForwardWindow = (t_float *)t_resizebytes(x->x_fftwInForwardWindow, x->x_window*sizeof(t_float), window*sizeof(t_float));
	x->x_fftwInBackWindow = (t_float *)t_resizebytes(x->x_fftwInBackWindow, x->x_window*sizeof(t_float), window*sizeof(t_float));
	
	x->x_blackman = (t_float *)t_resizebytes(x->x_blackman, x->x_window*sizeof(t_float), window*sizeof(t_float));
	x->x_cosine = (t_float *)t_resizebytes(x->x_cosine, x->x_window*sizeof(t_float), window*sizeof(t_float));
	x->x_hamming = (t_float *)t_resizebytes(x->x_hamming, x->x_window*sizeof(t_float), window*sizeof(t_float));
	x->x_hann = (t_float *)t_resizebytes(x->x_hann, x->x_window*sizeof(t_float), window*sizeof(t_float));

	x->x_window = window;
	x->x_windowHalf = windowHalf;

	// free the FFTW output buffer, and re-malloc according to new window
	fftwf_free(x->x_fftwOutForwardWindow);
	fftwf_free(x->x_fftwOutBackWindow);
	
	// destroy old plan, which depended on x->x_window
	fftwf_destroy_plan(x->x_fftwPlanForwardWindow); 
	fftwf_destroy_plan(x->x_fftwPlanBackWindow); 
	
	// allocate new FFTW output buffer memory
	x->x_fftwOutForwardWindow = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf+1);
	x->x_fftwOutBackWindow = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf+1);

	// create a new plan
	x->x_fftwPlanForwardWindow = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwInForwardWindow, x->x_fftwOutForwardWindow, FFTWPLANNERFLAG);
	x->x_fftwPlanBackWindow = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwInBackWindow, x->x_fftwOutBackWindow, FFTWPLANNERFLAG);

 	// we're supposed to initialize the input array after we create the plan
	for(i=0; i<x->x_window; i++)
 	{
		x->x_fftwInForwardWindow[i] = 0.0;
		x->x_fftwInBackWindow[i] = 0.0;
	}
	
	// initialize signal buffer
	for(i=0; i<(x->x_window*2); i++)
		x->x_signalBuffer[i] = 0.0;

	// re-init window functions
	tIDLib_blackmanWindow(x->x_blackman, x->x_window);
	tIDLib_cosineWindow(x->x_cosine, x->x_window);
	tIDLib_hammingWindow(x->x_hamming, x->x_window);
	tIDLib_hannWindow(x->x_hann, x->x_window);

	post("%s window size: %i", x->x_objSymbol->s_name, x->x_window);
}


static void tempo_tilde_overlap(t_tempo_tilde *x, t_floatarg o)
{
	// this change will be picked up the next time _dsp is called, where the samplerate will be updated to sp[0]->s_sr/x->x_overlap;
	x->x_overlap = (o<1)?1:o;

    post("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
}


static void tempo_tilde_x_windowFunction(t_tempo_tilde *x, t_floatarg f)
{
    f = (f<0)?0:f;
    f = (f>4)?4:f;
	x->x_windowFunction = f;

	switch(x->x_windowFunction)
	{
		case rectangular:
			post("%s window function: rectangular.", x->x_objSymbol->s_name);
			break;
		case blackman:
			post("%s window function: blackman.", x->x_objSymbol->s_name);
			break;
		case cosine:
			post("%s window function: cosine.", x->x_objSymbol->s_name);
			break;
		case hamming:
			post("%s window function: hamming.", x->x_objSymbol->s_name);
			break;
		case hann:
			post("%s window function: hann.", x->x_objSymbol->s_name);
			break;
		default:
			break;
	};
}


static void tempo_tilde_powerSpectrum(t_tempo_tilde *x, t_floatarg spec)
{
    spec = (spec<0)?0:spec;
    spec = (spec>1)?1:spec;
	x->x_powerSpectrum = spec;

	if(x->x_powerSpectrum)
		post("%s using power spectrum", x->x_objSymbol->s_name);
	else
		post("%s using magnitude spectrum", x->x_objSymbol->s_name);
}


static void tempo_tilde_squaredDiff(t_tempo_tilde *x, t_floatarg sd)
{
    sd = (sd<0)?0:sd;
    sd = (sd>1)?1:sd;
	x->x_squaredDiff = sd;
	
	if(x->x_squaredDiff)
		post("%s using squared difference", x->x_objSymbol->s_name);
	else
		post("%s using absolute value of difference", x->x_objSymbol->s_name);
}


static void tempo_tilde_normalize(t_tempo_tilde *x, t_floatarg norm)
{
    norm = (norm<0)?0:norm;
    norm = (norm>1)?1:norm;
	x->x_normalize = norm;

	if(x->x_normalize)
		post("%s spectrum normalization ON.", x->x_objSymbol->s_name);
	else
		post("%s spectrum normalization OFF.", x->x_objSymbol->s_name);
}


static void tempo_tilde_alignFirstPeak(t_tempo_tilde *x, t_floatarg align)
{
    align = (align<0)?0:align;
    align = (align)?1:align;
	x->x_alignFirstPeak = align;
	
	if(x->x_alignFirstPeak)
		post("%s align first peak: ON.", x->x_objSymbol->s_name);
	else
		post("%s align first peak: OFF.", x->x_objSymbol->s_name);
}


static void tempo_tilde_thresh(t_tempo_tilde *x, t_floatarg thresh)
{
    thresh = (thresh<0)?0:thresh;
	x->x_growthThresh = thresh;

	post("%s growth thresh: %f.", x->x_objSymbol->s_name, x->x_growthThresh);
}


static void tempo_tilde_freqRange(t_tempo_tilde *x, t_floatarg loFreq, t_floatarg hiFreq)
{
    loFreq = (loFreq<0)?0:loFreq;
    loFreq = (loFreq>(x->x_sr*0.5))?x->x_sr*0.5:loFreq;
    
	x->x_loBin = tIDLib_freq2bin(loFreq, x->x_window, x->x_sr);

    hiFreq = (hiFreq<0)?0:hiFreq;
    hiFreq = (hiFreq>(x->x_sr*0.5))?x->x_sr*0.5:hiFreq;
    
	x->x_hiBin = tIDLib_freq2bin(hiFreq, x->x_window, x->x_sr);
	
	post("%s frequency range: %fHz through %fHz (bin %i through %i).", x->x_objSymbol->s_name, loFreq, hiFreq, x->x_loBin, x->x_hiBin);
}


static void tempo_tilde_tempoRange(t_tempo_tilde *x, t_floatarg loTempo, t_floatarg hiTempo)
{	
    loTempo = (loTempo<10)?10:loTempo;
    loTempo = (loTempo>500)?500:loTempo;
	x->x_loTempo = loTempo;
	    
	loTempo = 60.0f/loTempo; // BPM to seconds
	loTempo *= x->x_sr; // seconds to samples
	loTempo /= x->x_hop;//samples to frames
	x->x_loTempoIdx = roundf(loTempo);

    hiTempo = (hiTempo<10)?10:hiTempo;
    hiTempo = (hiTempo>500)?500:hiTempo;
   	x->x_hiTempo = hiTempo;
 
	hiTempo = 60.0f/hiTempo; // BPM to seconds
	hiTempo *= x->x_sr; // seconds to samples
	hiTempo /= x->x_hop;//samples to frames
	x->x_hiTempoIdx = roundf(hiTempo);

// we could automatically resize the onsetsBuffer here based on the new loTempoIdx and current numHarm
//	x->x_onsetsBufSize = x->x_loTempoIdx*(x->x_numHarm+1);

	post("%s tempo range: %fBPM through %fBPM (frameIdx %i through %i).", x->x_objSymbol->s_name, x->x_loTempo, x->x_hiTempo, x->x_loTempoIdx, x->x_hiTempoIdx);
}


static void tempo_tilde_numHarm(t_tempo_tilde *x, t_floatarg h)
{
	x->x_numHarm = (h<2)?2:h;

// we could automatically resize the onsetsBuffer here based on the current loTempoIdx and new numHarm
//	x->x_onsetsBufSize = x->x_loTempoIdx*(x->x_numHarm+1);

    post("%s harmonics: %i", x->x_objSymbol->s_name, x->x_numHarm);
}


static void tempo_tilde_onsetsBufDur(t_tempo_tilde *x, t_floatarg n)
{
	t_float newSizeMs;
	t_uInt i, newSize, minSize;
	
    n = (n<100)?100:n;
	
	// ms to samples
	newSizeMs = (n*x->x_sr/1000.0f);
	
	// samples to frames
	newSize = newSizeMs/(t_float)x->x_hop;
	
	minSize = x->x_loTempoIdx*(x->x_numHarm+1);

	if(newSize<minSize)
		pd_error(x, "%s: requested size is too short based on tempo and harmonic setetings.", x->x_objSymbol->s_name);
	else
	{
		x->x_onsetsBuffer = (t_float *)t_resizebytes(x->x_onsetsBuffer, x->x_onsetsBufSize * sizeof(t_float), newSize * sizeof(t_float));
		x->x_listOut = (t_atom *)t_resizebytes(x->x_listOut, x->x_onsetsBufSize * sizeof(t_atom), newSize * sizeof(t_atom));
	
		x->x_onsetsBufSize = newSize;

		for(i=0; i<x->x_onsetsBufSize; i++)
			x->x_onsetsBuffer[i] = 0.0;

		post("%s onset buffer size: %f ms, %i frames.", x->x_objSymbol->s_name, n, x->x_onsetsBufSize);
	}
}


static void tempo_tilde_belowThreshDefault(t_tempo_tilde *x, t_floatarg val)
{
    x->x_belowThreshDefault = val;

	post("%s below threshold default value: %f.", x->x_objSymbol->s_name, x->x_belowThreshDefault);
}


static void tempo_tilde_hop(t_tempo_tilde *x, t_floatarg h)
{
	t_float loTempo, hiTempo;
	
	if(h > x->x_window)
	{
		post("%s frame hop cannot be more than current window size. Using half of current window size instead.", x->x_objSymbol->s_name);
		x->x_hop = x->x_windowHalf;
	}
	else if(h < 0)
	{
		post("%s frame hop must be > 0. Using half of current window size instead.", x->x_objSymbol->s_name);
		x->x_hop = x->x_windowHalf;
	}
	else
		x->x_hop = h;

	// recalc tempo indices based on new hop
	loTempo = x->x_loTempo;
	loTempo = 60.0f/loTempo; // BPM to seconds
	loTempo *= x->x_sr; // seconds to samples
	loTempo /= x->x_hop;//samples to frames
	x->x_loTempoIdx = roundf(loTempo);

   	hiTempo = x->x_hiTempo; 
	hiTempo = 60.0f/hiTempo; // BPM to seconds
	hiTempo *= x->x_sr; // seconds to samples
	hiTempo /= x->x_hop;//samples to frames
	x->x_hiTempoIdx = roundf(hiTempo);
	
    post("%s frame hop: %i", x->x_objSymbol->s_name, x->x_hop);
}


static void *tempo_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    t_tempo_tilde *x = (t_tempo_tilde *)pd_new(tempo_tilde_class);
	t_float sepFloat, loTempo, hiTempo;
	t_sampIdx i;
	
	x->x_tempo = outlet_new(&x->x_obj, &s_float);
	x->x_onsetList = outlet_new(&x->x_obj, gensym("list"));

	// store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
	x->x_objSymbol = s;

	switch(argc)
	{
		case 2:
			x->x_window = atom_getfloat(argv);
			if(x->x_window<MINWINDOWSIZE)
			{
				x->x_window = WINDOWSIZEDEFAULT;
				post("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, MINWINDOWSIZE, WINDOWSIZEDEFAULT);
			}
			
			sepFloat = atom_getfloat(argv+1);
			if(sepFloat > x->x_window)
			{
				post("%s WARNING: frame hop cannot be more than current window size. Using half of current window size instead.", x->x_objSymbol->s_name);
				x->x_hop = x->x_window*0.5;
			}
			else if(sepFloat < 0)
			{
				post("%s WARNING: frame hop must be > 0. Using half of current window size instead.", x->x_objSymbol->s_name);
				x->x_hop = x->x_window*0.5;
			}
			else
				x->x_hop = sepFloat;
			break;
	
		case 1:
			x->x_window = atom_getfloat(argv);
			if(x->x_window<MINWINDOWSIZE)
			{
				x->x_window = WINDOWSIZEDEFAULT;
				post("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, MINWINDOWSIZE, WINDOWSIZEDEFAULT);
			}

			x->x_hop = x->x_window*0.5;
			break;

		case 0:
			x->x_window = WINDOWSIZEDEFAULT;
			x->x_hop = x->x_window*0.5;
			break;
						
		default:
			post("%s WARNING: Too many arguments supplied. Using default window size of %i, and frame hop of %i.", x->x_objSymbol->s_name, WINDOWSIZEDEFAULT, (t_sampIdx)(WINDOWSIZEDEFAULT*0.5));
			x->x_window = WINDOWSIZEDEFAULT;
			x->x_hop = x->x_window*0.5;
			break;
	}

	x->x_windowHalf = x->x_window*0.5;
	x->x_loBin = 0;
	x->x_hiBin = x->x_windowHalf;
	x->x_growthThresh = 0.0;
	x->x_belowThreshDefault = 0.01;
	x->x_alignFirstPeak = true;
	x->x_numHarm = 6;
	
	x->x_sr = SAMPLERATEDEFAULT;
	x->x_n = BLOCKSIZEDEFAULT;
	x->x_overlap = 1;

	x->x_loTempo = 40;
	x->x_hiTempo = 160;

	loTempo = 60.0f/x->x_loTempo; // BPM to seconds
	loTempo *= x->x_sr; // seconds to samples
	loTempo /= x->x_hop;//samples to frames
	x->x_loTempoIdx = roundf(loTempo);

	hiTempo = 60.0f/x->x_hiTempo; // BPM to seconds
	hiTempo *= x->x_sr; // seconds to samples
	hiTempo /= x->x_hop;//samples to frames
	x->x_hiTempoIdx = roundf(hiTempo);
		
	x->x_windowFunction = blackman;
	x->x_normalize = false;
	x->x_powerSpectrum = false;
	x->x_lastDspTime = clock_getlogicaltime();
	x->x_squaredDiff = true;
	x->x_dspTicks = 0;
	
	x->x_onsetsBufSize = x->x_loTempoIdx*(x->x_numHarm+1);
	
	x->x_signalBuffer = (t_sample *)t_getbytes((x->x_window*2) * sizeof(t_sample));
	x->x_fftwInForwardWindow = (t_float *)t_getbytes(x->x_window * sizeof(t_float));
	x->x_fftwInBackWindow = (t_float *)t_getbytes(x->x_window * sizeof(t_float));
	x->x_onsetsBuffer = (t_float *)t_getbytes(x->x_onsetsBufSize * sizeof(t_float));
	
 	for(i=0; i<(x->x_window*2); i++)
		x->x_signalBuffer[i] = 0.0;

 	for(i=0; i<x->x_onsetsBufSize; i++)
		x->x_onsetsBuffer[i] = 0.0;
		
  	x->x_blackman = (t_float *)t_getbytes(x->x_window*sizeof(t_float));
  	x->x_cosine = (t_float *)t_getbytes(x->x_window*sizeof(t_float));
  	x->x_hamming = (t_float *)t_getbytes(x->x_window*sizeof(t_float));
  	x->x_hann = (t_float *)t_getbytes(x->x_window*sizeof(t_float));

 	// initialize signal windowing functions
	tIDLib_blackmanWindow(x->x_blackman, x->x_window);
	tIDLib_cosineWindow(x->x_cosine, x->x_window);
	tIDLib_hammingWindow(x->x_hamming, x->x_window);
	tIDLib_hannWindow(x->x_hann, x->x_window);

	// set up the FFTW output buffer. Is there no function to initialize it?
	x->x_fftwOutForwardWindow = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf+1);
	x->x_fftwOutBackWindow = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf+1);

	// DFT plan
	x->x_fftwPlanForwardWindow = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwInForwardWindow, x->x_fftwOutForwardWindow, FFTWPLANNERFLAG);
	x->x_fftwPlanBackWindow = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwInBackWindow, x->x_fftwOutBackWindow, FFTWPLANNERFLAG);

	// we're supposed to initialize the input array after we create the plan
 	for(i=0; i<x->x_window; i++)
 	{
		x->x_fftwInForwardWindow[i] = 0.0;
		x->x_fftwInBackWindow[i] = 0.0;
	}

	// create listOut memory
	x->x_listOut = (t_atom *)t_getbytes(x->x_onsetsBufSize*sizeof(t_atom));

    return (x);
}


static t_int *tempo_tilde_perform(t_int *w)
{
    t_uShortInt n;
	t_sampIdx i;
	
    t_tempo_tilde *x = (t_tempo_tilde *)(w[1]);

    t_sample *in = (t_float *)(w[2]);
    n = w[3];
 			
 	// shift signal buffer contents back.
	for(i=0; i<(x->x_window*2-n); i++)
		x->x_signalBuffer[i] = x->x_signalBuffer[i+n];
	
	// write new block to end of signal buffer.
	for(i=0; i<n; i++)
		x->x_signalBuffer[(x->x_window*2-n)+i] = in[i];

	x->x_dspTicks++;
	
 	if(x->x_dspTicks*n >= x->x_hop)
 	{
 		x->x_dspTicks = 0;
 		tempo_tilde_analyze(x);
 	}
 	
	x->x_lastDspTime = clock_getlogicaltime();

    return (w+4);
}


static void tempo_tilde_dsp(t_tempo_tilde *x, t_signal **sp)
{
	dsp_add(
		tempo_tilde_perform,
		3,
		x,
		sp[0]->s_vec,
		sp[0]->s_n
	);

// compare sr to stored sr and update if different
	if( sp[0]->s_sr != (x->x_sr*x->x_overlap) )
	{
		x->x_sr = sp[0]->s_sr/x->x_overlap;
	};

// compare n to stored n and update if different
	if( sp[0]->s_n != x->x_n )
	{
		x->x_n = sp[0]->s_n;
		x->x_lastDspTime = clock_getlogicaltime();
	}
};

static void tempo_tilde_free(t_tempo_tilde *x)
{	
	// free the input buffer memory
    t_freebytes(x->x_signalBuffer, (x->x_window*2)*sizeof(t_sample));
    t_freebytes(x->x_onsetsBuffer, x->x_onsetsBufSize*sizeof(t_float));

	// free FFTW stuff
    t_freebytes(x->x_fftwInForwardWindow, (x->x_window)*sizeof(t_float));
    t_freebytes(x->x_fftwInBackWindow, (x->x_window)*sizeof(t_float));
	fftwf_free(x->x_fftwOutForwardWindow);
	fftwf_free(x->x_fftwOutBackWindow);
	fftwf_destroy_plan(x->x_fftwPlanForwardWindow); 
	fftwf_destroy_plan(x->x_fftwPlanBackWindow); 
	
	// free the window memory
	t_freebytes(x->x_blackman, x->x_window*sizeof(t_float));
	t_freebytes(x->x_cosine, x->x_window*sizeof(t_float));
	t_freebytes(x->x_hamming, x->x_window*sizeof(t_float));
	t_freebytes(x->x_hann, x->x_window*sizeof(t_float));

	t_freebytes(x->x_listOut, x->x_onsetsBufSize*sizeof(t_atom));
}

void tempo_tilde_setup(void)
{
    tempo_tilde_class = 
    class_new(
    	gensym("tempo~"),
    	(t_newmethod)tempo_tilde_new,
    	(t_method)tempo_tilde_free,
        sizeof(t_tempo_tilde),
        CLASS_DEFAULT, 
        A_GIMME,
		0
    );

	class_addcreator(
		(t_newmethod)tempo_tilde_new,
		gensym("timbreIDLib/tempo~"),
		A_GIMME,
		0
	);

    CLASS_MAINSIGNALIN(tempo_tilde_class, t_tempo_tilde, x_f);

	class_addmethod(
		tempo_tilde_class,
		(t_method)tempo_tilde_print,
		gensym("print"),
		0
	);
	
	class_addmethod(
		tempo_tilde_class, 
        (t_method)tempo_tilde_window,
		gensym("window"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		tempo_tilde_class, 
        (t_method)tempo_tilde_hop,
		gensym("hop"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		tempo_tilde_class, 
        (t_method)tempo_tilde_overlap,
		gensym("overlap"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tempo_tilde_class, 
        (t_method)tempo_tilde_squaredDiff, 
		gensym("squared_diff"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tempo_tilde_class, 
        (t_method)tempo_tilde_normalize, 
		gensym("normalize"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tempo_tilde_class,
        (t_method)tempo_tilde_x_windowFunction,
		gensym("window_function"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tempo_tilde_class, 
        (t_method)tempo_tilde_powerSpectrum,
		gensym("power_spectrum"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tempo_tilde_class, 
        (t_method)tempo_tilde_alignFirstPeak,
		gensym("align_first_peak"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		tempo_tilde_class, 
        (t_method)tempo_tilde_thresh,
		gensym("thresh"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tempo_tilde_class, 
        (t_method)tempo_tilde_freqRange,
		gensym("freq_range"),
		A_DEFFLOAT,
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tempo_tilde_class, 
        (t_method)tempo_tilde_numHarm,
		gensym("harmonics"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		tempo_tilde_class, 
        (t_method)tempo_tilde_tempoRange,
		gensym("tempo_range"),
		A_DEFFLOAT,
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tempo_tilde_class, 
        (t_method)tempo_tilde_onsetsBufDur,
		gensym("onset_buf_duration"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tempo_tilde_class, 
        (t_method)tempo_tilde_belowThreshDefault,
		gensym("below_thresh_default"),
		A_DEFFLOAT,
		0
	);

    class_addmethod(
    	tempo_tilde_class,
    	(t_method)tempo_tilde_dsp,
    	gensym("dsp"),
    	0
    );
}
