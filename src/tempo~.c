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
    t_sampIdx x_window;
    t_sampIdx x_windowHalf;
    t_windowFunction x_windowFunction;
    t_uShortInt x_overlap;
    t_bool x_normalize;
    t_bool x_powerSpectrum;
    t_float x_growthThresh;
    t_binIdx x_loBin;
    t_binIdx x_hiBin;
    t_uInt x_onsetsBufSize;
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
	t_uInt x_hop;
	t_atom *x_listOut;
    t_outlet *x_onsetList;
    t_outlet *x_tempo;
    t_float x_f;

} t_tempo_tilde;


/* ------------------------ tempo~ -------------------------------- */

static void tempo_tilde_bang(t_tempo_tilde *x)
{
    t_sampIdx i, j, startIdx, window, windowHalf;
    t_float growth, *windowFuncPtr;
	double currentTime;

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

	if(x->x_normalize)
		tIDLib_normal(windowHalf+1, x->x_fftwInBackWindow);
		
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
		x->x_onsetsBuffer[x->x_onsetsBufSize-1] = 1.0; //set below thresh growth to 1.0 and NOT ZERO, which would ruin the running product in the HPS algo
	
	// align first peak to beginning?
	startIdx=0;
	while(x->x_onsetsBuffer[startIdx]==1.0)
	{
		startIdx++;

		if(startIdx>=x->x_onsetsBufSize-2)
			break;
	}
	
	for(i=startIdx, j=0; i<x->x_onsetsBufSize; i++, j++)
		SETFLOAT(x->x_listOut+j, x->x_onsetsBuffer[i]);

	for(; j<x->x_onsetsBufSize; j++)
		SETFLOAT(x->x_listOut+j, 1.0);
	

 	outlet_list(x->x_onsetList, 0, x->x_onsetsBufSize, x->x_listOut);
	outlet_float(x->x_tempo, growth);
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


static void tempo_tilde_hop(t_tempo_tilde *x, t_floatarg s)
{
	if(s > x->x_window)
	{
		post("%s frame hop cannot be more than current window size. Using half of current window size instead.", x->x_objSymbol->s_name);
		x->x_hop = x->x_windowHalf;
	}
	else if(s < 0)
	{
		post("%s frame hop must be > 0. Using half of current window size instead.", x->x_objSymbol->s_name);
		x->x_hop = x->x_windowHalf;
	}
	else
		x->x_hop = s;
		
    post("%s frame hop: %i", x->x_objSymbol->s_name, x->x_hop);
}


static void *tempo_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    t_tempo_tilde *x = (t_tempo_tilde *)pd_new(tempo_tilde_class);
	t_float sepFloat;
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
	x->x_onsetsBufSize = 600;
	x->x_sr = SAMPLERATEDEFAULT;
	x->x_n = BLOCKSIZEDEFAULT;
	x->x_overlap = 1;
	x->x_windowFunction = blackman;
	x->x_normalize = false;
	x->x_powerSpectrum = false;
	x->x_lastDspTime = clock_getlogicaltime();
	x->x_squaredDiff = true;
	x->x_dspTicks = 0;
	
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
 		tempo_tilde_bang(x);
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

	class_addbang(tempo_tilde_class, tempo_tilde_bang);

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

// TODO
/*
	class_addmethod(
		tempo_tilde_class, 
        (t_method)tempo_tilde_onsetsBufSize,
		gensym("onset_buf_size"),
		A_DEFFLOAT,
		0
	);
*/
	
    class_addmethod(
    	tempo_tilde_class,
    	(t_method)tempo_tilde_dsp,
    	gensym("dsp"),
    	0
    );
}
