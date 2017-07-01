/*

specFlux~

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.


version 0.7, September 26, 2016

- using FFTW now

*/

#include "tIDLib.h"

static t_class *specFlux_tilde_class;

typedef struct _specFlux_tilde
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
    double x_lastDspTime;
    t_sample *x_signalBuffer;
    t_float *x_fftwInForwardWindow;
    t_float *x_fftwInBackWindow;
    fftwf_complex *x_fftwOutForwardWindow;
    fftwf_complex *x_fftwOutBackWindow;
	fftwf_plan x_fftwPlanForwardWindow;
	fftwf_plan x_fftwPlanBackWindow;
    t_float *x_blackman;
    t_float *x_cosine;
    t_float *x_hamming;
    t_float *x_hann;
    t_bool x_squaredDiff;
	t_uInt x_separation;
	t_atom *x_listOut;
    t_outlet *x_fluxList;
    t_outlet *x_flux;
    t_float x_f;

} t_specFlux_tilde;


/* ------------------------ specFlux~ -------------------------------- */

static void specFlux_tilde_bang(t_specFlux_tilde *x)
{
    t_sampIdx i, j, window, windowHalf, bangSample;
    t_float flux, *windowFuncPtr;
	double currentTime;

    window = x->x_window;
    windowHalf = x->x_windowHalf;

	currentTime = clock_gettimesince(x->x_lastDspTime);
	bangSample = roundf((currentTime/1000.0)*x->x_sr);

	if(bangSample >= x->x_n)
        bangSample = x->x_n-1;
            
	// construct forward analysis window
	for(i=0, j=bangSample; i<window; i++, j++)
		x->x_fftwInForwardWindow[i] = x->x_signalBuffer[window + j];

	// construct back analysis window x->x_separation frames earlier
	for(i=0, j=bangSample; i<window; i++, j++)
		x->x_fftwInBackWindow[i] = x->x_signalBuffer[window - x->x_separation + j];
	
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
		
	flux=0;

	for(i=0; i<=windowHalf; i++)
	{
		t_float diff, val;

		diff = x->x_fftwInForwardWindow[i] - x->x_fftwInBackWindow[i];
		
		if(x->x_squaredDiff)
			val = diff*diff;
		else
			val = fabs(diff);
			
		SETFLOAT(x->x_listOut+i, diff);
		flux += val;
	}

 	outlet_list(x->x_fluxList, 0, windowHalf+1, x->x_listOut);
	outlet_float(x->x_flux, flux);
}


static void specFlux_tilde_print(t_specFlux_tilde *x)
{
	post("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr/x->x_overlap));
	post("%s block size: %i", x->x_objSymbol->s_name, (t_uShortInt)x->x_n);
	post("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
	post("%s window: %i", x->x_objSymbol->s_name, x->x_window);
	post("%s normalize: %i", x->x_objSymbol->s_name, x->x_normalize);
	post("%s power spectrum: %i", x->x_objSymbol->s_name, x->x_powerSpectrum);
	post("%s window function: %i", x->x_objSymbol->s_name, x->x_windowFunction);
	post("%s separation: %i", x->x_objSymbol->s_name, x->x_separation);
	post("%s squared difference: %i", x->x_objSymbol->s_name, x->x_squaredDiff);
}


static void specFlux_tilde_window(t_specFlux_tilde *x, t_floatarg w)
{
	t_sampIdx i, window, windowHalf;
	
	window = w;
	
	if(window<MINWINDOWSIZE)
	{
		window = WINDOWSIZEDEFAULT;
		post("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, MINWINDOWSIZE, WINDOWSIZEDEFAULT);
	}

	if(x->x_separation > window)
	{
		post("%s change in window size caused frame separation to be less than %i samples apart. Setting frame separation to half of current window size instead.", x->x_objSymbol->s_name, window);
		x->x_separation = window*0.5;
	}
	
	windowHalf = window*0.5;
	
	x->x_signalBuffer = (t_sample *)t_resizebytes(x->x_signalBuffer, (x->x_window*2+x->x_n) * sizeof(t_sample), (window*2+x->x_n) * sizeof(t_sample));
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
	x->x_fftwPlanForwardWindow = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwInForwardWindow, x->x_fftwOutForwardWindow, FFTW_MEASURE);
	x->x_fftwPlanBackWindow = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwInBackWindow, x->x_fftwOutBackWindow, FFTW_MEASURE);

 	// we're supposed to initialize the input array after we create the plan
	for(i=0; i<x->x_window; i++)
 	{
		x->x_fftwInForwardWindow[i] = 0.0;
		x->x_fftwInBackWindow[i] = 0.0;
	}
	
	// initialize signal buffer
	for(i=0; i<(x->x_window*2+x->x_n); i++)
		x->x_signalBuffer[i] = 0.0;

	// re-init window functions
	tIDLib_blackmanWindow(x->x_blackman, x->x_window);
	tIDLib_cosineWindow(x->x_cosine, x->x_window);
	tIDLib_hammingWindow(x->x_hamming, x->x_window);
	tIDLib_hannWindow(x->x_hann, x->x_window);

	post("%s window size: %i", x->x_objSymbol->s_name, x->x_window);
}


static void specFlux_tilde_overlap(t_specFlux_tilde *x, t_floatarg o)
{
	// this change will be picked up the next time _dsp is called, where the samplerate will be updated to sp[0]->s_sr/x->x_overlap;
	x->x_overlap = (o<1)?1:o;

    post("%s overlap: %i", x->x_objSymbol->s_name, x->x_overlap);
}


static void specFlux_tilde_x_windowFunction(t_specFlux_tilde *x, t_floatarg f)
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


static void specFlux_tilde_powerSpectrum(t_specFlux_tilde *x, t_floatarg spec)
{
    spec = (spec<0)?0:spec;
    spec = (spec>1)?1:spec;
	x->x_powerSpectrum = spec;

	if(x->x_powerSpectrum)
		post("%s using power spectrum", x->x_objSymbol->s_name);
	else
		post("%s using magnitude spectrum", x->x_objSymbol->s_name);
}


static void specFlux_tilde_squaredDiff(t_specFlux_tilde *x, t_floatarg sd)
{
    sd = (sd<0)?0:sd;
    sd = (sd>1)?1:sd;
	x->x_squaredDiff = sd;
	
	if(x->x_squaredDiff)
		post("%s using squared difference", x->x_objSymbol->s_name);
	else
		post("%s using absolute value of difference", x->x_objSymbol->s_name);
}


static void specFlux_tilde_normalize(t_specFlux_tilde *x, t_floatarg norm)
{
    norm = (norm<0)?0:norm;
    norm = (norm>1)?1:norm;
	x->x_normalize = norm;

	if(x->x_normalize)
		post("%s spectrum normalization ON.", x->x_objSymbol->s_name);
	else
		post("%s spectrum normalization OFF.", x->x_objSymbol->s_name);
}


static void specFlux_tilde_separation(t_specFlux_tilde *x, t_floatarg s)
{
	if(s > x->x_window)
	{
		post("%s frame separation cannot be more than current window size. Using half of current window size instead.", x->x_objSymbol->s_name);
		x->x_separation = x->x_windowHalf;
	}
	else if(s < 0)
	{
		post("%s frame separation must be > 0. Using half of current window size instead.", x->x_objSymbol->s_name);
		x->x_separation = x->x_windowHalf;
	}
	else
		x->x_separation = s;
		
    post("%s frame separation: %i", x->x_objSymbol->s_name, x->x_separation);
}


static void *specFlux_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
    t_specFlux_tilde *x = (t_specFlux_tilde *)pd_new(specFlux_tilde_class);
	t_float sepFloat;
	t_sampIdx i;
	
	x->x_flux = outlet_new(&x->x_obj, &s_float);
	x->x_fluxList = outlet_new(&x->x_obj, gensym("list"));

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
				post("%s WARNING: frame separation cannot be more than current window size. Using half of current window size instead.", x->x_objSymbol->s_name);
				x->x_separation = x->x_window*0.5;
			}
			else if(sepFloat < 0)
			{
				post("%s WARNING: frame separation must be > 0. Using half of current window size instead.", x->x_objSymbol->s_name);
				x->x_separation = x->x_window*0.5;
			}
			else
				x->x_separation = sepFloat;
			break;
	
		case 1:
			x->x_window = atom_getfloat(argv);
			if(x->x_window<MINWINDOWSIZE)
			{
				x->x_window = WINDOWSIZEDEFAULT;
				post("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, MINWINDOWSIZE, WINDOWSIZEDEFAULT);
			}

			x->x_separation = x->x_window*0.5;
			break;

		case 0:
			x->x_window = WINDOWSIZEDEFAULT;
			x->x_separation = x->x_window*0.5;
			break;
						
		default:
			post("%s WARNING: Too many arguments supplied. Using default window size of %i, and frame separation of %i.", x->x_objSymbol->s_name, WINDOWSIZEDEFAULT, (t_sampIdx)(WINDOWSIZEDEFAULT*0.5));
			x->x_window = WINDOWSIZEDEFAULT;
			x->x_separation = x->x_window*0.5;
			break;
	}

	x->x_windowHalf = x->x_window*0.5;
	x->x_sr = SAMPLERATEDEFAULT;
	x->x_n = BLOCKSIZEDEFAULT;
	x->x_overlap = 1;
	x->x_windowFunction = blackman;
	x->x_normalize = true;
	x->x_powerSpectrum = false;
	x->x_lastDspTime = clock_getlogicaltime();
	x->x_squaredDiff = false; // absolute value by default

	x->x_signalBuffer = (t_sample *)t_getbytes((x->x_window*2+x->x_n) * sizeof(t_sample));
	x->x_fftwInForwardWindow = (t_float *)t_getbytes(x->x_window * sizeof(t_float));
	x->x_fftwInBackWindow = (t_float *)t_getbytes(x->x_window * sizeof(t_float));
	
 	for(i=0; i<(x->x_window*2+x->x_n); i++)
		x->x_signalBuffer[i] = 0.0;

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
	x->x_fftwPlanForwardWindow = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwInForwardWindow, x->x_fftwOutForwardWindow, FFTW_MEASURE);
	x->x_fftwPlanBackWindow = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwInBackWindow, x->x_fftwOutBackWindow, FFTW_MEASURE);

	// we're supposed to initialize the input array after we create the plan
 	for(i=0; i<x->x_window; i++)
 	{
		x->x_fftwInForwardWindow[i] = 0.0;
		x->x_fftwInBackWindow[i] = 0.0;
	}

	// create listOut memory
	x->x_listOut = (t_atom *)t_getbytes((x->x_windowHalf+1)*sizeof(t_atom));

    return (x);
}


static t_int *specFlux_tilde_perform(t_int *w)
{
    t_uShortInt n;
	t_sampIdx i;
	
    t_specFlux_tilde *x = (t_specFlux_tilde *)(w[1]);

    t_sample *in = (t_float *)(w[2]);
    n = w[3];
 			
 	// shift signal buffer contents back.
	for(i=0; i<(x->x_window*2); i++)
		x->x_signalBuffer[i] = x->x_signalBuffer[i+n];
	
	// write new block to end of signal buffer.
	for(i=0; i<n; i++)
		x->x_signalBuffer[x->x_window*2+i] = in[i];
		
	x->x_lastDspTime = clock_getlogicaltime();

    return (w+4);
}


static void specFlux_tilde_dsp(t_specFlux_tilde *x, t_signal **sp)
{
	dsp_add(
		specFlux_tilde_perform,
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

// compare n to stored n and update/resize buffer if different
	if( sp[0]->s_n != x->x_n )
	{
		t_sampIdx i;

		x->x_signalBuffer = (t_sample *)t_resizebytes(x->x_signalBuffer, (x->x_window*2+x->x_n) * sizeof(t_sample), (x->x_window*2+sp[0]->s_n) * sizeof(t_sample));

		x->x_n = sp[0]->s_n;
		x->x_lastDspTime = clock_getlogicaltime();

		// init signal buffer
		for(i=0; i<(x->x_window*2+x->x_n); i++)
			x->x_signalBuffer[i] = 0.0;
	}
};

static void specFlux_tilde_free(t_specFlux_tilde *x)
{	
	// free the input buffer memory
    t_freebytes(x->x_signalBuffer, (x->x_window*2+x->x_n)*sizeof(t_sample));

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

	t_freebytes(x->x_listOut, (x->x_windowHalf+1)*sizeof(t_atom));
}

void specFlux_tilde_setup(void)
{
    specFlux_tilde_class = 
    class_new(
    	gensym("specFlux~"),
    	(t_newmethod)specFlux_tilde_new,
    	(t_method)specFlux_tilde_free,
        sizeof(t_specFlux_tilde),
        CLASS_DEFAULT, 
        A_GIMME,
		0
    );

    CLASS_MAINSIGNALIN(specFlux_tilde_class, t_specFlux_tilde, x_f);

	class_addbang(specFlux_tilde_class, specFlux_tilde_bang);

	class_addmethod(
		specFlux_tilde_class,
		(t_method)specFlux_tilde_print,
		gensym("print"),
		0
	);
	
	class_addmethod(
		specFlux_tilde_class, 
        (t_method)specFlux_tilde_window,
		gensym("window"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		specFlux_tilde_class, 
        (t_method)specFlux_tilde_separation,
		gensym("separation"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		specFlux_tilde_class, 
        (t_method)specFlux_tilde_overlap,
		gensym("overlap"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		specFlux_tilde_class, 
        (t_method)specFlux_tilde_squaredDiff, 
		gensym("squared_diff"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		specFlux_tilde_class, 
        (t_method)specFlux_tilde_normalize, 
		gensym("normalize"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		specFlux_tilde_class,
        (t_method)specFlux_tilde_x_windowFunction,
		gensym("window_function"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		specFlux_tilde_class, 
        (t_method)specFlux_tilde_powerSpectrum,
		gensym("power_spectrum"),
		A_DEFFLOAT,
		0
	);
	
    class_addmethod(
    	specFlux_tilde_class,
    	(t_method)specFlux_tilde_dsp,
    	gensym("dsp"),
    	0
    );
}