/*

specSlope

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.


version 0.7, September 26, 2016

- using FFTW now

*/

#include "tIDLib.h"

static t_class *specSlope_class;

typedef struct _specSlope
{
    t_object x_obj;
    t_symbol *x_objSymbol;
    t_float x_sr;
	t_sampIdx x_window;
	t_sampIdx x_windowHalf;
	t_windowFunction x_windowFunction;
	t_bool x_normalize;
	t_bool x_powerSpectrum;
	t_float *x_fftwIn;
    fftwf_complex *x_fftwOut;
	fftwf_plan x_fftwPlan;
	t_float *x_blackman;
	t_float *x_cosine;
	t_float *x_hamming;
	t_float *x_hann;
	t_word *x_vec;
	t_symbol *x_arrayName;
	t_sampIdx x_arrayPoints;
    t_outlet *x_slope;
} t_specSlope;


/* ------------------------ specSlope -------------------------------- */
static void specSlope_resizeWindow(t_specSlope *x, t_sampIdx oldWindow, t_sampIdx window, t_sampIdx startSamp, t_sampIdx *endSamp)
{			
	t_sampIdx windowHalf;
	
	windowHalf = window * 0.5;

	if(window<MINWINDOWSIZE)
	{
		window = WINDOWSIZEDEFAULT;
		windowHalf = window * 0.5;
		post("%s WARNING: window size must be %i or greater. Using default size of %i instead.", x->x_objSymbol->s_name, MINWINDOWSIZE, WINDOWSIZEDEFAULT);
		
		*endSamp = startSamp + window-1;
		if(*endSamp > x->x_arrayPoints)
			*endSamp = x->x_arrayPoints-1;
	}

	// hang on to these values for next time
	x->x_window = window;
	x->x_windowHalf = windowHalf;

	x->x_fftwIn = (t_float *)t_resizebytes(x->x_fftwIn, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));

	fftwf_free(x->x_fftwOut);
	fftwf_destroy_plan(x->x_fftwPlan); 
	// set up a new FFTW output buffer
	x->x_fftwOut = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf+1);
	// FFTW plan
	x->x_fftwPlan = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwIn, x->x_fftwOut, FFTW_MEASURE);

	x->x_blackman = (t_float *)t_resizebytes(x->x_blackman, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));
	x->x_cosine = (t_float *)t_resizebytes(x->x_cosine, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));
	x->x_hamming = (t_float *)t_resizebytes(x->x_hamming, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));
	x->x_hann = (t_float *)t_resizebytes(x->x_hann, oldWindow*sizeof(t_float), x->x_window*sizeof(t_float));

	tIDLib_blackmanWindow(x->x_blackman, x->x_window);
	tIDLib_cosineWindow(x->x_cosine, x->x_window);
	tIDLib_hammingWindow(x->x_hamming, x->x_window);
	tIDLib_hannWindow(x->x_hann, x->x_window);
}
		

static void specSlope_analyze(t_specSlope *x, t_floatarg start, t_floatarg n)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, j, window, startSamp, endSamp;
		t_float slope, *windowFuncPtr;

		startSamp = start;
		startSamp = (start<0)?0:start;

		if(n)
			endSamp = startSamp + n-1;
		else
			endSamp = startSamp + x->x_window-1;

		if(endSamp > x->x_arrayPoints)
			endSamp = x->x_arrayPoints-1;

		window = endSamp-startSamp+1;

		if(endSamp <= startSamp)
		{
			post("%s: bad range of samples.", x->x_objSymbol->s_name);
			return;
		}

		if(x->x_window != window)
			specSlope_resizeWindow(x, x->x_window, window, startSamp, &endSamp);
		
		// construct analysis window
		for(i=0, j=startSamp; j<=endSamp; i++, j++)
			x->x_fftwIn[i] = x->x_vec[j].w_float;

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

		// if windowFunction == 0, skip the windowing (rectangular)
		if(x->x_windowFunction!=rectangular)
			for(i=0; i<x->x_window; i++, windowFuncPtr++)
				x->x_fftwIn[i] *= *windowFuncPtr;

		fftwf_execute(x->x_fftwPlan);

		tIDLib_power(x->x_windowHalf+1, x->x_fftwOut, x->x_fftwIn);

		if(!x->x_powerSpectrum)
			tIDLib_mag(x->x_windowHalf+1, x->x_fftwIn);

		if(x->x_normalize)
			tIDLib_normal(x->x_windowHalf+1, x->x_fftwIn);
		
		slope = tIDLib_fitLineSlope(x->x_windowHalf+1, x->x_fftwIn);

		outlet_float(x->x_slope, slope);
	}
}


// analyze the whole damn array
static void specSlope_bang(t_specSlope *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx window, startSamp;
		startSamp = 0;
		window = x->x_arrayPoints;
		specSlope_analyze(x, startSamp, window);
	}
}


static void specSlope_set(t_specSlope *x, t_symbol *s)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(s, garray_class)))
		pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, s->s_name);
	else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
		pd_error(x, "%s: bad template for %s", s->s_name, x->x_objSymbol->s_name);
	else
	    x->x_arrayName = s;
}


static void specSlope_print(t_specSlope *x)
{
	post("%s array: %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
	post("%s samplerate: %i", x->x_objSymbol->s_name, (t_sampIdx)(x->x_sr));
	post("%s window: %i", x->x_objSymbol->s_name, x->x_window);
	post("%s normalize: %i", x->x_objSymbol->s_name, x->x_normalize);
	post("%s power spectrum: %i", x->x_objSymbol->s_name, x->x_powerSpectrum);
	post("%s window function: %i", x->x_objSymbol->s_name, x->x_windowFunction);
}


static void specSlope_samplerate(t_specSlope *x, t_floatarg sr)
{
	if(sr<MINSAMPLERATE)
		x->x_sr = MINSAMPLERATE;
	else
		x->x_sr = sr;
}


static void specSlope_windowFunction(t_specSlope *x, t_floatarg f)
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


static void specSlope_normalize(t_specSlope *x, t_floatarg norm)
{
    norm = (norm<0)?0:norm;
    norm = (norm>1)?1:norm;
	x->x_normalize = norm;

	if(x->x_normalize)
		post("%s spectrum normalization ON.", x->x_objSymbol->s_name);
	else
		post("%s spectrum normalization OFF.", x->x_objSymbol->s_name);
}


static void specSlope_powerSpectrum(t_specSlope *x, t_floatarg spec)
{
    spec = (spec<0)?0:spec;
    spec = (spec>1)?1:spec;
	x->x_powerSpectrum = spec;

	if(x->x_powerSpectrum)
		post("%s using power spectrum", x->x_objSymbol->s_name);
	else
		post("%s using magnitude spectrum", x->x_objSymbol->s_name);
}


static void *specSlope_new(t_symbol *s, int argc, t_atom *argv)
{
    t_specSlope *x = (t_specSlope *)pd_new(specSlope_class);
	t_sampIdx i;
//	t_garray *a;

	x->x_slope = outlet_new(&x->x_obj, gensym("list"));

	// store the pointer to the symbol containing the object name. Can access it for error and post functions via s->s_name
	x->x_objSymbol = s;
	
	switch(argc)
	{
		case 1:
			x->x_arrayName = atom_getsymbol(argv);
			/*
			if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
				pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
			else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
				pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
			*/
			break;
	
		case 0:
			post("%s: no array specified.", x->x_objSymbol->s_name);
			// a bogus array name to trigger the safety check in _analyze()
			x->x_arrayName = gensym("NOARRAYSPECIFIED");
			break;
		
		default:
			x->x_arrayName = atom_getsymbol(argv);
			/*
			if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
				pd_error(x, "%s: no array called %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
			else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
				pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
			*/
			post("%s WARNING: extra arguments ignored.", x->x_objSymbol->s_name);
			break;
	}
	
	x->x_sr = SAMPLERATEDEFAULT;
	x->x_window = WINDOWSIZEDEFAULT;
	x->x_windowHalf = x->x_window*0.5;
	x->x_windowFunction = blackman;
	x->x_normalize = true;
	x->x_powerSpectrum = false;
	
	x->x_fftwIn = (t_sample *)t_getbytes(x->x_window*sizeof(t_sample));

	// set up the FFTW output buffer. Is there no function to initialize it?
	x->x_fftwOut = (fftwf_complex *)fftwf_alloc_complex(x->x_windowHalf+1);

	// FFTW plan
	x->x_fftwPlan = fftwf_plan_dft_r2c_1d(x->x_window, x->x_fftwIn, x->x_fftwOut, FFTW_MEASURE); // FFTW_MEASURE may be slower than FFTW_ESTIMATE but more efficient after the first run?
	
	for(i=0; i<x->x_window; i++)
		x->x_fftwIn[i] = 0.0;

  	x->x_blackman = (t_float *)t_getbytes(x->x_window*sizeof(t_float));
  	x->x_cosine = (t_float *)t_getbytes(x->x_window*sizeof(t_float));
  	x->x_hamming = (t_float *)t_getbytes(x->x_window*sizeof(t_float));
  	x->x_hann = (t_float *)t_getbytes(x->x_window*sizeof(t_float));

 	// initialize signal windowing functions
	tIDLib_blackmanWindow(x->x_blackman, x->x_window);
	tIDLib_cosineWindow(x->x_cosine, x->x_window);
	tIDLib_hammingWindow(x->x_hamming, x->x_window);
	tIDLib_hannWindow(x->x_hann, x->x_window);

    return (x);
}


static void specSlope_free(t_specSlope *x)
{
	// free FFTW stuff
    t_freebytes(x->x_fftwIn, (x->x_window)*sizeof(t_float));
	fftwf_free(x->x_fftwOut);
	fftwf_destroy_plan(x->x_fftwPlan);
	
	// free the window memory
    t_freebytes(x->x_blackman, x->x_window*sizeof(t_float));
    t_freebytes(x->x_cosine, x->x_window*sizeof(t_float));
    t_freebytes(x->x_hamming, x->x_window*sizeof(t_float));
    t_freebytes(x->x_hann, x->x_window*sizeof(t_float));
}


void specSlope_setup(void)
{
    specSlope_class =
    class_new(
    	gensym("specSlope"),
    	(t_newmethod)specSlope_new,
    	(t_method)specSlope_free,
        sizeof(t_specSlope),
        CLASS_DEFAULT,
        A_GIMME,
		0
    );

	class_addbang(specSlope_class, specSlope_bang);

	class_addmethod(
		specSlope_class,
        (t_method)specSlope_analyze,
		gensym("analyze"),
        A_DEFFLOAT,
        A_DEFFLOAT,
		0
	);

	class_addmethod(
		specSlope_class,
		(t_method)specSlope_set,
		gensym("set"),
		A_SYMBOL,
		0
	);

	class_addmethod(
		specSlope_class,
		(t_method)specSlope_print,
		gensym("print"),
		0
	);

	class_addmethod(
		specSlope_class,
        (t_method)specSlope_samplerate,
		gensym("samplerate"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		specSlope_class,
        (t_method)specSlope_windowFunction,
		gensym("window_function"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		specSlope_class,
        (t_method)specSlope_normalize,
		gensym("normalize"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		specSlope_class,
        (t_method)specSlope_powerSpectrum,
		gensym("power_spectrum"),
		A_DEFFLOAT,
		0
	);
}
