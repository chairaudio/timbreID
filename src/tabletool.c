/*

tabletool - An array manipulation external.

Copyright 2010 William Brent

tabletool is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

tabletool is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.


version 0.7, September 26, 2016

*/


#include "tIDLib.h"
#include "g_canvas.h"

static t_class *tabletool_class;

typedef struct _tabletool
{
    t_object x_obj;
    t_symbol *x_objSymbol;
	t_word *x_vec;
    t_uInt x_randState;
	t_symbol *x_arrayName;
	t_bool x_storedFlag;
	t_float *x_originalData;
	t_sampIdx x_arrayPoints;
    t_outlet *x_info;
    t_outlet *x_list;

} t_tabletool;


/* ------------------------ tabletool -------------------------------- */
static void tabletool_dump(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		t_atom *listOut;

		listOut = (t_atom *)t_getbytes(x->x_arrayPoints*sizeof(t_atom));

		for(i=0; i<x->x_arrayPoints; i++)
			SETFLOAT(listOut+i, x->x_vec[i].w_float);

		outlet_list(x->x_list, 0, x->x_arrayPoints, listOut);

		// free local memory
		t_freebytes(listOut, x->x_arrayPoints * sizeof(t_atom));
 	}
}


static void tabletool_dumpRange(t_tabletool *x, t_float start, t_float finish)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_atom *listOut;
		t_sampIdx startIdx, finishIdx, range;
		
		start = (start<0)?0:start;
		start = (start>=x->x_arrayPoints)?x->x_arrayPoints-1:start;
		startIdx = start;

		finish = (finish<0)?0:finish;
		finish = (finish>=x->x_arrayPoints)?x->x_arrayPoints-1:finish;
		finishIdx = finish;
		
		if(startIdx>finishIdx)
		{
			t_sampIdx tmp;
			
			tmp = finishIdx;
			finishIdx = startIdx;
			startIdx = tmp;
		}
		
		range = finishIdx-startIdx+1;
		
		listOut = (t_atom *)t_getbytes(range*sizeof(t_atom));

		t_sampIdx i, j;
		
		for(i=0, j=startIdx; j<=finishIdx; i++, j++)
			SETFLOAT(listOut+i, x->x_vec[j].w_float);

		outlet_list(x->x_list, 0, range, listOut);

		// free local memory
		t_freebytes(listOut, range*sizeof(t_atom));
 	}
}


static void tabletool_size(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
		outlet_float(x->x_info, x->x_arrayPoints);
}


static void tabletool_range(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		t_float min, max, range;

		min = FLT_MAX;

		for(i=0; i<x->x_arrayPoints; i++)
			if(x->x_vec[i].w_float < min)
				min = x->x_vec[i].w_float;

		max = -FLT_MAX;

		for(i=0; i<x->x_arrayPoints; i++)
			if(x->x_vec[i].w_float > max)
				max = x->x_vec[i].w_float;
			
		range = max-min;
			
		outlet_float(x->x_info, range);
	}
}


static void tabletool_min(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		t_float min;
		t_atom indexOut;

		min = FLT_MAX;

		for(i=0; i<x->x_arrayPoints; i++)
			if(x->x_vec[i].w_float < min)
			{
				min = x->x_vec[i].w_float;
				SETFLOAT(&indexOut, i);
			}

		outlet_list(x->x_list, 0, 1, &indexOut);
		outlet_float(x->x_info, min);
	}
}


static void tabletool_max(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		t_float max;
		t_atom indexOut;

		max = -FLT_MAX;

		for(i=0; i<x->x_arrayPoints; i++)
			if(x->x_vec[i].w_float > max)
			{
				max = x->x_vec[i].w_float;
				SETFLOAT(&indexOut, i);
			}

		outlet_list(x->x_list, 0, 1, &indexOut);
		outlet_float(x->x_info, max);
	}
}


static void tabletool_maxMag(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		t_float max, maxVal;
		t_atom indexOut;

		max = -FLT_MAX;
		maxVal = 0;
		
		for(i=0; i<x->x_arrayPoints; i++)
			if(fabs(x->x_vec[i].w_float) > max)
			{
				max = fabs(x->x_vec[i].w_float);
				maxVal = x->x_vec[i].w_float;
				SETFLOAT(&indexOut, i);
			}

		outlet_list(x->x_list, 0, 1, &indexOut);
		outlet_float(x->x_info, maxVal);
	}
}


static void tabletool_nearest(t_tabletool *x, t_float val)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_atom indexOut;
		t_sampIdx i, idx;
		t_float dist;
	
		dist = FLT_MAX;
		idx = 0;

		for(i=0; i<x->x_arrayPoints; i++)
			if(fabs(x->x_vec[i].w_float - val) < dist)
			{
				dist = fabs(x->x_vec[i].w_float - val);
				SETFLOAT(&indexOut, i);
				idx = i;
			}

		outlet_list(x->x_list, 0, 1, &indexOut);
		outlet_float(x->x_info, x->x_vec[idx].w_float);
	}
}


static void tabletool_equals(t_tabletool *x, t_float val)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, numMatches, *indices;
			
		numMatches = 0;

		indices = (t_sampIdx *)t_getbytes(0);

		for(i=0; i<x->x_arrayPoints; i++)
			if(x->x_vec[i].w_float == val)
			{
				numMatches++;

				indices = (t_sampIdx *)t_resizebytes(indices, (numMatches-1)*sizeof(t_sampIdx), numMatches*sizeof(t_sampIdx));

				indices[numMatches-1] = i;
			}

		if(numMatches)
		{
			t_atom *listOut;
			
			listOut = (t_atom *)t_getbytes(numMatches*sizeof(t_atom));

			for(i=0; i<numMatches; i++)
				SETFLOAT(listOut+i, indices[i]);

			outlet_list(x->x_list, 0, numMatches, listOut);
			outlet_float(x->x_info, numMatches);

			// free local memory
			t_freebytes(listOut, numMatches * sizeof(t_atom));
		}
		else
			outlet_float(x->x_info, -1);

		// free local memory
		t_freebytes(indices, numMatches * sizeof(t_sampIdx));
 	}
}


static void tabletool_greater(t_tabletool *x, t_float val)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, numMatches, *indices;

		numMatches = 0;

		indices = (t_sampIdx *)t_getbytes(0);

		for(i=0; i<x->x_arrayPoints; i++)
			if(x->x_vec[i].w_float > val)
			{
				numMatches++;

				indices = (t_sampIdx *)t_resizebytes(indices, (numMatches-1)*sizeof(t_sampIdx), numMatches*sizeof(t_sampIdx));

				indices[numMatches-1] = i;
			}

		if(numMatches)
		{
			t_atom *listOut;

			listOut = (t_atom *)t_getbytes(numMatches*sizeof(t_atom));

			for(i=0; i<numMatches; i++)
				SETFLOAT(listOut+i, indices[i]);

			outlet_list(x->x_list, 0, numMatches, listOut);
			outlet_float(x->x_info, numMatches);

			// free local memory
			t_freebytes(listOut, numMatches * sizeof(t_atom));
		}
		else
			outlet_float(x->x_info, -1);

		// free local memory
		t_freebytes(indices, numMatches * sizeof(t_sampIdx));
 	}
}


static void tabletool_less(t_tabletool *x, t_float val)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, numMatches, *indices;

		numMatches = 0;

		indices = (t_sampIdx *)t_getbytes(0);

		for(i=0; i<x->x_arrayPoints; i++)
			if(x->x_vec[i].w_float < val)
			{
				numMatches++;

				indices = (t_sampIdx *)t_resizebytes(indices, (numMatches-1)*sizeof(t_sampIdx), numMatches*sizeof(t_sampIdx));

				indices[numMatches-1] = i;
			}

		if(numMatches)
		{
			t_atom *listOut;
		
			listOut = (t_atom *)t_getbytes(numMatches*sizeof(t_atom));

			for(i=0; i<numMatches; i++)
				SETFLOAT(listOut+i, indices[i]);

			outlet_list(x->x_list, 0, numMatches, listOut);
			outlet_float(x->x_info, numMatches);

			// free local memory
			t_freebytes(listOut, numMatches * sizeof(t_atom));
		}
		else
			outlet_float(x->x_info, -1);

		// free local memory
		t_freebytes(indices, numMatches * sizeof(t_sampIdx));
 	}
}


static void tabletool_between(t_tabletool *x, t_float lowBound, t_float hiBound)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, numMatches, *indices;
	
		numMatches = 0;

		indices = (t_sampIdx *)t_getbytes(0);

		for(i=0; i<x->x_arrayPoints; i++)
			if(x->x_vec[i].w_float > lowBound && x->x_vec[i].w_float < hiBound)
			{
				numMatches++;

				indices = (t_sampIdx *)t_resizebytes(indices, (numMatches-1)*sizeof(t_sampIdx), numMatches*sizeof(t_sampIdx));

				indices[numMatches-1] = i;
			}

		if(numMatches)
		{
			t_atom *listOut;

			listOut = (t_atom *)t_getbytes(numMatches*sizeof(t_atom));

			for(i=0; i<numMatches; i++)
				SETFLOAT(listOut+i, indices[i]);

			outlet_list(x->x_list, 0, numMatches, listOut);
			outlet_float(x->x_info, numMatches);

			// free local memory
			t_freebytes(listOut, numMatches * sizeof(t_atom));
		}
		else
			outlet_float(x->x_info, -1);

		// free local memory
		t_freebytes(indices, numMatches * sizeof(t_sampIdx));
 	}
}


static void tabletool_findZeroCrossings(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, numCrossings;
		t_float crossTotal;
		t_atom *listOut;
		
		listOut = (t_atom *)t_getbytes(0);
		
		numCrossings = 0;
		crossTotal = 0.0;
		
		for(i=1; i<x->x_arrayPoints; i++)
		{
			t_uChar crossing;
			crossing = abs(tIDLib_signum(x->x_vec[i].w_float) - tIDLib_signum(x->x_vec[i-1].w_float));
			
			crossTotal += crossing*0.5;
			
			if(crossing>0)
			{
				listOut = (t_atom *)t_resizebytes(listOut, numCrossings*(sizeof(t_atom)), (numCrossings+1)*(sizeof(t_atom)));
				SETFLOAT(listOut+numCrossings, i);
				numCrossings++;
			}
		}
		
		outlet_list(x->x_list, 0, numCrossings, listOut);
		outlet_float(x->x_info, crossTotal);
		
		// free the memory
		t_freebytes(listOut, numCrossings*sizeof(t_atom));
	}
}


static void tabletool_choose(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx nVal;
		t_uInt randVal;
		t_float choice;
		t_atom indexOut;

		randVal = x->x_randState;
		x->x_randState = randVal = randVal * 472940017 + 832416023; // from random in x_misc.c
		nVal = ((double)x->x_arrayPoints) * ((double)randVal)
			* (1./4294967296.);

		nVal = nVal % x->x_arrayPoints;
		SETFLOAT(&indexOut, nVal);

		choice = x->x_vec[nVal].w_float;

		outlet_list(x->x_list, 0, 1, &indexOut);
		outlet_float(x->x_info, choice);
	}
}


static void tabletool_const(t_tabletool *x, t_float val, t_float idx1, t_float idx2)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, intIdx1, intIdx2;

		intIdx1=(idx1<0)?0:idx1;
		intIdx1=(idx1>(x->x_arrayPoints-1))?(x->x_arrayPoints-1):idx1;
		
		intIdx2=(idx2<0)?0:idx2;
		intIdx2=(idx2>(x->x_arrayPoints-1))?(x->x_arrayPoints-1):idx2;

		if(intIdx1>intIdx2)
		{
			t_sampIdx tmp;
			tmp = intIdx2;
			intIdx2 = intIdx1;
			intIdx1 = tmp;
		}

		for(i=intIdx1; i<=intIdx2; i++)
			x->x_vec[i].w_float = val;

		garray_redraw(a);
	}
}


static void tabletool_series(t_tabletool *x, t_float startval, t_float endval, t_float exponent)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		float start, end, span, inc;
	
		exponent = (exponent<0)? 0 : exponent;

		start = startval;
		end = endval;
		span = end-start;
		inc = 1.0/(x->x_arrayPoints-1);

		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float = pow(inc*i, exponent) * span + start;

		garray_redraw(a);
	}
}


static void tabletool_randomWalk(t_tabletool *x, t_float start, t_float step, t_float lowlim, t_float uplim)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		double split;

		split = RAND_MAX * 0.5;

		x->x_vec[0].w_float = start;

		for(i=1; i<x->x_arrayPoints; i++)
		{
			if(rand()>split)
				x->x_vec[i].w_float = x->x_vec[i-1].w_float + step;
			else
				x->x_vec[i].w_float = x->x_vec[i-1].w_float - step;

			if(x->x_vec[i].w_float > uplim)
				x->x_vec[i].w_float = uplim - (x->x_vec[i].w_float-uplim);
			else if(x->x_vec[i].w_float < lowlim)
				x->x_vec[i].w_float = lowlim + (lowlim-x->x_vec[i].w_float);
		}

		garray_redraw(a);
	}
}


static void tabletool_offset(t_tabletool *x, t_float offset)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;

		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float += offset;

		garray_redraw(a);
 	}
}


static void tabletool_shift(t_tabletool *x, t_float s)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, j, lenTable;
		t_sLongInt shift;
		t_float *tableVals;

		tableVals = (t_float *)t_getbytes(x->x_arrayPoints*sizeof(t_float));

		shift = s;
		lenTable = x->x_arrayPoints;
		shift = shift % lenTable;

		for(i=0; i<lenTable; i++)
			tableVals[i] = x->x_vec[i].w_float;

		if(shift>0)
			for(i=shift, j=0; j<lenTable; i++, j++)
				x->x_vec[i%lenTable].w_float = tableVals[j];
		else // shift==0 and shift<0 handled by this case
			for(i=0, j=labs(shift); i<lenTable; i++, j++)
				x->x_vec[i].w_float = tableVals[j%lenTable];

		garray_redraw(a);

		// free local memory
		t_freebytes(tableVals, lenTable * sizeof(t_float));
 	}
}


static void tabletool_shift0(t_tabletool *x, t_float s)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, j, lenTable;
		t_sLongInt shift;
		t_float *tableVals;
	
		tableVals = (t_float *)t_getbytes(x->x_arrayPoints*sizeof(t_float));

		shift = s;
		lenTable = x->x_arrayPoints;
		shift = shift % lenTable;

		for(i=0; i<lenTable; i++)
			tableVals[i] = x->x_vec[i].w_float;

		if(shift>0)
		{
			// need to cast the signed long int shift variable to t_sampIdx to silence warning, because an unsigned long int can never be less than the lowest value of signed long int shift
			for(i=0; i<(t_sampIdx)shift; i++)
				x->x_vec[i].w_float = 0.0;

			for(i=shift, j=0; i<lenTable; i++, j++)
				x->x_vec[i].w_float = tableVals[j];
		}
		else
		{
			// shift==0 and shift<0 is handled by this case
			for(i=0, j=labs(shift); j<lenTable; i++, j++)
				x->x_vec[i].w_float = tableVals[j];

			for(; i<lenTable; i++)
				x->x_vec[i].w_float = 0.0;
		}		

		garray_redraw(a);

		// free local memory
		t_freebytes(tableVals, lenTable * sizeof(t_float));
 	}
}


static void tabletool_scale(t_tabletool *x, t_float scalar)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;

		// rescale
		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float *= scalar;

		garray_redraw(a);
	}
}


static void tabletool_invert(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		t_float max, min;

		max = -FLT_MAX;
		min = FLT_MAX;

		for(i=0; i<x->x_arrayPoints; i++)
		{
			if(x->x_vec[i].w_float < min)
				min = x->x_vec[i].w_float;

			if(x->x_vec[i].w_float > max)
				max = x->x_vec[i].w_float;
		}

		// rescale
		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float = max - x->x_vec[i].w_float + min;

		garray_redraw(a);
	}
}


static void tabletool_reverse(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, j;
		t_float *tableVals;

		tableVals = (t_float *)t_getbytes(x->x_arrayPoints*sizeof(t_float));

		// load table data
		for(i=0; i<x->x_arrayPoints; i++)
			tableVals[i] = x->x_vec[i].w_float;

		for(i=x->x_arrayPoints-1, j=0; j<x->x_arrayPoints; i--, j++)
			x->x_vec[j].w_float = tableVals[i];

		garray_redraw(a);

		// free local memory
		t_freebytes(tableVals, x->x_arrayPoints * sizeof(t_float));
    }
}


static void tabletool_remove(t_tabletool *x, t_float idx)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, idxInt, lastIdx;
		t_float val;
	
		lastIdx = x->x_arrayPoints-1;

		if(idx>lastIdx)
		{
	    	pd_error(x, "%s: index %lu out of range", x->x_objSymbol->s_name, (t_sampIdx)idx);
			return;
		}
		else if(idx<0)
		{
	    	pd_error(x, "%s: index %lu out of range", x->x_objSymbol->s_name, (t_sampIdx)idx);
			return;
		}
		else
			idxInt = idx;

		val = x->x_vec[idxInt].w_float;

		outlet_float(x->x_info, val);

		for(i=idxInt; i<lastIdx; i++)
			x->x_vec[i].w_float = x->x_vec[i+1].w_float;

		x->x_vec[lastIdx].w_float = 0.0;

		garray_redraw(a);
	}
}


static void tabletool_insert(t_tabletool *x, t_float idx, t_float val)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, idxInt, lastIdx;
		t_float *tableVals;
	
		tableVals = (t_float *)t_getbytes(x->x_arrayPoints*sizeof(t_float));

		// load table data
		for(i=0; i<x->x_arrayPoints; i++)
			tableVals[i] = x->x_vec[i].w_float;

		lastIdx = x->x_arrayPoints-1;

		if(idx>lastIdx)
		{
	    	pd_error(x, "%s: index %lu out of range", x->x_objSymbol->s_name, (t_sampIdx)idx);
			return;
		}
		else if(idx<0)
		{
	    	pd_error(x, "%s: index %lu out of range", x->x_objSymbol->s_name, (t_sampIdx)idx);
			return;
		}
		else
			idxInt = idx;

		outlet_float(x->x_info, x->x_vec[lastIdx].w_float);

		for(i=idxInt; i<lastIdx; i++)
			x->x_vec[i+1].w_float = tableVals[i];

		x->x_vec[idxInt].w_float = val;

		garray_redraw(a);

		// free local memory
		t_freebytes(tableVals, x->x_arrayPoints * sizeof(t_float));
	}
}


static void tabletool_round(t_tabletool *x, t_float res)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
		pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
	else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
		pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		t_float resRecip;

		if(res<=0)
			return;
	
		resRecip = 1.0/res;
		
		for(i=0; i<x->x_arrayPoints; i++)
		{
			x->x_vec[i].w_float *= resRecip;
			x->x_vec[i].w_float += res;
			x->x_vec[i].w_float = floor(x->x_vec[i].w_float);
			x->x_vec[i].w_float /= resRecip;
		}

		garray_redraw(a);
	}
}


static void tabletool_floor(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
		pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
	else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
		pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;

		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float = (t_sLongInt)x->x_vec[i].w_float;

		garray_redraw(a);
	}
}


static void tabletool_ceil(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
		pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
	else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
		pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;

		for(i=0; i<x->x_arrayPoints; i++)
			if( fabs(x->x_vec[i].w_float - (t_sLongInt)x->x_vec[i].w_float) > 0 )
				x->x_vec[i].w_float = 1 + (t_sLongInt)x->x_vec[i].w_float;

		garray_redraw(a);
	}
}


static void tabletool_mtof(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
		pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
	else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
		pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;

		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float = mtof(x->x_vec[i].w_float);

		garray_redraw(a);
	}
}


static void tabletool_ftom(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
		pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
	else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
		pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;

		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float = ftom(x->x_vec[i].w_float);

		garray_redraw(a);
	}
}


static void tabletool_dbtorms(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
		pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
	else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
		pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;

		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float = dbtorms(x->x_vec[i].w_float);

		garray_redraw(a);
	}
}


static void tabletool_rmstodb(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
		pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
	else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
		pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;

		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float = rmstodb(x->x_vec[i].w_float);

		garray_redraw(a);
	}
}


static void tabletool_bin2freq(t_tabletool *x, t_float n, t_float sr)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
		pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
	else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
		pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;

		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float = tIDLib_bin2freq(roundf(x->x_vec[i].w_float), n, sr);

		garray_redraw(a);
	}
}


static void tabletool_freq2bin(t_tabletool *x, t_float n, t_float sr)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
		pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
	else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
		pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;

		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float = tIDLib_freq2bin(x->x_vec[i].w_float, n, sr);

		garray_redraw(a);
	}
}


static void tabletool_bark2freq(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
		pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
	else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
		pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;

		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float = tIDLib_bark2freq(x->x_vec[i].w_float);

		garray_redraw(a);
	}
}


static void tabletool_freq2bark(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
		pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
	else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
		pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;

		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float = tIDLib_freq2bark(x->x_vec[i].w_float);

		garray_redraw(a);
	}
}


static void tabletool_clip(t_tabletool *x, t_float lower, t_float upper)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
		pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
	else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
		pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;

		for(i=0; i<x->x_arrayPoints; i++)
			if(x->x_vec[i].w_float > upper)
				x->x_vec[i].w_float = upper;
			else if(x->x_vec[i].w_float < lower)
				x->x_vec[i].w_float = lower;

		garray_redraw(a);
	}
}


static void tabletool_smooth(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
		pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
	else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
		pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		if(x->x_arrayPoints < 3)
			pd_error(x, "%s: table smoothing cannot be performed on tables with fewer than 3 elements", x->x_objSymbol->s_name);
		else
		{
			t_sampIdx i;
			t_float *tableVals, threeRecip;
	
			tableVals = (t_float *)t_getbytes(x->x_arrayPoints*sizeof(t_float));

			threeRecip = 1.0/3.0;

			// load table data
			for(i=0; i<x->x_arrayPoints; i++)
				tableVals[i] = x->x_vec[i].w_float;

			x->x_vec[0].w_float = (tableVals[0] + tableVals[1] + tableVals[2])*threeRecip;

			for(i=1; i<x->x_arrayPoints-1; i++)
				x->x_vec[i].w_float = (tableVals[i-1] + tableVals[i] + tableVals[i+1])*threeRecip;

			x->x_vec[x->x_arrayPoints-1].w_float = (tableVals[x->x_arrayPoints-1] + tableVals[x->x_arrayPoints-2] + tableVals[x->x_arrayPoints-3])*threeRecip;

			garray_redraw(a);

			// free local memory
			t_freebytes(tableVals, x->x_arrayPoints * sizeof(t_float));
		}
	}
}


static void tabletool_lace(t_tabletool *x, t_symbol *array2)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, j, array2pts, totalPoints;
		t_garray *b;
		t_atom *listOut;
		t_word *vec2;

		if(!(b = (t_garray *)pd_findbyclass(array2, garray_class)))
		{
			pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, array2->s_name);
			return;
		}
		else if(!garray_getfloatwords(b, (int *)&array2pts, &vec2))
		{
			pd_error(x, "%s: bad template for %s", array2->s_name, x->x_objSymbol->s_name);
			return;
		}

		totalPoints = x->x_arrayPoints + array2pts;
		listOut = (t_atom *)t_getbytes(totalPoints*sizeof(t_atom));

		for(i=0, j=0; i<(x->x_arrayPoints*2)-1; i+=2, j++)
		{
			if(j>(array2pts-1))
				break;

			SETFLOAT(listOut+i, x->x_vec[j].w_float);
			SETFLOAT(listOut+i+1, vec2[j].w_float);
		}

		for(; j<x->x_arrayPoints; i++, j++)
			SETFLOAT(listOut+i, x->x_vec[j].w_float);

		outlet_list(x->x_list, 0, totalPoints, listOut);

		// free local memory
		t_freebytes(listOut, totalPoints*sizeof(t_atom));
	}
}


static void tabletool_concatenate(t_tabletool *x, t_symbol *array2)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, j, array2pts, totalPoints;
		t_garray *b;
		t_word *vec2;
		t_atom *listOut;

		if(!(b = (t_garray *)pd_findbyclass(array2, garray_class)))
		{
			pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, array2->s_name);
			return;
		}
		else if(!garray_getfloatwords(b, (int *)&array2pts, &vec2))
		{
			pd_error(x, "%s: bad template for %s", array2->s_name, x->x_objSymbol->s_name);
			return;
		}

		totalPoints = x->x_arrayPoints + array2pts;
		listOut = (t_atom *)t_getbytes(totalPoints*sizeof(t_atom));

		for(i=0, j=0; i<x->x_arrayPoints; i++, j++)
			SETFLOAT(listOut+i, x->x_vec[j].w_float);


		for(j=0; i<totalPoints; i++, j++)
			SETFLOAT(listOut+i, vec2[j].w_float);

		outlet_list(x->x_list, 0, totalPoints, listOut);

		// free local memory
		t_freebytes(listOut, totalPoints*sizeof(t_atom));
	}
}


static void tabletool_scramble(t_tabletool *x)
{
    t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, length, range, nVal, *indices, *randIndices;
		t_float *tableVals;
		t_uInt randVal;

		tableVals = (t_float *)t_getbytes(x->x_arrayPoints*sizeof(t_float));

		indices = (t_sampIdx *)t_getbytes(x->x_arrayPoints*sizeof(t_sampIdx));

		randIndices = (t_sampIdx *)t_getbytes(x->x_arrayPoints*sizeof(t_sampIdx));

		for(i=0; i<x->x_arrayPoints; i++)
			indices[i] = i;

		for(i=0; i<x->x_arrayPoints; i++)
			tableVals[i] = x->x_vec[i].w_float;

		length = x->x_arrayPoints;

		for(i=0; i<length-1; i++)
		{
			t_float tmp;
			
			range = length-i;
			range = (range<1)?1:range;

			randVal = x->x_randState;
			x->x_randState = randVal = randVal * 472940017 + 832416023; // from random in x_misc.c
			nVal = ((double)range) * ((double)randVal)
				* (1.0/4294967296.0);

			nVal = nVal % range;
			nVal += i;

			tmp = tableVals[nVal];
			tableVals[nVal] = tableVals[i];
			tableVals[i] = tmp;
		}

		for(i=0; i<x->x_arrayPoints; i++)
			 x->x_vec[i].w_float = tableVals[i];

		garray_redraw(a);

		// free local memory
		t_freebytes(tableVals, x->x_arrayPoints * sizeof(t_float));
		t_freebytes(indices, x->x_arrayPoints * sizeof(t_sampIdx));
		t_freebytes(randIndices, x->x_arrayPoints * sizeof(t_sampIdx));
	}
}


static void tabletool_sort(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		t_float *tableVals;

		tableVals = (t_float *)t_getbytes(x->x_arrayPoints*sizeof(t_float));

		for(i=0; i<x->x_arrayPoints; i++)
			tableVals[i] = x->x_vec[i].w_float;

		tIDLib_bubbleSort(x->x_arrayPoints, tableVals);

		for(i=0; i<x->x_arrayPoints; i++)
			 x->x_vec[i].w_float = tableVals[i];

		garray_redraw(a);

		// free local memory
		t_freebytes(tableVals, x->x_arrayPoints * sizeof(t_float));
	}
}


static void tabletool_sort_range(t_tabletool *x, t_float start, t_float end)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, j, intStart, intEnd, length;
		t_float *tableVals;

		intStart = (start<0)?0:start;
		intStart = (start>=x->x_arrayPoints)?x->x_arrayPoints-1:start;
		intEnd = (end<0)?0:end;
		intEnd = (end>=x->x_arrayPoints)?x->x_arrayPoints-1:end;
		
		if(intStart>intEnd)
		{
			t_sampIdx tmp;
			tmp = intEnd;
			intEnd = intStart;
			intStart = tmp;
		}

		length = intEnd-intStart+1;

		tableVals = (t_float *)t_getbytes(length*sizeof(t_float));

		for(i=intStart, j=0; i<=intEnd; i++, j++)
			tableVals[j] = x->x_vec[i].w_float;

		tIDLib_bubbleSort(length, tableVals);

		for(i=intStart, j=0; i<=intEnd; i++, j++)
			 x->x_vec[i].w_float = tableVals[j];

		garray_redraw(a);

		// free local memory
		t_freebytes(tableVals, length * sizeof(t_float));
	}
}


static void tabletool_swap(t_tabletool *x, t_float idx1, t_float idx2)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i1, i2;
		t_float tmp;

		if(idx1<0 || idx1>x->x_arrayPoints-1 || idx2<0 || idx2>x->x_arrayPoints-1)
		{
			pd_error(x, "%s: index out of bounds", x->x_objSymbol->s_name);
			return;
		}
		else
		{
			i1 = idx1;
			i2 = idx2;
		}

		tmp = x->x_vec[i1].w_float;

		x->x_vec[i1].w_float = x->x_vec[i2].w_float;
		x->x_vec[i2].w_float = tmp;

		garray_redraw(a);
	}
}


static void tabletool_replace(t_tabletool *x, t_float searchVal, t_float newVal)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		
		for(i=0; i<x->x_arrayPoints; i++)
			if(x->x_vec[i].w_float==searchVal)
				x->x_vec[i].w_float = newVal;

		garray_redraw(a);
	}
}


static void tabletool_integrate(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		t_float sum;
		t_atom *listOut;

		listOut = (t_atom *)t_getbytes(x->x_arrayPoints*sizeof(t_atom));

		SETFLOAT(listOut, x->x_vec[0].w_float);

		sum = x->x_vec[0].w_float;

		for(i=1; i<x->x_arrayPoints; i++)
		{
			SETFLOAT(listOut+i, x->x_vec[i].w_float + sum);
			sum += x->x_vec[i].w_float;
		}

		outlet_list(x->x_list, 0, x->x_arrayPoints, listOut);

		// free local memory
		t_freebytes(listOut, x->x_arrayPoints * sizeof(t_atom));
 	}
}


static void tabletool_differentiate(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		t_atom *listOut;

		listOut = (t_atom *)t_getbytes(x->x_arrayPoints*sizeof(t_atom));

		SETFLOAT(listOut, x->x_vec[0].w_float);

		for(i=1; i<x->x_arrayPoints; i++)
			SETFLOAT(listOut+i, x->x_vec[i].w_float - x->x_vec[i-1].w_float);

		outlet_list(x->x_list, 0, x->x_arrayPoints, listOut);

		// free local memory
		t_freebytes(listOut, x->x_arrayPoints * sizeof(t_atom));
 	}
}


static void tabletool_add(t_tabletool *x, t_symbol *array2)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_garray *b;
		t_word *vec2;
		t_sampIdx i, array2pts, totalPoints;
		t_atom *listOut;

		if(!(b = (t_garray *)pd_findbyclass(array2, garray_class)))
		{
			pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, array2->s_name);
			return;
		}
		else if(!garray_getfloatwords(b, (int *)&array2pts, &vec2))
		{
			pd_error(x, "%s: bad template for %s", array2->s_name, x->x_objSymbol->s_name);
			return;
		}

		totalPoints = x->x_arrayPoints + array2pts;

		listOut = (t_atom *)t_getbytes(totalPoints*sizeof(t_atom));

		for(i=0; i<x->x_arrayPoints; i++)
			if(i>array2pts-1)
				break;
			else
				SETFLOAT(listOut+i, x->x_vec[i].w_float + vec2[i].w_float);

		outlet_list(x->x_list, 0, totalPoints, listOut);

		// free local memory
		t_freebytes(listOut, totalPoints * sizeof(t_atom));
	}
}


// undocumented. Also - shouldn't it just add directly to the table? And not output as a list at the right outlet? NOTE: not updated for timbreID 0.7 with new typedefs
/*
static void tabletool_add_range(t_tabletool *x, t_symbol *s, int argc, t_atom *argv)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		int i, j, sourcePts, ts, ss, se, length;
		t_garray *b;
		t_word *sourceVec;
		t_symbol *source;
		t_atom *listOut;

		ts = atom_getfloat(argv+0);
		source = atom_getsymbol(argv+1);
		ss = atom_getfloat(argv+2);
		se = atom_getfloat(argv+3);

		if(!(b = (t_garray *)pd_findbyclass(source, garray_class)))
		{
			pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, source->s_name);
			return;
		}
		else if(!garray_getfloatwords(b, (int *)&sourcePts, &sourceVec))
		{
			pd_error(x, "%s: bad template for %s", source->s_name, x->x_objSymbol->s_name);
			return;
		}

		if(ts<0)
			ts = 0;

		if(ss<0)
			ss = 0;

		if(se>sourcePts-1)
			se = sourcePts-1;

		length = se-ss+1;

		listOut = (t_atom *)t_getbytes(length*sizeof(t_atom));

		for(i=0, j=ss; j<=se; i++, j++)
			if(i > x->x_arrayPoints-1)
				break;
			else
				SETFLOAT(listOut+i, sourceVec[j].w_float + x->x_vec[ts+i].w_float);

		outlet_list(x->x_list, 0, length, listOut);

		// free local memory
		t_freebytes(listOut, length * sizeof(t_atom));
	}
}
*/

static void tabletool_subtract(t_tabletool *x, t_symbol *array2)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, array2pts, totalPoints;
		t_garray *b;
		t_word *vec2;
		t_atom *listOut;

		if(!(b = (t_garray *)pd_findbyclass(array2, garray_class)))
		{
			pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, array2->s_name);
			return;
		}
		else if(!garray_getfloatwords(b, (int *)&array2pts, &vec2))
		{
			pd_error(x, "%s: bad template for %s", array2->s_name, x->x_objSymbol->s_name);
			return;
		}

		totalPoints = x->x_arrayPoints + array2pts;

		listOut = (t_atom *)t_getbytes(totalPoints*sizeof(t_atom));

		for(i=0; i<x->x_arrayPoints; i++)
			if(i>array2pts-1)
				break;
			else
				SETFLOAT(listOut+i, x->x_vec[i].w_float - vec2[i].w_float);

		outlet_list(x->x_list, 0, totalPoints, listOut);

		// free local memory
		t_freebytes(listOut, totalPoints * sizeof(t_atom));
	}
}


static void tabletool_multiply(t_tabletool *x, t_symbol *array2)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_garray *b;
		t_word *vec2;
		t_atom *listOut;
		t_sampIdx i, array2pts, totalPoints;

		if(!(b = (t_garray *)pd_findbyclass(array2, garray_class)))
		{
			pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, array2->s_name);
			return;
		}
		else if(!garray_getfloatwords(b, (int *)&array2pts, &vec2))
		{
			pd_error(x, "%s: bad template for %s", array2->s_name, x->x_objSymbol->s_name);
			return;
		}

		totalPoints = x->x_arrayPoints + array2pts;

		listOut = (t_atom *)t_getbytes(totalPoints*sizeof(t_atom));

		for(i=0; i<x->x_arrayPoints; i++)
			if(i>array2pts-1)
				break;
			else
				SETFLOAT(listOut+i, x->x_vec[i].w_float * vec2[i].w_float);

		outlet_list(x->x_list, 0, totalPoints, listOut);

		// free local memory
		t_freebytes(listOut, totalPoints * sizeof(t_atom));
	}
}


static void tabletool_divide(t_tabletool *x, t_symbol *array2)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, array2pts, totalPoints;
		t_garray *b;
		t_word *vec2;
		t_atom *listOut;

		if(!(b = (t_garray *)pd_findbyclass(array2, garray_class)))
		{
			pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, array2->s_name);
			return;
		}
		else if(!garray_getfloatwords(b, (int *)&array2pts, &vec2))
		{
			pd_error(x, "%s: bad template for %s", array2->s_name, x->x_objSymbol->s_name);
			return;
		}

		totalPoints = x->x_arrayPoints + array2pts;

		listOut = (t_atom *)t_getbytes(totalPoints*sizeof(t_atom));

		for(i=0; i<x->x_arrayPoints; i++)
			if(i>array2pts-1)
				break;
			else
				if(vec2[i].w_float == 0.0)
				{
					pd_error(x, "%s: divide by zero detected", x->x_objSymbol->s_name);
					SETFLOAT(listOut+i, x->x_vec[i].w_float);
				}
				else
					SETFLOAT(listOut+i, x->x_vec[i].w_float / vec2[i].w_float);


		outlet_list(x->x_list, 0, totalPoints, listOut);

		// free local memory
		t_freebytes(listOut, totalPoints * sizeof(t_atom));
	}
}


static void tabletool_dot(t_tabletool *x, t_symbol *array1)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, array1pts, biggestN;
		t_garray *b;
		t_word *vec1;
		t_float dot, *vecBuffer, *vec1Buffer;

		if(!(b = (t_garray *)pd_findbyclass(array1, garray_class)))
		{
			pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, array1->s_name);
			return;
		}
		else if(!garray_getfloatwords(b, (int *)&array1pts, &vec1))
		{
			pd_error(x, "%s: bad template for %s", array1->s_name, x->x_objSymbol->s_name);
			return;
		}

		biggestN = 0;
		
		if(array1pts > x->x_arrayPoints)
			biggestN = array1pts;
		else
			biggestN = x->x_arrayPoints;

		vecBuffer = (t_float *)t_getbytes(biggestN*sizeof(t_float));
		vec1Buffer = (t_float *)t_getbytes(biggestN*sizeof(t_float));

		for(i=0; i<biggestN; i++)
		{
			vecBuffer[i] = 0.0;
			vec1Buffer[i] = 0.0;
		}

		for(i=0; i<biggestN; i++)
		{
			if(i<x->x_arrayPoints)
				vecBuffer[i] = x->x_vec[i].w_float;
			
			if(i<array1pts)
				vec1Buffer[i] = vec1[i].w_float;
		}
			
		dot = tIDLib_dotProd(biggestN, vecBuffer, vec1Buffer);

		// free local memory
		t_freebytes(vecBuffer, biggestN * sizeof(t_float));
		t_freebytes(vec1Buffer, biggestN * sizeof(t_float));

		outlet_float(x->x_info, dot);
	}
}


static void tabletool_euclid(t_tabletool *x, t_symbol *array1)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, array1pts, biggestN;
		t_garray *b;
		t_word *vec1;
		t_float dist, *vecBuffer, *vec1Buffer, *vecWeights;

		if(!(b = (t_garray *)pd_findbyclass(array1, garray_class)))
		{
			pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, array1->s_name);
			return;
		}
		else if(!garray_getfloatwords(b, (int *)&array1pts, &vec1))
		{
			pd_error(x, "%s: bad template for %s", array1->s_name, x->x_objSymbol->s_name);
			return;
		}

		biggestN = 0;
		
		if(array1pts > x->x_arrayPoints)
			biggestN = array1pts;
		else
			biggestN = x->x_arrayPoints;

		vecBuffer = (t_float *)t_getbytes(biggestN*sizeof(t_float));
		vec1Buffer = (t_float *)t_getbytes(biggestN*sizeof(t_float));
		vecWeights = (t_float *)t_getbytes(biggestN*sizeof(t_float));

		for(i=0; i<biggestN; i++)
		{
			vecBuffer[i] = 0.0;
			vec1Buffer[i] = 0.0;
			// these weights aren't used in tabletool really. Just need to set all weights to 1.0 in order to use tIDLib_ distance functions
			vecWeights[i] = 1.0;
		}

		for(i=0; i<biggestN; i++)
		{
			if(i<x->x_arrayPoints)
				vecBuffer[i] = x->x_vec[i].w_float;
			
			if(i<array1pts)
				vec1Buffer[i] = vec1[i].w_float;
		}
			
		dist = tIDLib_euclidDist(biggestN, vecBuffer, vec1Buffer, vecWeights, true);

		// free local memory
		t_freebytes(vecBuffer, biggestN * sizeof(t_float));
		t_freebytes(vec1Buffer, biggestN * sizeof(t_float));
		t_freebytes(vecWeights, biggestN * sizeof(t_float));

		outlet_float(x->x_info, dist);
	}
}


static void tabletool_taxi(t_tabletool *x, t_symbol *array1)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, array1pts, biggestN;
		t_garray *b;
		t_word *vec1;
		t_float dist, *vecBuffer, *vec1Buffer, *vecWeights;

		if(!(b = (t_garray *)pd_findbyclass(array1, garray_class)))
		{
			pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, array1->s_name);
			return;
		}
		else if(!garray_getfloatwords(b, (int *)&array1pts, &vec1))
		{
			pd_error(x, "%s: bad template for %s", array1->s_name, x->x_objSymbol->s_name);
			return;
		}

		biggestN = 0;
		
		if(array1pts > x->x_arrayPoints)
			biggestN = array1pts;
		else
			biggestN = x->x_arrayPoints;

		vecBuffer = (t_float *)t_getbytes(biggestN*sizeof(t_float));
		vec1Buffer = (t_float *)t_getbytes(biggestN*sizeof(t_float));
		vecWeights = (t_float *)t_getbytes(biggestN*sizeof(t_float));

		for(i=0; i<biggestN; i++)
		{
			vecBuffer[i] = 0.0;
			vec1Buffer[i] = 0.0;
			// these weights aren't used in tabletool really. Just need to set all weights to 1.0 in order to use tIDLib_ distance functions
			vecWeights[i] = 1.0;
		}

		for(i=0; i<biggestN; i++)
		{
			if(i<x->x_arrayPoints)
				vecBuffer[i] = x->x_vec[i].w_float;
			
			if(i<array1pts)
				vec1Buffer[i] = vec1[i].w_float;
		}
			
		dist = tIDLib_taxiDist(biggestN, vecBuffer, vec1Buffer, vecWeights);

		// free local memory
		t_freebytes(vecBuffer, biggestN * sizeof(t_float));
		t_freebytes(vec1Buffer, biggestN * sizeof(t_float));
		t_freebytes(vecWeights, biggestN * sizeof(t_float));

		outlet_float(x->x_info, dist);
	}
}


static void tabletool_corr(t_tabletool *x, t_symbol *array1)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{ 
		t_sampIdx i, array1pts, biggestN;
		t_garray *b;
		t_word *vec1;
		t_float corr, *vecBuffer, *vec1Buffer;
	
		if(!(b = (t_garray *)pd_findbyclass(array1, garray_class)))
		{
			pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, array1->s_name);
			return;
		}
		else if(!garray_getfloatwords(b, (int *)&array1pts, &vec1))
		{
			pd_error(x, "%s: bad template for %s", array1->s_name, x->x_objSymbol->s_name);
			return;
		}

		biggestN = 0;
		
		if(array1pts > x->x_arrayPoints)
			biggestN = array1pts;
		else
			biggestN = x->x_arrayPoints;

		vecBuffer = (t_float *)t_getbytes(biggestN*sizeof(t_float));
		vec1Buffer = (t_float *)t_getbytes(biggestN*sizeof(t_float));

		for(i=0; i<biggestN; i++)
		{
			vecBuffer[i] = 0.0;
			vec1Buffer[i] = 0.0;
		}

		for(i=0; i<biggestN; i++)
		{
			if(i<x->x_arrayPoints)
				vecBuffer[i] = x->x_vec[i].w_float;
			
			if(i<array1pts)
				vec1Buffer[i] = vec1[i].w_float;
		}
			
		corr = tIDLib_corr(biggestN, vecBuffer, vec1Buffer);

		// free local memory
		t_freebytes(vecBuffer, biggestN * sizeof(t_float));
		t_freebytes(vec1Buffer, biggestN * sizeof(t_float));

		outlet_float(x->x_info, corr);
	}
}


static void tabletool_abs(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;

		for(i=0; i<x->x_arrayPoints; i++)
			if(x->x_vec[i].w_float < 0.0)
				x->x_vec[i].w_float *= -1.0;

		garray_redraw(a);
 	}
}


static void tabletool_reciprocal(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;

		for(i=0; i<x->x_arrayPoints; i++)
		{
			t_float thisVal;
			
			thisVal = x->x_vec[i].w_float;
			
			if(thisVal!=0.0)
				x->x_vec[i].w_float = 1.0/thisVal;
			else
				pd_error(x, "%s: divide by zero detected", x->x_objSymbol->s_name);
		}

		garray_redraw(a);
 	}
}


static void tabletool_sum(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		t_float sum;

		sum = 0.0;

		for(i=0; i<x->x_arrayPoints; i++)
			sum += x->x_vec[i].w_float;

		outlet_float(x->x_info, sum);
 	}
}


static void tabletool_mean(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		t_float mean;

		mean = 0.0;

		for(i=0; i<x->x_arrayPoints; i++)
			mean += x->x_vec[i].w_float;

		mean /= x->x_arrayPoints;

		outlet_float(x->x_info, mean);
 	}
}


static void tabletool_geomean(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		double numPointsRecip, mean, *nthRoots;

		nthRoots = (double *)t_getbytes(x->x_arrayPoints*sizeof(double));

		mean = 1.0;
		numPointsRecip = 1.0/x->x_arrayPoints;

		// take the nth roots first so as not to lose data.
		for(i=0; i<x->x_arrayPoints; i++)
			nthRoots[i] = pow(x->x_vec[i].w_float, numPointsRecip);

		// take the product of nth roots
		for(i=0; i<x->x_arrayPoints; i++)
			mean *= nthRoots[i];

		outlet_float(x->x_info, mean);

		// free local memory
		t_freebytes(nthRoots, x->x_arrayPoints*sizeof(double));
 	}
}


static void tabletool_std(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		t_float *tableVals, sum, mean, std;

		tableVals = (t_float *)t_getbytes(x->x_arrayPoints*sizeof(t_float));

		for(i=0; i<x->x_arrayPoints; i++)
			tableVals[i] = x->x_vec[i].w_float;

		sum = 0.0;

		for(i=0; i<x->x_arrayPoints; i++)
			sum += tableVals[i];

		mean = sum/x->x_arrayPoints;

		sum = 0.0;

		// center & square the data
		for(i=0; i<x->x_arrayPoints; i++)
		{
			tableVals[i] -= mean;
			tableVals[i] *= tableVals[i];
			sum += tableVals[i];
		}

		std = sum/(x->x_arrayPoints-1.0);
		std = sqrt(std);

		outlet_float(x->x_info, std);

		// free local memory
		t_freebytes(tableVals, x->x_arrayPoints * sizeof(t_float));
 	}
}


static void tabletool_bestFitLine(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
/*
		t_sampIdx i;
		t_float slope, *dataBuf;
		
		dataBuf = (t_float *)t_getbytes(x->x_arrayPoints*sizeof(t_float));

		for(i=0; i<x->x_arrayPoints; i++)
			dataBuf[i] = x->x_vec[i].w_float;
			
		slope = tIDLib_fitLineSlope(x->x_arrayPoints, dataBuf);

		t_freebytes(dataBuf, x->x_arrayPoints*sizeof(t_float));

		outlet_float(x->x_info, slope);
*/
		pd_error(x, "%s: method not implemented", x->x_objSymbol->s_name);
 	}
}


static void tabletool_normalize(t_tabletool *x, t_float min, t_float max)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		t_float range_data, rangeFinal, smallest, largest;

		smallest = FLT_MAX;

		for(i=0; i<x->x_arrayPoints; i++)
			if(x->x_vec[i].w_float < smallest)
				smallest = x->x_vec[i].w_float;

		largest = smallest;

		for(i=0; i<x->x_arrayPoints; i++)
			if(x->x_vec[i].w_float > largest)
				largest = x->x_vec[i].w_float;

		rangeFinal = max - min;
		range_data = largest - smallest;

		if(rangeFinal <= 0)
			rangeFinal = 1;

		if(range_data <= 0)
			range_data = 1;

		range_data = 1.0/range_data;

		// shift everything >= 0
		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float -= smallest;

		// scale to 0-1 range
		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float *= range_data;

		// scale to requested range
		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float *= rangeFinal;

		// offset downward according to min
		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float += min;

		garray_redraw(a);
 	}
}


static void tabletool_normalize_sum(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		t_float sum, smallest;

		smallest = FLT_MAX;

		for(i=0; i<x->x_arrayPoints; i++)
			if(x->x_vec[i].w_float < smallest)
				smallest = x->x_vec[i].w_float;

		// shift everything >= 0
		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float -= smallest;

		sum = 0.0;

		for(i=0; i<x->x_arrayPoints; i++)
			sum += x->x_vec[i].w_float;

		sum = 1.0/sum;

		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float *= sum;

		garray_redraw(a);
 	}
}


static void tabletool_set(t_tabletool *x, t_symbol *s)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(s, garray_class)))
		pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, s->s_name);
	else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
		pd_error(x, "%s: bad template for %s", s->s_name, x->x_objSymbol->s_name);
	else
	{
	    x->x_arrayName = s;
	}
}


static void tabletool_fitBounds(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
		pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
	else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
		pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_glist *thisGlist;
		t_sampIdx i, x1, x2;
		t_float max, min, y1, y2;

		thisGlist = garray_getglist(a);

		max = -FLT_MAX;
		// find max
		for(i=0; i<x->x_arrayPoints; i++)
			if(x->x_vec[i].w_float > max)
				max = x->x_vec[i].w_float;

		min = FLT_MAX;
		// find max
		for(i=0; i<x->x_arrayPoints; i++)
			if(x->x_vec[i].w_float < min)
				min = x->x_vec[i].w_float;

		x1 = 0;
		x2 = x->x_arrayPoints-1;
		y1 = max;
		y2 = min;

		thisGlist->gl_x1 = x1;
		thisGlist->gl_x2 = x2;
		thisGlist->gl_y1 = y1;
		thisGlist->gl_y2 = y2;

		if (thisGlist->gl_x2 == thisGlist->gl_x1 ||
			thisGlist->gl_y2 == thisGlist->gl_y1)
		{
			pd_error(x, "%s: empty bounds rectangle", x->x_objSymbol->s_name);
			x1 = y1 = 0;
			x2 = y2 = 1;
		}

		glist_redraw(thisGlist);
		garray_redraw(a);
	}
}


static void tabletool_copy(t_tabletool *x, t_symbol *source)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, sourcePts;
		t_garray *b;
		t_word *sourceVec;
	
		if(!(b = (t_garray *)pd_findbyclass(source, garray_class)))
		{
			pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, source->s_name);
			return;
		}
		else if(!garray_getfloatwords(b, (int *)&sourcePts, &sourceVec))
		{
			pd_error(x, "%s: bad template for %s", source->s_name, x->x_objSymbol->s_name);
			return;
		}

		for(i=0; i<sourcePts; i++)
			if(i > x->x_arrayPoints-1)
				break;
			else
				x->x_vec[i].w_float = sourceVec[i].w_float;

		for(; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float = 0.0;

		garray_redraw(a);
	}
}


static void tabletool_copy_range(t_tabletool *x, t_symbol *s, int argc, t_atom *argv)
{
	t_garray *a;

	if(argc<4)
	{
		pd_error(x, "%s: too few arguments for %s method", x->x_objSymbol->s_name, s->s_name);
		return;
	}

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, j, sourcePts, targetStart, sourceStart, sourceEnd;
		t_float targetStartFloat, sourceStartFloat;
		t_garray *b;
		t_symbol *source;
		t_word *sourceVec;

		targetStartFloat = atom_getfloat(argv+0);
		source = atom_getsymbol(argv+1);
		sourceStartFloat = atom_getfloat(argv+2);
		sourceEnd = atom_getfloat(argv+3);

		if(!(b = (t_garray *)pd_findbyclass(source, garray_class)))
		{
			pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, source->s_name);
			return;
		}
		else if(!garray_getfloatwords(b, (int *)&sourcePts, &sourceVec))
		{
			pd_error(x, "%s: bad template for %s", source->s_name, x->x_objSymbol->s_name);
			return;
		}

		targetStart = (targetStartFloat<0)?0:targetStartFloat;
		sourceStart = (sourceStartFloat<0)?0:sourceStartFloat;
		sourceEnd = (sourceEnd>sourcePts-1)?sourcePts-1:sourceEnd;

		if(sourceStart>sourceEnd)
		{
			t_sampIdx tmp;
			tmp = sourceEnd;
			sourceEnd = sourceStart;
			sourceStart = tmp;
		}
		
		for(i=targetStart, j=sourceStart; j<=sourceEnd; i++, j++)
			if(i>x->x_arrayPoints-1)
				break;
			else
				x->x_vec[i].w_float = sourceVec[j].w_float;

		garray_redraw(a);
	}
}


static void tabletool_blackman(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		t_float *wPtr;
		
		wPtr = (t_float *)t_getbytes(x->x_arrayPoints*sizeof(t_float));

		tIDLib_blackmanWindow(wPtr, x->x_arrayPoints);
		
		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float = *(wPtr+i);

		garray_redraw(a);
		
		t_freebytes(wPtr, x->x_arrayPoints*sizeof(t_float));
	}
}


static void tabletool_cosine(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		t_float *wPtr;
		
		wPtr = (t_float *)t_getbytes(x->x_arrayPoints*sizeof(t_float));

		tIDLib_cosineWindow(wPtr, x->x_arrayPoints);
		
		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float = *(wPtr+i);

		garray_redraw(a);
		
		t_freebytes(wPtr, x->x_arrayPoints*sizeof(t_float));
	}
}


static void tabletool_hamming(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		t_float *wPtr;
		
		wPtr = (t_float *)t_getbytes(x->x_arrayPoints*sizeof(t_float));

		tIDLib_hammingWindow(wPtr, x->x_arrayPoints);
		
		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float = *(wPtr+i);

		garray_redraw(a);
		
		t_freebytes(wPtr, x->x_arrayPoints*sizeof(t_float));
	}
}


static void tabletool_hann(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		t_float *wPtr;
		
		wPtr = (t_float *)t_getbytes(x->x_arrayPoints*sizeof(t_float));

		tIDLib_hannWindow(wPtr, x->x_arrayPoints);
		
		for(i=0; i<x->x_arrayPoints; i++)
			x->x_vec[i].w_float = *(wPtr+i);

		garray_redraw(a);
		
		t_freebytes(wPtr, x->x_arrayPoints*sizeof(t_float));
	}
}


static void tabletool_randFill(t_tabletool *x, t_float min, t_float max)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		t_float range;

		range = max - min;

		for(i=0; i<x->x_arrayPoints; i++)
		{
			x->x_vec[i].w_float = rand();
			x->x_vec[i].w_float /= RAND_MAX;
			x->x_vec[i].w_float *= range;
			x->x_vec[i].w_float += min;
		}

		garray_redraw(a);
 	}
}


static void tabletool_peaks(t_tabletool *x, t_float threshPct)
{
	t_garray *a;

	threshPct = (threshPct<0.0)?0.0:threshPct;
	threshPct = (threshPct>100.0)?100.0:threshPct;
	
	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		t_float *flagsBuf, *dataBuf, maxPeakVal, minPeakVal, minVal, maxPeakRange, thresh;
		flagsBuf = (t_float *)t_getbytes(x->x_arrayPoints*sizeof(t_float));
		dataBuf = (t_float *)t_getbytes(x->x_arrayPoints*sizeof(t_float));


		minVal = FLT_MAX;
		maxPeakVal = -FLT_MAX;
		minPeakVal = FLT_MAX;

		for(i=0; i<x->x_arrayPoints; i++)
			if(x->x_vec[i].w_float<minVal)
				minVal = x->x_vec[i].w_float;

		for(i=0; i<x->x_arrayPoints; i++)
			dataBuf[i] = x->x_vec[i].w_float;
			
		tIDLib_peaksValleys(x->x_arrayPoints, dataBuf, flagsBuf, &minPeakVal, &maxPeakVal);
		
		maxPeakRange = maxPeakVal - minVal;
	
		thresh = maxPeakRange * (threshPct/100.0);

		for(i=0; i<x->x_arrayPoints; i++)
		{
			// 0.5 in the flagsBuf means a half peak, which we'll ignore
			if(flagsBuf[i]>0.5)
			{
				if(dataBuf[i]-minVal>=thresh)
				{
					t_atom indexOut;

					SETFLOAT(&indexOut, i);
					outlet_list(x->x_list, 0, 1, &indexOut);
					outlet_float(x->x_info, x->x_vec[i].w_float);
				}
			}
		}
		
		// free local memory
		t_freebytes(flagsBuf, x->x_arrayPoints*sizeof(t_float));
		t_freebytes(dataBuf, x->x_arrayPoints*sizeof(t_float));
	}	
}


static void tabletool_valleys(t_tabletool *x, t_float threshPct)
{
	t_garray *a;

	threshPct = (threshPct<0.0)?0.0:threshPct;
	threshPct = (threshPct>100.0)?100.0:threshPct;
	
	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;
		t_float *flagsBuf, *dataBuf, maxPeakVal, minPeakVal, maxVal, minVal, maxPeakRange, thresh;
		flagsBuf = (t_float *)t_getbytes(x->x_arrayPoints*sizeof(t_float));
		dataBuf = (t_float *)t_getbytes(x->x_arrayPoints*sizeof(t_float));

		maxVal = -FLT_MAX;
		minVal = FLT_MAX;
		maxPeakVal = -FLT_MAX;
		minPeakVal = FLT_MAX;

		for(i=0; i<x->x_arrayPoints; i++)
			if(x->x_vec[i].w_float>maxVal)
				maxVal = x->x_vec[i].w_float;
				
		for(i=0; i<x->x_arrayPoints; i++)
			if(x->x_vec[i].w_float<minVal)
				minVal = x->x_vec[i].w_float;

		// invert the data so that valleys become peaks
		for(i=0; i<x->x_arrayPoints; i++)
			dataBuf[i] = maxVal - x->x_vec[i].w_float + minVal;

		// get new minVal after inversion and do the same as in the peaks function
		minVal = FLT_MAX;
		for(i=0; i<x->x_arrayPoints; i++)
			if(x->x_vec[i].w_float<minVal)
				minVal = x->x_vec[i].w_float;

		tIDLib_peaksValleys(x->x_arrayPoints, dataBuf, flagsBuf, &minPeakVal, &maxPeakVal);
		
		maxPeakRange = maxPeakVal - minVal;
	
		thresh = maxPeakRange * (threshPct/100.0);

		for(i=0; i<x->x_arrayPoints; i++)
		{
			// 0.5 in the flagsBuf means a half peak, which we'll ignore
			if(flagsBuf[i]>0.5)
			{
				if(dataBuf[i]-minVal>=thresh)
				{
					t_atom indexOut;

					SETFLOAT(&indexOut, i);
					outlet_list(x->x_list, 0, 1, &indexOut);
					outlet_float(x->x_info, x->x_vec[i].w_float);
				}
			}
		}
		
		// free local memory
		t_freebytes(flagsBuf, x->x_arrayPoints*sizeof(t_float));
		t_freebytes(dataBuf, x->x_arrayPoints*sizeof(t_float));
	}	
}


static void tabletool_peaksThresh(t_tabletool *x, t_float min, t_float max)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i, peakMaxIdx;
		t_bool peakFlag;
		t_float peakMax, tmp;

		if(min > max)
		{
			tmp = min;
			min = max;
			max = tmp;
		}

		peakFlag = false;
		peakMax = 0;
		peakMaxIdx = 0;

		for(i=0; i<x->x_arrayPoints-1; i++)
		{
			if(x->x_vec[i].w_float > max)
				peakFlag = true;

			if(peakFlag && (x->x_vec[i].w_float > peakMax))
			{
				peakMax = x->x_vec[i].w_float;
				peakMaxIdx = i;
			}
			else if(peakFlag && (x->x_vec[i].w_float < min))
			{
				t_atom indexOut;

				SETFLOAT(&indexOut, peakMaxIdx);
				outlet_list(x->x_list, 0, 1, &indexOut);
				outlet_float(x->x_info, peakMax);

				peakFlag = false;
				peakMax = 0;
				peakMaxIdx = 0;
			}
		}
	}
}


static void tabletool_store(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		t_sampIdx i;

		x->x_originalData = (t_float *)t_getbytes(x->x_arrayPoints*sizeof(t_float));

		// load table data
		for(i=0; i<x->x_arrayPoints; i++)
			x->x_originalData[i] = x->x_vec[i].w_float;

		x->x_storedFlag = true;
	}
}


static void tabletool_restore(t_tabletool *x)
{
	t_garray *a;

	if(x->x_storedFlag)
	{
		if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
			pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
		else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
			pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
		else
		{
			t_sampIdx i;

			// load original data
			for(i=0; i<x->x_arrayPoints; i++)
				x->x_vec[i].w_float = x->x_originalData[i];

			garray_redraw(a);
		}
	}
	else
		pd_error(x, "%s: table not previously stored", x->x_objSymbol->s_name);
}


static void tabletool_wipe(t_tabletool *x)
{
	if(x->x_storedFlag)
	{
		t_freebytes(x->x_originalData, x->x_arrayPoints*sizeof(t_float));
		x->x_storedFlag = false;
	}
}


static void tabletool_change(t_tabletool *x)
{
	t_garray *a;

	if(!(a = (t_garray *)pd_findbyclass(x->x_arrayName, garray_class)))
        pd_error(x, "%s: no array named %s", x->x_objSymbol->s_name, x->x_arrayName->s_name);
    else if(!garray_getfloatwords(a, (int *)&x->x_arrayPoints, &x->x_vec))
    	pd_error(x, "%s: bad template for %s", x->x_arrayName->s_name, x->x_objSymbol->s_name);
	else
	{
		if(x->x_storedFlag)
		{
			t_sampIdx i;
	
			for(i=0; i<x->x_arrayPoints; i++)
			{
				if(x->x_vec[i].w_float != x->x_originalData[i])
				{
					t_atom indexOut;
					
					SETFLOAT(&indexOut, i);
					outlet_list(x->x_list, 0, 1, &indexOut);
					outlet_float(x->x_info, x->x_vec[i].w_float);
					x->x_originalData[i] = x->x_vec[i].w_float;
				}
			}
		}
		else
			tabletool_store(x);
 	}
}


static void *tabletool_new(t_symbol *s)
{
    t_tabletool *x = (t_tabletool *)pd_new(tabletool_class);

	x->x_info = outlet_new(&x->x_obj, &s_float);
	x->x_list = outlet_new(&x->x_obj, gensym("list"));

	x->x_objSymbol = gensym("tabletool");

	x->x_arrayName = s;

	x->x_storedFlag = false;
	x->x_randState = (t_uInt)clock_getlogicaltime(); // seed with (t_uInt) logical time

    return (x);
}


static void tabletool_free(t_tabletool *x)
{
	if(x->x_storedFlag)
	    t_freebytes(x->x_originalData, x->x_arrayPoints*sizeof(t_float));
}


void tabletool_setup(void)
{
    tabletool_class =
    class_new(
    	gensym("tabletool"),
    	(t_newmethod)tabletool_new,
    	(t_method)tabletool_free,
        sizeof(t_tabletool),
        CLASS_DEFAULT,
        A_DEFSYMBOL,
		0
    );

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_dump,
		gensym("dump"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_dumpRange,
		gensym("dump_range"),
		A_DEFFLOAT,
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		tabletool_class,
		(t_method)tabletool_size,
		gensym("size"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_range,
		gensym("range"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_min,
		gensym("min"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_max,
		gensym("max"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_maxMag,
		gensym("max_mag"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_nearest,
		gensym("nearest"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_equals,
		gensym("equals"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_greater,
		gensym("greater"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_less,
		gensym("less"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_between,
		gensym("between"),
		A_DEFFLOAT,
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_findZeroCrossings,
		gensym("find_zero_crossings"),
		0
	);
	
	class_addmethod(
		tabletool_class,
		(t_method)tabletool_choose,
		gensym("choose"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_const,
		gensym("const"),
		A_DEFFLOAT,
		A_DEFFLOAT,
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_series,
		gensym("series"),
		A_DEFFLOAT,
		A_DEFFLOAT,
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_randomWalk,
		gensym("random_walk"),
		A_DEFFLOAT,
		A_DEFFLOAT,
		A_DEFFLOAT,
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_offset,
		gensym("offset"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_shift,
		gensym("shift"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_shift0,
		gensym("shift0"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_invert,
		gensym("invert"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_reverse,
		gensym("reverse"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_remove,
		gensym("remove"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_insert,
		gensym("insert"),
		A_DEFFLOAT,
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_round,
		gensym("round"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		tabletool_class,
		(t_method)tabletool_floor,
		gensym("floor"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_ceil,
		gensym("ceil"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_mtof,
		gensym("mtof"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_ftom,
		gensym("ftom"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_dbtorms,
		gensym("dbtorms"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_rmstodb,
		gensym("rmstodb"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_bin2freq,
		gensym("bin2freq"),
		A_DEFFLOAT,
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_freq2bin,
		gensym("freq2bin"),
		A_DEFFLOAT,
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_bark2freq,
		gensym("bark2freq"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_freq2bark,
		gensym("freq2bark"),
		0
	);
	
	class_addmethod(
		tabletool_class,
		(t_method)tabletool_clip,
		gensym("clip"),
		A_DEFFLOAT,
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		tabletool_class,
		(t_method)tabletool_smooth,
		gensym("smooth"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_lace,
		gensym("lace"),
		A_DEFSYMBOL,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_concatenate,
		gensym("concatenate"),
		A_DEFSYMBOL,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_scramble,
		gensym("scramble"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_sort,
		gensym("sort"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_sort_range,
		gensym("sort_range"),
		A_DEFFLOAT,
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_swap,
		gensym("swap"),
		A_DEFFLOAT,
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_replace,
		gensym("replace"),
		A_DEFFLOAT,
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		tabletool_class,
		(t_method)tabletool_integrate,
		gensym("integrate"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_differentiate,
		gensym("differentiate"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_add,
		gensym("add"),
		A_DEFSYMBOL,
		0
	);

/*
	class_addmethod(
		tabletool_class,
		(t_method)tabletool_add_range,
		gensym("add_range"),
		A_GIMME,
		0
	);
*/

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_subtract,
		gensym("subtract"),
		A_DEFSYMBOL,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_multiply,
		gensym("multiply"),
		A_DEFSYMBOL,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_divide,
		gensym("divide"),
		A_DEFSYMBOL,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_dot,
		gensym("dot"),
		A_DEFSYMBOL,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_euclid,
		gensym("euclid"),
		A_DEFSYMBOL,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_taxi,
		gensym("taxi"),
		A_DEFSYMBOL,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_corr,
		gensym("corr"),
		A_DEFSYMBOL,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_abs,
		gensym("abs"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_reciprocal,
		gensym("reciprocal"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_sum,
		gensym("sum"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_mean,
		gensym("mean"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_geomean,
		gensym("geomean"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_std,
		gensym("std"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_bestFitLine,
		gensym("best_fit_line"),
		0
	);
	
	class_addmethod(
		tabletool_class,
		(t_method)tabletool_scale,
		gensym("scale"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_normalize,
		gensym("normalize"),
		A_DEFFLOAT,
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_normalize_sum,
		gensym("normalize_sum"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_set,
		gensym("set"),
		A_SYMBOL,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_fitBounds,
		gensym("fit_bounds"),
		0
	);
	
	class_addmethod(
		tabletool_class,
		(t_method)tabletool_copy,
		gensym("copy"),
		A_DEFSYMBOL,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_copy_range,
		gensym("copy_range"),
		A_GIMME,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_blackman,
		gensym("blackman"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_cosine,
		gensym("cosine"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_hamming,
		gensym("hamming"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_hann,
		gensym("hann"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_randFill,
		gensym("randFill"),
		A_DEFFLOAT,
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_peaks,
		gensym("peaks"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_valleys,
		gensym("valleys"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_peaksThresh,
		gensym("peaks_thresh"),
		A_DEFFLOAT,
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_store,
		gensym("store"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_restore,
		gensym("restore"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_wipe,
		gensym("wipe"),
		0
	);

	class_addmethod(
		tabletool_class,
		(t_method)tabletool_change,
		gensym("change"),
		0
	);
}