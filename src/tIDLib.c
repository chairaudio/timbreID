/*

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "tIDLib.h"

/* ---------------- conversion functions ---------------------- */
t_binIdx tIDLib_freq2bin(t_float freq, t_float n, t_float sr)
{
	t_binIdx bin;

	bin = roundf(freq*n/sr);

	return(bin);
}

t_float tIDLib_bin2freq(t_binIdx bin, t_float n, t_float sr)
{
	t_float freq;

	freq = bin*sr/n;

	return(freq);
}

// can have a flag for the other way to calculate it: 13*arctan(0.00076*freq) + 3.5*arctan((freq/7500)^2)
t_float tIDLib_freq2bark(t_float freq)
{
	t_float barkFreq;
	t_freq2barkFormula formula; // TODO: this should be an enum
	
	formula = 0;

	switch(formula)
	{
		case freq2barkFormula0:
			barkFreq = 6.0*asinh(freq/600.0);
			break;
		case freq2barkFormula1:
			barkFreq = 13*atan(0.00076*freq) + 3.5*atan(powf((freq/7500), 2));
			break;
		case freq2barkFormula2:
			barkFreq = ((26.81*freq)/(1960+freq))-0.53;
			if(barkFreq<2)
				barkFreq += 0.15*(2-barkFreq);
			else if(barkFreq>20.1)
				barkFreq += 0.22*(barkFreq-20.1);
			break;
		default:
			barkFreq = 0;
			break;
	}
	
	barkFreq = (barkFreq<0)?0:barkFreq;
	
	return(barkFreq);
}

t_float tIDLib_bark2freq(t_float bark)
{
	t_float freq;
	t_bark2freqFormula formula; // TODO: this should be an enum
	
	// formula 2 was used in timbreID 0.6
	formula = 0;

	switch(formula)
	{
		case bark2freqFormula0:
			freq = 600.0 * sinh(bark/6.0);
			break;
		case bark2freqFormula1:
			freq = 53548.0/(bark*bark - 52.56*bark + 690.39);
			break;
		case bark2freqFormula2:
			freq = 1960/(26.81/(bark+0.53) - 1);
			freq = (freq<0)?0:freq;
			break;
		default:
			freq = 0;
	}
	
	return(freq);
}

t_float tIDLib_freq2mel(t_float freq)
{
	t_float mel;

	mel = 1127*log(1+(freq/700));
	mel = (mel<0)?0:mel;
	return(mel);
}

t_float tIDLib_mel2freq(t_float mel)
{
	t_float freq;

	freq = 700 * (exp(mel/1127) - 1);
//	freq = 700 * (exp(mel/1127.01048) - 1);
	freq = (freq<0)?0:freq;
	return(freq);
}
/* ---------------- END conversion functions ---------------------- */




/* ---------------- utility functions ---------------------- */
void tIDLib_linspace(t_float *ramp, t_float start, t_float finish, t_binIdx n)
{
	t_float diffInc;
	t_binIdx i;
	
	diffInc = (finish-start)/(t_float)(n-1);
	
	ramp[0] = start;
	for(i=1; i<n; i++)
		ramp[i] = ramp[i-1]+diffInc;
}


t_sChar tIDLib_signum(t_float input)
{
	t_sChar sign;
	
	sign=0;
	
	if(input>0)
		sign = 1;
	else if(input<0)
		sign = -1;
	else
		sign = 0;
	
	return(sign);
}


t_float tIDLib_fitLineSlope(t_sampIdx n, t_float *input)
{
	t_float x[n], y[n], diffsX[n], diffsY[n];	
	t_float sumX, sumY, meanX, meanY;
	t_float dividend, divisor, slope;
	t_sampIdx i;
	
	sumX = 0;
	sumY = 0;
	
	for(i=0; i<n; i++, input++)
	{
		x[i] = i;
		y[i] = fabs(*input);
		
		sumX += x[i];
		sumY += y[i];
	}
	
	meanX = sumX/n;
	meanY = sumY/n;
	dividend = 0;
	divisor = 0;
	
	for(i=0; i<n; i++, input++)
	{
		diffsX[i] = x[i] - meanX;
		diffsY[i] = y[i] - meanY;
		dividend += diffsX[i] * diffsY[i];
		divisor += diffsX[i] * diffsX[i];
	}
	
	// scale slope according to window size. This would give a perfect 0-1 ramp over N samples a slope of 1.0
	slope = (dividend/divisor) * n;
	return(slope);
}


void tIDLib_bubbleSort(t_sampIdx n, t_float *list)
{	
	t_sampIdx i;
	
	for(i=0; i<n; i++)
	{
		t_bool flag = false;
		t_sampIdx j;
		
		for(j=0; j<(n-1); j++)
		{
			if(list[j] > list[j+1])
			{
				t_float tmp;

				flag = true;

				tmp = list[j+1];
				list[j+1] = list[j];
				list[j] = tmp;
			}
		}

		if(!flag)
			break;
	}
}


void tIDLib_knnInfoBubbleSort(t_uShortInt n, t_instance *instances)
{
	t_uShortInt i;

	for(i=0; i<n; i++)
	{
		t_bool flag = false;
		t_uShortInt j;
		
		for(j=0; j<(n-1); j++)
		{
			if(instances[j].knnInfo.safeDist > instances[j+1].knnInfo.safeDist)
			{
				t_knnInfo tmp;

				flag = true;

				tmp = instances[j+1].knnInfo;
				instances[j+1].knnInfo = instances[j].knnInfo;
				instances[j].knnInfo = tmp;
			}
		}

		if(!flag)
			break;
	}
}


void tIDLib_sortKnnInfo(t_uShortInt k, t_instanceIdx numInstances, t_instanceIdx prevMatch, t_instance *instances)
{
	t_instanceIdx *topMatches;
	t_uShortInt i;

	topMatches = (t_instanceIdx *)t_getbytes(k * sizeof(t_instanceIdx));

	for(i=0; i<k; i++)
	{
		t_float maxBest;
		t_instanceIdx j, topIdx;
		
		maxBest = FLT_MAX;
		topIdx = 0;
			
		for(j=0; j<numInstances; j++)
			if(instances[j].knnInfo.dist < maxBest)
				if(instances[j].knnInfo.idx != prevMatch) // doesn't include previous match - this is good
				{
					maxBest = instances[j].knnInfo.dist;
					topIdx = j;
				};
		
		instances[topIdx].knnInfo.dist = FLT_MAX;
	
		topMatches[i] = instances[topIdx].knnInfo.idx;
	}

	for(i=0; i<k; i++)
	{
		t_knnInfo tmp;

		tmp = instances[i].knnInfo;
		instances[i].knnInfo = instances[topMatches[i]].knnInfo;
		instances[topMatches[i]].knnInfo = tmp;
	}

	// sort the top K matches, because they have the lowest distances in the whole list, but they're not sorted
	tIDLib_knnInfoBubbleSort(k, instances);

	// free memory
	t_freebytes(topMatches, k*sizeof(t_instanceIdx));

	// now, the list passed to the function will have the first k elements in order, 
	// and these elements have the lowest distances in the whole list.
}


t_float tIDLib_dotProd(t_attributeIdx n, t_float *v1, t_float *v2)
{
	t_float dot;
	t_attributeIdx i;

	dot=0.0;

	for(i=0; i<n; i++)
		dot += v1[i]*v2[i];

	return(dot);
}


t_float tIDLib_euclidDist(t_attributeIdx n, t_float *v1, t_float *v2, t_float *weights, t_bool sqroot)
{
	t_float diff, dist;
	t_attributeIdx i;

	diff=dist=0.0;

	for(i=0; i<n; i++)
	{
		diff = v1[i] - v2[i];
		dist += diff*diff*weights[i];
	}

	// if the square root flag is on, take the square root
	if(sqroot)
		dist = sqrt(dist);
	
	return(dist);
}


t_float tIDLib_taxiDist(t_attributeIdx n, t_float *v1, t_float *v2, t_float *weights)
{
	t_float diff, dist;
	t_attributeIdx i;

	diff=dist=0.0;

	for(i=0; i<n; i++)
	{
		diff = v1[i] - v2[i];
		dist += fabs(diff)*weights[i];
	}
	
	return(dist);
}


t_float tIDLib_corr(t_attributeIdx n, t_float *v1, t_float *v2)
{
	t_float corr, sum1, mean1, std1, sum2, mean2, std2, *v1centered, *v2centered;
	t_attributeIdx i;

	v1centered = (t_float *)t_getbytes(n*sizeof(t_float));
	v2centered = (t_float *)t_getbytes(n*sizeof(t_float));
	
	sum1 = sum2 = mean1 = mean2 = std1 = std2 = corr = 0.0;

	for(i=0; i<n; i++)
	{
		sum1 += v1[i];
		sum2 += v2[i];
	}

	mean1 = sum1/n;
	mean2 = sum2/n;

	for(i=0; i<n; i++)
	{
		v1centered[i] = v1[i]-mean1;
		v2centered[i] = v2[i]-mean2;
	}

	for(i=0; i<n; i++)
	{
		std1 += v1centered[i]*v1centered[i];
		std2 += v2centered[i]*v2centered[i];
	}

	std1 = sqrt(std1/n);
	std2 = sqrt(std2/n);

	for(i=0; i<n; i++)
		corr += v1centered[i]*v2centered[i];

	corr /= n;
	corr = corr/(std1*std2);


	// free local memory
	t_freebytes(v1centered, n * sizeof(t_float));
	t_freebytes(v2centered, n * sizeof(t_float));

	return(corr);
}


void tIDLib_peaksValleys(t_sampIdx n, t_float *data, t_float *flags, t_float *minVal, t_float *maxVal)
{
	t_float slopeBuf[n];
	t_float localMaxVal, localMinVal;
	t_sampIdx i;

	for(i=0; i<n; i++)
	{
		t_float slope;
		
		if(i>0)
			slope = data[i] - data[i-1];
		else
			slope = data[i] - 0;
	
		slopeBuf[i] = slope;
	}

	localMaxVal = -FLT_MAX;
	localMinVal = FLT_MAX;
	
	flags[0] = tIDLib_signum(slopeBuf[0]);
	for(i=1; i<n; i++)
	{
		t_sChar change;

		change = tIDLib_signum(slopeBuf[i]) - tIDLib_signum(slopeBuf[i-1]);

		switch(change)
		{
			// peak case
			case -2:
				// write to i-1, because the peak moment occurred just before this observed change in slope direction
				flags[i-1] = 1;
				if(data[i-1]>localMaxVal)
					localMaxVal = data[i-1];
				break;
			// plateau case
			case 0:
				flags[i-1] = 0;
				break;
			// half-peak case
			case 1:
				flags[i-1] = 0.5;
				if(data[i-1]>localMaxVal)
					localMaxVal = data[i-1];
				break;
			// valley case
			case 2:
				flags[i-1] = -1;
				if(data[i-1]<localMinVal)
					localMinVal = data[i-1];
				break;
			default:
				flags[i-1] = 0;
				break;
		}
	}

	*minVal = localMinVal;
	*maxVal = localMaxVal;
	// after this point, all flags indices are at i. This is because we wrote the peak/valley/plateau flags to flags[i-1] in the switch above, so they're written at the maximum point of the peak.
	return;
}
/* ---------------- END utility functions ---------------------- */




/* ---------------- filterbank functions ---------------------- */
t_binIdx tIDLib_nearestBinIndex(t_float target, t_float *binFreqs, t_binIdx n)
{
	t_binIdx i, bin;
	t_float dist, bestDist;

	bestDist = FLT_MAX;
	dist = 0.0;
	bin = 0;

	for(i=0; i<n; i++)
		if((dist=fabs(binFreqs[i]-target)) < bestDist)
		{	
			bestDist = dist;
			bin = i;
		}

	return(bin);
}


// TODO: this should probably take in a pointer to x->x_filterbank, resize it to sizeFilterFreqs-2, then write the filter bound freqs in x_filterbank.filterFreqs[0] and [1]. Then x_filterbank can be passed to _createFilterbank with the filter bound freqs known, and no need for a separate x_filterFreqs buffer
// resizes the filterFreqs array and fills it with the Hz values for the Bark filter boundaries. Reports the new number of Bark frequency band boundaries based on the desired spacing. The size of the corresponding filterbank would be sizeFilterFreqs-2
t_filterIdx tIDLib_getBarkBoundFreqs(t_float **filterFreqs, t_filterIdx oldSizeFilterFreqs, t_float spacing, t_float sr)
{
	t_filterIdx i, sizeFilterFreqs;
	t_float sumBark;

	if(spacing<0.1 || spacing>6.0)
	{
		spacing = BARKSPACINGDEFAULT;
		post("Bark spacing must be between 0.1 and 6.0 Barks. Using default spacing of %f Barks instead.", BARKSPACINGDEFAULT);
	}

	sumBark = 0.0;
	sizeFilterFreqs = 0;

	while( (tIDLib_bark2freq(sumBark)<=(sr*0.5)) && (sumBark<=MAXBARKS) )
	{
		sizeFilterFreqs++;
		sumBark += spacing;
	}

	(*filterFreqs) = (t_float *)t_resizebytes((*filterFreqs), oldSizeFilterFreqs, sizeFilterFreqs*sizeof(t_float));

	// First filter boundary should be at 0Hz
	(*filterFreqs)[0] = 0.0;

	// reset the running Bark sum to the first increment past 0 Barks
	sumBark = spacing;
	
	// fill up filterFreqs with the Hz values of all Bark range boundaries
	for(i=1; i<sizeFilterFreqs; i++)
	{
		(*filterFreqs)[i] = tIDLib_bark2freq(sumBark);
		sumBark += spacing;
	}

	return(sizeFilterFreqs);
}


// TODO: this should probably take in a pointer to x->x_filterbank, resize it to sizeFilterFreqs-2, then write the filter bound freqs in x_filterbank.filterFreqs[0] and [1]. Then x_filterbank can be passed to _createFilterbank with the filter bound freqs known, and no need for a separate x_filterFreqs buffer
t_filterIdx tIDLib_getMelBoundFreqs(t_float **filterFreqs, t_filterIdx oldSizeFilterFreqs, t_float spacing, t_float sr)
{
	t_filterIdx i, sizeFilterFreqs;
	t_float sumMel;

	if(spacing<5 || spacing >1000)
	{
		spacing = MELSPACINGDEFAULT;
		post("mel spacing must be between 5 and 1000 mels. Using default spacing of %f mels instead.", MELSPACINGDEFAULT);
	}

	sumMel = 0.0;
	sizeFilterFreqs = 0;

	while( (tIDLib_mel2freq(sumMel)<=(sr*0.5)) && (sumMel<=MAXMELS) )
	{
		sizeFilterFreqs++;
		sumMel += spacing;
	}

	(*filterFreqs) = (t_float *)t_resizebytes((*filterFreqs), oldSizeFilterFreqs, sizeFilterFreqs*sizeof(t_float));

	// First filter boundary should be at 0Hz
	(*filterFreqs)[0] = 0.0;

	// reset the running Bark sum to the first increment past 0 mels
	sumMel = spacing;

	// fill up filterFreqs with the Hz values of all mel range boundaries
	for(i=1; i<sizeFilterFreqs; i++)
	{
		(*filterFreqs)[i] = tIDLib_mel2freq(sumMel);
		sumMel += spacing;
	}
	
	return(sizeFilterFreqs);
}


void tIDLib_createFilterbank(t_float *filterFreqs, t_filter **filterbank, t_filterIdx oldNumFilters, t_filterIdx newNumFilters, t_float window, t_float sr)
{
	t_binIdx windowHalf, windowHalfPlus1;
	t_float *binFreqs;
	t_filterIdx i;

	windowHalf = window*0.5;	
	windowHalfPlus1 = windowHalf+1;

	// create local memory
	binFreqs = (t_float *)t_getbytes(windowHalfPlus1*sizeof(t_float));

	// free memory for each filter
	for(i=0; i<oldNumFilters; i++)
		t_freebytes((*filterbank)[i].filter, (*filterbank)[i].size*sizeof(t_float));

	(*filterbank) = (t_filter *)t_resizebytes((*filterbank), oldNumFilters*sizeof(t_filter), newNumFilters*sizeof(t_filter));

	// initialize indices
	for(i=0; i<newNumFilters; i++)
	{
		t_uChar j;
		
		for(j=0; j<2; j++)
			(*filterbank)[i].indices[j] = 0;
	}

	// initialize filterbank sizes
	for(i=0; i<newNumFilters; i++)
		(*filterbank)[i].size = 0;

	{
		t_binIdx i;
		// first, find the actual freq for each bin based on current window size
		for(i=0; i<windowHalfPlus1; i++)
			binFreqs[i] = tIDLib_bin2freq(i, window, sr);
	}
	
	if(newNumFilters >= windowHalfPlus1)
	{
		post("timbreID WARNING: current filterbank is invalid. For window size N, filterbank size must be less than N/2+1. Change filter spacing or window size accordingly.");
		
		t_binIdx i;
		
		// fill the first filter indices with bins 0 - N/2
		for(i=0; i<windowHalfPlus1; i++)
		{
			t_binIdx j;
	
			(*filterbank)[i].size = 1;
			(*filterbank)[i].indices[0] = i;
			(*filterbank)[i].indices[1] = i;
			(*filterbank)[i].filterFreqs[0] = binFreqs[i];
			(*filterbank)[i].filterFreqs[1] = binFreqs[i];
			
			(*filterbank)[i].filter = (t_float *)t_getbytes((*filterbank)[i].size*sizeof(t_float));

			// initialize this filter
			for(j=0; j<(*filterbank)[i].size; j++)
				(*filterbank)[i].filter[j] = 0.0;
		}

		// fill out additional filter indices with the Nyquist bin, just so we avoid crashes
		for(; i<newNumFilters; i++)
		{
			t_binIdx j;

			(*filterbank)[i].size = 1;
			(*filterbank)[i].indices[0] = windowHalf;
			(*filterbank)[i].indices[1] = windowHalf;
			(*filterbank)[i].filterFreqs[0] = binFreqs[windowHalf];
			(*filterbank)[i].filterFreqs[1] = binFreqs[windowHalf];
						
			(*filterbank)[i].filter = (t_float *)t_getbytes((*filterbank)[i].size*sizeof(t_float));

			// initialize this filter
			for(j=0; j<(*filterbank)[i].size; j++)
				(*filterbank)[i].filter[j] = 0.0;
		}
	}
	else
	{
		t_filterIdx i;
		
		// finally, build the filterbank
		for(i=1; i<=newNumFilters; i++)
		{
			t_binIdx j, startIdx, peakIdx, finishIdx, filterWidth, upN, downN;
			t_float *upRamp, *downRamp;
			
			startIdx = peakIdx = finishIdx = 0;

			startIdx = tIDLib_nearestBinIndex(filterFreqs[i-1], binFreqs, windowHalfPlus1);
			peakIdx = tIDLib_nearestBinIndex(filterFreqs[i], binFreqs, windowHalfPlus1);
			finishIdx = tIDLib_nearestBinIndex(filterFreqs[i+1], binFreqs, windowHalfPlus1);
		
			// TODO: should check that startIdx<finishIdx, and swap if not. and peakIdx should be between those bounds
			
			// grab memory for this filter
			filterWidth = finishIdx-startIdx + 1;
			filterWidth = (filterWidth<1)?1:filterWidth;
				
			(*filterbank)[i-1].size = filterWidth; // store the sizes for freeing memory later
			(*filterbank)[i-1].filter = (t_float *)t_getbytes(filterWidth*sizeof(t_float));

			// initialize this filter
			for(j=0; j<filterWidth; j++)
				(*filterbank)[i-1].filter[j] = 0.0;

			// some special cases for very narrow filter widths
			switch(filterWidth)
			{
				t_binIdx j, k;

				case 1:
					(*filterbank)[i-1].filter[0] = 1.0;
					break;
				
				// no great way to do a triangle with a filter width of 2, so might as well average
				case 2:
					(*filterbank)[i-1].filter[0] = 0.5;
					(*filterbank)[i-1].filter[1] = 0.5;
					break;

				// with 3 and greater, we can use our ramps
				default:
					upN = peakIdx-startIdx+1;
					upRamp = (t_float *)t_getbytes(upN*sizeof(t_float));
					tIDLib_linspace(upRamp, 0.0, 1.0, upN);

					downN = finishIdx-peakIdx+1;
					downRamp = (t_float *)t_getbytes(downN*sizeof(t_float));
					tIDLib_linspace(downRamp, 1.0, 0.0, downN);
					
					// copy into (*filterbank)[i-1].filter
					for(j=0; j<upN; j++)
						(*filterbank)[i-1].filter[j] = upRamp[j];

					// start at k=1 because k=0 will be the peak (i.e., 1.0)
					for(k=1; k<downN; j++, k++)
						(*filterbank)[i-1].filter[j] = downRamp[k];

					// clip the triangle within 0 and 1, just in case
					for(j=0; j<filterWidth; j++)
					{
						if((*filterbank)[i-1].filter[j] < 0.0)
							(*filterbank)[i-1].filter[j] = 0.0;

						if((*filterbank)[i-1].filter[j] > 1.0)
							(*filterbank)[i-1].filter[j] = 1.0;
					}
						
					t_freebytes(upRamp, upN*sizeof(t_float));
					t_freebytes(downRamp, downN*sizeof(t_float));
					break;
			};

			(*filterbank)[i-1].indices[0] = startIdx;
			(*filterbank)[i-1].indices[1] = finishIdx;
			(*filterbank)[i-1].filterFreqs[0] = binFreqs[startIdx];
			(*filterbank)[i-1].filterFreqs[1] = binFreqs[finishIdx];
		}
	}
	
	// free local memory
	t_freebytes(binFreqs, windowHalfPlus1*sizeof(t_float));
}


void tIDLib_specFilterBands(t_binIdx n, t_filterIdx numFilters, t_float *spectrum, t_filter *filterbank, t_bool normalize)
{
	t_float smoothedSpec[numFilters], totalEnergy;
	t_filterIdx i;

	totalEnergy = 0;
	
	for(i=0; i<numFilters; i++)
	{
		t_binIdx j;
	   	smoothedSpec[i] = 0.0;

		for(j=filterbank[i].indices[0]; j<=filterbank[i].indices[1]; j++)
	    		smoothedSpec[i] += spectrum[j];

		smoothedSpec[i] /= filterbank[i].size;
		
		totalEnergy += smoothedSpec[i];
	};
	
	// Check that the spectrum window size N is larger than the number of filters, otherwise we'll be writing to invalid memory indices
	if(n>=numFilters)
	{
		t_filterIdx i;

		for(i=0; i<numFilters; i++)
			if(normalize)
				spectrum[i] = smoothedSpec[i]/totalEnergy;
			else
				spectrum[i] = smoothedSpec[i];
	}
	else
	{
		t_filterIdx i;

		for(i=0; i<n; i++)
			if(normalize)
				spectrum[i] = smoothedSpec[i]/totalEnergy;
			else
				spectrum[i] = smoothedSpec[i];
	}
}


void tIDLib_filterbankMultiply(t_float *spectrum, t_bool normalize, t_bool filterAvg, t_filter *filterbank, t_filterIdx numFilters)
{
	t_float sumSum, *filterPower;
	t_filterIdx i;

	// create local memory
	filterPower = (t_float *)t_getbytes(numFilters*sizeof(t_float));

 	sumSum = 0;

	for(i=0; i<numFilters; i++)
	{
	   	t_float sum = 0.0;
		t_binIdx j, k;
		
		for(j=filterbank[i].indices[0], k=0; j<=filterbank[i].indices[1]; j++, k++)
	    	sum += spectrum[j] * filterbank[i].filter[k];

		k = filterbank[i].size;
		if(filterAvg)
			sum /= k;

		filterPower[i] = sum;  // get the total power.  another weighting might be better.

 		sumSum += sum;  // normalize so power in all bands sums to 1
	};

	if(normalize)
	{
		// prevent divide by 0
		if(sumSum==0)
			sumSum=1;
		else
			sumSum = 1.0/sumSum; // take the reciprocal here to save a divide below
	}
	else
		sumSum=1.0;
	
	{
		t_binIdx i;
		for(i=0; i<numFilters; i++)
			spectrum[i] = filterPower[i]*sumSum;
	}
	
	// free local memory
	t_freebytes(filterPower, numFilters*sizeof(t_float));
}


void tIDLib_cosineTransform(t_float *output, t_sample *input, t_filterIdx numFilters)
{
	t_float piOverNfilters;
	t_filterIdx i;

	piOverNfilters = M_PI/numFilters; // save multiple divides below

	for(i=0; i<numFilters; i++)
    {
		t_filterIdx k;
	   	output[i] = 0.0;

		for(k=0; k<numFilters; k++)
 	    	output[i] += input[k] * cos(i * (k+0.5) * piOverNfilters);  // DCT-II
	};
}
/* ---------------- END filterbank functions ---------------------- */




/* ---------------- stat computation functions ---------------------- */
t_float tIDLib_computeCentroid(t_binIdx n, t_float *spectrum, t_float *freqList, t_float energySum)
{
	t_float dividend, divisor, centroid;
	t_binIdx i;

	dividend = centroid = 0.0;
	divisor = energySum;

	for(i=0; i<n; i++)
		dividend += spectrum[i]*freqList[i];  // weight by bin freq

	if(divisor<=0.0)
		return(-1.0);
	else
	{
		centroid = dividend/divisor;
		return(centroid);
	}
}

t_float tIDLib_computeSpread(t_binIdx n, t_float *spectrum, t_float *freqList, t_float energySum, t_float centroid)
{
	t_float dividend, divisor, spread;
	t_binIdx i;

	dividend = spread = 0.0;
 	divisor = energySum;

	for(i=0; i<n; i++)
		dividend += powf((freqList[i] - centroid), 2) * spectrum[i];

	if(divisor<=0.0)
		return(-1.0);
	else
	{
		spread = sqrt(dividend/divisor);
		return(spread);
	}
}

t_float tIDLib_computeSkewness(t_binIdx n, t_float *spectrum, t_float *freqList, t_float energySum, t_float centroid, t_float spread)
{
	t_float dividend, divisor, skewness;
	t_binIdx i;

	dividend = skewness = 0.0;
	divisor = energySum;

	for(i=0; i<n; i++)
		dividend += powf((freqList[i] - centroid), 3) * spectrum[i];

	if(divisor<=0.0 || spread<=0.0)
		return(-1.0);
	else
	{
		spread = powf(spread, 3);	
		skewness = (dividend/divisor)/spread;
		return(skewness);
	}
}

t_float tIDLib_computeKurtosis(t_binIdx n, t_float *spectrum, t_float *freqList, t_float energySum, t_float centroid, t_float spread)
{
	t_float dividend, divisor, kurtosis;
	t_binIdx i;

	dividend = kurtosis = 0.0;
	divisor = energySum;

	for(i=0; i<n; i++)
		dividend += powf((freqList[i] - centroid), 4) * spectrum[i];


	if(divisor<=0.0 || spread<=0.0)
		return(-1.0);
	else
	{
		spread = powf(spread, 4);	
		kurtosis = (dividend/divisor)/spread;
		return(kurtosis);
	}
}
/* ---------------- END stat computation functions ---------------------- */




/* ---------------- windowing buffer functions ---------------------- */

void tIDLib_blackmanWindow(t_float *wPtr, t_sampIdx n)
{
	t_sampIdx i;
	for(i=0; i<n; i++, wPtr++)
    	*wPtr = 0.42 - (0.5 * cos(2*M_PI*i/n)) + (0.08 * cos(4*M_PI*i/n));
}

void tIDLib_cosineWindow(t_float *wPtr, t_sampIdx n)
{
	t_sampIdx i;
	for(i=0; i<n; i++, wPtr++)
    	*wPtr = sin(M_PI*i/n);
}

void tIDLib_hammingWindow(t_float *wPtr, t_sampIdx n)
{
	t_sampIdx i;
	for(i=0; i<n; i++, wPtr++)
    	*wPtr = 0.5 - (0.46 * cos(2*M_PI*i/n));
}

void tIDLib_hannWindow(t_float *wPtr, t_sampIdx n)
{
	t_sampIdx i;
	for(i=0; i<n; i++, wPtr++)
    	*wPtr = 0.5 * (1 - cos(2*M_PI*i/n));
}
/* ---------------- END windowing buffer functions ---------------------- */




/* ---------------- dsp utility functions ---------------------- */
t_float tIDLib_ampDB(t_sampIdx n, t_sample *input)
{
	t_float rms, dB;
	t_sampIdx i;	

	rms = 0.0;
	
	for(i=0; i<n; i++, input++)
		rms += *input * *input;
	
	rms /= n;
	rms = sqrt(rms);
	
	// from Pd source: d_arithmetic.c
	dB = rmstodb(rms);
	
	return(dB);
}

void tIDLib_peakSample(t_sampIdx n, t_float *input, t_sampIdx *peakIdx, t_float *peakVal)
{
	t_sampIdx i;
	
	*peakVal = -FLT_MAX;
	*peakIdx = ULONG_MAX;
	
	for(i=0; i<n; i++, input++)
	{
		if(fabs(*input) > *peakVal)
		{
			*peakIdx = i;
			*peakVal = fabs(*input);
		}
	}
}

t_sampIdx tIDLib_findAttackStartSamp(t_sampIdx n, t_float *input, t_float sampMagThresh, t_uShortInt numSampsThresh)
{
	t_sampIdx i, j, startSamp;
		
	startSamp = ULONG_MAX;
	
	i=n;
	
	while(i--)
	{
		if(fabs(input[i]) <= sampMagThresh)
		{
			t_uShortInt sampCount;

			sampCount = 1;
			j=i;

			while(j--)
			{
				if(fabs(input[j]) <= sampMagThresh)
				{
					sampCount++;
				
					if(sampCount>=numSampsThresh)
					{
						startSamp = j;
						return(startSamp);
					}
				}
				else
					break;
			}
		}
	}

	return(startSamp);
}

// this could also return the location of the zero crossing
t_float tIDLib_zeroCrossingRate(t_sampIdx n, t_sample *input)
{
	t_float crossings;
	t_sampIdx i;

	crossings = 0.0;
	
	for(i=1; i<n; i++)
		crossings += abs(tIDLib_signum(input[i]) - tIDLib_signum(input[i-1]));
	
	crossings *= 0.5;
		
	return(crossings);
}

void tIDLib_power(t_binIdx n, void *fftw_out, t_float *powBuf)
{
	fftwf_complex *fftw_out_local = (fftwf_complex *)fftw_out;

	while(n--)
    	powBuf[n] = (fftw_out_local[n][0] * fftw_out_local[n][0]) + (fftw_out_local[n][1] * fftw_out_local[n][1]);

}

void tIDLib_mag(t_binIdx n, t_float *input)
{
	while(n--)
	{
	    *input = sqrt(*input);
	    input++;
	}
}

// for normalizing spectra so sum of energy==1
void tIDLib_normal(t_binIdx n, t_float *input)
{
	t_float sum, normScalar;
	t_binIdx i;

	sum = normScalar = 0.0;
	
	for(i=0; i<n; i++)
		sum += input[i];

	sum = (sum==0.0)?1.0:sum;
			
	normScalar = 1.0/sum;
	
	while(n--)
		*input++ *= normScalar;
}

// for normalizing bipolar waveforms so that peak amplitude==1
void tIDLib_normalPeak(t_binIdx n, t_float *input)
{
	t_float max, normScalar;
	t_binIdx i;

	max = normScalar = 0.0;
	
	for(i=0; i<n; i++)
		if(fabs(input[i])>max)
			max = fabs(input[i]);

	max = (max==0.0)?1.0:max;
			
	normScalar = 1.0/max;
	
	while(n--)
		*input++ *= normScalar;
}

void tIDLib_log(t_binIdx n, t_float *input)
{
	while (n--)
    {
		// if to protect against log(0)
    	if(*input==0.0)
    		*input = 0.0;
    	else
	        *input = log(*input);
	        
        input++;
    };
}
/* ---------------- END dsp utility functions ---------------------- */

