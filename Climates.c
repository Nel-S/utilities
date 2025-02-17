#include "Climates.h"

/*
	TODO: Check if all climate samples are truly truncated after four decimal places, as Cubiomes suggests.
	If so, I have a lot of fixing to do.
*/

// The maximum value that a Perlin sample can return.
const double U_MAX_PERLIN_VALUE = 1.0363538112118025;
// The minimum value that a Perlin sample can return.
const double U_MIN_PERLIN_VALUE = -U_MAX_PERLIN_VALUE;
// The values that, alongside their negative inverses, encompass the middle N% likelihood of an average Perlin sample.
const double U_PERLIN_BENCHMARKS[NUMBER_OF_PERCENTILES] = {0.0295836, 0.0611534, 0.0945022, 0.1298043, 0.167999, 0.2099861, 0.257027, 0.317638, 0.4037141, U_MAX_PERLIN_VALUE};

// The maximum value each climate sample can return. 
const double U_MAX_CLIMATE_AMPLITUDES[NP_MAX] = {20./9  * U_MAX_PERLIN_VALUE, 320./189 * U_MAX_PERLIN_VALUE, 267./73 * U_MAX_PERLIN_VALUE,
										   75./31 * U_MAX_PERLIN_VALUE,  28./3   * U_MAX_PERLIN_VALUE,  20./7  * U_MAX_PERLIN_VALUE};
// The minimum value each climate sample can return.
const double U_MIN_CLIMATE_AMPLITUDES[NP_MAX] = {20./9  * U_MIN_PERLIN_VALUE, 320./189 * U_MIN_PERLIN_VALUE, 267./73 * U_MIN_PERLIN_VALUE,
										   75./31 * U_MIN_PERLIN_VALUE,  28./3   * U_MIN_PERLIN_VALUE,  20./7  * U_MIN_PERLIN_VALUE};

// The total number of Perlin octaves across all climates.
const int U_NUMBER_OF_OCTAVES = 46;
// The number of Perlin octaves each climate has.
// If a constant compile-time variant is needed instead, `sizeof(U_<MAX/MIN>_<climate>_OCTAVE_AMPLITUDE_SUMS)/sizeof(*U_<MAX/MIN>_<climate>_OCTAVE_AMPLITUDE_SUMS)`
// will return the same value.
const int U_CLIMATE_NUMBER_OF_OCTAVES[NP_MAX] = {4, 4, 18, 8, 6, 6};
// The BiomeNoise-set amplitude for each climate.
const double U_SET_CLIMATE_AMPLITUDES[NP_MAX] = {5./4, 10./9, 3./2, 25./18, 5./4, 5./4};

// The maximum value that all temperature Perlin octaves after the `i`th one can return. (This is with the climate amplitude factored in.)
const double U_MAX_TEMP_OCTAVE_AMPLITUDE_SUMS[4] = {80./63 * U_MAX_PERLIN_VALUE, 20./63 * U_MAX_PERLIN_VALUE, 10./63 * U_MAX_PERLIN_VALUE, 0.};

// The minimum value that all temperature Perlin octaves after the `i`th one can return. (This is with the climate amplitude factored in.)
const double U_MIN_TEMP_OCTAVE_AMPLITUDE_SUMS[4] = {80./63 * U_MIN_PERLIN_VALUE, 20./63 * U_MIN_PERLIN_VALUE, 10./63 * U_MIN_PERLIN_VALUE, 0.};

// The maximum value that all temperature Perlin octaves after the `i`th one can return. (This is with the climate amplitude factored in.)
const double U_MAX_HUMID_OCTAVE_AMPLITUDE_SUMS[4] = {640./567 * U_MAX_PERLIN_VALUE, 320./567 * U_MAX_PERLIN_VALUE, 160./567 * U_MAX_PERLIN_VALUE, 0.};

// The minimum value that all temperature Perlin octaves after the `i`th one can return. (This is with the climate amplitude factored in.)
const double U_MIN_HUMID_OCTAVE_AMPLITUDE_SUMS[4] = {640./567 * U_MIN_PERLIN_VALUE, 320./567 * U_MIN_PERLIN_VALUE, 160./567 * U_MIN_PERLIN_VALUE, 0.};

// The maximum value that all continentalness Perlin octaves after the `i`th one can return. (This is with the climate amplitude factored in.)
const double U_MAX_CONT_OCTAVE_AMPLITUDE_SUMS[18] = {1485./511  * U_MAX_PERLIN_VALUE, 1101./511 * U_MAX_PERLIN_VALUE, 909./511 * U_MAX_PERLIN_VALUE, 717./511 * U_MAX_PERLIN_VALUE,
													   75./73   * U_MAX_PERLIN_VALUE,  333./511 * U_MAX_PERLIN_VALUE, 237./511 * U_MAX_PERLIN_VALUE, 141./511 * U_MAX_PERLIN_VALUE,
													   93./511  * U_MAX_PERLIN_VALUE,   45./511 * U_MAX_PERLIN_VALUE,  33./511 * U_MAX_PERLIN_VALUE,   3./73  * U_MAX_PERLIN_VALUE,
													   15./511  * U_MAX_PERLIN_VALUE,    9./511 * U_MAX_PERLIN_VALUE,   6./511 * U_MAX_PERLIN_VALUE,   3./511 * U_MAX_PERLIN_VALUE,
														3./1022 * U_MAX_PERLIN_VALUE,    0.};
// The minimum value that all continentalness Perlin octaves after the `i`th one can return. (This is with the climate amplitude factored in.)
const double U_MIN_CONT_OCTAVE_AMPLITUDE_SUMS[18] = {1485./511  * U_MIN_PERLIN_VALUE, 1101./511 * U_MIN_PERLIN_VALUE, 909./511 * U_MIN_PERLIN_VALUE, 717./511 * U_MIN_PERLIN_VALUE,
													   75./73   * U_MIN_PERLIN_VALUE,  333./511 * U_MIN_PERLIN_VALUE, 237./511 * U_MIN_PERLIN_VALUE, 141./511 * U_MIN_PERLIN_VALUE,
													   93./511  * U_MIN_PERLIN_VALUE,   45./511 * U_MIN_PERLIN_VALUE,  33./511 * U_MIN_PERLIN_VALUE,   3./73  * U_MIN_PERLIN_VALUE,
													   15./511  * U_MIN_PERLIN_VALUE,    9./511 * U_MIN_PERLIN_VALUE,   6./511 * U_MIN_PERLIN_VALUE,   3./511 * U_MIN_PERLIN_VALUE,
														3./1022 * U_MIN_PERLIN_VALUE,    0.};
// The maximum value that all erosion Perlin octaves after the `i`th one can return. (This is with the climate amplitude factored in.)
const double U_MAX_EROS_OCTAVE_AMPLITUDE_SUMS[8] = {475./279 * U_MAX_PERLIN_VALUE, 275./279 * U_MAX_PERLIN_VALUE, 175./279  * U_MAX_PERLIN_VALUE, 25./93  * U_MAX_PERLIN_VALUE,
													 50./279 * U_MAX_PERLIN_VALUE,  25./279 * U_MAX_PERLIN_VALUE,  25./558 * U_MAX_PERLIN_VALUE,   0.};
// The minimum value that all erosion Perlin octaves after the `i`th one can return. (This is with the climate amplitude factored in.)
const double U_MIN_EROS_OCTAVE_AMPLITUDE_SUMS[8] = {475./279 * U_MIN_PERLIN_VALUE, 275./279 * U_MIN_PERLIN_VALUE, 175./279  * U_MIN_PERLIN_VALUE, 25./93  * U_MIN_PERLIN_VALUE,
													 50./279 * U_MIN_PERLIN_VALUE,  25./279 * U_MIN_PERLIN_VALUE,  25./558 * U_MIN_PERLIN_VALUE,   0.};
// The maximum value that all shift Perlin octaves after the `i`th one can return. (This is with the climate amplitude factored in.)
// Note that this is not the same as the Depth climate, which is not included here.
const double U_MAX_SHIFT_OCTAVE_AMPLITUDE_SUMS[6] = {20./3 * U_MAX_PERLIN_VALUE, 4 * U_MAX_PERLIN_VALUE, 8./3 * U_MAX_PERLIN_VALUE, 4./3 * U_MAX_PERLIN_VALUE,
													  2./3 * U_MAX_PERLIN_VALUE, 0.};
// The minimum value that all shift Perlin octaves after the `i`th one can return. (This is with the climate amplitude factored in.)
// Note that this is not the same as the Depth climate, which is not included here.
const double U_MIN_SHIFT_OCTAVE_AMPLITUDE_SUMS[6] = {20./3 * U_MIN_PERLIN_VALUE, 4 * U_MIN_PERLIN_VALUE, 8./3 * U_MIN_PERLIN_VALUE, 4./3 * U_MIN_PERLIN_VALUE,
													  2./3 * U_MIN_PERLIN_VALUE, 0.};
// The maximum value that all weirdness Perlin octaves after the `i`th one can return. (This is with the climate amplitude factored in.)
const double U_MAX_WEIRD_OCTAVE_AMPLITUDE_SUMS[6] = {20./9  * U_MAX_PERLIN_VALUE, 100./63 * U_MAX_PERLIN_VALUE, 20./21 * U_MAX_PERLIN_VALUE, 20./63 * U_MAX_PERLIN_VALUE,
													 10./63 * U_MAX_PERLIN_VALUE, 0.};
// The minimum value that all weirdness Perlin octaves after the `i`th one can return. (This is with the climate amplitude factored in.)
const double U_MIN_WEIRD_OCTAVE_AMPLITUDE_SUMS[6] = {20./9  * U_MIN_PERLIN_VALUE, 100./63 * U_MIN_PERLIN_VALUE, 20./21 * U_MIN_PERLIN_VALUE, 20./63 * U_MIN_PERLIN_VALUE,
													 10./63 * U_MIN_PERLIN_VALUE, 0.};



int U_initClimateBoundsArray(const int climate, const double desiredClimateSample, const int percentile, double *array, const size_t arraySize) {
	static const double *CLIMATE_ARRAYS[] = {U_MAX_TEMP_OCTAVE_AMPLITUDE_SUMS,  U_MIN_TEMP_OCTAVE_AMPLITUDE_SUMS,
											 U_MAX_HUMID_OCTAVE_AMPLITUDE_SUMS, U_MIN_HUMID_OCTAVE_AMPLITUDE_SUMS,
											 U_MAX_CONT_OCTAVE_AMPLITUDE_SUMS,  U_MIN_CONT_OCTAVE_AMPLITUDE_SUMS,
											 U_MAX_EROS_OCTAVE_AMPLITUDE_SUMS,  U_MIN_EROS_OCTAVE_AMPLITUDE_SUMS,
											 U_MAX_SHIFT_OCTAVE_AMPLITUDE_SUMS, U_MIN_SHIFT_OCTAVE_AMPLITUDE_SUMS,
											 U_MAX_WEIRD_OCTAVE_AMPLITUDE_SUMS, U_MIN_WEIRD_OCTAVE_AMPLITUDE_SUMS};
	const double *OCTAVE_RANGES_ARRAY = CLIMATE_ARRAYS[2*climate + (desiredClimateSample < 0.)];
	int count = 0;
	for (size_t i = 0; i < min(arraySize, (size_t)(U_CLIMATE_NUMBER_OF_OCTAVES[climate])); ++i, ++count) {
		array[i] = desiredClimateSample - OCTAVE_RANGES_ARRAY[i]*U_PERLIN_BENCHMARKS[percentile]/U_MAX_PERLIN_VALUE;
		if (climate == NP_SHIFT && U_CLIMATE_NUMBER_OF_OCTAVES[NP_SHIFT] + i < arraySize) {
			array[U_CLIMATE_NUMBER_OF_OCTAVES[NP_SHIFT] + i] = array[i];
			++count;
		}
	}
	return count;
}

void U_manualBNinit(BiomeNoise *biomeNoise, const bool largeBiomesFlag) {
	if (biomeNoise->mc < MC_1_18) return;
	const int LB_MULT = 3*largeBiomesFlag + 1;

	biomeNoise->climate[NP_TEMPERATURE].amplitude = biomeNoise->climate[NP_SHIFT].amplitude = biomeNoise->climate[NP_WEIRDNESS].amplitude = 5./4;
	biomeNoise->climate[NP_HUMIDITY].amplitude = 10./9;
	biomeNoise->climate[NP_CONTINENTALNESS].amplitude = 3./2;
	biomeNoise->climate[NP_EROSION].amplitude = 25./18;
	biomeNoise->climate[NP_TEMPERATURE].octA.octaves = &biomeNoise->oct[0];
	biomeNoise->climate[NP_TEMPERATURE].octB.octaves = &biomeNoise->oct[2];
	biomeNoise->climate[NP_HUMIDITY].octA.octaves = &biomeNoise->oct[4];
	biomeNoise->climate[NP_HUMIDITY].octB.octaves = &biomeNoise->oct[6];
	biomeNoise->climate[NP_CONTINENTALNESS].octA.octaves = &biomeNoise->oct[8];
	biomeNoise->climate[NP_CONTINENTALNESS].octB.octaves = &biomeNoise->oct[17];
	biomeNoise->climate[NP_EROSION].octA.octaves = &biomeNoise->oct[26];
	biomeNoise->climate[NP_EROSION].octB.octaves = &biomeNoise->oct[30];
	biomeNoise->climate[NP_SHIFT].octA.octaves = &biomeNoise->oct[34];
	biomeNoise->climate[NP_SHIFT].octB.octaves = &biomeNoise->oct[37];
	biomeNoise->climate[NP_WEIRDNESS].octA.octaves = &biomeNoise->oct[40];
	biomeNoise->climate[NP_WEIRDNESS].octB.octaves = &biomeNoise->oct[43];
	biomeNoise->climate[NP_TEMPERATURE].octA.octcnt = biomeNoise->climate[NP_TEMPERATURE].octB.octcnt = biomeNoise->climate[NP_HUMIDITY].octA.octcnt = biomeNoise->climate[NP_HUMIDITY].octB.octcnt = 2;
	biomeNoise->climate[NP_CONTINENTALNESS].octA.octcnt = biomeNoise->climate[NP_CONTINENTALNESS].octB.octcnt = 9;
	biomeNoise->climate[NP_EROSION].octA.octcnt = biomeNoise->climate[NP_EROSION].octB.octcnt = 4;
	biomeNoise->climate[NP_SHIFT].octA.octcnt = biomeNoise->climate[NP_SHIFT].octB.octcnt = biomeNoise->climate[NP_WEIRDNESS].octA.octcnt = biomeNoise->climate[NP_WEIRDNESS].octB.octcnt = 3;
	// Temperature
	biomeNoise->oct[0].amplitude   = biomeNoise->oct[2].amplitude   = 16./21;
	biomeNoise->oct[0].lacunarity  = biomeNoise->oct[2].lacunarity  = 1./(1024*LB_MULT);
	biomeNoise->oct[1].amplitude   = biomeNoise->oct[3].amplitude   = 8./63;
	biomeNoise->oct[1].lacunarity  = biomeNoise->oct[3].lacunarity  = 1./(256*LB_MULT);
	// Humidity
	biomeNoise->oct[4].amplitude   = biomeNoise->oct[6].amplitude   = 32./63;
	biomeNoise->oct[4].lacunarity  = biomeNoise->oct[6].lacunarity  = 1./(256*LB_MULT);
	biomeNoise->oct[5].amplitude   = biomeNoise->oct[7].amplitude   = 16./63;
	biomeNoise->oct[5].lacunarity  = biomeNoise->oct[7].lacunarity  = 1./(128*LB_MULT);
	// Continentalness
	biomeNoise->oct[8].amplitude   = biomeNoise->oct[17].amplitude  = 256./511;
	biomeNoise->oct[8].lacunarity  = biomeNoise->oct[17].lacunarity = 1./(512*LB_MULT);
	biomeNoise->oct[9].amplitude   = biomeNoise->oct[18].amplitude  = biomeNoise->oct[10].amplitude = biomeNoise->oct[19].amplitude  = 128./511;
	biomeNoise->oct[9].lacunarity  = biomeNoise->oct[18].lacunarity = 1./(256*LB_MULT);
	biomeNoise->oct[10].lacunarity = biomeNoise->oct[19].lacunarity = 1./(128*LB_MULT);
	biomeNoise->oct[11].amplitude  = biomeNoise->oct[20].amplitude  = 64./511;
	biomeNoise->oct[11].lacunarity = biomeNoise->oct[20].lacunarity = 1./(64*LB_MULT);
	biomeNoise->oct[12].amplitude  = biomeNoise->oct[21].amplitude  = 32./511;
	biomeNoise->oct[12].lacunarity = biomeNoise->oct[21].lacunarity = 1./(32*LB_MULT);
	biomeNoise->oct[13].amplitude  = biomeNoise->oct[22].amplitude  = 8./511;
	biomeNoise->oct[13].lacunarity = biomeNoise->oct[22].lacunarity = 1./(16*LB_MULT);
	biomeNoise->oct[14].amplitude  = biomeNoise->oct[23].amplitude  = 4./511;
	biomeNoise->oct[14].lacunarity = biomeNoise->oct[23].lacunarity = 1./(8*LB_MULT);
	biomeNoise->oct[15].amplitude  = biomeNoise->oct[24].amplitude  = 2./511;
	biomeNoise->oct[15].lacunarity = biomeNoise->oct[24].lacunarity = 1./(4*LB_MULT);
	biomeNoise->oct[16].amplitude  = biomeNoise->oct[25].amplitude  = 1./511;
	biomeNoise->oct[16].lacunarity = biomeNoise->oct[25].lacunarity = 1./(2*LB_MULT);
	// Erosion
	biomeNoise->oct[26].amplitude  = biomeNoise->oct[30].amplitude  = 16./31;
	biomeNoise->oct[26].lacunarity = biomeNoise->oct[30].lacunarity = 1./(512*LB_MULT);
	biomeNoise->oct[27].amplitude  = biomeNoise->oct[31].amplitude  = 8./31;
	biomeNoise->oct[27].lacunarity = biomeNoise->oct[31].lacunarity = 1./(256*LB_MULT);
	biomeNoise->oct[28].amplitude  = biomeNoise->oct[32].amplitude  = 2./31;
	biomeNoise->oct[28].lacunarity = biomeNoise->oct[32].lacunarity = 1./(64*LB_MULT);
	biomeNoise->oct[29].amplitude  = biomeNoise->oct[33].amplitude  = 1./31;
	biomeNoise->oct[29].lacunarity = biomeNoise->oct[33].lacunarity = 1./(32*LB_MULT);
	// Shift
	biomeNoise->oct[34].amplitude  = biomeNoise->oct[37].amplitude  = 8./3;
	biomeNoise->oct[34].lacunarity = biomeNoise->oct[37].lacunarity = 1./8;
	biomeNoise->oct[35].amplitude  = biomeNoise->oct[38].amplitude  = 4./3;
	biomeNoise->oct[35].lacunarity = biomeNoise->oct[38].lacunarity = 1./4;
	biomeNoise->oct[36].amplitude  = biomeNoise->oct[39].amplitude  = 2./3;
	biomeNoise->oct[36].lacunarity = biomeNoise->oct[39].lacunarity = 1./2;
	// Weirdness
	biomeNoise->oct[40].amplitude  = biomeNoise->oct[43].amplitude  = biomeNoise->oct[41].amplitude = biomeNoise->oct[44].amplitude = 32./63;
	biomeNoise->oct[40].lacunarity = biomeNoise->oct[43].lacunarity = 1./128;
	biomeNoise->oct[41].lacunarity = biomeNoise->oct[44].lacunarity = 1./64;
	biomeNoise->oct[42].amplitude  = biomeNoise->oct[45].amplitude  = 8./63;
	biomeNoise->oct[42].lacunarity = biomeNoise->oct[45].lacunarity = 1./32;
}

void U_initPerlin(PerlinNoise *octave, Xoroshiro *xoroshiro) {
	octave->a = xNextDouble(xoroshiro) * 256.;
	octave->b = xNextDouble(xoroshiro) * 256.;
	octave->c = xNextDouble(xoroshiro) * 256.;
	uint8_t *idx = octave->d;
	for (int i = 0; i < 256; ++i) idx[i] = i;
	for (int i = 0; i < 255; ++i) {
		int j = xNextInt(xoroshiro, 256 - i) + i;
		uint_fast8_t n = idx[i];
		idx[i] = idx[j];
		idx[j] = n;
	}
	idx[256] = idx[0];
	double i2 = floor(octave->b);
	double d2 = octave->b - i2;
	octave->h2 = (int)i2;
	octave->d2 = d2;
	octave->t2 = d2*d2*d2 * (d2 * (d2*6. - 15.) + 10.);
}

void U_initClimate(const int climate, PerlinNoise *octaves, const uint64_t seed, const bool largeBiomesFlag) {
	double px = 0., pz = 0., climateSample;
	U_initAndSampleClimateBounded(climate, octaves, &px, &pz, NULL, NULL, &seed, largeBiomesFlag, &climateSample);
}

void U_initAllClimates(PerlinNoise *octaves, const uint64_t seed, const bool largeBiomesFlag) {
	for (int climate = 0; climate < NP_MAX; ++climate) U_initClimate(climate, octaves, seed, largeBiomesFlag);
}

// TODO: See if any way exists to specify `octaves` as `const` while preserving deferral to `U_initAndSampleClimateBounded`, or see if it doesn't make a difference performance-wise
double U_sampleClimate(const int climate, const PerlinNoise *octaves, double *px, double *pz, const bool largeBiomesFlag) {
	double climateSample = -INFINITY;
	// U_initAndSampleClimateBounded(climate, octaves, px, pz, NULL, NULL, NULL, largeBiomesFlag, &climateSample);
	if (climate == NP_SHIFT) U_sampleShiftBounded(octaves, px, pz, NULL, NULL);
	else U_sampleClimateBounded(climate, octaves, *px, *pz, NULL, NULL, largeBiomesFlag, &climateSample);
	return climateSample;
}

int U_sampleShiftBounded(const PerlinNoise *octaves, double *px, double *pz, const double *lowerBounds, const double *upperBounds) {
	// return U_initAndSampleClimateBounded(climate, octaves, px, pz, lowerBounds, upperBounds, NULL, largeBiomesFlag, climateSample);
	const double OFF = 337./331;

	*px = floor(*px/4.), *pz = floor(*pz/4.);
	const double origPx = *px, origPz = *pz;
	*px += 8./3 * samplePerlin(&octaves[34], origPx/8., 0, origPz/8., 0., 0.);
	if ((lowerBounds && *px < lowerBounds[0]) || (upperBounds && *px > upperBounds[0])) return 0;
	*pz += 8./3 * samplePerlin(&octaves[34], origPz/8., origPx/8., 0, 0., 0.);
	if ((lowerBounds && *pz < lowerBounds[6]) || (upperBounds && *pz > upperBounds[6])) return 1;
	*px += 8./3 * samplePerlin(&octaves[37], origPx/8.*OFF, 0, origPz/8.*OFF, 0, 0);
	if ((lowerBounds && *px < lowerBounds[1]) || (upperBounds && *px > upperBounds[1])) return 2;
	*pz += 8./3 * samplePerlin(&octaves[37], origPz/8.*OFF, origPx/8.*OFF, 0, 0, 0);
	if ((lowerBounds && *pz < lowerBounds[7]) || (upperBounds && *pz > upperBounds[7])) return 3;
	*px += 4./3 * samplePerlin(&octaves[35], origPx/4., 0, origPz/4., 0, 0);
	if ((lowerBounds && *px < lowerBounds[2]) || (upperBounds && *px > upperBounds[2])) return 4;
	*pz += 4./3 * samplePerlin(&octaves[35], origPz/4., origPx/4., 0, 0, 0);
	if ((lowerBounds && *pz < lowerBounds[8]) || (upperBounds && *pz > upperBounds[8])) return 5;
	*px += 4./3 * samplePerlin(&octaves[38], origPx/4.*OFF, 0, origPz/4.*OFF, 0, 0);
	if ((lowerBounds && *px < lowerBounds[3]) || (upperBounds && *px > upperBounds[3])) return 6;
	*pz += 4./3 * samplePerlin(&octaves[38], origPz/4.*OFF, origPx/4.*OFF, 0, 0, 0);
	if ((lowerBounds && *pz < lowerBounds[9]) || (upperBounds && *pz > upperBounds[9])) return 7;
	*px += 2./3 * samplePerlin(&octaves[36], origPx/2., 0, origPz/2., 0, 0);
	if ((lowerBounds && *px < lowerBounds[4]) || (upperBounds && *px > upperBounds[4])) return 8;
	*pz += 2./3 * samplePerlin(&octaves[36], origPz/2., origPx/2., 0, 0, 0);
	if ((lowerBounds && *pz < lowerBounds[10]) || (upperBounds && *pz > upperBounds[10])) return 9;
	*px += 2./3 * samplePerlin(&octaves[39], origPx/2.*OFF, 0, origPz/2.*OFF, 0, 0);
	if ((lowerBounds && *px < lowerBounds[5]) || (upperBounds && *px > upperBounds[5])) return 10;
	*pz += 2./3 * samplePerlin(&octaves[39], origPz/2.*OFF, origPx/2.*OFF, 0, 0, 0);
	return 11 + ((!lowerBounds || *pz > lowerBounds[11]) && (!upperBounds || *pz < upperBounds[11]));
}

// TODO: See if any way exists to specify `octaves` as `const` while preserving deferral to `U_initAndSampleClimateBounded`, or see if it doesn't make a difference performance-wise
int U_sampleClimateBounded(const int climate, const PerlinNoise *octaves, double px, double pz, const double *lowerBounds, const double *upperBounds, const bool largeBiomesFlag, double *climateSample) {
	// return U_initAndSampleClimateBounded(climate, octaves, px, pz, lowerBounds, upperBounds, NULL, largeBiomesFlag, climateSample);
	if (!climateSample && climate != NP_SHIFT) return 0;
	const int LB_MULT = 3*largeBiomesFlag + 1;
	const double OFF = 337./331;

	switch (climate) {
		case NP_TEMPERATURE:
		*climateSample = 20./21 * samplePerlin(&octaves[0], px/(1024.*LB_MULT), 0., pz/(1024.*LB_MULT), 0., 0.);
		if ((lowerBounds && *climateSample < lowerBounds[0]) || (upperBounds && *climateSample > upperBounds[0])) return 0;
		*climateSample += 20./21 * samplePerlin(&octaves[2], px/(1024*LB_MULT)*OFF, 0, pz/(1024*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[1]) || (upperBounds && *climateSample > upperBounds[1])) return 1;
		*climateSample += 10./63 * samplePerlin(&octaves[1], px/(256*LB_MULT), 0, pz/(256*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[2]) || (upperBounds && *climateSample > upperBounds[2])) return 2;
		*climateSample += 10./63 * samplePerlin(&octaves[3], px/(256*LB_MULT)*OFF, 0, pz/(256*LB_MULT)*OFF, 0, 0);
		return 3 + ((!lowerBounds || *climateSample > lowerBounds[3]) && (!upperBounds || *climateSample < upperBounds[3]));

		case NP_HUMIDITY:
		*climateSample = 320./567 * samplePerlin(&octaves[4], px/(256*LB_MULT), 0, pz/(256*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[0]) || (upperBounds && *climateSample > upperBounds[0])) return 0;
		*climateSample += 320./567 * samplePerlin(&octaves[6], px/(256*LB_MULT)*OFF, 0, pz/(256*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[1]) || (upperBounds && *climateSample > upperBounds[1])) return 
		*climateSample += 160./567 * samplePerlin(&octaves[5], px/(128*LB_MULT), 0, pz/(128*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[2]) || (upperBounds && *climateSample > upperBounds[2])) return 2;
		*climateSample += 160./567 * samplePerlin(&octaves[7], px/(128*LB_MULT)*OFF, 0, pz/(128*LB_MULT)*OFF, 0, 0);
		return 3 + ((!lowerBounds || *climateSample > lowerBounds[3]) && (!upperBounds || *climateSample < upperBounds[3]));

		case NP_CONTINENTALNESS:
		*climateSample = 384./511 * samplePerlin(&octaves[8], px/(512*LB_MULT), 0, pz/(512*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[0]) || (upperBounds && *climateSample > upperBounds[0])) return 0;
		*climateSample += 384./511 * samplePerlin(&octaves[17], px/(512*LB_MULT)*OFF, 0, pz/(512*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[1]) || (upperBounds && *climateSample > upperBounds[1])) return 1;
		*climateSample += 192./511 * samplePerlin(&octaves[9], px/(256*LB_MULT), 0, pz/(256*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[2]) || (upperBounds && *climateSample > upperBounds[2])) return 2;
		*climateSample += 192./511 * samplePerlin(&octaves[18], px/(256*LB_MULT)*OFF, 0, pz/(256*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[3]) || (upperBounds && *climateSample > upperBounds[3])) return 3;
		*climateSample += 192./511 * samplePerlin(&octaves[10], px/(128*LB_MULT), 0, pz/(128*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[4]) || (upperBounds && *climateSample > upperBounds[4])) return 4;
		*climateSample += 192./511 * samplePerlin(&octaves[19], px/(128*LB_MULT)*OFF, 0, pz/(128*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[5]) || (upperBounds && *climateSample > upperBounds[5])) return 5;
		*climateSample += 96./511 * samplePerlin(&octaves[11], px/(64*LB_MULT), 0, pz/(64*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[6]) || (upperBounds && *climateSample > upperBounds[6])) return 6;
		*climateSample += 96./511 * samplePerlin(&octaves[20], px/(64*LB_MULT)*OFF, 0, pz/(64*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[7]) || (upperBounds && *climateSample > upperBounds[7])) return 7;
		*climateSample += 48./511 * samplePerlin(&octaves[12], px/(32*LB_MULT), 0, pz/(32*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[8]) || (upperBounds && *climateSample > upperBounds[8])) return 8;
		*climateSample += 48./511 * samplePerlin(&octaves[21], px/(32*LB_MULT)*OFF, 0, pz/(32*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[9]) || (upperBounds && *climateSample > upperBounds[9])) return 9;
		*climateSample += 12./511 * samplePerlin(&octaves[13], px/(16*LB_MULT), 0, pz/(16*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[10]) || (upperBounds && *climateSample > upperBounds[10])) return 10;
		*climateSample += 12./511 * samplePerlin(&octaves[22], px/(16*LB_MULT)*OFF, 0, pz/(16*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[11]) || (upperBounds && *climateSample > upperBounds[11])) return 11;
		*climateSample += 6./511 * samplePerlin(&octaves[14], px/(8*LB_MULT), 0, pz/(8*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[12]) || (upperBounds && *climateSample > upperBounds[12])) return 12;
		*climateSample += 6./511 * samplePerlin(&octaves[23], px/(8*LB_MULT)*OFF, 0, pz/(8*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[13]) || (upperBounds && *climateSample > upperBounds[13])) return 13;
		*climateSample += 3./511 * samplePerlin(&octaves[15], px/(4*LB_MULT), 0, pz/(4*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[14]) || (upperBounds && *climateSample > upperBounds[14])) return 14;
		*climateSample += 3./511 * samplePerlin(&octaves[24], px/(4*LB_MULT)*OFF, 0, pz/(4*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[15]) || (upperBounds && *climateSample > upperBounds[15])) return 15;
		*climateSample += 3./1022 * samplePerlin(&octaves[16], px/(2*LB_MULT), 0, pz/(2*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[16]) || (upperBounds && *climateSample > upperBounds[16])) return 16;
		*climateSample += 3./1022 * samplePerlin(&octaves[25], px/(2*LB_MULT)*OFF, 0, pz/(2*LB_MULT)*OFF, 0, 0);
		return 17 + ((!lowerBounds || *climateSample > lowerBounds[17]) && (!upperBounds || *climateSample < upperBounds[17]));

		case NP_EROSION:
		*climateSample = 200./279 * samplePerlin(&octaves[26], px/(512*LB_MULT), 0, pz/(512*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[0]) || (upperBounds && *climateSample > upperBounds[0])) return 0;
		*climateSample += 200./279 * samplePerlin(&octaves[30], px/(512*LB_MULT)*OFF, 0, pz/(512*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[1]) || (upperBounds && *climateSample > upperBounds[1])) return 1;
		*climateSample += 100./279 * samplePerlin(&octaves[27], px/(256*LB_MULT), 0, pz/(256*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[2]) || (upperBounds && *climateSample > upperBounds[2])) return 2;
		*climateSample += 100./279 * samplePerlin(&octaves[31], px/(256*LB_MULT)*OFF, 0, pz/(256*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[3]) || (upperBounds && *climateSample > upperBounds[3])) return 3;
		*climateSample += 25./279 * samplePerlin(&octaves[28], px/(64*LB_MULT), 0, pz/(64*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[4]) || (upperBounds && *climateSample > upperBounds[4])) return 4;
		*climateSample += 25./279 * samplePerlin(&octaves[32], px/(64*LB_MULT)*OFF, 0, pz/(64*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[5]) || (upperBounds && *climateSample > upperBounds[5])) return 5;
		*climateSample += 25./558 * samplePerlin(&octaves[29], px/(32*LB_MULT), 0, pz/(32*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[6]) || (upperBounds && *climateSample > upperBounds[6])) return 6;
		*climateSample += 25./558 * samplePerlin(&octaves[33], px/(32*LB_MULT)*OFF, 0, pz/(32*LB_MULT)*OFF, 0, 0);
		return 7 + ((!lowerBounds || *climateSample > lowerBounds[7]) && (!upperBounds || *climateSample < upperBounds[7]));

		// case NP_SHIFT: *px = floor(*px/4.), *pz = floor(*pz/4.);
		// const double origPx = *px, origPz = *pz;
		// *px += 8./3 * samplePerlin(&octaves[34], origPx/8., 0, origPz/8., 0., 0.);
		// if ((lowerBounds && *px < lowerBounds[0]) || (upperBounds && *px > upperBounds[0])) return 0;
		// *pz += 8./3 * samplePerlin(&octaves[34], origPz/8., origPx/8., 0, 0., 0.);
		// if ((lowerBounds && *pz < lowerBounds[6]) || (upperBounds && *pz > upperBounds[6])) return 1;
		// *px += 8./3 * samplePerlin(&octaves[37], origPx/8.*OFF, 0, origPz/8.*OFF, 0, 0);
		// if ((lowerBounds && *px < lowerBounds[1]) || (upperBounds && *px > upperBounds[1])) return 2;
		// *pz += 8./3 * samplePerlin(&octaves[37], origPz/8.*OFF, origPx/8.*OFF, 0, 0, 0);
		// if ((lowerBounds && *pz < lowerBounds[7]) || (upperBounds && *pz > upperBounds[7])) return 3;
		// *px += 4./3 * samplePerlin(&octaves[35], origPx/4., 0, origPz/4., 0, 0);
		// if ((lowerBounds && *px < lowerBounds[2]) || (upperBounds && *px > upperBounds[2])) return 4;
		// *pz += 4./3 * samplePerlin(&octaves[35], origPz/4., origPx/4., 0, 0, 0);
		// if ((lowerBounds && *pz < lowerBounds[8]) || (upperBounds && *pz > upperBounds[8])) return 5;
		// *px += 4./3 * samplePerlin(&octaves[38], origPx/4.*OFF, 0, origPz/4.*OFF, 0, 0);
		// if ((lowerBounds && *px < lowerBounds[3]) || (upperBounds && *px > upperBounds[3])) return 6;
		// *pz += 4./3 * samplePerlin(&octaves[38], origPz/4.*OFF, origPx/4.*OFF, 0, 0, 0);
		// if ((lowerBounds && *pz < lowerBounds[9]) || (upperBounds && *pz > upperBounds[9])) return 7;
		// *px += 2./3 * samplePerlin(&octaves[36], origPx/2., 0, origPz/2., 0, 0);
		// if ((lowerBounds && *px < lowerBounds[4]) || (upperBounds && *px > upperBounds[4])) return 8;
		// *pz += 2./3 * samplePerlin(&octaves[36], origPz/2., origPx/2., 0, 0, 0);
		// if ((lowerBounds && *pz < lowerBounds[10]) || (upperBounds && *pz > upperBounds[10])) return 9;
		// *px += 2./3 * samplePerlin(&octaves[39], origPx/2.*OFF, 0, origPz/2.*OFF, 0, 0);
		// if ((lowerBounds && *px < lowerBounds[5]) || (upperBounds && *px > upperBounds[5])) return 10;
		// *pz += 2./3 * samplePerlin(&octaves[39], origPz/2.*OFF, origPx/2.*OFF, 0, 0, 0);
		// return 11 + ((!lowerBounds || *pz > lowerBounds[11]) && (!upperBounds || *pz < upperBounds[11]));

		case NP_WEIRDNESS:
		*climateSample = 40./63 * samplePerlin(&octaves[40], px/128., 0, pz/128., 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[0]) || (upperBounds && *climateSample > upperBounds[0])) return 0;
		*climateSample += 40./63 * samplePerlin(&octaves[43], px/128.*OFF, 0, pz/128.*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[1]) || (upperBounds && *climateSample > upperBounds[1])) return 1;
		*climateSample += 40./63 * samplePerlin(&octaves[41], px/64., 0, pz/64., 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[2]) || (upperBounds && *climateSample > upperBounds[2])) return 2;
		*climateSample += 40./63 * samplePerlin(&octaves[44], px/64.*OFF, 0, pz/64.*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[3]) || (upperBounds && *climateSample > upperBounds[3])) return 3;
		*climateSample += 10./63 * samplePerlin(&octaves[42], px/32., 0, pz/32., 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[4]) || (upperBounds && *climateSample > upperBounds[4])) return 4;
		*climateSample += 10./63 * samplePerlin(&octaves[45], px/32.*OFF, 0, pz/32.*OFF, 0, 0);
		return 5 + ((!lowerBounds || *climateSample > lowerBounds[5]) && (!upperBounds || *climateSample < upperBounds[5]));
	}
	return 0;
}

int U_initAndSampleClimateBounded(const int climate, PerlinNoise *octaves, double *px, double *pz, const double *lowerBounds, const double *upperBounds, const uint64_t *seedIfInitializingClimate, const bool largeBiomesFlag, double *climateSample) {
	if (!climateSample && climate != NP_SHIFT) return 0;
	const int LB_MULT = 3*largeBiomesFlag + 1;
	const double OFF = 337./331;

	Xoroshiro pxr = {0, 0}, pxr2;
	uint64_t xloo, xhio, xlo = 0, xhi = 0, xlo2 = 0, xhi2 = 0;
	if (seedIfInitializingClimate) {
		xSetSeed(&pxr, *seedIfInitializingClimate);
		xloo = xNextLong(&pxr);
		xhio = xNextLong(&pxr);
	}
	switch (climate) {
		case NP_TEMPERATURE:
		if (seedIfInitializingClimate) {
			pxr.lo = xloo ^ (largeBiomesFlag ? 0x944b0073edf549db : 0x5c7e6b29735f0d7f);
			pxr.hi = xhio ^ (largeBiomesFlag ? 0x4ff44347e9d22b96 : 0xf7d86f1bbc734988); // md5 "minecraft:temperature_large" or "minecraft:temperature"
			xlo = xNextLong(&pxr);
			xhi = xNextLong(&pxr);
			pxr2.lo = xlo ^ (largeBiomesFlag ? 0xb198de63a8012672 : 0x36d326eed40efeb2);
			pxr2.hi = xhi ^ (largeBiomesFlag ? 0x7b84cad43ef7b5a8 : 0x5be9ce18223c636a); // md5 "octave_-12" or "octave_-10"
			U_initPerlin(&octaves[0], &pxr2);
		}
		*climateSample = 20./21 * samplePerlin(&octaves[0], *px/(1024.*LB_MULT), 0., *pz/(1024.*LB_MULT), 0., 0.);
		if ((lowerBounds && *climateSample < lowerBounds[0]) || (upperBounds && *climateSample > upperBounds[0])) return 0;

		if (seedIfInitializingClimate) {
			xlo2 = xNextLong(&pxr);
			xhi2 = xNextLong(&pxr);
			pxr2.lo = xlo2 ^ (largeBiomesFlag ? 0xb198de63a8012672 : 0x36d326eed40efeb2);
			pxr2.hi = xhi2 ^ (largeBiomesFlag ? 0x7b84cad43ef7b5a8 : 0x5be9ce18223c636a); // md5 "octave_-12" or "octave_-10"
			U_initPerlin(&octaves[2], &pxr2);
		}
		*climateSample += 20./21 * samplePerlin(&octaves[2], *px/(1024*LB_MULT)*OFF, 0, *pz/(1024*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[1]) || (upperBounds && *climateSample > upperBounds[1])) return 1;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo ^ (largeBiomesFlag ? 0x36d326eed40efeb2 : 0x0ef68ec68504005e);
			pxr2.hi = xhi ^ (largeBiomesFlag ? 0x5be9ce18223c636a : 0x48b6bf93a2789640); // md5 "octave_-10" or "octave_-8"
			U_initPerlin(&octaves[1], &pxr2);
		}
		*climateSample += 10./63 * samplePerlin(&octaves[1], *px/(256*LB_MULT), 0, *pz/(256*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[2]) || (upperBounds && *climateSample > upperBounds[2])) return 2;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo2 ^ (largeBiomesFlag ? 0x36d326eed40efeb2 : 0x0ef68ec68504005e);
			pxr2.hi = xhi2 ^ (largeBiomesFlag ? 0x5be9ce18223c636a : 0x48b6bf93a2789640); // md5 "octave_-10" or "octave_-8"
			U_initPerlin(&octaves[3], &pxr2);
		}
		*climateSample += 10./63 * samplePerlin(&octaves[3], *px/(256*LB_MULT)*OFF, 0, *pz/(256*LB_MULT)*OFF, 0, 0);
		return 3 + ((!lowerBounds || *climateSample > lowerBounds[3]) && (!upperBounds || *climateSample < upperBounds[3]));

		case NP_HUMIDITY:
		if (seedIfInitializingClimate) {
			pxr.lo = xloo ^ (largeBiomesFlag ? 0x71b8ab943dbd5301 : 0x81bb4d22e8dc168e);
			pxr.hi = xhio ^ (largeBiomesFlag ? 0xbb63ddcf39ff7a2b : 0xf1c8b4bea16303cd); // md5 "minecraft:vegetation_large" or "minecraft:vegetation"
			xlo = xNextLong(&pxr);
			xhi = xNextLong(&pxr);
			pxr2.lo = xlo ^ (largeBiomesFlag ? 0x36d326eed40efeb2 : 0x0ef68ec68504005e);
			pxr2.hi = xhi ^ (largeBiomesFlag ? 0x5be9ce18223c636a : 0x48b6bf93a2789640); // md5 "octave_-10" or "octave_-8"
			U_initPerlin(&octaves[4], &pxr2);
		}
		*climateSample = 320./567 * samplePerlin(&octaves[4], *px/(256*LB_MULT), 0, *pz/(256*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[0]) || (upperBounds && *climateSample > upperBounds[0])) return 0;

		if (seedIfInitializingClimate) {
			xlo2 = xNextLong(&pxr);
			xhi2 = xNextLong(&pxr);
			pxr2.lo = xlo2 ^ (largeBiomesFlag ? 0x36d326eed40efeb2 : 0x0ef68ec68504005e);
			pxr2.hi = xhi2 ^ (largeBiomesFlag ? 0x5be9ce18223c636a : 0x48b6bf93a2789640); // md5 "octave_-10" or "octave_-8"
			U_initPerlin(&octaves[6], &pxr2);
		}
		*climateSample += 320./567 * samplePerlin(&octaves[6], *px/(256*LB_MULT)*OFF, 0, *pz/(256*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[1]) || (upperBounds && *climateSample > upperBounds[1])) return 1;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo ^ (largeBiomesFlag ? 0x082fe255f8be6631 : 0xf11268128982754f);
			pxr2.hi = xhi ^ (largeBiomesFlag ? 0x4e96119e22dedc81 : 0x257a1d670430b0aa); // md5 "octave_-9" or "octave_-7"
			U_initPerlin(&octaves[5], &pxr2);
		}
		*climateSample += 160./567 * samplePerlin(&octaves[5], *px/(128*LB_MULT), 0, *pz/(128*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[2]) || (upperBounds && *climateSample > upperBounds[2])) return 2;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo2 ^ (largeBiomesFlag ? 0x082fe255f8be6631 : 0xf11268128982754f);
			pxr2.hi = xhi2 ^ (largeBiomesFlag ? 0x4e96119e22dedc81 : 0x257a1d670430b0aa); // md5 "octave_-9" or "octave_-7"
			U_initPerlin(&octaves[7], &pxr2);
		}
		*climateSample += 160./567 * samplePerlin(&octaves[7], *px/(128*LB_MULT)*OFF, 0, *pz/(128*LB_MULT)*OFF, 0, 0);
		return 3 + ((!lowerBounds || *climateSample > lowerBounds[3]) && (!upperBounds || *climateSample < upperBounds[3]));

		case NP_CONTINENTALNESS:
		if (seedIfInitializingClimate) {
			pxr.lo = xloo ^ (largeBiomesFlag ? 0x9a3f51a113fce8dc : 0x83886c9d0ae3a662);
			pxr.hi = xhio ^ (largeBiomesFlag ? 0xee2dbd157e5dcdad : 0xafa638a61b42e8ad); // md5 "minecraft:continentalness_large" or "minecraft:continentalness"
			xlo = xNextLong(&pxr);
			xhi = xNextLong(&pxr);
			pxr2.lo = xlo ^ (largeBiomesFlag ? 0x0fd787bfbc403ec3 : 0x082fe255f8be6631);
			pxr2.hi = xhi ^ (largeBiomesFlag ? 0x74a4a31ca21b48b8 : 0x4e96119e22dedc81); // md5 "octave_-11" or "octave_-9"
			U_initPerlin(&octaves[8], &pxr2);
		}
		*climateSample = 384./511 * samplePerlin(&octaves[8], *px/(512*LB_MULT), 0, *pz/(512*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[0]) || (upperBounds && *climateSample > upperBounds[0])) return 0;

		if (seedIfInitializingClimate) {
			xlo2 = xNextLong(&pxr);
			xhi2 = xNextLong(&pxr);
			pxr2.lo = xlo2 ^ (largeBiomesFlag ? 0x0fd787bfbc403ec3 : 0x082fe255f8be6631);
			pxr2.hi = xhi2 ^ (largeBiomesFlag ? 0x74a4a31ca21b48b8 : 0x4e96119e22dedc81); // md5 "octave_-11" or "octave_-9"
			U_initPerlin(&octaves[17], &pxr2);
		}
		*climateSample += 384./511 * samplePerlin(&octaves[17], *px/(512*LB_MULT)*OFF, 0, *pz/(512*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[1]) || (upperBounds && *climateSample > upperBounds[1])) return 1;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo ^ (largeBiomesFlag ? 0x36d326eed40efeb2 : 0x0ef68ec68504005e);
			pxr2.hi = xhi ^ (largeBiomesFlag ? 0x5be9ce18223c636a : 0x48b6bf93a2789640); // md5 "octave_-10" or "octave_-8"
			U_initPerlin(&octaves[9], &pxr2);
		}
		*climateSample += 192./511 * samplePerlin(&octaves[9], *px/(256*LB_MULT), 0, *pz/(256*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[2]) || (upperBounds && *climateSample > upperBounds[2])) return 2;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo2 ^ (largeBiomesFlag ? 0x36d326eed40efeb2 : 0x0ef68ec68504005e);
			pxr2.hi = xhi2 ^ (largeBiomesFlag ? 0x5be9ce18223c636a : 0x48b6bf93a2789640); // md5 "octave_-10" or "octave_-8"
			U_initPerlin(&octaves[18], &pxr2);
		}
		*climateSample += 192./511 * samplePerlin(&octaves[18], *px/(256*LB_MULT)*OFF, 0, *pz/(256*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[3]) || (upperBounds && *climateSample > upperBounds[3])) return 3;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo ^ (largeBiomesFlag ? 0x082fe255f8be6631 : 0xf11268128982754f);
			pxr2.hi = xhi ^ (largeBiomesFlag ? 0x4e96119e22dedc81 : 0x257a1d670430b0aa); // md5 "octave_-9" or "octave_-7"
			U_initPerlin(&octaves[10], &pxr2);
		}
		*climateSample += 192./511 * samplePerlin(&octaves[10], *px/(128*LB_MULT), 0, *pz/(128*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[4]) || (upperBounds && *climateSample > upperBounds[4])) return 4;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo2 ^ (largeBiomesFlag ? 0x082fe255f8be6631 : 0xf11268128982754f);
			pxr2.hi = xhi2 ^ (largeBiomesFlag ? 0x4e96119e22dedc81 : 0x257a1d670430b0aa); // md5 "octave_-9" or "octave_-7"
			U_initPerlin(&octaves[19], &pxr2);
		}
		*climateSample += 192./511 * samplePerlin(&octaves[19], *px/(128*LB_MULT)*OFF, 0, *pz/(128*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[5]) || (upperBounds && *climateSample > upperBounds[5])) return 5;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo ^ (largeBiomesFlag ? 0x0ef68ec68504005e : 0xe51c98ce7d1de664);
			pxr2.hi = xhi ^ (largeBiomesFlag ? 0x48b6bf93a2789640 : 0x5f9478a733040c45); // md5 "octave_-8" or "octave_-6"
			U_initPerlin(&octaves[11], &pxr2);
		}
		*climateSample += 96./511 * samplePerlin(&octaves[11], *px/(64*LB_MULT), 0, *pz/(64*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[6]) || (upperBounds && *climateSample > upperBounds[6])) return 6;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo2 ^ (largeBiomesFlag ? 0x0ef68ec68504005e : 0xe51c98ce7d1de664);
			pxr2.hi = xhi2 ^ (largeBiomesFlag ? 0x48b6bf93a2789640 : 0x5f9478a733040c45); // md5 "octave_-8" or "octave_-6"
			U_initPerlin(&octaves[20], &pxr2);
		}
		*climateSample += 96./511 * samplePerlin(&octaves[20], *px/(64*LB_MULT)*OFF, 0, *pz/(64*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[7]) || (upperBounds && *climateSample > upperBounds[7])) return 7;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo ^ (largeBiomesFlag ? 0xf11268128982754f : 0x6d7b49e7e429850a);
			pxr2.hi = xhi ^ (largeBiomesFlag ? 0x257a1d670430b0aa : 0x2e3063c622a24777); // md5 "octave_-7" or "octave_-5"
			U_initPerlin(&octaves[12], &pxr2);
		}
		*climateSample += 48./511 * samplePerlin(&octaves[12], *px/(32*LB_MULT), 0, *pz/(32*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[8]) || (upperBounds && *climateSample > upperBounds[8])) return 8;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo2 ^ (largeBiomesFlag ? 0xf11268128982754f : 0x6d7b49e7e429850a);
			pxr2.hi = xhi2 ^ (largeBiomesFlag ? 0x257a1d670430b0aa : 0x2e3063c622a24777); // md5 "octave_-7" or "octave_-5"
			U_initPerlin(&octaves[21], &pxr2);
		}
		*climateSample += 48./511 * samplePerlin(&octaves[21], *px/(32*LB_MULT)*OFF, 0, *pz/(32*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[9]) || (upperBounds && *climateSample > upperBounds[9])) return 9;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo ^ (largeBiomesFlag ? 0xe51c98ce7d1de664 : 0xbd90d5377ba1b762);
			pxr2.hi = xhi ^ (largeBiomesFlag ? 0x5f9478a733040c45 : 0xc07317d419a7548d); // md5 "octave_-6" or "octave_-4"
			U_initPerlin(&octaves[13], &pxr2);
		}
		*climateSample += 12./511 * samplePerlin(&octaves[13], *px/(16*LB_MULT), 0, *pz/(16*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[10]) || (upperBounds && *climateSample > upperBounds[10])) return 10;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo2 ^ (largeBiomesFlag ? 0xe51c98ce7d1de664 : 0xbd90d5377ba1b762);
			pxr2.hi = xhi2 ^ (largeBiomesFlag ? 0x5f9478a733040c45 : 0xc07317d419a7548d); // md5 "octave_-6" or "octave_-4"
			U_initPerlin(&octaves[22], &pxr2);
		}
		*climateSample += 12./511 * samplePerlin(&octaves[22], *px/(16*LB_MULT)*OFF, 0, *pz/(16*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[11]) || (upperBounds && *climateSample > upperBounds[11])) return 11;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo ^ (largeBiomesFlag ? 0x6d7b49e7e429850a : 0x53d39c6752dac858);
			pxr2.hi = xhi ^ (largeBiomesFlag ? 0x2e3063c622a24777 : 0xbcd1c5a80ab65b3e); // md5 "octave_-5" or "octave_-3"
			U_initPerlin(&octaves[14], &pxr2);
		}
		*climateSample += 6./511 * samplePerlin(&octaves[14], *px/(8*LB_MULT), 0, *pz/(8*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[12]) || (upperBounds && *climateSample > upperBounds[12])) return 12;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo2 ^ (largeBiomesFlag ? 0x6d7b49e7e429850a : 0x53d39c6752dac858);
			pxr2.hi = xhi2 ^ (largeBiomesFlag ? 0x2e3063c622a24777 : 0xbcd1c5a80ab65b3e); // md5 "octave_-5" or "octave_-3"
			U_initPerlin(&octaves[23], &pxr2);
		}
		*climateSample += 6./511 * samplePerlin(&octaves[23], *px/(8*LB_MULT)*OFF, 0, *pz/(8*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[13]) || (upperBounds && *climateSample > upperBounds[13])) return 13;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo ^ (largeBiomesFlag ? 0xbd90d5377ba1b762 : 0xb4a24d7a84e7677b);
			pxr2.hi = xhi ^ (largeBiomesFlag ? 0xc07317d419a7548d : 0x023ff9668e89b5c4); // md5 "octave_-4" or "octave_-2"
			U_initPerlin(&octaves[15], &pxr2);
		}
		*climateSample += 3./511 * samplePerlin(&octaves[15], *px/(4*LB_MULT), 0, *pz/(4*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[14]) || (upperBounds && *climateSample > upperBounds[14])) return 14;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo2 ^ (largeBiomesFlag ? 0xbd90d5377ba1b762 : 0xb4a24d7a84e7677b);
			pxr2.hi = xhi2 ^ (largeBiomesFlag ? 0xc07317d419a7548d : 0x023ff9668e89b5c4); // md5 "octave_-4" or "octave_-2"
			U_initPerlin(&octaves[24], &pxr2);
		}
		*climateSample += 3./511 * samplePerlin(&octaves[24], *px/(4*LB_MULT)*OFF, 0, *pz/(4*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[15]) || (upperBounds && *climateSample > upperBounds[15])) return 15;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo ^ (largeBiomesFlag ? 0x53d39c6752dac858 : 0xdffa22b534c5f608);
			pxr2.hi = xhi ^ (largeBiomesFlag ? 0xbcd1c5a80ab65b3e : 0xb9b67517d3665ca9); // md5 "octave_-3" or "octave_-1"
			U_initPerlin(&octaves[16], &pxr2);
		}
		*climateSample += 3./1022 * samplePerlin(&octaves[16], *px/(2*LB_MULT), 0, *pz/(2*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[16]) || (upperBounds && *climateSample > upperBounds[16])) return 16;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo2 ^ (largeBiomesFlag ? 0x53d39c6752dac858 : 0xdffa22b534c5f608);
			pxr2.hi = xhi2 ^ (largeBiomesFlag ? 0xbcd1c5a80ab65b3e : 0xb9b67517d3665ca9); // md5 "octave_-3" or "octave_-1"
			U_initPerlin(&octaves[25], &pxr2);
		}
		*climateSample += 3./1022 * samplePerlin(&octaves[25], *px/(2*LB_MULT)*OFF, 0, *pz/(2*LB_MULT)*OFF, 0, 0);
		return 17 + ((!lowerBounds || *climateSample > lowerBounds[17]) && (!upperBounds || *climateSample < upperBounds[17]));

		case NP_EROSION:
		if (seedIfInitializingClimate) {
			pxr.lo = xloo ^ (largeBiomesFlag ? 0x8c984b1f8702a951 : 0xd02491e6058f6fd8);
			pxr.hi = xhio ^ (largeBiomesFlag ? 0xead7b1f92bae535f : 0x4792512c94c17a80); // md5 "minecraft:erosion_large" or "minecraft:erosion"
			xlo = xNextLong(&pxr);
			xhi = xNextLong(&pxr);
			pxr2.lo = xlo ^ (largeBiomesFlag ? 0x0fd787bfbc403ec3 : 0x082fe255f8be6631);
			pxr2.hi = xhi ^ (largeBiomesFlag ? 0x74a4a31ca21b48b8 : 0x4e96119e22dedc81); // md5 "octave_-11" or "octave_-9"
			U_initPerlin(&octaves[26], &pxr2);
		}
		*climateSample = 200./279 * samplePerlin(&octaves[26], *px/(512*LB_MULT), 0, *pz/(512*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[0]) || (upperBounds && *climateSample > upperBounds[0])) return 0;

		if (seedIfInitializingClimate) {
			xlo2 = xNextLong(&pxr);
			xhi2 = xNextLong(&pxr);
			pxr2.lo = xlo2 ^ (largeBiomesFlag ? 0x0fd787bfbc403ec3 : 0x082fe255f8be6631);
			pxr2.hi = xhi2 ^ (largeBiomesFlag ? 0x74a4a31ca21b48b8 : 0x4e96119e22dedc81); // md5 "octave_-11" or "octave_-9"
			U_initPerlin(&octaves[30], &pxr2);
		}
		*climateSample += 200./279 * samplePerlin(&octaves[30], *px/(512*LB_MULT)*OFF, 0, *pz/(512*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[1]) || (upperBounds && *climateSample > upperBounds[1])) return 1;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo ^ (largeBiomesFlag ? 0x36d326eed40efeb2 : 0x0ef68ec68504005e);
			pxr2.hi = xhi ^ (largeBiomesFlag ? 0x5be9ce18223c636a : 0x48b6bf93a2789640); // md5 "octave_-10" or "octave_-8"
			U_initPerlin(&octaves[27], &pxr2);
		}
		*climateSample += 100./279 * samplePerlin(&octaves[27], *px/(256*LB_MULT), 0, *pz/(256*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[2]) || (upperBounds && *climateSample > upperBounds[2])) return 2;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo2 ^ (largeBiomesFlag ? 0x36d326eed40efeb2 : 0x0ef68ec68504005e);
			pxr2.hi = xhi2 ^ (largeBiomesFlag ? 0x5be9ce18223c636a : 0x48b6bf93a2789640); // md5 "octave_-10" or "octave_-8"
			U_initPerlin(&octaves[31], &pxr2);
		}
		*climateSample += 100./279 * samplePerlin(&octaves[31], *px/(256*LB_MULT)*OFF, 0, *pz/(256*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[3]) || (upperBounds && *climateSample > upperBounds[3])) return 3;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo ^ (largeBiomesFlag ? 0x0ef68ec68504005e : 0xe51c98ce7d1de664);
			pxr2.hi = xhi ^ (largeBiomesFlag ? 0x48b6bf93a2789640 : 0x5f9478a733040c45); // md5 "octave_-8" or "octave_-6"
			U_initPerlin(&octaves[28], &pxr2);
		}
		*climateSample += 25./279 * samplePerlin(&octaves[28], *px/(64*LB_MULT), 0, *pz/(64*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[4]) || (upperBounds && *climateSample > upperBounds[4])) return 4;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo2 ^ (largeBiomesFlag ? 0x0ef68ec68504005e : 0xe51c98ce7d1de664);
			pxr2.hi = xhi2 ^ (largeBiomesFlag ? 0x48b6bf93a2789640 : 0x5f9478a733040c45); // md5 "octave_-8" or "octave_-6"
			U_initPerlin(&octaves[32], &pxr2);
		}
		*climateSample += 25./279 * samplePerlin(&octaves[32], *px/(64*LB_MULT)*OFF, 0, *pz/(64*LB_MULT)*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[5]) || (upperBounds && *climateSample > upperBounds[5])) return 5;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo ^ (largeBiomesFlag ? 0xf11268128982754f : 0x6d7b49e7e429850a);
			pxr2.hi = xhi ^ (largeBiomesFlag ? 0x257a1d670430b0aa : 0x2e3063c622a24777); // md5 "octave_-7" or "octave_-5"
			U_initPerlin(&octaves[29], &pxr2);
		}
		*climateSample += 25./558 * samplePerlin(&octaves[29], *px/(32*LB_MULT), 0, *pz/(32*LB_MULT), 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[6]) || (upperBounds && *climateSample > upperBounds[6])) return 6;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo2 ^ (largeBiomesFlag ? 0xf11268128982754f : 0x6d7b49e7e429850a);
			pxr2.hi = xhi2 ^ (largeBiomesFlag ? 0x257a1d670430b0aa : 0x2e3063c622a24777); // md5 "octave_-7" or "octave_-5"
			U_initPerlin(&octaves[33], &pxr2);
		}
		*climateSample += 25./558 * samplePerlin(&octaves[33], *px/(32*LB_MULT)*OFF, 0, *pz/(32*LB_MULT)*OFF, 0, 0);
		return 7 + ((!lowerBounds || *climateSample > lowerBounds[7]) && (!upperBounds || *climateSample < upperBounds[7]));

		case NP_SHIFT: *px = floor(*px/4.), *pz = floor(*pz/4.);
		const double origPx = *px, origPz = *pz;
		if (seedIfInitializingClimate) {
			pxr.lo = xloo ^ 0x080518cf6af25384;
			pxr.hi = xhio ^ 0x3f3dfb40a54febd5; // md5 "minecraft:offset"
			xlo = xNextLong(&pxr);
			xhi = xNextLong(&pxr);
			pxr2.lo = xlo ^ 0x53d39c6752dac858;
			pxr2.hi = xhi ^ 0xbcd1c5a80ab65b3e; // md5 "octave_-3"
			U_initPerlin(&octaves[34], &pxr2);
		}
		*px += 8./3 * samplePerlin(&octaves[34], origPx/8., 0, origPz/8., 0., 0.);
		if ((lowerBounds && *px < lowerBounds[0]) || (upperBounds && *px > upperBounds[0])) return 0;
		*pz += 8./3 * samplePerlin(&octaves[34], origPz/8., origPx/8., 0, 0., 0.);
		if ((lowerBounds && *pz < lowerBounds[6]) || (upperBounds && *pz > upperBounds[6])) return 1;

		if (seedIfInitializingClimate) {
			xlo2 = xNextLong(&pxr);
			xhi2 = xNextLong(&pxr);
			pxr2.lo = xlo2 ^ 0x53d39c6752dac858;
			pxr2.hi = xhi2 ^ 0xbcd1c5a80ab65b3e; // md5 "octave_-3"
			U_initPerlin(&octaves[37], &pxr2);
		}
		*px += 8./3 * samplePerlin(&octaves[37], origPx/8.*OFF, 0, origPz/8.*OFF, 0, 0);
		if ((lowerBounds && *px < lowerBounds[1]) || (upperBounds && *px > upperBounds[1])) return 2;
		*pz += 8./3 * samplePerlin(&octaves[37], origPz/8.*OFF, origPx/8.*OFF, 0, 0, 0);
		if ((lowerBounds && *pz < lowerBounds[7]) || (upperBounds && *pz > upperBounds[7])) return 3;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo ^ 0xb4a24d7a84e7677b;
			pxr2.hi = xhi ^ 0x023ff9668e89b5c4; // md5 "octave_-2"
			U_initPerlin(&octaves[35], &pxr2);
		}
		*px += 4./3 * samplePerlin(&octaves[35], origPx/4., 0, origPz/4., 0, 0);
		if ((lowerBounds && *px < lowerBounds[2]) || (upperBounds && *px > upperBounds[2])) return 4;
		*pz += 4./3 * samplePerlin(&octaves[35], origPz/4., origPx/4., 0, 0, 0);
		if ((lowerBounds && *pz < lowerBounds[8]) || (upperBounds && *pz > upperBounds[8])) return 5;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo2 ^ 0xb4a24d7a84e7677b;
			pxr2.hi = xhi2 ^ 0x023ff9668e89b5c4; // md5 "octave_-2"
			U_initPerlin(&octaves[38], &pxr2);
		}
		*px += 4./3 * samplePerlin(&octaves[38], origPx/4.*OFF, 0, origPz/4.*OFF, 0, 0);
		if ((lowerBounds && *px < lowerBounds[3]) || (upperBounds && *px > upperBounds[3])) return 6;
		*pz += 4./3 * samplePerlin(&octaves[38], origPz/4.*OFF, origPx/4.*OFF, 0, 0, 0);
		if ((lowerBounds && *pz < lowerBounds[9]) || (upperBounds && *pz > upperBounds[9])) return 7;
		
		if (seedIfInitializingClimate) {
			pxr2.lo = xlo ^ 0xdffa22b534c5f608;
			pxr2.hi = xhi ^ 0xb9b67517d3665ca9; // md5 "octave_-1"
			U_initPerlin(&octaves[36], &pxr2);
		}
		*px += 2./3 * samplePerlin(&octaves[36], origPx/2., 0, origPz/2., 0, 0);
		if ((lowerBounds && *px < lowerBounds[4]) || (upperBounds && *px > upperBounds[4])) return 8;
		*pz += 2./3 * samplePerlin(&octaves[36], origPz/2., origPx/2., 0, 0, 0);
		if ((lowerBounds && *pz < lowerBounds[10]) || (upperBounds && *pz > upperBounds[10])) return 9;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo2 ^ 0xdffa22b534c5f608;
			pxr2.hi = xhi2 ^ 0xb9b67517d3665ca9; // md5 "octave_-1"
			U_initPerlin(&octaves[39], &pxr2);
		}
		*px += 2./3 * samplePerlin(&octaves[39], origPx/2.*OFF, 0, origPz/2.*OFF, 0, 0);
		if ((lowerBounds && *px < lowerBounds[5]) || (upperBounds && *px > upperBounds[5])) return 10;
		*pz += 2./3 * samplePerlin(&octaves[39], origPz/2.*OFF, origPx/2.*OFF, 0, 0, 0);
		return 11 + ((!lowerBounds || *pz > lowerBounds[11]) && (!upperBounds || *pz < upperBounds[11]));

		case NP_WEIRDNESS:
		if (seedIfInitializingClimate) {
			pxr.lo = xloo ^ 0xefc8ef4d36102b34;
			pxr.hi = xhio ^ 0x1beeeb324a0f24ea; // md5 "minecraft:ridge"
			xlo = xNextLong(&pxr);
			xhi = xNextLong(&pxr);
			pxr2.lo = xlo ^ 0xf11268128982754f;
			pxr2.hi = xhi ^ 0x257a1d670430b0aa; // md5 "octave_-7"
			U_initPerlin(&octaves[40], &pxr2);
		}
		*climateSample = 40./63 * samplePerlin(&octaves[40], *px/128., 0, *pz/128., 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[0]) || (upperBounds && *climateSample > upperBounds[0])) return 0;

		if (seedIfInitializingClimate) {
			xlo2 = xNextLong(&pxr);
			xhi2 = xNextLong(&pxr);
			pxr2.lo = xlo2 ^ 0xf11268128982754f;
			pxr2.hi = xhi2 ^ 0x257a1d670430b0aa; // md5 "octave_-7"
			U_initPerlin(&octaves[43], &pxr2);
		}
		*climateSample += 40./63 * samplePerlin(&octaves[43], *px/128.*OFF, 0, *pz/128.*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[1]) || (upperBounds && *climateSample > upperBounds[1])) return 1;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo ^ 0xe51c98ce7d1de664;
			pxr2.hi = xhi ^ 0x5f9478a733040c45; // md5 "octave_-6"
			U_initPerlin(&octaves[41], &pxr2);
		}
		*climateSample += 40./63 * samplePerlin(&octaves[41], *px/64., 0, *pz/64., 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[2]) || (upperBounds && *climateSample > upperBounds[2])) return 2;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo2 ^ 0xe51c98ce7d1de664;
			pxr2.hi = xhi2 ^ 0x5f9478a733040c45; // md5 "octave_-6"
			U_initPerlin(&octaves[44], &pxr2);
		}
		*climateSample += 40./63 * samplePerlin(&octaves[44], *px/64.*OFF, 0, *pz/64.*OFF, 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[3]) || (upperBounds && *climateSample > upperBounds[3])) return 3;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo ^ 0x6d7b49e7e429850a;
			pxr2.hi = xhi ^ 0x2e3063c622a24777; // md5 "octave_-5"
			U_initPerlin(&octaves[42], &pxr2);
		}
		*climateSample += 10./63 * samplePerlin(&octaves[42], *px/32., 0, *pz/32., 0, 0);
		if ((lowerBounds && *climateSample < lowerBounds[4]) || (upperBounds && *climateSample > upperBounds[4])) return 4;

		if (seedIfInitializingClimate) {
			pxr2.lo = xlo2 ^ 0x6d7b49e7e429850a;
			pxr2.hi = xhi2 ^ 0x2e3063c622a24777; // md5 "octave_-5"
			U_initPerlin(&octaves[45], &pxr2);
		}
		*climateSample += 10./63 * samplePerlin(&octaves[45], *px/32.*OFF, 0, *pz/32.*OFF, 0, 0);
		return 5 + ((!lowerBounds || *climateSample > lowerBounds[5]) && (!upperBounds || *climateSample < upperBounds[5]));
	}
	return 0;
}