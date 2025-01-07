#include "Spawn.h"

// The maximum fitness that each climate by itself can equate to.
const double U_MAX_CLIMATE_FITNESSES[NP_MAX] = {169783107.122, 56953151.7314, 1521389678.6949422, 227197622.444, 0., 384556370.747};

double U_getFitness(const Pos *coord, const double *temperature, const double *humidity, const double *continentalness, const double *erosion, const double *weirdness, const bool post1_21_1) {
	double offset, fitness = 0.;
	// U_getFitnessBounded(coord, temperature, humidity, continentalness, erosion, weirdness, NULL, &fitness);
	// Distance
	if (coord) {
		uint64_t squaredEuclid = (uint64_t)(coord->x) * coord->x + (uint64_t)(coord->z) * coord->z;
		fitness = post1_21_1 ? squaredEuclid : (squaredEuclid*squaredEuclid)/390625.;
	}

	const double MULT = post1_21_1 ? 419430400000000. : 100000000.;
	// Continentalness
	if (continentalness) {
		offset = *continentalness < -0.11 ? *continentalness + 0.11 : *continentalness > 1. ? *continentalness - 1. : 0.;
		fitness += MULT*offset*offset;
	}

	// Weirdness
	if (weirdness) {
		double absWeirdness = fabs(*weirdness);
		offset = absWeirdness < 0.16 ? 0.16 - absWeirdness : absWeirdness > 1 ? absWeirdness - 1. : 0.;
		fitness += MULT*offset*offset;
	}

	// Erosion
	if (erosion) {
		double absErosion = fabs(*erosion);
		offset = absErosion > 1. ? absErosion - 1. : 0.;
		fitness += MULT*offset*offset;
	}

	// Temperature
	if (temperature) {
		double absTemperature = fabs(*temperature);
		offset = absTemperature > 1. ? absTemperature - 1. : 0.;
		fitness += MULT*offset*offset;
	}

	// Humidity
	if (humidity) {
		double absHumidity = fabs(*humidity);
		offset = absHumidity > 1. ? absHumidity - 1. : 0.;
		fitness += MULT*offset*offset;
	}

	return fitness;
}

bool U_getFitnessBounded(const Pos *coord, const double *temperature, const double *humidity, const double *continentalness, const double *erosion, const double *weirdness, const double *upperBound, const bool post1_21_1, double *fitness) {
	if (!fitness) return false;
	// Distance
	if (coord) {
		uint64_t squaredEuclid = (uint64_t)(coord->x) * coord->x + (uint64_t)(coord->z) * coord->z;
		*fitness = post1_21_1 ? squaredEuclid : (squaredEuclid*squaredEuclid)/390625.;
	}
	if (upperBound && *fitness >= *upperBound) return false;

	const double MULT = post1_21_1 ? 419430400000000. : 100000000.;
	double offset;
	// Continentalness
	if (continentalness) {
		offset = *continentalness < -0.11 ? *continentalness + 0.11 : *continentalness > 1. ? *continentalness - 1. : 0.;
		*fitness += MULT*offset*offset;
	}
	if (upperBound && *fitness >= *upperBound) return false;

	// Weirdness
	if (weirdness) {
		double absWeirdness = fabs(*weirdness);
		offset = absWeirdness < 0.16 ? 0.16 - absWeirdness : absWeirdness > 1. ? absWeirdness - 1. : 0.;
		*fitness += MULT*offset*offset;
	}
	if (upperBound && *fitness >= *upperBound) return false;

	// Erosion
	if (erosion) {
		double absErosion = fabs(*erosion);
		offset = absErosion > 1. ? absErosion - 1. : 0.;
		*fitness += MULT*offset*offset;
	}
	if (upperBound && *fitness >= *upperBound) return false;

	// Temperature
	if (temperature) {
		double absTemperature = fabs(*temperature);
		offset = absTemperature > 1. ? absTemperature - 1. : 0.;
		*fitness += MULT*offset*offset;
	}
	if (upperBound && *fitness >= *upperBound) return false;

	// Humidity
	if (humidity) {
		double absHumidity = fabs(*humidity);
		offset = absHumidity > 1. ? absHumidity - 1. : 0.;
		*fitness += MULT*offset*offset;
	}
	return !upperBound || *fitness < *upperBound;
}

double U_sampleAndGetFitness(const Pos *coord, const PerlinNoise *oct, const bool post1_21_1, const bool largeBiomesFlag) {
	// U_sampleAndGetFitnessBounded(coord, oct, NULL, largeBiomesFlag, &fitness);
	if (!coord) return -INFINITY;
	// Distance
	uint64_t squaredEuclid = (uint64_t)(coord->x) * coord->x + (uint64_t)(coord->z) * coord->z;
	double fitness = post1_21_1 ? squaredEuclid : (squaredEuclid*squaredEuclid)/390625.;
	const double MULT = post1_21_1 ? 419430400000000. : 100000000.;

	double px = coord->x, pz = coord->z;
	U_sampleClimate(NP_SHIFT, oct, &px, &pz, largeBiomesFlag);
	// Continentalness
	double sample = U_sampleClimate(NP_CONTINENTALNESS, oct, &px, &pz, largeBiomesFlag);
	double offset = sample < -0.11 ? sample + 0.11 : sample > 1. ? sample - 1. : 0.;
	fitness += MULT*offset*offset;

	// Weirdness
	sample = fabs(U_sampleClimate(NP_WEIRDNESS, oct, &px, &pz, largeBiomesFlag));
	offset = sample < 0.16 ? 0.16 - sample : sample > 1. ? sample - 1. : 0.;
	fitness += MULT*offset*offset;

	// Erosion
	sample = fabs(U_sampleClimate(NP_EROSION, oct, &px, &pz, largeBiomesFlag));
	offset = sample > 1. ? sample - 1. : 0.;
	fitness += MULT*offset*offset;

	// Temperature
	sample = fabs(U_sampleClimate(NP_TEMPERATURE, oct, &px, &pz, largeBiomesFlag));
	offset = sample > 1. ? sample - 1. : 0.;
	fitness += MULT*offset*offset;

	// Humidity
	sample = fabs(U_sampleClimate(NP_HUMIDITY, oct, &px, &pz, largeBiomesFlag));
	offset = sample > 1. ? sample - 1. : 0.;
	fitness += MULT*offset*offset;

	return fitness;
}

bool U_sampleAndGetFitnessBounded(const Pos *coord, const PerlinNoise *oct, const double *upperBound, const bool post1_21_1, const bool largeBiomesFlag, double *fitness) {
	if (!coord || !fitness) return false;
	// Distance
	uint64_t squaredEuclid = (uint64_t)(coord->x) * coord->x + (uint64_t)(coord->z) * coord->z;
	*fitness = post1_21_1 ? squaredEuclid : (squaredEuclid*squaredEuclid)/390625.;
	if (upperBound && *fitness >= *upperBound) return false;

	const double MULT = post1_21_1 ? 419430400000000. : 100000000.;
	double px = coord->x, pz = coord->z;
	U_sampleClimate(NP_SHIFT, oct, &px, &pz, largeBiomesFlag);

	// Continentalness
	double sample = U_sampleClimate(NP_CONTINENTALNESS, oct, &px, &pz, largeBiomesFlag);
	double offset = sample < -0.11 ? sample + 0.11 : sample > 1. ? sample - 1. : 0.;
	*fitness += MULT*offset*offset;
	if (upperBound && *fitness >= *upperBound) return false;

	// Weirdness
	sample = fabs(U_sampleClimate(NP_WEIRDNESS, oct, &px, &pz, largeBiomesFlag));
	offset = sample < 0.16 ? 0.16 - sample : sample > 1. ? sample - 1. : 0.;
	*fitness += MULT*offset*offset;
	if (upperBound && *fitness >= *upperBound) return false;

	// Erosion
	sample = fabs(U_sampleClimate(NP_EROSION, oct, &px, &pz, largeBiomesFlag));
	offset = sample > 1. ? sample - 1. : 0.;
	*fitness += MULT*offset*offset;
	if (upperBound && *fitness >= *upperBound) return false;

	// Temperature
	sample = fabs(U_sampleClimate(NP_TEMPERATURE, oct, &px, &pz, largeBiomesFlag));
	offset = sample > 1. ? sample - 1. : 0.;
	*fitness += MULT*offset*offset;
	if (upperBound && *fitness >= *upperBound) return false;

	// Humidity
	sample = fabs(U_sampleClimate(NP_HUMIDITY, oct, &px, &pz, largeBiomesFlag));
	offset = sample > 1. ? sample - 1. : 0.;
	*fitness += MULT*offset*offset;
	return !upperBound || *fitness < *upperBound;
}

double U_getEffectiveDistance(const double fitness, const bool post1_21_1) {
	if (fitness < 0.) return -INFINITY;
	return post1_21_1 ? sqrt(fitness) : 25.*pow(fitness, 1./4);
}

double U_getEffectiveTemperature(const double fitness, const bool post1_21_1) {
	if (fitness < 0.) return -INFINITY;
	return 1. + sqrt(fitness)/(post1_21_1 ? 20480000. : 10000.);
}

double U_getEffectiveHumidity(const double fitness, const bool post1_21_1) {
	// return U_getEffectiveTemperature(fitness);
	if (fitness < 0.) return -INFINITY;
	return 1. + sqrt(fitness)/(post1_21_1 ? 20480000. : 10000.);
}

double U_getEffectiveContinentalness(const double fitness, const bool post1_21_1) {
	if (fitness < 0) return -INFINITY;
	return -0.11 - sqrt(fitness)/(post1_21_1 ? 20480000. : 10000.);
}

double U_getEffectiveErosion(const double fitness, const bool post1_21_1) {
	// return U_getEffectiveTemperature(fitness);
	if (fitness < 0.) return -INFINITY;
	return 1. + sqrt(fitness)/(post1_21_1 ? 20480000. : 10000.);
}

double U_getEffectiveWeirdnessOuter(const double fitness, const bool post1_21_1) {
	// return U_getEffectiveTemperature(fitness);
	if (fitness < 0.) return -INFINITY;
	return 1. + sqrt(fitness)/(post1_21_1 ? 20480000. : 10000.);
}

double U_getEffectiveWeirdnessInner(const double fitness, const bool post1_21_1) {
	if (fitness < 0.) return -INFINITY;
	if (fitness > (post1_21_1 ? 67108864000000. : 2560000.)) return INFINITY;
	return 0.16 - sqrt(fitness)/(post1_21_1 ? 20480000. : 10000.);
}

bool U_singleStageSpawnBounded(PerlinNoise *oct, const SpawnResult *firstStageChosenResult, const double fitnessLowerBound, const bool post1_21_1, const bool largeBiomesFlag, SpawnResult *chosenResult, double *chosenFitness) {
	#ifndef _SPAWN_TABLES_ARE_PRESENT
		double fitness, bestFitness = INFINITY;
		const double maxRad = 512.*(1 + 3*!firstStageChosenResult);
		const double radInc = 32.*(1 + 3*!firstStageChosenResult);
		for (double rad = 0.; rad <= maxRad; rad += radInc) {
			for (double ang = 0.; ang <= U_TWO_PI; ang += rad ? radInc/rad : INFINITY) {
				Pos pos = {sin(ang) * rad, cos(ang) * rad};
				if (!U_sampleAndGetFitnessBounded(&pos, oct, &bestFitness, post1_21_1, largeBiomesFlag, &fitness)) continue;
				if (chosenResult) {
					chosenResult->pos.x = pos.x;
					chosenResult->pos.z = pos.z;
				}
				bestFitness = fitness;
				if (*chosenFitness) *chosenFitness = fitness;
				if (bestFitness < fitnessLowerBound) return false;
			}
		}
		return true;
	#else
	static const int *TABLES[] = {U_SPAWN_FIRST_STAGE_VALS, U_SPAWN_SECOND_STAGE_VALS};
	const int **CHOSEN_TABLE = firstStageChosenResult ? U_SPAWN_SECOND_STAGE_VALS[firstStageChosenResult->index] : U_SPAWN_FIRST_STAGE_VALS;
	double bestFitness = INFINITY;
	// TODO: Continue as soon as an individual samplePerlin pushes fitness over fitness?
	for (size_t i = 0; i < sizeof(CHOSEN_TABLE)/sizeof(*CHOSEN_TABLE); ++i) {
		double fitness = CHOSEN_TABLE[i][U_spawn_table_fitness];
		if (fitness >= bestFitness) continue;
		Pos pos = {CHOSEN_TABLE[i][U_spawn_table_x], CHOSEN_TABLE[i][U_spawn_table_z]};
		if (!U_sampleAndGetFitnessBounded(&pos, oct, &bestFitness, post1_21_1, largeBiomesFlag, &fitness)) continue;
		if (chosenResult) chosenResult = i;
		if (*chosenFitness) *chosenFitness = fitness;
		if (bestFitness < fitnessLowerBound) return false;
	}
	return true;
	#endif
}

bool U_firstStageSpawnBounded(PerlinNoise *oct, const double fitnessLowerBound, const bool post1_21_1, const bool largeBiomesFlag, SpawnResult *chosenResult, double *chosenFitness) {
	#ifndef _SPAWN_TABLES_ARE_PRESENT
		// Fallback in case tables were not linked
		double fitness, bestFitness = INFINITY;
		for (double rad = 0.; rad <= 2048.; rad += 512.) {
			for (double ang = 0.; ang <= U_TWO_PI; ang += rad ? 512./rad : INFINITY) {
				Pos pos = {sin(ang) * rad, cos(ang) * rad};
				if (!U_sampleAndGetFitnessBounded(&pos, oct, &bestFitness, post1_21_1, largeBiomesFlag, &fitness)) continue;
				if (chosenResult) {
					chosenResult->pos.x = pos.x;
					chosenResult->pos.z = pos.z;
				}
				bestFitness = fitness;
				if (*chosenFitness) *chosenFitness = fitness;
				if (bestFitness < fitnessLowerBound) return false;
			}
		}
		return true;
	#else
		double bestFitness = INFINITY;
		// TODO: Continue as soon as an individual samplePerlin pushes fitness over fitness?
		for (size_t i = 0; i < sizeof(U_SPAWN_FIRST_STAGE_VALS)/sizeof(*U_SPAWN_FIRST_STAGE_VALS); ++i) {
			double fitness = U_SPAWN_FIRST_STAGE_VALS[i][U_spawn_table_fitness];
			if (fitness >= bestFitness) continue;
			Pos pos = {U_SPAWN_FIRST_STAGE_VALS[i][U_spawn_table_x], U_SPAWN_FIRST_STAGE_VALS[i][U_spawn_table_z]};
			if (!U_sampleAndGetFitnessBounded(&pos, oct, &bestFitness, largeBiomesFlag, &fitness)) continue;
			if (chosenResult) chosenResult = i;
			if (*chosenFitness) *chosenFitness = fitness;
			bestFitness = fitness;
			if (bestFitness < fitnessLowerBound) return false;
		}
		return true;
	#endif
}

bool U_secondStageSpawnBounded(PerlinNoise *oct, const SpawnResult *firstStageChosenResult, const double firstStageChosenFitness, const double fitnessLowerBound, const bool post1_21_1, const bool largeBiomesFlag, SpawnResult *chosenResult, double *chosenFitness) {
	#ifndef _SPAWN_TABLES_ARE_PRESENT
		// Fallback in case tables were not linked
		double fitness, bestFitness = firstStageChosenFitness;
		for (double rad = 32.; rad <= 512.; rad += 32.) {
			for (double ang = 0.; ang <= U_TWO_PI; ang += 32./rad) {
				Pos pos = {firstStageChosenResult->pos.x + (int)(sin(ang) * rad), firstStageChosenResult->pos.z + (int)(cos(ang) * rad)};
				if (!U_sampleAndGetFitnessBounded(&pos, oct, &bestFitness, post1_21_1, largeBiomesFlag, &fitness)) continue;
				if (chosenResult) {
					chosenResult->pos.x = pos.x;
					chosenResult->pos.z = pos.z;
				}
				bestFitness = fitness;
				if (*chosenFitness) *chosenFitness = fitness;
				if (bestFitness < fitnessLowerBound) return false;
			}
		}
		return true;
	#else
		double bestFitness = firstStageChosenFitness;
		// TODO: Continue as soon as an individual samplePerlin pushes fitness over fitness?
		for (size_t i = 0; i < sizeof(U_SPAWN_SECOND_STAGE_VALS[firstStageChosenResult->index])/sizeof(*U_SPAWN_SECOND_STAGE_VALS[firstStageChosenResult->index]); ++i) {
			double fitness = U_SPAWN_SECOND_STAGE_VALS[firstStageChosenResult->index][i][U_spawn_table_fitness];
			if (fitness >= bestFitness) continue;
			Pos pos = {U_SPAWN_SECOND_STAGE_VALS[firstStageChosenResult->index][i][U_spawn_table_x], U_SPAWN_SECOND_STAGE_VALS[firstStageChosenResult->index][i][U_spawn_table_z]};
			if (!U_sampleAndGetFitnessBounded(&pos, oct, &bestFitness, post1_21_1, largeBiomesFlag, &fitness)) continue;
			if (chosenResult) chosenResult = i;
			if (*chosenFitness) *chosenFitness = fitness;
			bestFitness = fitness;
			if (bestFitness < fitnessLowerBound) return false;
		}
		return true;
	#endif
}