#ifndef _SPAWN_H
#define _SPAWN_H

#include "Climates.h"

enum {U_spawn_table_x, U_spawn_table_z, U_spawn_table_fitness, U_spawn_table_entrySize};
extern const int U_SPAWN_FIRST_STAGE_VALS[66][U_spawn_table_entrySize];
extern const double U_MAX_CLIMATE_FITNESSES[NP_MAX], U_SPAWN_SECOND_STAGE_VALS[66][862][U_spawn_table_entrySize];

/* Since the second-stage spawn lookup table is so drastically large (170000+ elements), the spawn functions are configured to return their outputs as
   - indices of the spawn tables, if they were linked; or
   - raw coordinates, if they were not.
   WARNING: 90% sure Spawn Lookup Tables.c needs to be linked before Spawn.h/Spawn.c when compiling*/
typedef union SpawnResult SpawnResult;
union SpawnResult {
	int index;
	Pos pos;
};

#ifdef __cplusplus
extern "C" {
#endif

// #define _SPAWN_TABLES_ARE_PRESENT

/* To standardize fetching values from the SpawnResults, this helper function returns
   - `true`, if the index attribute should be invoked, or
   - `false`, if the pos attribute should be invoked.*/
static inline bool areSpawnResultsIndices() {
	#ifdef _SPAWN_TABLES_ARE_PRESENT
		return true;
	#else
		return false;
	#endif
}

// Returns the fitness of a point, given its coordinate and its climate values.
// All parameters are nullable if desired.
double U_getFitness(const Pos *coord, const double *temperature, const double *humidity, const double *continentalness, const double *erosion, const double *weirdness, const bool post1_21_1);
// Calculates the fitness of a point, aborting early (and returning false) if the fitness ever rises above `*upperBound`.
// Otherwise returns true and stores the fitness value in `*fitness`.
// All parameters except `*fitness` are nullable if desired.
bool U_getFitnessBounded(const Pos *coord, const double *temperature, const double *humidity, const double *continentalness, const double *erosion, const double *weirdness, const double *upperBound, const bool post1_21_1, double *fitness);

// Samples each climate to return the fitness of a point.
// All climates must have been initialized beforehand (with e.g. `U_initAllClimates()`.)
double U_sampleAndGetFitness(const Pos *coord, const PerlinNoise *oct, const bool post1_21_1, const bool largeBiomesFlag);

// Samples each climate to calculate the fitness of a point, aborting early (and returning false) if the fitness ever rises above `*upperBound`.
// Otherwise returns true and stores the fitness value in `*fitness`.
// All climates must have been initialized beforehand (with e.g. `U_initAllClimates()`.)
bool U_sampleAndGetFitnessBounded(const Pos *coord, const PerlinNoise *oct, const double *upperBound, const bool post1_21_1, const bool largeBiomesFlag, double *fitness);

// Returns the Euclidean distance that would equate to the provided fitness value.
double U_getEffectiveDistance(const double fitness, const bool post1_21_1);

// Returns the positive temperature that would equate to the provided fitness value.
// The equivalent negative temperature can be obtained by negating the output.
double U_getEffectiveTemperature(const double fitness, const bool post1_21_1);

// Returns the positive humidity that would equate to the provided fitness value.
// The equivalent negative humidity can be obtained by negating the output.
double U_getEffectiveHumidity(const double fitness, const bool post1_21_1);

// Returns the **negative** continentalness that would equate to the provided fitness value.
// The equivalent positive continentalness would be the same value as `U_getEffectiveTemperature(fitness)`.
double U_getEffectiveContinentalness(const double fitness, const bool post1_21_1);

// Returns the positive erosion that would equate to the provided fitness value.
// The equivalent negative erosion can be obtained by negating the output.
double U_getEffectiveErosion(const double fitness, const bool post1_21_1);

// Returns the positive weirdness > 1 that would equate to the provided fitness value.
// The equivalent negative weirdness < -1 can be obtained by negating the output.
double U_getEffectiveWeirdnessOuter(const double fitness, const bool post1_21_1);

// Returns the positive weirdness < 0.16 that would equate to the provided fitness value, or `INFINITY` if none exists.
// The equivalent negative weirdness > -0.16 can be obtained by negating the output.
double U_getEffectiveWeirdnessInner(const double fitness, const bool post1_21_1);

// Emulates either one of the first two stages of the spawn algorithm (first stage if `firstStageChosenResult` is NULL, second stage if it is provided).
// Aborts early and returns false if the current chosen point's fitness ever drops below `fitnessLowerBound`.
// Otherwise returns true, and stores the chosen result in `*chosenResult` and the final fitness value in `*chosenFitness`.
bool U_singleStageSpawnBounded(PerlinNoise *oct, const SpawnResult *firstStageChosenResult, const double fitnessLowerBound, const bool post1_21_1, const bool largeBiomesFlag, SpawnResult *chosenResult, double *chosenFitness);

// Emulates the first stage of the spawn algorithm, aborting early (and returning false) if the current chosen point's fitness ever drops below `fitnessLowerBound`.
// Otherwise returns true, and stores the chosen result in `*chosenResult` and the final fitness value in `*chosenFitness`.
bool U_firstStageSpawnBounded(PerlinNoise *oct, const double fitnessLowerBound, const bool post1_21_1, const bool largeBiomesFlag, SpawnResult *chosenResult, double *chosenFitness);

// Emulates the second stage of the spawn algorithm, aborting early (and returning false) if the current chosen point's fitness ever drops below `fitnessLowerBound`.
// Otherwise returns true and stores the chosen coordinate index fitness value in `*fitness`.
bool U_secondStageSpawnBounded(PerlinNoise *oct, const SpawnResult *firstStageChosenResult, const double firstStageChosenFitness, const double fitnessLowerBound, const bool post1_21_1, const bool largeBiomesFlag, SpawnResult *chosenResult, double *chosenFitness);

#ifdef __cplusplus
}
#endif

#endif