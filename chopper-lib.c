//
// Created by Gregory Tucker, ESS ERIC on 2023-06-01.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef CHOPPER_LIB_CHOPPER_LIB_H
#include "chopper-lib.h"
#endif

/******************************** range functions ************************************/
void range_sort(range a){
  if (a.maximum < a.minimum){
    double tmp = a.minimum;
    a.minimum = a.maximum;
    a.maximum = tmp;
  }
}
int classify_range_overlap(const range * a, const range * b){
  //  |A    A|   |B   B| ... or ... |B   B|   |A   A|
  if (a->maximum < b->minimum || b->maximum < a->minimum) return 0;
  // there should be *some* overlap at this point:
  if (a->minimum == b->minimum && a->maximum == b->maximum) return 1; // identical
  // |2| -> B on right end; |3| -> A on right end; (+) -> one inside the other, (-) -> overlapping subregion
  if (a->minimum <= b->minimum && a->maximum >= b->maximum) return 3; // |A|BB|A|
  if (b->minimum <= a->minimum && b->maximum >= a->maximum) return 2; // |B|AA|B|
  if (a->minimum < b->minimum && a->maximum < b->maximum) return -2; // |A|BA|B|
  if (b->minimum < a->minimum && b->maximum < a->maximum) return -3; // |B|AB|A|
  return 0;
}
int compare_ranges(const range * a, const range * b){
  if (a->minimum > b->minimum) return 1;
  if (a->minimum < b->minimum) return -1;
  return 0;
}
// This gateway function is used along with qsort, which handles only void pointers
int compare_sorted_ranges(const void * ptr_a, const void * ptr_b){
  return compare_ranges((range *) ptr_a, (range *) ptr_b);
}
/******************************** range_set functions ************************************/
range_set range_set_sort(range_set s){
  // sort all sub-ranges:
  for (unsigned i=0; i<s.count; ++i) range_sort(s.ranges[i]);
  // sort the sub-ranges by minimum
  qsort(s.ranges, s.count, sizeof(*(s.ranges)), compare_sorted_ranges);
  // combine overlapping ranges:
  unsigned overlapping = 0;
  for (unsigned i=1; i<s.count; ++i) if (s.ranges[i-1].maximum >= s.ranges[i].minimum) ++overlapping;
  if (overlapping){
    range_set new_s;
    new_s.count = s.count - overlapping;
    new_s.ranges = calloc(new_s.count, sizeof(range));
    // copy the first element
    new_s.ranges[0].minimum = s.ranges[0].minimum;
    new_s.ranges[0].maximum = s.ranges[0].maximum;

    unsigned copied = 1;
    for (unsigned i=1; i<s.count; ++i) if (s.ranges[i-1].maximum >= s.ranges[i].minimum) {
        // combine the lower bound of the end of the new ranges and the i_th range upper bound over the last new range:
        new_s.ranges[copied-1].maximum = s.ranges[i].maximum;
      } else {
        new_s.ranges[copied].minimum = s.ranges[i].minimum;
        new_s.ranges[copied++].maximum = s.ranges[i].maximum;
      }
    if (copied != s.count - overlapping) printf("Expected to copy %u but copied %u ranges!\n", s.count - overlapping, copied);
    // // free the now-old range_set before we lose its handle ... this is dangerous
    // if (s.ranges) free(s.ranges);
    // recursively re-sort in case we missed overlapping ranges
    return range_set_sort(new_s);
  } else {
    return s;
  }
}

range_set range_intersection(range_set ain, range_set bin){
  range_set a = range_set_sort(ain);
  range_set b = range_set_sort(bin);
  range_set out;
  out.count = 0;
  unsigned i=0, j=0;
  while (i < a.count && j < b.count){
    switch (abs(classify_range_overlap(a.ranges + i, b.ranges + j))) {
      case 3: out.count++; ++j; break; // -3: |B|AB|A|, 3: |A|BB|A|, both increment B
      case 2: out.count++; ++i; break; // -2: |A|BA|B|, 2: |B|AA|B|, both increment A
      case 1: out.count++; ++i; ++j; break; // identical, increment both
      default: {
        // no overlap, so increment the range with the smaller minimum edge
        int comp = compare_ranges(a.ranges + i, b.ranges + j);
        if (comp < 0) ++i;
        if (comp > 0) ++j;
        if (comp == 0) {
          printf("This should not be possible");
          exit(-1);
        }
      }
    }
  }
  // we now know *how many* output sub-ranges there will be, and can allocate the output structure
  out.ranges = calloc(out.count, sizeof(range));
  // and can go through again actually assigning outputs:
  i = 0;
  j = 0;
  unsigned k=0;
  range tmp;
  while (i < a.count && j < b.count){
    switch (classify_range_overlap(a.ranges + i, b.ranges + j)) {
      case -3: { // |B|AB|A|; keep (A[i]min, B[j]max) increment j since A extends to higher value
        out.ranges[k].minimum = a.ranges[i].minimum;
        out.ranges[k].maximum = b.ranges[j].maximum;
        ++j; ++k;
        break;
      }
      case -2: { //|A|BA|B|; keep (B[j]min, A[i]max) increment i since B extends to higher value
        out.ranges[k].minimum = b.ranges[j].minimum;
        out.ranges[k].maximum = a.ranges[i].maximum;
        ++i; ++k;
        break;
      }
      case 1: { // identical, keep either a[i] or b[j] and increment both i & j
        out.ranges[k].minimum = a.ranges[i].minimum;
        out.ranges[k].maximum = a.ranges[i].maximum;
        ++i; ++j; ++k;
        break;
      }
      case 2: { // |B|AA|B|; keep a[i] and increment i
        out.ranges[k].minimum = a.ranges[i].minimum;
        out.ranges[k].maximum = a.ranges[i].maximum;
        ++i; ++k;
        break;
      }
      case 3: { // |A|BB|A|; keep b[j] and increment j
        out.ranges[k].minimum = b.ranges[j].minimum;
        out.ranges[k].maximum = b.ranges[j].maximum;
        ++j; ++k;
        break;
      }
      default: {
        if (compare_ranges(a.ranges + i, b.ranges + j) < 0) ++i; else ++j;
      }
    }
  }
  // memory management: (did the sort function make a new range_set?)
  if (a.ranges != NULL && a.ranges != ain.ranges) free(a.ranges);
  if (b.ranges != NULL && b.ranges != bin.ranges) free(b.ranges);
  return out;
}

/********************** local helper functions *************************/
void print_help(const char * program_name) {
  printf("Usage: %s [parameter name]=value ...\n", program_name);
  printf("\tValid parameter names: xxNtype for xx in (ps, fo, bw), N in (1, 2), and type in (speed, phase)\n");
  exit(0);
}

void print_one(const char * name, double value) {
  printf("  %s = % 10.4f\n", name, value);
}
/*
 * toff=fabs(t-atan2(x,yprime)/omega - delay - (jitter ? jitter*randnorm():0));
   // does neutron hit outside slit? *
   if (fmod(toff+To/2.0,Tg)>To) ABSORB;
 */

/****************** chopper train functionality ************************/
range_set chopper_inverse_velocity_windows(unsigned count, const chopper_parameters * choppers,
                                           double inv_v_min, double inv_v_max, double latest_emission){
  range_set limits;
  limits.count = 1;
  limits.ranges = calloc(1, sizeof(range));
  limits.ranges[0].minimum = inv_v_min;
  limits.ranges[0].maximum = inv_v_max;
  for (unsigned i=0; i<count && limits.count; ++i) if (choppers[i].speed) {
    // the period of the chopper is a positive time
    double tau = 1.0 / fabs(choppers[i].speed);
    // the delta time of the chopper is half the time it takes to rotate through the angle of the slit
    double dt = choppers[i].angle / 360.0 / 2.0 * tau;
    // to match McStas DiskChopper, the delay time depends on the absolute value of the speed?
    double t0 = choppers[i].phase / 360.0 / fabs(choppers[i].speed);
    // find the smallest n for which (t0 + dt + n * tau) / d >= inv_v_min
    int n_min = (int) floor((choppers[i].path * inv_v_min - t0 - dt) / tau);
    // find the largest n for which (t0 - dt + n * tau) / d <= inv_v_max
    int n_max = (int) ceil((choppers[i].path * inv_v_max - t0 + dt) / tau);
    // collect the ranges for each of the n in (n_min, n_max) into a set:
    range_set ith;
    ith.count = (unsigned)(n_max - n_min + 1);
    ith.ranges = calloc(ith.count, sizeof(range));
    if (ith.ranges == NULL) {
      printf("Out of memory\n");
      exit(-1);
    }
    unsigned c=0;
    for (unsigned j=0; j < ith.count; ++j) {
      double n_tau = tau * (double) (n_min + (int) j);
      // let the minimum 1/v come from the *end* of the pulse:
      double j_min = (t0 - dt + n_tau - latest_emission) / choppers[i].path;
      double j_max = (t0 + dt + n_tau) / choppers[i].path;
      j_min = j_min < inv_v_min ? inv_v_min : j_min > inv_v_max ? inv_v_max : j_min;
      j_max = j_max < inv_v_min ? inv_v_min : j_max > inv_v_max ? inv_v_max : j_max;
      if (j_min < j_max) {
        ith.ranges[c].minimum = j_min;
        ith.ranges[c++].maximum = j_max;
      }
    }
    ith.count = c;
    // find the intersection of this chopper with the running set
    range_set new_limits = range_intersection(limits, ith);
    // clean-up allocated memory, ensuring limits or ith is not erased if transferred to new_limits:
    if ((!new_limits.count || new_limits.ranges != limits.ranges) && limits.ranges)  free(limits.ranges);
    if ((!new_limits.count || new_limits.ranges != ith.ranges) && ith.ranges) free(ith.ranges);
    // rename the new limits in preparation of returning or the next loop
    limits = new_limits;
  }
  return limits;
}


unsigned chopper_inverse_velocity_limits(double * lower, double * upper,
                                         unsigned count, const chopper_parameters * choppers,
                                         double inv_v_min, double inv_v_max, double latest_emission){
  range_set limits = chopper_inverse_velocity_windows(count, choppers, inv_v_min, inv_v_max, latest_emission);
  if (limits.count){
    *lower = limits.ranges[0].minimum;
    *upper = limits.ranges[limits.count-1].maximum;
  }
  if (limits.ranges) free(limits.ranges);
  return limits.count;
}

// Use the McStas defines if possible, or define them ourselves
#ifndef V2K
#define V2K 1.58825361e-3     /* Convert v[m/s] to k[1/AA] */
#endif
#ifndef K2V
#define K2V 629.622368        /* Convert k[1/AA] to v[m/s] */
#endif
#ifndef PI
#define PI 3.14159265358979323846
#endif
unsigned chopper_wavelength_limits(double * lower, double * upper, unsigned count, const chopper_parameters * choppers,
                                   double lambda_min, double lambda_max, double latest_emission){
  unsigned windows = chopper_inverse_velocity_limits(lower, upper, count, choppers,
                                                     lambda_min * V2K / 2 / PI, lambda_max * V2K / 2 / PI,
                                                     latest_emission);
  if (windows) {
    *lower *= K2V * 2 * PI;
    *upper *= K2V * 2 * PI;
  }
  return windows;
}
