#include <stdio.h>
#include <stdlib.h>

#include <egads.h>

#include "mpi.h"

void func() {}

enum StatForFunc {
  EG_open_failed = -5,
  EG_loadModel_failed = -4,
  EG_makeTransform_failed = -3,
  EG_copyObject_failed = -2,
  EG_saveModel_failed = -1,
  successful = 0
};

int ConvertStepToEGADS(const char* inputstep,
                       const char* outputegads,
                       double scaling_factor) {
  ego context, model, scaled_model;
  int stat;

  // Initialize EGADS context
  stat = EG_open(&context);
  if ( stat != EGADS_SUCCESS ) {
    printf("EG_open failed with status %d\n", stat);
    return EG_open_failed;
  }

  // Load the STEP file
  stat = EG_loadModel(context, 0, inputstep, &model);
  if ( stat != EGADS_SUCCESS ) {
    printf("EG_loadModel failed with status %d\n", stat);
    EG_close(context);
    return EG_loadModel_failed;
  }

  // Assign scaling factor
  double scale = scaling_factor;

  // clang-format off

  // Construct 4x4 scaling matrix (row-major order)
  double xform[16] = { 
    scale, 0.0, 0.0, 0.0,
    0.0, scale, 0.0, 0.0,
    0.0, 0.0, scale, 0.0,
    0.0, 0.0, 0.0,   1.0 };

  // clang-format on

  // Create transform object
  ego xformObj;
  stat = EG_makeTransform(context, xform, &xformObj);
  if ( stat != EGADS_SUCCESS ) {
    printf("EG_makeTransform failed with status %d\n", stat);
    EG_deleteObject(model);
    EG_close(context);
    return EG_makeTransform_failed;
  }

  // Apply the transform
  stat = EG_copyObject(model, xformObj, &scaled_model);
  if ( stat != EGADS_SUCCESS ) {
    printf("EG_copyObject failed with status %d\n", stat);
    EG_deleteObject(xformObj);
    EG_deleteObject(model);
    EG_close(context);
    return EG_copyObject_failed;
  }

  // Save the scaled model
  stat = EG_saveModel(scaled_model, outputegads);
  if ( stat != EGADS_SUCCESS ) {
    printf("EG_saveModel failed with status %d\n", stat);
    EG_deleteObject(scaled_model);
    EG_deleteObject(xformObj);
    EG_deleteObject(model);
    EG_close(context);
    return EG_saveModel_failed;
  }

  printf(
      "Successfully converted and scaled %s to %s\n", inputstep, outputegads);

  // Clean up
  EG_deleteObject(scaled_model);
  EG_deleteObject(xformObj);
  EG_deleteObject(model);
  EG_close(context);

  return successful;
}

int main(int argc, char* argv[]) {
   if ( argc < 3 ) {
     printf("Usage: %s input.step output.egads [scaling_factor]\n", argv[0]);
     return 1;
   }

	/*
	const char* in =
      "/home/ubuntu/Desktop/RefineWB/ProjectRefine/Assets/oneraM6.step";
  const char* out =
      "/home/ubuntu/Desktop/RefineWB/ProjectRefine/Assets/oneraM6.egads";
	*/

	const char* in = argv[1];
	const char* out = argv[2];

	printf("  Input Step File: %s\n", in);
	printf("Output Egads File: %s\n", out);

  double scaling_factor;
  if ( argc > 3 ) {
    scaling_factor = atof(argv[3]);
  }
  else {
    // scaling_factor = 1.0 / (2.54 * 100.0);
    scaling_factor = 1.0;
  }

  // int stat = ConvertStepToEGADS(argv[1], argv[2], scaling_factor);
  int stat = ConvertStepToEGADS(in, out, scaling_factor);

  if ( stat != successful ) {
    return 1;
  }

  return 0;
}