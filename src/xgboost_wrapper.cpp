#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>

#include <xgboost/data.h>
#include <xgboost/c_api.h>

extern "C" {
  int predict(char *model_name, float *feature, int nrow, int nfea, float *output){
    BoosterHandle booster;
    bst_ulong out_len;

    XGBoosterCreate(NULL, 0, &booster) ;
    if (int err = XGBoosterLoadModel(booster,model_name)) {
      std::cerr << "load model error : " << model_name << std::endl ;
      return err;
    }

    DMatrixHandle input;
    if(int err = XGDMatrixCreateFromMat(feature, nrow, nfea, -1, &input)){
      std::cerr << "XGDMatrixCreateFromMat error" << std::endl;
      return err;
    }

    const float *outXGB;

    if(int err = XGBoosterPredict(booster, input, 0, 0, 0, &out_len, &outXGB)){
      std::cerr << "xgb predict error" << std::endl;
      return err;
    }

    for(int i=0; i<int(out_len); i++)
      output[i] = outXGB[i];

    XGDMatrixFree(input);
    XGBoosterFree(booster);
    return 0;
  }
}
