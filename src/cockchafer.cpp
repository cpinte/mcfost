#include "cockchafer.hpp"

int predict(char *model, float *feature, int nrow, int nfea, const float **output, long unsigned int *out_len){
  BoosterHandle booster;
  
  XGBoosterCreate(NULL, 0, &booster);
  if(int err = XGBoosterLoadModel(booster,model)){
    std::cerr << "load model error" << std::endl;
    return err;
  }
    
  DMatrixHandle input;
  if(int err = XGDMatrixCreateFromMat(feature, nrow, nfea, -1, &input)){
    std::cerr << "XGDMatrixCreateFromMat error" << std::endl;
    return err;
  }
  
  if(int err = XGBoosterPredict(booster, input, 0, 0, out_len, output)){
    std::cerr << "xgb predict error" << std::endl;
    return err;
  }
  
  XGDMatrixFree(input);
  XGBoosterFree(booster);
  
  return 0;
}
