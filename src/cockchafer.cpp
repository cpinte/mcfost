#include "cockchafer.hpp"
#include "iostream"

int predict(char *model_name, float *feature, int nrow, int nfea, const float *output){
  BoosterHandle booster;
  bst_ulong out_len;

  XGBoosterCreate(NULL, 0, &booster) ;
  std::cout << " Trying to read " << model_name << std::endl ;
  if (int err = XGBoosterLoadModel(booster,model_name)) {
    std::cerr << "load model error : " << model_name << std::endl ;
    return err;
  }

  DMatrixHandle input;
  if(int err = XGDMatrixCreateFromMat(feature, nrow, nfea, -1, &input)){
    std::cerr << "XGDMatrixCreateFromMat error" << std::endl;
    return err;
  }

  if(int err = XGBoosterPredict(booster, input, 0, 0, &out_len, &output)){
    std::cerr << "xgb predict error" << std::endl;
    return err;
  }

  XGDMatrixFree(input);
  XGBoosterFree(booster);
  return 0;
}
