#include "cockchafer.hpp"
#include "iostream"

int predict(int *model, float *feature, int *nrow, int *nfea, const float *output){
  BoosterHandle booster;
  bst_ulong out_len;

  std::string model_str ;

  if (*model == 1)
    model_str = "model_Tgas.raw" ;

  XGBoosterCreate(NULL, 0, &booster);
  if(int err = XGBoosterLoadModel(booster,model_str.c_str())){
    std::cerr << "load model error : " << model_str << " !" << std::endl;
    return err;
  }

  DMatrixHandle input;
  if(int err = XGDMatrixCreateFromMat(feature, *nrow, *nfea, -1, &input)){
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
