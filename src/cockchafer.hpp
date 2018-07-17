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
    int predict(char *model, float *feature, int nrow, int nfea, const float **output, long unsigned int *out_len);
}
